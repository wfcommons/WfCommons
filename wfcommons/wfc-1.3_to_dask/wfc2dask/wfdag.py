"""
The DAG corresponding to a workflow

Vertices are instances of the inner class WFDAGTask. A WFDAGTask instance is an enriched WFCTask
"""
from __future__ import annotations
import logging
from wfc2dask.wfctask import WFCTask
from code_templates.workflow_task import WorkflowTask
# Logging setup
logger = logging.getLogger(__name__)


class WFDAG:
    class WFDAGTask:
        """
        A vertex of the DAG, i.e. a WFCTask with a private id (dag_id) and the set of its parent private ids (a set
        of dag_id)
        """
        def __init__(self, wfctask: WFCTask, index: int):
            self.wfctask = wfctask  # The underlying WFCTask
            self.dag_id = 'dv_%d' % index  # dv stands for DAG Vertex
            self.dag_parents = set()

    def __init__(self, workflow_name):
        # This is just a placeholder to store each wfctask
        self.workflow_name = workflow_name
        self.wftasks = {}
        # The vertices in the DAG, i.e. WFDAGTask instances
        self.dag_tasks = {}  # Indexed by the WFDAGTask id
        self.ordered_tasks = []

    def add_task(self, task: WFCTask) -> None:
        if task.name in self.wftasks:
            raise Exception("Duplicate task named '%s'" % task.name)
        self.wftasks[task.name] = task

    @staticmethod
    def from_tasks(tasks: list[WFCTask], wfname: str) -> WFDAG:
        wfdag = WFDAG(wfname)
        for wftask in tasks:
            wfdag.add_task(wftask)
        return wfdag

    @staticmethod
    def _resolve_reference(references: dict[str], key: str) -> str:
        """
        Resolve "references"

        :param references: a dictionary of "references", each key is
        + either a user reference (task name, input name or output name) associated to a dag_id value
        + or a dag_id associated to itself
        :param key: the key to resolve
        :return: the resolved reference, i.e. a dag_id

        Notes:
        + This method is expected to raise a KeyError exception if the key is unknown, that is, in the case of
          "external" inputs (an input which is not the output of a task)
        + It does not guarantee the consistency of the DAG, e.g. if a task has a parent task which has not been defined
        """
        _key = key
        already_seen = set()
        while references[_key] != _key:
            if _key in already_seen:
                # Let's prevent cycles (shouldn't happen if the wf is correctly specified)
                raise Exception("Cycle detected for key '%s': %s", key, already_seen)
            already_seen.add(_key)
            _key = references[_key]
        return _key

    def _build_dag_first_pass(self):
        """
        First pass: Collect all dependencies
        Note: The dependencies between two tasks can be defined in parents, in childrens, as a files.input,
        and/or as a files.output. However, there is no guaranteed completeness, e.g. if task1 is in the parents
        of task2, task1 children do not have to mention task2
        The purpose of the references dictionary is to collect the dependencies
        + task.id -> task.id (this is our stop condition when resolving a dependency)
        + task.name -> task.id

        :return: a dictionary of "references", see _resolve_reference for a detailed explanation
        """
        references = {}
        for task_index, wftask in enumerate(self.wftasks.values()):
            task = WFDAG.WFDAGTask(wftask, task_index)
            self.dag_tasks[task.dag_id] = task
            references[task.dag_id] = task.dag_id        # Necessary stop condition when resolving a reference
            references[task.wfctask.name] = task.dag_id  #
            # wfctask.parents is a list of task.names (str). What the task depends on
            for parent in task.wfctask.parents:
                task.dag_parents.add(parent)
            # wfctask inputs. What the task depends on (bis repetita)
            for inp in task.wfctask.inputs:
                # It is not known (yet) if inp is the output of another task or an "external" input
                # For the moment, we just add it to the parents. We will resolve that reference in the second pass
                task.dag_parents.add(inp)
            # wfctask outputs . What the task creates (and other tasks possibly depend on)
            for out in task.wfctask.outputs:
                if out in references:
                    # TODO Make this a specific Exception
                    raise Exception(
                        "'%s' in '%s' is also an output of '%s?!" % (out, task.wfctask.name, references[out]))
                references[out] = task.dag_id  # This output is a direct dependendy on this task
            pass
        return references

    def _build_dag_second_pass(self, references: dict[str]) -> None:
        """
        :param references: the dictionary of references built in the first pass
        :return: None

        The parent dependencies (parent tasks and inputs) of each WFDAGTask have been built in the first pass
        Now we need to consolidate them, that is, when a parent reference is encountered, we need to:
        + Replace it with a dag_id if the reference can be resolved
        + Remove it if the reference cannot be resolved, i.e. in the case of "external" inputs),
        We then need to ensure that the resulting set of dag_ids contains each dag_id only once (this is
        ensured by the use of a Python set)
        """
        for task in self.dag_tasks.values():
            consolidated_parents = set()
            for parent in task.dag_parents:
                try:
                    # Try to find the dag_id associated to parent
                    consolidated_parents.add(WFDAG._resolve_reference(references, parent))
                except KeyError:
                    # Tell the user that the parent is an external input
                    logger.info("In task '%s', '%s' seems to be an external input" % (task.wfctask.name, parent))
            task.dag_parents = consolidated_parents
        pass

    def build_dag(self) -> None:
        """
        This method is called once all tasks have been added through add_task

        TODO Need to check idempotence
        """
        references = self._build_dag_first_pass()
        self._build_dag_second_pass(references)
        self.order_tasks()

    def order_tasks(self) -> None:
        """
        :return: None

        The tasks just need to be ordered in the logical workflow sense:
        + Tasks without parents first,
        + Then tasks that are children of those (and that won't have any other parent)
        + The children of those children,
        + ... and so on

        We iterate over the list of tasks to be ordered. If there is no cycle in the graph,
        each iteration will find one task that can be executed before the others

        Complexity is \thera(n^2) where n = number of tasks.
        TODO Provide faster implementation (but dask is limited in size anyway...)
        """
        tasks_to_order = list(self.dag_tasks.keys())  # That's what we need to order
        for task in self.dag_tasks.values():
            # Backup the set of parent dag_ids. We will need them later
            task.dag_parents_backup = set(task.dag_parents)
        while len(tasks_to_order) != 0:
            logger.debug("Bef: %d tasks with ordered size %d" % (len(tasks_to_order), len(self.ordered_tasks)))
            tasks_without_parents = []
            tasks_with_parents = []
            # Iterate over the set of tasks to find those without parents (i.e. that can be executed at this stage)
            for task_id in tasks_to_order:
                if len(self.dag_tasks[task_id].dag_parents) == 0:
                    tasks_without_parents.append(task_id)
                else:
                    tasks_with_parents.append(task_id)
            # Check if there is a cycle
            # No need to check if len(tasks_with_parents) != 0 since tasks from tasks_to_order are
            # either with or without parents and len(tasks_to_order) is != 0
            if len(tasks_without_parents) == 0:
                raise Exception("Cycle detected in workflow")
            self.ordered_tasks.append(tasks_without_parents)
            # Now remove the dag_ids of tasks_without_parents from the set of parents of tasks of tasks_with_parents
            # (since they will be executed before)
            for task_id in tasks_with_parents:
                for _id in tasks_without_parents:
                    if _id in self.dag_tasks[task_id].dag_parents:
                        self.dag_tasks[task_id].dag_parents.remove(_id)
            # What remains to order in the next iteration is the set of tasks_with_parents
            tasks_to_order = tasks_with_parents
            logger.debug("Aft: %d tasks with ordered size %d" % (len(tasks_to_order), len(self.ordered_tasks)))
        # Restore the deleted parents and delete the backup
        for task in self.dag_tasks.values():
            task.dag_parents = set(task.dag_parents_backup)
            del task.dag_parents_backup

    def __repr__(self) -> str:
        """
        :return: a representation of the DAG by level
        """
        if len(self.ordered_tasks) == 0:
            self.build_dag()
        rep = '\n'
        for level, tasks in enumerate(self.ordered_tasks):
            rep += "Level %d: " % level
            rep += "; ".join(['%s (%s)' % (_id, self.dag_tasks[_id].wfctask.name) for _id in tasks])
            rep += "\n"
        return rep

    def dask_wftasks_codelines(self, randomizer_varname: str) -> list[str]:
        """
        Build the code definining all tasks in the workflow, i.e. WorkflowTask instances
        :param randomizer_varname: the name of the randomizer
        :return: the non-indented Python lines of code used to instantiate the WorkflowTask instances
        """
        codelines = ["randomizer = random.Random(seed)",
                     "TASKS = {}"]
        for dag_task in self.dag_tasks.values():
            _workflow_task = WorkflowTask(dag_id=dag_task.dag_id,
                                          name=dag_task.wfctask.name,
                                          command_arguments=dag_task.wfctask.command,
                                          inputs=list(dag_task.wfctask.inputs),
                                          outputs=list(dag_task.wfctask.outputs)
                                          )
            code = _workflow_task.pythonize(randomizer_varname)
            codelines.append("TASKS['%s'] = %s" % (dag_task.dag_id, code[0]))
            codelines.extend([codeline for codeline in code[1:]])
        return codelines

    def dask_codelines(self) -> list[str]:
        self.build_dag()
        logger.debug('%s' % self)  # Display the DAG
        noindent_python_codelines = self.dask_wftasks_codelines("randomizer")
        # client.submit() lines
        for level, task_ids in enumerate(self.ordered_tasks):
            noindent_python_codelines.append("# Level %d (%d tasks)" % (level + 1, len(task_ids)))
            for task_id in task_ids:
                task = self.dag_tasks[task_id]
                noindent_python_codelines.append("# Task %s (%s)" % (task_id, task.wfctask.name))
                fut_task_varname = "fut_%s" % task_id
                fut_inputs_list = ", ".join(["fut_%s" % parent_id for parent_id in task.dag_parents])
                codeline = "%s = client.submit(execute_task, TASKS['%s'], [%s])" % (fut_task_varname,
                                                                                    task_id,
                                                                                    fut_inputs_list)
                noindent_python_codelines.append(codeline)
        # future.result() lines in reverse order
        # We call result() for each future: This guarantees that all tasks will be executed (and we
        # delegate to dask the task to find leaves in the DAG at various levels)
        # Notes:
        # + if the future is Finished, this has no overhead
        # + there cannot be deadlocks (or the previous client.submit() lines would be syntactically incorrect: a future
        # would be used before its creation)
        for level, task_ids in enumerate(reversed(self.ordered_tasks)):
            noindent_python_codelines.append("# Level %d (%d tasks)" % (len(self.ordered_tasks)-level, len(task_ids)))
            for task_id in reversed(task_ids):
                task = self.dag_tasks[task_id]
                noindent_python_codelines.append("# Task %s (%s)" % (task_id, task.wfctask.name))
                fut_task_varname = "fut_%s" % task_id
                codeline = "TASKS['%s'] = %s.result()" % (task_id, fut_task_varname)
                noindent_python_codelines.append(codeline)
        return noindent_python_codelines
