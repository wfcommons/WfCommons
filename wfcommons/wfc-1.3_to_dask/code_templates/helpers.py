"""
This is the method that should allow the execution of a task

"""
# Those imports are required when pretending to run the commands
import os
import pathlib
import time
import logging
from workflow_task import WorkflowTask


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def execute_task(task: WorkflowTask, fut_inputs_list) -> WorkflowTask:
    """
    :param task: The task to be executed (it holds all relevant information)
    :param fut_inputs_list: Unused here but necessary for dask to build its own DAG
    :return:
    """
    logger.info("Executing task %s/%s: %s / in=%s / out=%s" % (task.name, task.dag_id, task.command_arguments, task.inputs, task.outputs))
    start = time.time()
    if task.simulate or task.command_arguments is None or len(task.command_arguments) == 0:
        logger.info("Simulating execution of task %s" % task.name)
        # Pretend we do something/Wait some time
        task.simulate_execution()
        for output in task.outputs:
            logger.debug("Simulating %s => %s" % (task.command_arguments, output))
            pathlib.Path(output).touch()
    else:
        command = " ".join(task.command_arguments)
        logger.info("Running command for task %s/%s: %s" % (task.name, task.dag_id, command))
        os.system(command)  # TODO Use subprocess?
    task.execution_time = time.time()-start
    logger.info("End of task %s/%s (%f)" % (task.name, task.dag_id, task.execution_time))
    return task
