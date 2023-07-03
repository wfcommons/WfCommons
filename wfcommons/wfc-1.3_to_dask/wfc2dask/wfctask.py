"""
A minimal implementation  of a task according to https://wfcommons.org/format
Hence the name: WFCTask = WorkFlowCommonTask

We don't care about the task execution features here, just the minimum minimorum to build up
the DAG representing the overall workflow and how to execute the task if it's not simulated.

Namely, that means that we only support the following elements in the JSON document:
- name
- command.arguments
- parents
- tasks.[link=input|output,name]
Any other elements can be present but will be ignored

Important note:
- childrens (sic!) element contents are IGNORED
Defining "children" of a task is a violation of the Chain-of-responsibility pattern: When defining a task, I know
what my task depends on, but I cannot anticipate which tasks will depend on it
"""
from __future__ import annotations
import logging


# Logging setup
logger = logging.getLogger(__name__)


class WFCTask:
    def __init__(self):
        self.name = None  # Ignoring type
        self.command = None
        self.parents = set()
        self.inputs = set()
        self.outputs = set()
        pass

    @staticmethod
    def from_json(o_task: dict):
        """
        :param o_task: a task represented by a JSON object/dictionary
        :return: a WFCTask instance
        """
        logger.debug("'%s'" % o_task)
        task = WFCTask()
        task.name = o_task["name"]
        try:
            for parent in o_task["parents"]:
                task.parents.add(parent)
        except KeyError:
            pass  # No parents element for the task. Ignore
        if "childrens" in o_task and len(o_task["childrens"]):
            logger.warning("Non-empty 'childrens' element found in task '%s' is IGNORED" % task.name)
        if "files" in o_task:  # From the sepc, files does not have to be present
            for file in o_task["files"]:
                if file["link"] == "input":
                    task.inputs.add(file["name"])
                elif file["link"] == "output":
                    task.outputs.add(file["name"])
                else:
                    logger.debug(f"{file['link']} not supported: must be either input or output")
        try:
            task.command = o_task["command"]["arguments"]
        except KeyError:
            pass  # If the command is not specified, we will simulate it by creating empty output files
        return task

    @staticmethod
    def load(filename: str) -> tuple[list[WFCTask], str]:
        """
        :param filename: the name of a file containing a representation of a WorkFlowCommon workflow
        :return: the list of tasks in filename
        """
        import json
        with open(filename) as fp:
            return WFCTask.loads(json.load(fp))

    @staticmethod
    def loads(json_workflow: dict) -> tuple[list[WFCTask], str]:
        """
        :param json_workflow: A JSON document describing a workflow
        :return: the list of tasks in the JSON document
        """
        tasks = []
        for o_task in json_workflow["workflow"]["tasks"]:
            tasks.append(WFCTask.from_json(o_task))
        return tasks, json_workflow["name"]
