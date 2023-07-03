"""

A class used to define and help execute a task in a dask workflow

You don't need to modify this unless you know what you're doing

TODO The JSON should actually be a WFCTask
"""
import random
import time


class WorkflowTask:
    def __init__(self,
                 dag_id: str = None,
                 name: str = None,
                 command_arguments: list[str] = None,
                 inputs: list[str] = None,
                 outputs: list[str] = None,
                 simulate: bool = False,
                 randomizer: random.Random = random.Random(),
                 simulate_minimum_execution_time: float = 0.1,
                 simulate_maximum_execution_time: float = 1.1,
                 execution_time: float = None,  # This is an execution output
                 ):
        self.dag_id = dag_id
        self.name = name
        self.command_arguments = command_arguments
        self.inputs = inputs
        self.outputs = outputs
        self.simulate = simulate
        self.randomizer = randomizer
        self.simulate_minimum_execution_time = simulate_minimum_execution_time
        self.simulate_maximum_execution_time = simulate_maximum_execution_time
        self.execution_time = execution_time

    def simulate_execution(self):
        time.sleep(self.randomizer.uniform(self.simulate_minimum_execution_time,
                                           self.simulate_maximum_execution_time))

    def pythonize(self, randomizer_varname: str = "randomizer") -> list[str]:
        codelines = ["WorkflowTask(dag_id = '%s'," % self.dag_id,
                     "             name = '%s'," % self.name,
                     "             command_arguments = %s," % self.command_arguments,
                     "             inputs = %s," % self.inputs,
                     "             outputs = %s," % self.outputs,
                     "             simulate = simulate,",
                     "             randomizer = %s," % randomizer_varname,
                     "             simulate_minimum_execution_time = %s," % self.simulate_minimum_execution_time,
                     "             simulate_maximum_execution_time = %s," % self.simulate_maximum_execution_time,
                     "             )"]
        return codelines
