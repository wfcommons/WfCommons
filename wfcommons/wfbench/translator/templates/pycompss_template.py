# IMPORTS
import sys
import os
# PyCOMPSs imports
# All task arguments: https://compss-doc.readthedocs.io/en/stable/Sections/02_App_Development/02_Python/01_1_Task_definition/Sections/04_Task_parameters_summary.html
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.binary import binary
from pycompss.api.mpi import mpi
from pycompss.api.parameter import *
# All API functions: https://compss-doc.readthedocs.io/en/stable/Sections/02_App_Development/02_Python/01_2_Synchronization/01_API.html?highlight=compss_open#api-summary
from pycompss.api.api import compss_open
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier

# @binary
# @constraint(computing_units=24)
# @mpi(runner="mpirun", binary="gmx_mpi", computing_nodes=1)
# @task
# def task_id():
#     pass
#
# def main_program(arg1, arg2):
#     # Execute task
#     task_id()
#
# if __name__ == "__main__":
#     arg1 = sys.argv[1]
#     arg2 = sys.argv[2]
#     main_program(arg1, arg2)

# Generated code goes here
