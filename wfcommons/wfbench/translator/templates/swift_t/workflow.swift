import files;
import io;
import python;
import string;
import unix;

global const string flowcept_start =
"""
import time
workflow_id = "%s"
time.sleep(30)
""";

global const string flowcept =
"""
import logging
import pathlib
import time
from flowcept.flowcept_api.flowcept_controller import Flowcept

logging.basicConfig(
    level=logging.INFO,
    format="[WfBench][%%(asctime)s][%%(levelname)s] %%(message)s",
    datefmt="%%H:%%M:%%S",
    handlers=[logging.StreamHandler()]
)

workflow_id = "%s"
workflow_name = "%s"
out_files = [%s]

__import__("logging").info("Flowcept Starting")
flowcept_agent = Flowcept(workflow_id=workflow_id, workflow_name=workflow_name, bundle_exec_id=workflow_id, start_persistence=False, save_workflow=True)

try:
    flowcept_agent.start()
except Exception:
    import traceback
    traceback.print_exc()

remaining_files = set(out_files)

while remaining_files:
    found_files = set()
    for f in remaining_files:
        if pathlib.Path(f).exists():
            found_files.add(f)
    remaining_files -= found_files
    if not remaining_files:
        break
    time.sleep(1)
    
try:
    flowcept_agent.stop()
except Exception:
    import traceback
    traceback.print_exc()

__import__("logging").info("Flowcept Completed")
""";

string command = 
"""
import logging
import os
import sys
import pathlib
import signal
import socket
import subprocess
import time

__import__("logging").basicConfig(
    level=logging.INFO,
    format="[WfBench][%%(asctime)s][%%(levelname)s] %%(message)s",
    datefmt="%%H:%%M:%%S",
    handlers=[logging.StreamHandler()]
)

wfbench = "%s"
task_name = "%s"
input_file = ["%s"]
gpu_work = int(%i)
cpu_work = int(%i)
percent_cpu = %f
cpu_threads = int(10 * percent_cpu)
output_file = "%s"
output_file_size = int(%i)
dep = %i
workflow_id = "%s"
task_id = f"{workflow_id}_{task_name}"

if 'workflow_id':
    __import__("logging").info("Running with Flowcept.")
    from flowcept import Flowcept, FlowceptTask
    fc = Flowcept(workflow_id=workflow_id,
                bundle_exec_id=workflow_id,
                start_persistence=False, save_workflow=False)
    fc.start()
    fc_task = FlowceptTask(workflow_id=workflow_id, task_id=task_id, used={
      'workflow_id': workflow_id,
      'name': task_name,
      'percent-cpu': percent_cpu,
      'cpu-work': cpu_work,
      'gpu-work': gpu_work
    })

__import__("logging").info(f"Starting {task_name} Benchmark on {socket.gethostname()}")

from importlib.util import spec_from_file_location, module_from_spec
from importlib.machinery import SourceFileLoader
loader = SourceFileLoader("wfbench", wfbench)
spec = spec_from_file_location("wfbench", wfbench, loader=loader)
mod = module_from_spec(spec)
spec.loader.exec_module(mod)
mod.run(
    name=task_name,
    workflow_id=workflow_id,
    percent_cpu=percent_cpu,
    cpu_work=cpu_work,
    mem=None,
    gpu_work=gpu_work,
    output_files=f'{{"{output_file}": {output_file_size}}}',
    input_files=str(input_file).replace("'", '"'),
    with_flowcept=True,
    silent=False,
    debug=False,
    rundir=None,
    path_lock=None,
    path_cores=None,
    time_limit=None,
    num_chunks=10,
)

__import__("logging").info(f"Benchmark {task_name} completed!")

if 'workflow_id':
    fc_task.end()
    fc.stop()
""";

# Generated code goes here