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

logging.info("Flowcept Starting")
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

logging.info("Flowcept Completed")
""";

string command = 
"""
import logging
import os
import pathlib
import signal
import socket
import subprocess
import time
from pathos.helpers import mp as multiprocessing

logging.basicConfig(
    level=logging.INFO,
    format="[WfBench][%%(asctime)s][%%(levelname)s] %%(message)s",
    datefmt="%%H:%%M:%%S",
    handlers=[logging.StreamHandler()]
)

cpu_benchmark = "%s"
task_name = "%s"
files_list = ["%s"]
gpu_work = int(%i)
cpu_work = int(%i)
percent_cpu = %f
cpu_threads = int(10 * percent_cpu)
output_data = {"%s": int(%i)}
dep = %i
workflow_id = "%s"
task_id = f"{workflow_id}_{task_name}"

if 'workflow_id':
    logging.info("Running with Flowcept.")
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

logging.info(f"Starting {task_name} Benchmark on {socket.gethostname()}")

procs = []
cpu_queue = multiprocessing.Queue()
logging.debug(f"Working directory: {os.getcwd()}")

logging.debug("Starting IO benchmark...")
io_proc = None
termination_event = multiprocessing.Event()

io_proc = multiprocessing.Process(
    target=lambda inputs=files_list, outputs=output_data, cpu_queue=cpu_queue, 
           termination_event=termination_event: (
        memory_limit := 10 * 1024 * 1024,
        [open(name, "wb").close() for name in outputs],
        io_completed := 0,
        bytes_read := {name: 0 for name in inputs},
        bytes_written := {name: 0 for name in outputs},
        input_sizes := {name: __import__("os").path.getsize(name) for name in inputs},
        [
            (
                cpu_percent := cpu_queue.get(timeout=1.0),                
                should_exit := termination_event.is_set(),
                (
                    while_loop_var := True,
                    [
                        (
                            new_val := (
                                cpu_queue.get(timeout = 1.0)
                                if not cpu_queue.empty() else None
                            ),
                            cpu_percent := (
                                max(cpu_percent, new_val) 
                                if new_val is not None else cpu_percent
                            ),
                            while_loop_var := (
                                new_val is not None and not cpu_queue.empty()
                            )
                        )
                        for _ in range(100) if while_loop_var
                    ],
                    bytes_to_read := {
                        name: max(0, int(size * (cpu_percent / 100) - bytes_read[name]))
                        for name, size in input_sizes.items()
                    },
                    bytes_to_write := {
                        name: max(0, int(size * (cpu_percent / 100) - bytes_written[name]))
                        for name, size in outputs.items()
                    },
                    logging.debug("Starting IO Read Benchmark..."),
                    in_file := list(bytes_to_read.keys())[0],
                    in_size := list(bytes_to_read.values())[0],
                    open(in_file, "rb").read(int(in_size)),
                    logging.debug("Completed IO Read Benchmark!"),
                    out_file := list(output_data.keys())[0],
                    out_size := list(output_data.values())[0],
                    logging.debug(f"Writing output file '{out_file}'"),
                    open(out_file, "ab").write(__import__("os").urandom(int(out_size))),
                    bytes_read.update({
                        name: bytes_read[name] + bytes_to_read[name]
                        for name in bytes_to_read
                    }),
                    bytes_written.update({
                        name: bytes_written[name] + bytes_to_write[name]
                        for name in bytes_to_write
                    }),
                    
                    logging.debug(f"Bytes Read: {bytes_read}"),
                    logging.debug(f"Bytes Written: {bytes_written}"),
                    io_completed := cpu_percent,
                ) if cpu_percent is not None else time.sleep(0.1),
                not (should_exit or io_completed >= 100)
            )
            for _ in range(1000000)
            if not (io_completed >= 100 or termination_event.is_set())
        ],
        logging.info("IO benchmark completed")
    )
)
io_proc.start()
procs.append(io_proc)

if cpu_work > 0:
    logging.info(f"Starting CPU and Memory Benchmarks for {task_name}...")

    mem_threads = 10 - cpu_threads
    cpu_work_per_thread = int(cpu_work / cpu_threads)

    cpu_procs = []
    mem_procs = []
    cpu_prog = [f"{cpu_benchmark}", f"{cpu_work_per_thread}"]
    mem_prog = ["stress-ng", "--vm", f"{mem_threads}",
                "--vm-bytes", "0.05%%", "--vm-keep"]

    for i in range(cpu_threads):
        cpu_proc = subprocess.Popen(cpu_prog, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cpu_procs.append(cpu_proc)
        monitor_thread = multiprocessing.Process(
            target=lambda proc=cpu_proc, queue=cpu_queue: 
                [
                    queue.put(float(line.strip().split()[1].strip('%%')))
                    for line in iter(proc.stdout.readline, "") 
                    if line.strip() and line.strip().startswith("Progress:")
                ]
        )
        monitor_thread.start()

    if mem_threads > 0:
        mem_proc = subprocess.Popen(mem_prog, preexec_fn=os.setsid)
        mem_procs.append(mem_proc)

    procs.extend(cpu_procs)
    for proc in procs:
        if isinstance(proc, subprocess.Popen):
            proc.wait()
    if io_proc is not None and io_proc.is_alive():
        io_proc.join()

    for mem_proc in mem_procs:
        try:
            os.kill(mem_proc.pid, signal.SIGKILL)
        except subprocess.TimeoutExpired:
            logging.debug("Memory process did not terminate; force-killing.")
    subprocess.Popen(["pkill", "-f", "stress-ng"]).wait()

    logging.debug("Completed CPU and Memory Benchmarks!")
    
logging.info(f"Benchmark {task_name} completed!")

if 'workflow_id':
    fc_task.end()
    fc.stop()
""";

# Generated code goes here