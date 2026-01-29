import pathlib
import shutil
import tarfile
import os
import io
import sys
import docker
import networkx
from docker.errors import ImageNotFound

from wfcommons.common import Workflow


def _create_fresh_local_dir(path: str) -> pathlib.Path:
    dirpath = pathlib.Path(path)
    if dirpath.exists():
        shutil.rmtree(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    return dirpath

def _remove_local_dir_if_it_exists(path: str) -> None:
    dirpath = pathlib.Path(path)
    if dirpath.exists():
        shutil.rmtree(dirpath)


def _make_tarfile_of_wfcommons():
    source_dir = os.getcwd() # This assumes the testing is run from the root
    tar_stream = io.BytesIO()
    with tarfile.open(fileobj=tar_stream, mode='w') as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    tar_stream.seek(0)
    return tar_stream


def _install_WfCommons_on_container(container):
    # sys.stderr.write("Installing WfCommons on the container...\n")
    # Copy the WfCommons code to it (removing stuff that should be removed)
    target_path = '/tmp/'  # inside container
    tar_data = _make_tarfile_of_wfcommons()
    container.put_archive(target_path, tar_data)
    # Cleanup files from the host
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/build/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/*.egg-info/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark.o", stdout=True,
                                           stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark", stdout=True,
                                           stderr=True)

    # Install WfCommons on the container (to install wfbench and cpu-benchmark really)
    exit_code, output = container.exec_run("sudo python3 -m pip install . --break-system-packages",
                                           workdir="/tmp/WfCommons", stdout=True, stderr=True)
    if exit_code != 0:
        raise RuntimeError("Failed to install WfCommons on the container")

def _start_docker_container(backend, mounted_dir, working_dir, bin_dir, command=None):
    if command is None:
        command = ["sleep", "infinity"]
    # Pulling the Docker image
    client = docker.from_env()
    image_name = f"wfcommons/wfcommons-testing-{backend}"

    try:
        image = client.images.get(image_name)
        sys.stderr.write(f"[{backend}] Image '{image_name}' is available locally\n")
    except ImageNotFound:
        sys.stderr.write(f"[{backend}] Pulling image '{image_name}'...\n")
        client.images.pull(image_name)

    # Launch the docker container to actually run the translated workflow
    sys.stderr.write(f"[{backend}] Starting Docker container...\n")
    container = client.containers.run(
        image=image_name,
        command=command,
        volumes={mounted_dir: {'bind': mounted_dir, 'mode': 'rw'}},
        working_dir=working_dir,
        tty=True,
        detach=True
    )

    # Installing WfCommons on container
    _install_WfCommons_on_container(container)

    # Copy over the wfbench and cpu-benchmark executables to where they should go on the container
    if bin_dir:
        sys.stderr.write(f"[{backend}] Copying wfbench and cpu-benchmark...\n")
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which wfbench` " + bin_dir],
                                               stdout=True, stderr=True)
        if exit_code != 0:
            raise RuntimeError("Failed to copy wfbench script to the bin directory")
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which cpu-benchmark` " + bin_dir],
                                               stdout=True, stderr=True)
        if exit_code != 0:
            raise RuntimeError("Failed to copy cpu-benchmark executable to the bin directory")
    else:
        sys.stderr.write(f"[{backend}] Not Copying wfbench and cpu-benchmark...\n")

    container.backend = backend
    return container

def _shutdown_docker_container_and_remove_image(container):
    image = container.image
    sys.stderr.write(f"[{container.backend}] Terminating container if need be...\n")
    try:
        container.stop()
        container.remove()
    except Exception as e:
        pass

    # Remove the images as we go, if running on GitHub
    if os.getenv('CI') or os.getenv('GITHUB_ACTIONS'):
        sys.stderr.write(f"[{container.backend}] Removing Docker image...\n")
        try:
            image.remove(force=True)
        except Exception as e:
            sys.stderr.write(f"[{container.backend}] Warning: Error while removing image: {e}\n")

def _get_total_size_of_directory(directory_path: str):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            total_size += os.path.getsize(filepath)
    return total_size

def _compare_workflows(workflow_1: Workflow, workflow_2: Workflow):
    
    # Test the number of tasks
    assert (len(workflow_1.tasks) == len(workflow_2.tasks))
    # Test the task graph topology
    assert (networkx.is_isomorphic(workflow_1, workflow_2))
    # Test the total file size sum
    workflow1_input_bytes, workflow2_input_bytes = 0, 0
    workflow1_output_bytes, workflow2_output_bytes = 0, 0
    for workflow1_task, workflow2_task in zip(workflow_1.tasks.values(), workflow_2.tasks.values()):
        # sys.stderr.write(f"WORKFLOW1: {workflow1_task.task_id}  WORKFLOW2 TASK: {workflow2_task.task_id}\n")
        for input_file in workflow1_task.input_files:
            # sys.stderr.write(f"WORKFLOW1 INPUT FILE: {input_file.file_id} {input_file.size}\n")
            workflow1_input_bytes += input_file.size
        for input_file in workflow2_task.input_files:
            # sys.stderr.write(f"WORKFLOW2 INPUT FILE: {input_file.file_id} {input_file.size}\n")
            workflow2_input_bytes += input_file.size
        for output_file in workflow1_task.output_files:
            # sys.stderr.write(f"WORKFLOW1 OUTPUT FILE: {output_file.file_id} {output_file.size}\n")
            workflow1_output_bytes += output_file.size
        for output_file in workflow2_task.output_files:
            # sys.stderr.write(f"WORKFLOW2 OUTPUT FILE: {output_file.file_id} {output_file.size}\n")
            workflow2_output_bytes += output_file.size
    assert (workflow1_input_bytes == workflow2_input_bytes)
    assert (workflow1_output_bytes == workflow2_output_bytes)