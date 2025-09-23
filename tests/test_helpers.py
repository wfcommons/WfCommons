import pathlib
import shutil
import tarfile
import os
import io
import sys
import docker
from docker.errors import ImageNotFound

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
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/build/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/*.egg-info/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark.o", stdout=True,
                                           stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark", stdout=True,
                                           stderr=True)

    # Install WfCommons on the container (to install wfbench and cpu-benchmark really)
    exit_code, output = container.exec_run("sudo python3 -m pip install . --break-system-packages",
                                           workdir="/tmp/WfCommons", stdout=True, stderr=True)

def _start_docker_container(backend, mounted_dir, working_dir, bin_dir, command=None):
    if command is None:
        command = ["sleep", "infinity"]
    # Pulling the Docker image
    client = docker.from_env()
    image_name = f"wfcommons/wfcommons-testing-{backend}"

    try:
        image = client.images.get(image_name)
        sys.stderr.write(f"Image '{image_name}' is available locally\n")
    except ImageNotFound:
        sys.stderr.write(f"Pulling image '{image_name}'...\n")
        client.images.pull(image_name)

    # Launch the docker container to actually run the translated workflow
    sys.stderr.write("Starting Docker container...\n")
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
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which wfbench` " + bin_dir],
                                               stdout=True, stderr=True)
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which cpu-benchmark` " + bin_dir],
                                               stdout=True, stderr=True)

    return container

def _get_total_size_of_directory(directory_path: str):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            total_size += os.path.getsize(filepath)
    return total_size