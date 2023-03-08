"""
Feel free to modify this to target your local dask configuration

Lots of info there:
https://docs.dask.org/en/stable/configuration.html
https://dask.pydata.org/en/latest/scheduling.html
"""
from dask.distributed import Client


def build_dask_client():
    cpu_count = 8        # default is 4 for me
    threads_per_cpu = 4  # default value is 4 for me
    return Client(n_workers=cpu_count, threads_per_worker=threads_per_cpu)
