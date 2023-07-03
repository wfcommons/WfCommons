# Objective
wfc2dask.py is a tool that attempts to translate workflows specified using the WFCommons workflow format (version 1.3; 
see https://wfcommons.org/format or https://github.com/wfcommons/wfformat) to a workflow actionnable through dask

# TL;DR
* You need to install dask (see https://docs.dask.org/en/stable/install.html that will tell you to run `conda install
dask` or `pip install 'dask[distributed]'`).  //TODO Check if this works if dask is not installed//
* `python wfc2dask.py -h` to get help
* `python wfc2dask.py <workflow_filename>` creates a directory (named out by default) where the Python code needed
  to execute the workflow is stored
* `cd out` and then `python application.py` to run the workflow in a local dask

# The Gory Details
## Configuration
* The dask client is defined in `/out/dask_client.py`. Feel free to modify to your local configuration.
* The workflow tasks and dask DAG are defined in `/out/run_workflow.py`. Feel free to modify as well but take into 
  account the comments of the `__doc__` of that file
* That's it for the configuration. You don't need to change the other files

## Implementation
`/wfc2dask.py` contains the main. The WFCommons JSON document is digested using `/wfc2dask/wfctask.py` and its tasks 
are grouped in a `WFDag` defined in `/wfc2dask/wfdag.py`. Once all tasks have been ingested, the DAG is built and Python
code is written to the output directory.

## Contents of the 'samples' directory 
### Contents of the 'samples/unittests' directory

Files necessary for the proper execution of unit tests

* hello-world-sequence.json:
The first task (hello) creates a file whose contents are used in the second task (world). 
The DAG is the following:
```
(init) -> hello -> world -> (end)
```

* hello-world-join.json:
Three tasks in this workflow. The 'hello' and 'world' tasks each create a different file. Those files
are then used in the third 'join' task
The DAG is the following:
```
      /-> hello -\
(init)            -> join -> (end)
      \-> world -/
```

## Unit tests/Coverage
### Unit tests
```commandline
python -m unittest discover -s tests
````

## Coverage
```
conda install coverage  # or pip install coverage
coverage run --source=. --omit='code_templates/*' -m unittest discover -s tests/ && coverage html
```
Open htmlcov/index.html
