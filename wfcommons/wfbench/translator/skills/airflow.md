# Skill: Airflow

## Triggers
DAG(, BashOperator, airflow, task_id=, from airflow, PythonOperator, default_args

## Description
Apache Airflow: Python-based workflow orchestration using Directed Acyclic Graphs (DAGs).
DAGs define task dependencies using Python operators and the >> dependency syntax.

## Domain Knowledge
### Airflow Concepts
- **DAG**: A Directed Acyclic Graph defined with `DAG(dag_id, ...)` or `@dag` decorator. Each DAG is a workflow.
- **Operators**: Task definitions. Common types:
  - `BashOperator(task_id=..., bash_command=...)`: shell commands
  - `PythonOperator(task_id=..., python_callable=...)`: Python functions
  - `DockerOperator`, `KubernetesPodOperator`, etc.
- **task_id**: The unique identifier for each task within a DAG.
- **Dependencies**: Set via bitshift operators:
  - `task_a >> task_b` means task_a runs before task_b (task_a is parent)
  - `task_a >> [task_b, task_c]` means task_a is parent of both task_b and task_c
  - `[task_a, task_b] >> task_c` means both are parents of task_c
- **default_args**: Default parameters applied to all tasks in the DAG.

### Mapping to WfFormat
- DAG name (dag_id) -> workflow name
- Each operator instance -> one task in specification
- `task_id` parameter -> task name and id
- `>>` dependency chains -> parent/children relationships
- `bash_command` -> command.program and command.arguments
- `python_callable` -> command.program (function name)
- Tasks with no upstream dependencies -> root tasks (no parents)

### Dependency Inference
- `a >> b` -> a is parent of b, b is child of a
- `a >> [b, c]` -> a is parent of both b and c
- `[a, b] >> c` -> both a and b are parents of c
- Chain: `a >> b >> c` -> a is parent of b, b is parent of c
- Tasks not appearing on the right side of >> are root tasks

## Examples

### Input
```python
from airflow import DAG
from airflow.operators.bash import BashOperator

with DAG("etl_pipeline", schedule_interval="@daily") as dag:
    extract = BashOperator(task_id="extract", bash_command="python extract.py")
    transform = BashOperator(task_id="transform", bash_command="python transform.py")
    load = BashOperator(task_id="load", bash_command="python load.py")

    extract >> transform >> load
```

### Expected Output
```json
{
  "name": "etl_pipeline",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "extract",
          "id": "extract",
          "parents": [],
          "children": ["transform"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "transform",
          "id": "transform",
          "parents": ["extract"],
          "children": ["load"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "load",
          "id": "load",
          "parents": ["transform"],
          "children": [],
          "inputFiles": [],
          "outputFiles": []
        }
      ],
      "files": []
    },
    "execution": {
      "makespanInSeconds": 0,
      "executedAt": "1970-01-01T00:00:00Z",
      "tasks": [
        {
          "id": "extract",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "python", "arguments": ["extract.py"]},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "transform",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "python", "arguments": ["transform.py"]},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "load",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "python", "arguments": ["load.py"]},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        }
      ],
      "machines": [{"nodeName": "unknown"}]
    }
  }
}
```
