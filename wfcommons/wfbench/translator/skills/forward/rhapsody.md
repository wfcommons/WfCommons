# Skill: Rhapsody

## Triggers
rhapsody, rhapsody.api, rhapsody.backends, ComputeTask, ConcurrentExecutionBackend, DragonExecutionBackendV3, DaskExecutionBackend, Session, submit_tasks, wait_tasks

## Description
Generate Rhapsody (radical-cybertools/rhapsody) Python scripts from WfFormat JSON.
Rhapsody is the asyncio-based successor in the RADICAL ecosystem to RADICAL-Pilot's
imperative Session/PilotManager/TaskManager API. It has no built-in DAG/dependencies
field — task ordering is expressed by submitting one topological level at a time and
awaiting it before submitting the next.

## Domain Knowledge

### Rhapsody Concepts
- **Session**: Top-level execution context, constructed with a list of backends.
  Usable as `async with Session([backend]) as session:` or directly with manual
  `await backend.shutdown()`.
- **Backend**: Concrete execution engine. `ConcurrentExecutionBackend` for local
  thread/process pools (laptop-friendly, plain `python` launcher).
  `DragonExecutionBackendV3` for HPC clusters (requires `dragon` launcher and the
  Dragon runtime). `DaskExecutionBackend` for Dask clusters.
- **ComputeTask**: A single task description. Holds an `executable` (path) +
  `arguments` (list[str]), OR a Python `function` + `args`/`kwargs`. Optional
  `input_files` / `output_files` (list[str]), `capture_stdio` (bool), and
  `task_backend_specific_kwargs` for backend-specific config.
- **AITask**: For LLM inference workloads. Not used when translating ordinary
  HPC/scientific WfFormat workflows — skip it unless the WfFormat task names
  clearly describe inference.
- **No declarative DAG**: Rhapsody has no `parents`, `children`, `depends_on`
  field on tasks. Dependencies are expressed by **submitting topological levels
  sequentially**:
    ```
    await session.submit_tasks(level_N_tasks)
    await session.wait_tasks(level_N_tasks)
    # then build level N+1 tasks (whose inputs are now on disk), submit, wait...
    ```
- **No pilot-sandbox staging**: Unlike RADICAL-Pilot's `pilot:///` / `task:///`
  schemes, Rhapsody tasks share the session's working directory (cwd by default).
  File dependencies between tasks are expressed by writing outputs to and reading
  inputs from the **same shared `./data/` directory**. The translator's job is to
  ensure that filenames are consistent across the producer's `output_files` and
  the consumer's `input_files`.
- **All code is async**: Top-level entrypoint must be `async def main()` invoked
  with `asyncio.run(main())`. Backend construction is awaited:
  `backend = await ConcurrentExecutionBackend()`.

### Mapping from WfFormat to Rhapsody

1. **Boilerplate**:
   ```python
   import asyncio
   from rhapsody.api import ComputeTask, Session
   from rhapsody.backends import ConcurrentExecutionBackend
   # Or, for HPC:
   # from rhapsody.backends import DragonExecutionBackendV3
   ```

2. **Topological sort**: From WfFormat `parents`/`children`, compute dependency
   levels. Level 0 = root tasks (no parents). Level N = tasks whose parents are
   all in levels < N.

3. **Backend selection**: Look at the `TARGET OPTIONS` section of the prompt.
   - `backend=concurrent` (default): use `ConcurrentExecutionBackend()` and a
     plain `python script.py` launcher in the script's docstring/comment.
   - `backend=dragon`: use `DragonExecutionBackendV3()` and note the required
     `dragon script.py` launcher.

4. **ComputeTask per task**: For each WfFormat execution task:
   - `uid`: from the WfFormat `id` field (Rhapsody auto-generates if omitted but
     keep the WfFormat IDs for traceability).
   - `executable`: derive from `command.program`. If it's `wfbench`, use
     `"bin/wfbench"`. Otherwise prefer a path under `bin/` (e.g. `"bin/split_fasta"`).
   - `arguments`: from `command.arguments`. Pass as `list[str]`.
   - `input_files`: list of `"data/<file_id>"` strings, one per `inputFile`.
   - `output_files`: list of `"data/<file_id>"` strings, one per `outputFile`.
   - `capture_stdio=True`: include this so stdout/stderr are inspectable on each task.
   - Do NOT invent `cores_per_rank` / `gpus_per_rank` / `memory` top-level
     fields — they don't exist on `ComputeTask`. If `coreCount > 1` and the
     backend is Dragon, attach `task_backend_specific_kwargs={"process_template": {}}`
     (single multi-core process) or `process_templates=[(coreCount, {})]`. For
     the Concurrent backend, ignore `coreCount`.

5. **Submit-wait loop**: For each dependency level:
   ```python
   await session.submit_tasks(level_N_tasks)
   await session.wait_tasks(level_N_tasks)
   ```
   `wait_tasks` returns when all the level's tasks have reached a terminal state.

6. **Cleanup**: Use `async with session:` so backends shut down even on errors.
   Always call `await backend.shutdown()` (or rely on the context manager) before
   exiting `main()`.

### ComputeTask Key Attributes (post-execution, read-only)
- `state`: `"DONE"` / `"FAILED"` / etc.
- `stdout`, `stderr`: strings (populated when `capture_stdio=True`)
- `exit_code`: int (executable tasks)
- `return_value`: the Python return value (function tasks)
- `exception`: error info if `state == "FAILED"`

### What NOT to emit
- No `rp.Session`, `rp.PilotManager`, `rp.TaskManager`, `rp.TaskDescription`,
  `rp.PilotDescription`, `rp.LINK`, `rp.COPY`, `rp.TRANSFER`, `pilot:///`,
  `task:///`, `client:///`. Those are RADICAL-Pilot, not Rhapsody.
- No `depends_on`, `parents`, `children`, `dependencies=` fields on `ComputeTask`.
- No top-level `cores`, `runtime`, `resource` arguments — those belonged to
  `rp.PilotDescription`, not Rhapsody.

## Examples

### Input (WfFormat JSON)
```json
{
  "name": "blast-workflow",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "split",
          "id": "split_00",
          "parents": [],
          "children": ["blast_00", "blast_01"],
          "inputFiles": ["sequences.fasta"],
          "outputFiles": ["chunk_0.fasta", "chunk_1.fasta"]
        },
        {
          "name": "blast",
          "id": "blast_00",
          "parents": ["split_00"],
          "children": ["merge_00"],
          "inputFiles": ["chunk_0.fasta"],
          "outputFiles": ["result_0.xml"]
        },
        {
          "name": "blast",
          "id": "blast_01",
          "parents": ["split_00"],
          "children": ["merge_00"],
          "inputFiles": ["chunk_1.fasta"],
          "outputFiles": ["result_1.xml"]
        },
        {
          "name": "merge",
          "id": "merge_00",
          "parents": ["blast_00", "blast_01"],
          "children": [],
          "inputFiles": ["result_0.xml", "result_1.xml"],
          "outputFiles": ["final_results.xml"]
        }
      ]
    },
    "execution": {
      "tasks": [
        {"id": "split_00", "command": {"program": "wfbench", "arguments": ["--task-name", "split"]}, "coreCount": 1},
        {"id": "blast_00", "command": {"program": "wfbench", "arguments": ["--task-name", "blast"]}, "coreCount": 1},
        {"id": "blast_01", "command": {"program": "wfbench", "arguments": ["--task-name", "blast"]}, "coreCount": 1},
        {"id": "merge_00", "command": {"program": "wfbench", "arguments": ["--task-name", "merge"]}, "coreCount": 1}
      ]
    }
  }
}
```

### Expected Output (backend=concurrent, run with: `python run_workflow.py`)
```python
#!/usr/bin/env python
"""Rhapsody workflow translated from WfFormat. Run: python run_workflow.py"""
import asyncio

from rhapsody.api import ComputeTask, Session
from rhapsody.backends import ConcurrentExecutionBackend


async def main():
    backend = await ConcurrentExecutionBackend()
    session = Session([backend])

    async with session:
        # ---- Level 0: split_00 ----
        split_00 = ComputeTask(
            uid="split_00",
            executable="bin/wfbench",
            arguments=["--task-name", "split"],
            input_files=["data/sequences.fasta"],
            output_files=["data/chunk_0.fasta", "data/chunk_1.fasta"],
            capture_stdio=True,
        )
        level_0 = [split_00]
        await session.submit_tasks(level_0)
        await session.wait_tasks(level_0)

        # ---- Level 1: blast_00, blast_01 ----
        blast_00 = ComputeTask(
            uid="blast_00",
            executable="bin/wfbench",
            arguments=["--task-name", "blast"],
            input_files=["data/chunk_0.fasta"],
            output_files=["data/result_0.xml"],
            capture_stdio=True,
        )
        blast_01 = ComputeTask(
            uid="blast_01",
            executable="bin/wfbench",
            arguments=["--task-name", "blast"],
            input_files=["data/chunk_1.fasta"],
            output_files=["data/result_1.xml"],
            capture_stdio=True,
        )
        level_1 = [blast_00, blast_01]
        await session.submit_tasks(level_1)
        await session.wait_tasks(level_1)

        # ---- Level 2: merge_00 ----
        merge_00 = ComputeTask(
            uid="merge_00",
            executable="bin/wfbench",
            arguments=["--task-name", "merge"],
            input_files=["data/result_0.xml", "data/result_1.xml"],
            output_files=["data/final_results.xml"],
            capture_stdio=True,
        )
        level_2 = [merge_00]
        await session.submit_tasks(level_2)
        await session.wait_tasks(level_2)

        for t in (split_00, blast_00, blast_01, merge_00):
            print(f"{t.uid}: {t.state}")


if __name__ == "__main__":
    asyncio.run(main())
```

### Expected Output (backend=dragon, run with: `dragon run_workflow.py`)
Same structure as above with these substitutions:
```python
from rhapsody.backends import DragonExecutionBackendV3
# ...
backend = await DragonExecutionBackendV3()
# Each ComputeTask gets, when coreCount > 1:
#   task_backend_specific_kwargs={"process_template": {}}
```
And the module docstring should read: `"""... Run: dragon run_workflow.py"""`.

## Validation

### Required elements
- import asyncio
- from rhapsody.api
- ComputeTask
- Session
- submit_tasks
- wait_tasks
- asyncio.run

### Syntax check
command: python -m py_compile {file}
