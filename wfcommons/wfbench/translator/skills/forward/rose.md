# Skill: Rose

## Triggers
rose, rose.al, rose.metrics, SequentialActiveLearner, WorkflowEngine, radical.asyncflow, asyncflow, executable_task, function_task, RadicalExecutionBackend, LocalExecutionBackend

## Description
Generate ROSE (radical-cybertools/ROSE) Python scripts from WfFormat JSON.
ROSE ("RADICAL Optimal & Smart-Surrogate Explorer") is a RADICAL-ecosystem
toolkit for expressing and executing workflows on HPC. It is built on
`radical.asyncflow`'s `WorkflowEngine`, which runs tasks over a pluggable
execution backend. Unlike RADICAL-Pilot's imperative submit-wait loop or
Rhapsody's per-level `submit_tasks`/`wait_tasks`, ROSE/asyncflow expresses the
DAG **declaratively through the Python call graph**: a task is invoked to get a
future, and a dependency is declared by passing that future as an argument to a
downstream task. AsyncFlow schedules independent tasks concurrently and defers
each dependent task until its inputs are ready.

## Domain Knowledge

### ROSE / AsyncFlow Concepts
- **Execution backend**: The engine that actually runs tasks.
  - `RadicalExecutionBackend({'runtime': 30, 'resource': 'local.localhost'})`
    (imported from `rhapsody.backends`) runs on RADICAL resources — HPC clusters
    or `local.localhost`. This is ROSE's default backend.
  - `LocalExecutionBackend(ThreadPoolExecutor(max_workers=N))` (imported from
    `radical.asyncflow`) runs tasks in a local thread pool — laptop-friendly, no
    RADICAL stack required, and all tasks share the process working directory.
  - Backend construction is awaited: `backend = await RadicalExecutionBackend(...)`.
- **WorkflowEngine**: The orchestrator. Created with
  `flow = await WorkflowEngine.create(backend=backend)`. It owns the task
  decorators (`flow.executable_task`, `flow.function_task`) and the scheduler.
- **executable_task**: Decorator for a task that runs a shell command. The
  decorated `async` function **returns the command line as a single string**
  (e.g. `return "bin/wfbench --name blast_00 --cpu-work 500"`). Calling the
  decorated function returns a *future*, not the string.
- **function_task**: Decorator for a task that runs a Python callable and returns
  a Python value. Not needed when translating ordinary HPC/scientific WfFormat
  workflows whose tasks are executables — use `executable_task`. (function_task
  is what ROSE's active-learning learners use to evaluate models; skip it unless
  the WfFormat tasks clearly describe in-process Python computation.)
- **Dependencies are the call graph**: There is **no** `depends_on`, `parents`,
  or `children` field. To make task B depend on task A, pass A's future into B:
    ```python
    a_fut = task_a()            # future for A
    b_fut = task_b(a_fut)       # B waits for A; a_fut passed as a dependency
    await b_fut                 # awaiting a leaf drives the whole sub-graph
    ```
  Tasks whose futures are never passed to each other run concurrently.
- **Active-learning API is NOT the target here**: ROSE also ships
  `SequentialActiveLearner` with `@acl.simulation_task` / `@acl.training_task` /
  `@acl.active_learn_task` decorators and `await acl.start(max_iter=...)`. That
  is for iterative surrogate-training loops, **not** for a general WfFormat DAG.
  Do NOT emit `SequentialActiveLearner`, `acl.start`, or the learner decorators
  when translating an ordinary WfFormat workflow — use the plain
  `WorkflowEngine` + `executable_task` + future-passing pattern.
- **Shared working directory**: Tasks read/write files relative to the working
  directory. wfbench reads/writes the files named in its `--input-files` /
  `--output-files` arguments relative to its cwd, so use bare file IDs (no
  directory prefix). The producer's output filename and the consumer's input
  filename must match exactly. Run the generated script from a directory
  containing `bin/wfbench` (or a symlink to it).
- **All code is async**: Entry point is `async def main()` invoked with
  `asyncio.run(main())`.

### Mapping from WfFormat to ROSE

1. **Boilerplate** (backend=radical, the default):
   ```python
   import asyncio
   from radical.asyncflow import WorkflowEngine
   from rhapsody.backends import RadicalExecutionBackend
   ```
   For backend=local:
   ```python
   import asyncio
   from concurrent.futures import ThreadPoolExecutor
   from radical.asyncflow import WorkflowEngine, LocalExecutionBackend
   ```

2. **Backend selection**: Look at the `TARGET OPTIONS` section of the prompt.
   - `backend=radical` (default): construct
     `await RadicalExecutionBackend({'runtime': 30, 'resource': 'local.localhost'})`.
     Note the launcher in the docstring: `python run_workflow.py` (a RADICAL
     stack must be installed).
   - `backend=local`: construct
     `await LocalExecutionBackend(ThreadPoolExecutor(max_workers=4))`.

3. **One decorated task function per workflow task.** Emit an explicit
   `@flow.executable_task`-decorated `async def <task_id>(*deps):` for **every**
   WfFormat task. Do NOT collapse repeated tasks into a `for i in range(...)`
   loop, list comprehension, or generator — the workflow's task IDs, file IDs,
   and argument values are not always strictly sequential, and looping silently
   drops or duplicates tasks when the IDs skip or the arguments differ per task.
   If the source workflow has 195 blastall tasks, the generated script must
   contain 195 distinct decorated task functions. Long files are acceptable;
   correctness over brevity. Give each function a `*deps` parameter so it can
   absorb parent futures (which express ordering) even though the command string
   ignores them.

4. **Command string per task**: For each WfFormat execution task the decorated
   function returns a single shell command string built from `command.program`
   and `command.arguments`:
   - **program**: if it's `wfbench`, use `"bin/wfbench"`; otherwise prefer a
     path under `bin/` (e.g. `"bin/split_fasta"`).
   - **arguments**: join `command.arguments` after the program with single
     spaces. Because `executable_task` returns a *shell* command line (not a
     tokenized `argv` list), the already-joined `"--flag value"` strings from
     WfFormat are fine as-is — do NOT tokenize them the way the Rhapsody skill
     does. Just concatenate: `"bin/wfbench " + " ".join(arguments)`.
   - **JSON quoting in --output-files / --input-files**: WfFormat encodes these
     with Python repr (single quotes: `{'foo': 502513}`), but wfbench's
     argparser runs them through `json.loads()`, which needs strict JSON (double
     quotes) — and because this is a shell command line, the value must be a
     single shell token. So wrap the value in single quotes and use double
     quotes inside:
     - `--output-files {'foo': 502513}` → `--output-files '{"foo": 502513}'`
     - `--input-files ['bar']`          → `--input-files '["bar"]'`

5. **Wire the DAG with futures**: After defining the task functions, build the
   call graph. For each task, call its function passing the futures of all its
   parent tasks as positional arguments, and keep the returned future in a dict
   keyed by task id:
   ```python
   futures = {}
   futures["split_00"] = split_00()
   futures["blast_00"] = blast_00(futures["split_00"])
   futures["merge_00"] = merge_00(futures["blast_00"], futures["blast_01"])
   ```
   Emit the calls in topological order so every parent future exists before it
   is referenced.

6. **Drive and collect**: Await the leaf futures (tasks with no children) to run
   the whole graph — `await asyncio.gather(*leaf_futures)`. Awaiting a leaf
   pulls in all of its transitive dependencies.

7. **Shutdown**: Always `await flow.shutdown()` at the end of `main()` (in a
   `try/finally` so the backend is released even on error).

### What NOT to emit
- No `rp.Session`, `rp.PilotManager`, `rp.TaskManager`, `rp.TaskDescription`,
  `pilot:///`, `task:///`, `client:///` — that is RADICAL-Pilot.
- No `ComputeTask`, `Session([backend])`, `session.submit_tasks`,
  `session.wait_tasks` — that is Rhapsody.
- No `SequentialActiveLearner`, `acl.start`, `@acl.simulation_task`,
  `@acl.training_task`, `@acl.active_learn_task`, `rose.metrics` — that is
  ROSE's active-learning layer, not a general DAG.
- No `depends_on`, `parents`, `children`, or `dependencies=` fields. Dependencies
  are expressed only by passing futures as arguments.
- No `for i in range(...)` loops, list comprehensions, or generator expressions
  to define task functions in bulk. Each task is a distinct decorated function.
- No tokenized `argv` list for the command — `executable_task` returns one shell
  command string.

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
        {"id": "split_00", "command": {"program": "wfbench", "arguments": ["--name split_00", "--cpu-work 100"]}, "coreCount": 1},
        {"id": "blast_00", "command": {"program": "wfbench", "arguments": ["--name blast_00", "--cpu-work 500"]}, "coreCount": 1},
        {"id": "blast_01", "command": {"program": "wfbench", "arguments": ["--name blast_01", "--cpu-work 500"]}, "coreCount": 1},
        {"id": "merge_00", "command": {"program": "wfbench", "arguments": ["--name merge_00", "--cpu-work 50"]}, "coreCount": 1}
      ]
    }
  }
}
```

### Expected Output (backend=radical, run with: `python run_workflow.py`)
```python
#!/usr/bin/env python
"""ROSE / radical.asyncflow workflow translated from WfFormat. Run: python run_workflow.py"""
import asyncio

from radical.asyncflow import WorkflowEngine
from rhapsody.backends import RadicalExecutionBackend


async def main():
    backend = await RadicalExecutionBackend(
        {"runtime": 30, "resource": "local.localhost"}
    )
    flow = await WorkflowEngine.create(backend=backend)

    # ---- Task definitions (one executable_task per workflow task) ----
    @flow.executable_task
    async def split_00(*deps):
        return "bin/wfbench --name split_00 --cpu-work 100"

    @flow.executable_task
    async def blast_00(*deps):
        return "bin/wfbench --name blast_00 --cpu-work 500"

    @flow.executable_task
    async def blast_01(*deps):
        return "bin/wfbench --name blast_01 --cpu-work 500"

    @flow.executable_task
    async def merge_00(*deps):
        return "bin/wfbench --name merge_00 --cpu-work 50"

    try:
        # ---- DAG wiring: dependencies are futures passed as arguments ----
        futures = {}
        futures["split_00"] = split_00()
        futures["blast_00"] = blast_00(futures["split_00"])
        futures["blast_01"] = blast_01(futures["split_00"])
        futures["merge_00"] = merge_00(futures["blast_00"], futures["blast_01"])

        # ---- Drive the graph by awaiting the leaf tasks ----
        await asyncio.gather(futures["merge_00"])
    finally:
        await flow.shutdown()


if __name__ == "__main__":
    asyncio.run(main())
```

### Expected Output (backend=local, run with: `python run_workflow.py`)
Same structure with these substitutions:
```python
from concurrent.futures import ThreadPoolExecutor
from radical.asyncflow import WorkflowEngine, LocalExecutionBackend
# ...
backend = await LocalExecutionBackend(ThreadPoolExecutor(max_workers=4))
flow = await WorkflowEngine.create(backend=backend)
```

## Validation

### Required elements
- import asyncio
- from radical.asyncflow import
- WorkflowEngine
- executable_task
- flow.shutdown
- asyncio.run

### Syntax check
command: python -m py_compile {file}
