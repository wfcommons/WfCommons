# Skill: Radical Pilot

## Triggers
radical, radical.pilot, rp.Session, rp.TaskDescription, rp.PilotDescription, radical_pilot

## Description
Generate RADICAL-Pilot (RP) Python scripts from WfFormat JSON.
RP uses procedural submit-wait-submit for task dependencies (no declarative DAG).

## Domain Knowledge
### RADICAL-Pilot Concepts
- **Session**: Top-level RP context. Created with `rp.Session()`, must be closed at the end.
- **PilotManager / PilotDescription**: Manages compute pilots. A pilot reserves resources on a target machine.
- **TaskManager**: Submits tasks to pilots. Tasks are `rp.TaskDescription` dicts.
- **TaskDescription**: Defines a single task — executable, arguments, staging, resource requirements.
- **No declarative DAG**: RP has NO `depends_on`, `parents`, or `children` field. Dependencies are
  expressed by grouping tasks into topological levels and calling `tmgr.submit_tasks()` then
  `tmgr.wait_tasks()` for each level before submitting the next.
- **Staging directives**: Files are moved between tasks via the **pilot sandbox** (`pilot:///`).
  A task copies output to `pilot:///`, the next task links/copies from there.

### Mapping from WfFormat

1. **Boilerplate**: Always generate the session/pilot/task manager setup:
   ```python
   import radical.pilot as rp

   session = rp.Session()
   pmgr    = rp.PilotManager(session=session)
   tmgr    = rp.TaskManager(session=session)

   pilot = pmgr.submit_pilots(rp.PilotDescription({
       'resource': 'local.localhost',
       'cores'   : 4,
       'runtime' : 30,
       'exit_on_error': False
   }))
   tmgr.add_pilots(pilot)
   pilot.wait(rp.PMGR_ACTIVE)
   ```

2. **Topological sort**: From WfFormat `parents`/`children`, compute dependency levels.
   Level 0 = root tasks (no parents). Level N = tasks whose parents are all in levels < N.

3. **TaskDescription per task**:
   - `uid`: task ID from WfFormat `id` field
   - `executable`: from `command.program` in execution tasks. If the task uses `wfbench`,
     use the path `bin/wfbench` or `bin/<program>`.
   - `arguments`: from `command.arguments` in execution tasks
   - `cores_per_rank`: from `coreCount` (default 1)
   - `input_staging`: for each `inputFile`, if the file is produced by a parent task,
     stage from `pilot:///<parent_task_id>_<file_id>` to `task:///<file_id>` using `rp.LINK`.
     For workflow-level input files (no parent produces them), stage from `client:///data/<file_id>`.
   - `output_staging`: for each `outputFile`, stage from `task:///<file_id>` to
     `pilot:///<task_id>_<file_id>` using `rp.COPY`. For leaf tasks, also stage to
     `client:///data/<file_id>` using `rp.TRANSFER`.

4. **Submit-wait loop**: For each dependency level:
   ```python
   tasks_level_N = tmgr.submit_tasks(level_N_descriptions)
   tmgr.wait_tasks([t.uid for t in tasks_level_N])
   ```

5. **Cleanup**: Always end with `session.close()` in a `try/finally` block.

### TaskDescription Key Attributes
- `uid` (str): Unique task ID
- `executable` (str): Program to run
- `arguments` ([str]): Command-line arguments
- `cores_per_rank` (int): CPU cores per task (default 1)
- `gpus_per_rank` (float): GPUs per task (default 0)
- `mem_per_rank` (int): Memory in MB per task (default 0)
- `input_staging` ([dict]): Input file staging directives
- `output_staging` ([dict]): Output file staging directives
- `pre_exec` ([str]): Shell commands before executable
- `environment` (dict): Environment variables

### Staging Directive Format
```python
{
    'source': '<scheme>:///<path>',
    'target': '<scheme>:///<path>',
    'action': rp.LINK | rp.COPY | rp.TRANSFER,
    'flags' : rp.CREATE_PARENTS
}
```
Schemes: `client://` (launch machine), `pilot://` (pilot sandbox), `task://` (task sandbox)

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
          "inputFiles": [{"id": "sequences.fasta", "sizeInBytes": 50000}],
          "outputFiles": [
            {"id": "chunk_0.fasta", "sizeInBytes": 25000},
            {"id": "chunk_1.fasta", "sizeInBytes": 25000}
          ]
        },
        {
          "name": "blast",
          "id": "blast_00",
          "parents": ["split_00"],
          "children": ["merge_00"],
          "inputFiles": [{"id": "chunk_0.fasta", "sizeInBytes": 25000}],
          "outputFiles": [{"id": "result_0.xml", "sizeInBytes": 10000}]
        },
        {
          "name": "blast",
          "id": "blast_01",
          "parents": ["split_00"],
          "children": ["merge_00"],
          "inputFiles": [{"id": "chunk_1.fasta", "sizeInBytes": 25000}],
          "outputFiles": [{"id": "result_1.xml", "sizeInBytes": 10000}]
        },
        {
          "name": "merge",
          "id": "merge_00",
          "parents": ["blast_00", "blast_01"],
          "children": [],
          "inputFiles": [
            {"id": "result_0.xml", "sizeInBytes": 10000},
            {"id": "result_1.xml", "sizeInBytes": 10000}
          ],
          "outputFiles": [{"id": "final_results.xml", "sizeInBytes": 20000}]
        }
      ],
      "files": []
    },
    "execution": {
      "makespanInSeconds": 120,
      "executedAt": "1970-01-01T00:00:00Z",
      "tasks": [
        {
          "id": "split_00",
          "runtimeInSeconds": 5,
          "command": {"program": "wfbench", "arguments": ["--task-name", "split"]},
          "coreCount": 1
        },
        {
          "id": "blast_00",
          "runtimeInSeconds": 50,
          "command": {"program": "wfbench", "arguments": ["--task-name", "blast"]},
          "coreCount": 1
        },
        {
          "id": "blast_01",
          "runtimeInSeconds": 50,
          "command": {"program": "wfbench", "arguments": ["--task-name", "blast"]},
          "coreCount": 1
        },
        {
          "id": "merge_00",
          "runtimeInSeconds": 10,
          "command": {"program": "wfbench", "arguments": ["--task-name", "merge"]},
          "coreCount": 1
        }
      ],
      "machines": [{"nodeName": "unknown"}]
    }
  }
}
```

### Expected Output
```python
#!/usr/bin/env python
import radical.pilot as rp

session = rp.Session()
try:
    pmgr = rp.PilotManager(session=session)
    tmgr = rp.TaskManager(session=session)

    pilot = pmgr.submit_pilots(rp.PilotDescription({
        'resource'     : 'local.localhost',
        'cores'        : 4,
        'runtime'      : 30,
        'exit_on_error': False
    }))
    tmgr.add_pilots(pilot)
    pilot.wait(rp.PMGR_ACTIVE)

    # ---- Level 0: split_00 (no parents) ----
    td_split_00 = rp.TaskDescription({
        'uid'           : 'split_00',
        'executable'    : 'bin/wfbench',
        'arguments'     : ['--task-name', 'split'],
        'cores_per_rank': 1,
        'input_staging' : [
            {'source': 'client:///data/sequences.fasta',
             'target': 'task:///sequences.fasta',
             'action': rp.TRANSFER, 'flags': rp.CREATE_PARENTS}
        ],
        'output_staging': [
            {'source': 'task:///chunk_0.fasta',
             'target': 'pilot:///split_00_chunk_0.fasta',
             'action': rp.COPY},
            {'source': 'task:///chunk_1.fasta',
             'target': 'pilot:///split_00_chunk_1.fasta',
             'action': rp.COPY}
        ]
    })
    tasks_0 = tmgr.submit_tasks([td_split_00])
    tmgr.wait_tasks([t.uid for t in tasks_0])

    # ---- Level 1: blast_00, blast_01 ----
    td_blast_00 = rp.TaskDescription({
        'uid'           : 'blast_00',
        'executable'    : 'bin/wfbench',
        'arguments'     : ['--task-name', 'blast'],
        'cores_per_rank': 1,
        'input_staging' : [
            {'source': 'pilot:///split_00_chunk_0.fasta',
             'target': 'task:///chunk_0.fasta',
             'action': rp.LINK}
        ],
        'output_staging': [
            {'source': 'task:///result_0.xml',
             'target': 'pilot:///blast_00_result_0.xml',
             'action': rp.COPY}
        ]
    })
    td_blast_01 = rp.TaskDescription({
        'uid'           : 'blast_01',
        'executable'    : 'bin/wfbench',
        'arguments'     : ['--task-name', 'blast'],
        'cores_per_rank': 1,
        'input_staging' : [
            {'source': 'pilot:///split_00_chunk_1.fasta',
             'target': 'task:///chunk_1.fasta',
             'action': rp.LINK}
        ],
        'output_staging': [
            {'source': 'task:///result_1.xml',
             'target': 'pilot:///blast_01_result_1.xml',
             'action': rp.COPY}
        ]
    })
    tasks_1 = tmgr.submit_tasks([td_blast_00, td_blast_01])
    tmgr.wait_tasks([t.uid for t in tasks_1])

    # ---- Level 2: merge_00 ----
    td_merge_00 = rp.TaskDescription({
        'uid'           : 'merge_00',
        'executable'    : 'bin/wfbench',
        'arguments'     : ['--task-name', 'merge'],
        'cores_per_rank': 1,
        'input_staging' : [
            {'source': 'pilot:///blast_00_result_0.xml',
             'target': 'task:///result_0.xml',
             'action': rp.LINK},
            {'source': 'pilot:///blast_01_result_1.xml',
             'target': 'task:///result_1.xml',
             'action': rp.LINK}
        ],
        'output_staging': [
            {'source': 'task:///final_results.xml',
             'target': 'client:///data/final_results.xml',
             'action': rp.TRANSFER}
        ]
    })
    tasks_2 = tmgr.submit_tasks([td_merge_00])
    tmgr.wait_tasks([t.uid for t in tasks_2])

finally:
    session.close()
```
