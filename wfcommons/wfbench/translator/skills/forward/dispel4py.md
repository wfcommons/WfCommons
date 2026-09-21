# Skill: dispel4py

## Triggers
dispel4py, d4py, stream-d4py, StreamingFlow, GenericPE, IterativePE, ProducerPE, ConsumerPE, SimpleFunctionPE, WorkflowGraph, graph.connect, numprocesses, inputconnections, outputconnections, timed_multi, timed_simple, multi_process, simple_process

## Description
Generate dispel4py (StreamingFlow/d4py, PyPI `stream-d4py`) Python workflow scripts from
WfFormat JSON. dispel4py is a *streaming* dataflow system: a workflow is a graph of
Processing Elements (PEs) connected by named in-memory connections, and each PE is
replicated across processes at run time. It is NOT a task-graph scheduler -- there are no
per-task submissions, no files passed between steps, and no dependency list. The graph
describes PE *types* and how their ports are wired; the mapping decides how many processes
each type gets.

## Domain Knowledge

### dispel4py concepts
- **PE (Processing Element)**: a class that processes one data item at a time. Data
  arrives on named input connections and leaves on named output connections. Base
  classes (`dispel4py.base`):
  - `GenericPE`: full control. Declare ports in `__init__` with `self._add_input("name")`
    / `self._add_output("name")`, implement `process(self, inputs)`, and emit with
    `self.write("port", data)`. Use this whenever a PE has more than one input or
    output port.
  - `IterativePE`: exactly one input named `input` and one output named `output`.
    Implement `_process(self, data)` and `return` the result.
  - `ProducerPE`: output only (named `output`), for the PE that starts the stream.
    Implement `_process(self, inputs)` and `self.write("output", item)` per item.
  - `ConsumerPE`: input only (named `input`), for the PE that ends the stream.
    Implement `_process(self, data)` with no return.
- **WorkflowGraph**: `graph = WorkflowGraph()` then
  `graph.connect(src_pe, "src_port", dst_pe, "dst_port")`. Connections are between PE
  *objects*, never between process ranks.
- **numprocesses**: `pe.numprocesses = 4` asks for 4 parallel instances of that PE.
  Instances of the same PE share one stream: each data item goes to exactly ONE of
  them (see Groupings).
- **Module-level `graph`**: the script must expose the built graph as a module-level
  variable. dispel4py imports the module and looks for it (`-a NAME` selects a
  differently-named attribute). Do not run the graph yourself from `__main__` --
  the `dispel4py` command executes it.

### Mapping from WfFormat to dispel4py

**This target inverts the usual rule.** The base instructions say one WfFormat task
becomes one task in the target system. For dispel4py that is wrong and produces
unrunnable code. WfFormat records one task per PE *instance* (per rank), so a workflow
with 98 tasks is typically 7 PE types replicated across 14 ranks each. Emit **one PE
class per task category**, and express the instance count as `numprocesses`. Never emit
one class, one object, or one `graph.connect` call per WfFormat task.

1. **Group tasks into PE types.** Use the task's `category` field when present; it holds
   the dispel4py PE id. Otherwise strip the trailing id from `id`/`name`: both
   `NormalizeDataPE1_ID0000007` and `NormalizeDataPE1_00000094` yield
   `NormalizeDataPE1`. The PE type name becomes the Python class name; keep it verbatim
   so the workflow can be traced back to the source instance.

2. **numprocesses per type** = the number of WfFormat tasks in that group. Assign it
   explicitly (`normalize.numprocesses = 3`) rather than leaving the default.

3. **A category is `pe.name` + the graph index -- split it.** `WorkflowNode.__init__`
   assigns `pe.id = pe.name + str(node_counter)`, where the counter increments as nodes
   join the graph, and the monitoring trace records that id as `pe_id`. So the WfFormat
   category `NormalizeDataPE1` means the PE was named `NormalizeDataPE` and was the
   second node added; `read0` means a PE named `read`, added first.
   - Strip the trailing digits off each category to get the PE name, and read those
     digits as the node's position. Sanity-check it: the stripped numbers must be exactly
     `0..N-1`, each once. If they are not, the trailing digits are part of the real name
     -- use the category verbatim and do not try to reproduce indices.
   - Name the class after the stripped name when it is a valid identifier, and set
     `pe.name` explicitly anyway.
   - Pin the order with `graph.add(pe)` for every PE, in index order, *before* any
     `graph.connect(...)` call. Connecting alone would assign indices in the order PEs
     happen to appear in the connect calls, which silently shifts every id.

   Getting the name wrong does more than cosmetic damage: `-d` root-input keys are
   matched against `pe.name` and `pe.id` (`get_inputs` tries `inputs[pe]`, then
   `inputs[pe.name]`, then `inputs[pe.id]`). If neither matches, the source PE is invoked
   without its root input and dies with `KeyError: 'input'` inside a worker process --
   under the `multi` mappings the other ranks then block forever instead of reporting the
   error.

4. **Connections come from the file ids, not just from parents/children.** The converter
   that produces these instances names every in-memory stream
   `<pe_id>:<rank>.<port>` -- for example `DeterministicPrecheckPE2:4.resolved` is the
   stream written by rank 4 of `DeterministicPrecheckPE2` on its output port `resolved`.
   So:
   - A file id matching `<pe>:<rank>.<port>` is a **stream**, not a file. The text after
     the last `.` is the producing PE's **output port name** -- declare that port with
     `_add_output("resolved")` and write to it with `self.write("resolved", data)`.
   - Collapse to type level: every stream produced by any rank of PE `A` on port `p`,
     consumed by any rank of PE `B`, becomes exactly one
     `graph.connect(a, "p", b, "<in_port>")`. Several ranks producing the same port is
     one connection, not one per rank.
   - If the file ids carry no `<pe>:<rank>.<port>` shape (an older instance), fall back
     to the `parents`/`children` DAG and use the default port names `output` / `input`.

5. **Consumer port names are not recoverable -- name them deterministically.** WfFormat
   stores only the producing side of each stream. Apply this rule:
   - A PE type consuming exactly one upstream port: name its input `input`.
   - A PE type consuming several distinct upstream ports: declare one input per upstream
     port, named after that upstream port (a PE fed by `...:N.resolved` and
     `...:N.output` from the agent declares inputs `resolved` and `output`). Add a
     comment noting the names were reconstructed, because the original script may have
     used different ones.

6. **Real files stay files.** A file id that is *not* in `<pe>:<rank>.<port>` form (e.g.
   `sensor_data_parallel_100.json`, `agentic_parallel_results.jsonl`) is a genuine
   on-disk file:
   - On a **root** task, the source PE reads it. Make that PE a `ProducerPE` that takes
     the path from its root input (`file_path = inputs["input"]`), reads it, and writes
     one stream item per record. Pass the path at run time with
     `-d '{"<pe_name>": [{"input": "sensor_data_parallel_100.json"}]}'`.
   - On a **leaf** task, the sink PE appends to it -- a `ConsumerPE` opening the file in
     append mode per item.

7. **Runtime.** WfFormat's `runtimeInSeconds` is per instance, per whole run -- not per
   item. If `execution.tasks[].command.program` is `wfbench`, shell out to it with the
   recorded arguments. Otherwise emit a placeholder `_process` that does the PE's real
   work if it is evident from the category name, and busy-waits for the recorded
   per-item share otherwise (`total runtime / count of items`), with a comment saying
   the body is a stand-in. Never silently emit an empty `_process`.

8. **Ordering.** Emit PE classes in topological order (source first), then instantiate,
   then set `numprocesses`, then `graph.connect(...)` calls in the same order. Identical
   input must produce identical output.

### Groupings and data routing
The default routing between two connected PE types is `ShuffleCommunication`:
round-robin, so **each data item reaches exactly one instance** of the destination PE.
That is what you want for data-parallel fan-out; do not emit any broadcast logic for it.
Only set a grouping when the source workflow needs different semantics:
- `pe.inputconnections["input"]["grouping"] = "global"` -- all items to one instance
  (`AllToOneCommunication`); use for an aggregator that must see everything.
- `pe.inputconnections["input"]["grouping"] = "all"` -- every item copied to every
  instance (`OneToAllCommunication`); this is the only true broadcast.
- `pe.inputconnections["input"]["grouping"] = ["key1", "key2"]` -- hash-partition by
  those keys (`GroupByCommunication`); use when instances must own a key space.
Ports can also declare it at construction: `self._add_input("input", grouping="global")`.
WfFormat does not record groupings, so default to shuffle (emit nothing) unless the
instance shows an obvious all-to-one structure (a PE type with `numprocesses = 1`
receiving from a wider upstream type still works fine under shuffle -- do not add
`global` for that).

### Running the generated workflow
Put the exact command in the module docstring. The runner is:

```
dispel4py <mapping> <module> [-i ITERATIONS] [-n SIZE] [-d JSON | -f FILE] [-a ATTR]
```

- Mappings (from `dispel4py/new/mappings.py`): `simple`, `multi`, `mpi`, `redis`,
  `spark`, `dyn_multi`, `dyn_auto_multi`, `dyn_redis`, `dyn_auto_redis`,
  `hybrid_redis`, and the monitoring variants `timed_simple`, `timed_multi`,
  `timed_mpi`. The timed variants additionally write `monitor_*` trace files and accept
  `--timing-dir`, `--timing-prefix` and `--run-id`.
- `-i N` is the number of iterations (items the source produces), not the process count.
- `-n SIZE` is the total process budget and **must be at least**
  `1 per source PE + sum(numprocesses of every non-source PE)`; dispel4py refuses with
  "Graph is larger than job size" otherwise. State the computed minimum in the docstring.
- A PE with `numprocesses > 1` gets exactly that many instances. A PE left at
  `numprocesses = 1` is instead *scaled* by the process budget, and a source PE always
  gets exactly one process regardless of what it asks for.

### What NOT to emit
- One PE class, PE object, or `graph.connect` call per WfFormat task. Group by category.
- `wfbench`-style file staging between PEs, `inputFiles`/`outputFiles` handling that
  writes intermediate files, or any `open()` of a stream id. Streams are in memory.
- Rank or instance ids anywhere in the graph construction (`NormalizeDataPE1:2` is a
  trace artifact, never a Python identifier).
- A `if __name__ == "__main__":` block that builds or runs the graph itself, `mpi4py`
  imports, or a manual process pool. The mapping does all of that.
- `graph.connect()` with a port that was never declared with `_add_input` /
  `_add_output` on the PE -- dispel4py fails at run time, not at import.

## Examples

### Input (WfFormat JSON, abbreviated)
```json
{
  "schemaVersion": "1.6",
  "workflow": {
    "specification": {
      "tasks": [
        {"id": "read0_ID0000015", "name": "read0_ID0000015", "category": "read0",
         "parents": [], "children": ["NormalizeDataPE1_ID0000007"],
         "inputFiles": ["sensor_data_parallel_100.json"],
         "outputFiles": ["read0:0.output"]},
        {"id": "NormalizeDataPE1_ID0000007", "category": "NormalizeDataPE1",
         "parents": ["read0_ID0000015"], "children": ["DeterministicPrecheckPE2_ID0000004"],
         "inputFiles": ["read0:0.output"],
         "outputFiles": ["NormalizeDataPE1:1.output"]},
        {"id": "DeterministicPrecheckPE2_ID0000004", "category": "DeterministicPrecheckPE2",
         "inputFiles": ["NormalizeDataPE1:1.output"],
         "outputFiles": ["DeterministicPrecheckPE2:4.resolved",
                         "DeterministicPrecheckPE2:4.needs_agent"]},
        {"id": "DecisionMergePE3_ID0000002", "category": "DecisionMergePE3",
         "inputFiles": ["DeterministicPrecheckPE2:4.resolved",
                        "ParallelLLMSensorAgentPE4:9.output"],
         "outputFiles": ["DecisionMergePE3:7.output"]}
      ],
      "files": [{"id": "sensor_data_parallel_100.json", "sizeInBytes": 114841},
                {"id": "read0:0.output", "sizeInBytes": 0}]
    },
    "execution": {"tasks": [{"id": "read0_ID0000015", "runtimeInSeconds": 0.0004}]}
  }
}
```
Three `NormalizeDataPE1_*` tasks, three `DeterministicPrecheckPE2_*`, four
`ParallelLLMSensorAgentPE4_*` and two `DecisionMergePE3_*` tasks in the full instance
become `numprocesses` 3, 3, 4 and 2 -- not 12 classes.

### Expected Output (`run_workflow.py`)
```python
#!/usr/bin/env python
"""
dispel4py workflow generated from WfFormat.

Run (16 processes = 1 source + 3 + 3 + 4 + 2 + 2 + 1):

    dispel4py timed_multi run_workflow.py -i 100 -n 16 \
        -d '{"read0": [{"input": "sensor_data_parallel_100.json"}]}'
"""

import json
import time

from dispel4py.base import ConsumerPE, GenericPE, IterativePE, ProducerPE
from dispel4py.workflow_graph import WorkflowGraph

OUTPUT_FILE = "agentic_parallel_results.jsonl"


class read(ProducerPE):
    """Root PE: reads the real input file and emits one item per record."""

    def _process(self, inputs):
        file_path = inputs["input"]
        with open(file_path, "r", encoding="utf-8") as handle:
            for event in json.load(handle):
                self.write("output", event)


class NormalizeDataPE(IterativePE):
    """Single input/output: IterativePE gives ports "input" and "output"."""

    def _process(self, data):
        time.sleep(0.000006)  # stand-in for the recorded per-item runtime
        return data


class DeterministicPrecheckPE(GenericPE):
    """Two output ports, taken from the stream ids .resolved and .needs_agent."""

    def __init__(self):
        super().__init__()
        self._add_input("input")
        self._add_output("resolved")
        self._add_output("needs_agent")

    def process(self, inputs):
        data = inputs["input"]
        if data.get("confidence", 0.0) >= 0.9:
            self.write("resolved", data)
        else:
            self.write("needs_agent", data)


class ParallelLLMSensorAgentPE(IterativePE):

    def _process(self, data):
        time.sleep(0.02)  # stand-in for the recorded per-item runtime
        return data


class DecisionMergePE(GenericPE):
    """Input port names reconstructed from the upstream port names."""

    def __init__(self):
        super().__init__()
        self._add_input("resolved")
        self._add_input("output")
        self._add_output("output")

    def process(self, inputs):
        for port in ("resolved", "output"):
            if port in inputs:
                self.write("output", inputs[port])


class ActionExecutorPE(IterativePE):

    def _process(self, data):
        time.sleep(0.00006)  # stand-in for the recorded per-item runtime
        return data


class ResultWriterPE(ConsumerPE):
    """Leaf PE: appends to the real output file."""

    def _process(self, data):
        with open(OUTPUT_FILE, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(data) + "\n")


read0 = read()
normalize = NormalizeDataPE()
precheck = DeterministicPrecheckPE()
agent = ParallelLLMSensorAgentPE()
merge = DecisionMergePE()
execute = ActionExecutorPE()
write_results = ResultWriterPE()

# Name = category minus its trailing graph index, so that dispel4py rebuilds the
# original ids: "read" + 0 -> read0, "NormalizeDataPE" + 1 -> NormalizeDataPE1, ...
read0.name = "read"
normalize.name = "NormalizeDataPE"
precheck.name = "DeterministicPrecheckPE"
agent.name = "ParallelLLMSensorAgentPE"
merge.name = "DecisionMergePE"
execute.name = "ActionExecutorPE"
write_results.name = "ResultWriterPE"

read0.numprocesses = 1
normalize.numprocesses = 3
precheck.numprocesses = 3
agent.numprocesses = 4
merge.numprocesses = 2
execute.numprocesses = 2
write_results.numprocesses = 1

graph = WorkflowGraph()

# Added in category-index order, before connecting: this is what fixes each
# pe.id (name + index) to the category recorded in the source instance.
for pe in (read0, normalize, precheck, merge, agent, execute, write_results):
    graph.add(pe)

graph.connect(read0, "output", normalize, "input")
graph.connect(normalize, "output", precheck, "input")
graph.connect(precheck, "resolved", merge, "resolved")
graph.connect(precheck, "needs_agent", agent, "input")
graph.connect(agent, "output", merge, "output")
graph.connect(merge, "output", execute, "input")
graph.connect(execute, "output", write_results, "input")
```

### Re-tracing the generated workflow
When the goal is a round trip (WfFormat -> dispel4py -> WfFormat), run it under a timed
mapping so the monitor artifacts are written next to the script:

```
dispel4py timed_multi run_workflow.py -i 100 -n 16 \
    -d '{"read0": [{"input": "sensor_data_parallel_100.json"}]}' \
    --timing-dir monitoring_multi --timing-prefix monitor
```

## Validation

### Required elements
- from dispel4py.workflow_graph import WorkflowGraph
- WorkflowGraph()
- graph.connect
- numprocesses

### Syntax check
command: python -m py_compile {file}
