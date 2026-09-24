# wfstream

**Predict what a dispel4py workflow will cost before you run it at scale.**

Someone asks: *how long would this take with 500 processes and 10,000 sensor
readings, and how much CPU and memory would it need?*

Running it to find out is expensive — and if the workflow calls an LLM, slow and
billable too. So instead:

> Run it **small** a few times → learn what each part costs → build a synthetic
> workflow of the size you care about → predict from that.

---

## The idea in one picture

```
   a few small real runs              a workflow you never ran
  ┌──────────────────────┐          ┌──────────────────────────┐
  │ dispel4py traces     │          │ synthetic 500-task        │
  │ (monitor_*.csv/json) │          │ WfFormat instance         │
  └──────────┬───────────┘          └────────────┬─────────────┘
             │                                   │
      what shape is it?  ──── recipe ────────────┘
      what does it cost? ──── stats ──────────────┐
                                                  ▼
                                    runtime · CPU · memory
                                      per PE and overall
```

Two separate things are learned, because they answer different questions:

| learned | by | answers |
|---|---|---|
| **shape** | WfChef recipe | what does this workflow look like at 500 tasks? |
| **cost** | `resource_stats` | what does each PE cost per item? |

---

## Why not just use WfCommons as-is?

WfCommons was built for **task workflows**: a DAG of jobs where each runs once,
and the workflow finishes when the critical path finishes.

dispel4py is a **streaming** system, and three things differ:

**Scaling means more processes, not more work.** WfChef grows a workflow by
copying subgraphs. Copy a dispel4py pipeline and you get two pipelines with two
readers — nonsense. Scaling dispel4py means giving a PE `numprocesses > 1`. So
`streaming_recipe` replicates PE instances instead, never the source.

**Everything runs at once.** PEs are concurrent processes exchanging items in
memory. The workflow is not done when a path completes — it is done when the
**slowest stage** has drained the stream. Adding up stage times would be wrong.

**There are no files between stages.** PEs stream over named connections, so
WfFormat's file dependencies are recorded as zero-byte streams. Only the first
and last PE touch real files.

---

## What the prediction is

An analytical model, not a scheduler:

```
items reaching a PE  = items × selectivity          not every stage sees every item
how long it is busy  = ceil(items / instances)      its busiest instance's share
                       × seconds per item
total runtime        = slowest stage + fill/drain   everything else overlaps
```

Three details make it accurate, each forced by real measurements:

**Selectivity.** A precheck filters, so the expensive LLM stage sees 5% of the
stream. Assume otherwise and you are 20× wrong on the only stage that matters.

**The ceiling.** Items are whole things. Five items across four instances go
2/1/1/1, so the stage lasts as long as the instance holding two — and one
instance sits idle. That really happens in the traces.

**Fill follows the graph.** One item must reach the bottleneck before it can
start, and the last must leave it afterwards. That is the *slowest path* in and
out — branches beside the bottleneck overlap it and are not charged.

CPU is learned separately, never derived from runtime: a PE waiting on an LLM
call burns ~10% of a core for ten seconds, while a trivial PE burns ~100% for
microseconds.

Memory is a per-process interpreter baseline (~71 MB, mostly imports) plus what
each PE adds on top. Process RSS is never summed — under `timed_simple` every PE
shares one process, so summing would count the same memory seven times.

---

## Using it

Two entry points, for the two situations a registry runs into:

```python
from wfcommons.wfstream import on_new_workflow, on_new_size_run

# A workflow nobody has seen: convert, cook, learn, generate, predict
result = on_new_workflow(trace_dirs, num_tasks=500, items=10_000)
result["synthetic"]    # the generated WfFormat instance
result["simulation"]   # runtime, CPU, memory, per PE and overall

# A known workflow was run for real at a new size: fold it in, stop
on_new_size_run([new_trace_dir])
```

Case A writes three things beside the synthetic instance: the instance itself,
`*.prediction.json` (the full result, and what it was predicted from), and
`*.summary.txt` — the plain-language version:

```
climate: 500 processes, 10,000 items

  Runtime    70.3 seconds   (between 69.8 seconds and 70.9 seconds)
  CPU        5.5 cores on average, 388 core-seconds in total
  Memory     38.0 GB across 500 processes

  The time goes almost entirely to LLMSensorAgentPE4 (100% of it), which
  handles 500 of the 10,000 items across 83 processes.

  Worth knowing:
    - Processes are mostly idle -- 5.5 cores busy out of 500. Adding processes
      will not help unless LLMSensorAgentPE4 gets more of them.
```

### Making a cooked recipe discoverable

`create_recipe` writes a `pyproject.toml` at the build directory's root that
declares the recipe as a `workflow_recipes` entry point — but nothing installs
it, so the recipe stays invisible to `wfchef ls` and to `get_recipe`, which is
how the rest of WfCommons finds a recipe.

`register()` closes that gap:

```python
from wfcommons.wfstream import build_recipe

cooked = build_recipe.cook(corpus, build_dir, "climate")
build_recipe.register(build_dir)        # -> pip install; 'climate_recipe'
```

or `python -m wfcommons.wfstream.build_recipe --register`, or
`on_new_workflow(..., register=True)`.

Once registered, the recipe behaves like a built-in: `wfchef ls` lists it, and
`load_recipe` resolves it through its entry point rather than a hardcoded path —
so it works in a fresh environment that only pip-installed the package.

`install()` is the faster alternative that skips pip: it copies the data over
whichever copy of the recipe actually resolves. Use it to refresh a recipe that
is already in place, and `register()` for one that is not.

### Scoring

Case B scores itself. Before a new run is folded into the statistics, the model
has never seen it — so predicting it then is a genuine held-out test. Each result
is appended to `accuracy.jsonl`, and the accuracy record builds up on its own as
runs arrive.

Every step is also its own module and CLI:

```bash
python -m wfcommons.wfstream.convert_traces monitoring_*      # traces → WfFormat
python -m wfcommons.wfstream.build_recipe --install           # learn the shape
python -m wfcommons.wfstream.resource_stats monitoring_* -o stats.json
python -m wfcommons.wfstream.generate_workflows 500           # build one
python -m wfcommons.wfstream.simulate climate-500.json -s stats.json -i 10000
```

| module | role |
|---|---|
| `dispel_fwd_converter` | traces → WfFormat (the parser) |
| `convert_traces` | build the instance corpus |
| `update_traces` | add new runs to an existing corpus |
| `build_recipe` | cook and install a WfChef recipe |
| `streaming_recipe` | the scaling rule: replicate PEs, not pipelines |
| `generate_workflows` | write synthetic instances, validate, attach costs |
| `resource_stats` | per-PE time / CPU / memory from the CSVs |
| `simulate` | the prediction |
| `pipeline` | the two entry points above |
| `config` | defaults for one workflow |

---

## What you need

Traces from **at least two runs of different sizes**, collected with resource
sampling on. One run is never enough: WfChef finds repeated structure by
comparing instances of different sizes, and a single pipeline with one instance
per PE has nothing that repeats.

Ideally a `timed_simple` run (one process per PE — the shape generation grows
from) and `timed_multi` runs at a couple of widths (where memory can actually be
attributed, since each rank gets its own process).

dispel4py itself is **not** required. wfstream reads traces; it does not run
workflows.

---

## What it does not do

- **Selectivity is extrapolated as a ratio.** If a filter passes a roughly
  constant *number* of items rather than a constant fraction, large predictions
  will be too high. This is the least certain part of the model.
- **Memory does not grow with the stream.** It is modelled from shape alone,
  which matches traces where RSS plateaus. A PE that accumulates per item breaks
  that.
- **No end-to-end wall clock is recorded** by dispel4py's monitoring, only
  per-PE service time. The runtime model is validated stage by stage.
- **Variance bounds precision.** Two identical runs of the sample workflow
  differed by 1.7× on the dominant stage, because the LLM API does. Predictions
  carry a range for that reason.
- **One mapping for the whole workflow** — every PE is assumed to run the same
  way (all separate processes, or all shared).

See [`examples/wfstream_dispel4py.ipynb`](../../examples/wfstream_dispel4py.ipynb)
for the whole thing end to end, with real numbers and the model checked against a
discrete-event simulation.
