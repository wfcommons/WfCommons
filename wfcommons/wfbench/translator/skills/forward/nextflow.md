# Skill: Forward Nextflow

## Triggers
nextflow, .nf

## Description
Generate Nextflow workflow (.nf) files from WfFormat JSON.

## Domain Knowledge
### Nextflow Output Structure
Generate a complete Nextflow `.nf` file with these components:

1. **Parameter declarations**: `params.pwd`, `params.simulate`, `params.help`
2. **Processes**: One `process NAME { ... }` block per WfFormat task
3. **Workflow block**: Wires processes together respecting the dependency DAG

### Mapping from WfFormat
- Each task in `specification.tasks` -> one Nextflow `process` block
  - Process name = task ID (with `.` replaced by `_`)
  - `parents` -> input channel dependencies
  - `children` -> output channels
  - `inputFiles` -> `input:` declarations with `val` type (collected from parent outputs)
  - `outputFiles` -> `output:` declarations as `val` tuples
- Task script: `bash ${pwd}/bin/script_<task_id>.sh`
  - Include simulate mode: `${params.simulate ? 'sleep 1' : "bash ${pwd}/bin/script_<task_id>.sh"}`
- Root tasks (no parents): called directly in workflow block
- Tasks with parents: receive collected input channels via `.mix()` and `.collect()`

### Code Pattern
For each task, generate:
1. A `process` block with input/output declarations and script
2. A helper `function_<task_id>` that handles channel wiring

The workflow block calls `bootstrap()` then chains all task functions:
```
workflow {
    results = bootstrap()
    results = function_task1(results)
    results = function_task2(results)
    ...
}
```

### File Dependencies
- Each task should reference `bin/script_<task_id>.sh` for its execution
- Input/output files live under `data/` directory
- Binary files live under `bin/` directory

## Examples

### Input (WfFormat JSON)
```json
{
  "name": "example-workflow",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "preprocess",
          "id": "preprocess_00",
          "parents": [],
          "children": ["analyze_00"],
          "inputFiles": [{"id": "raw_data.csv", "sizeInBytes": 1000}],
          "outputFiles": [{"id": "clean_data.csv", "sizeInBytes": 500}]
        },
        {
          "name": "analyze",
          "id": "analyze_00",
          "parents": ["preprocess_00"],
          "children": [],
          "inputFiles": [{"id": "clean_data.csv", "sizeInBytes": 500}],
          "outputFiles": [{"id": "results.json", "sizeInBytes": 200}]
        }
      ],
      "files": []
    }
  }
}
```

### Expected Output
```nextflow
params.simulate = false
params.pwd = null
params.help = null
pwd = null

// ... parameter validation ...

process preprocess_00() {
    input:
    output:
        val([clean_data.csv])
    script:
        """
        ${params.simulate ? 'sleep 1' : "bash ${pwd}/bin/script_preprocess_00.sh"}
        """
}

process analyze_00() {
    input:
        val analyze_00_ifs
    output:
    script:
        def clean_data.csv = analyze_00_ifs[0]
        """
        ${params.simulate ? 'sleep 1' : "bash ${pwd}/bin/script_analyze_00.sh"}
        """
}

def bootstrap() {
    def outputs = [:]
    return outputs
}

// ... task functions and workflow block ...
```
