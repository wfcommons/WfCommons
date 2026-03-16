# Skill: CWL

## Triggers
cwlVersion, class: Workflow, class: CommandLineTool, .cwl, steps:, baseCommand, InlineJavascriptRequirement

## Description
Common Workflow Language (CWL): a YAML/JSON-based standard for describing analysis workflows.
Files use the .cwl extension. Workflows define steps that reference CommandLineTools;
inputs/outputs are typed and linked across steps.

## Domain Knowledge
### CWL Concepts
- **Workflow class**: Top-level document with `class: Workflow`. Contains `inputs`, `outputs`, and `steps`.
- **CommandLineTool class**: Defines a single executable step with `baseCommand`, `inputs`, `outputs`.
- **Steps**: Each entry in the `steps:` section becomes a WfFormat task.
- **Input/Output linking**: Steps reference other steps' outputs via `source:` fields (e.g., `step1/output_name`).
- **Types**: CWL uses typed inputs/outputs: `File`, `File[]`, `string`, `int`, `Directory`, etc.
- **Requirements**: `DockerRequirement`, `InlineJavascriptRequirement`, `ScatterFeatureRequirement`, etc.
- **Scatter**: Runs a step in parallel over an array input (similar to map).

### Mapping to WfFormat
- Each entry in `steps:` -> one task in specification
- Step name/id -> task name and id
- `source:` references pointing to other steps -> parent/children relationships
- Steps with no `source:` referencing other steps -> root tasks (no parents)
- `File` or `File[]` typed inputs -> inputFiles
- `File` or `File[]` typed outputs -> outputFiles
- `baseCommand` -> command.program
- `arguments` -> command.arguments

### Dependency Inference
- If step B has `in: [{id: x, source: stepA/output}]`, then stepA is parent of stepB
- Workflow-level inputs (not from other steps) indicate root task inputs
- `linkMerge: merge_flattened` indicates multiple parents feeding into one step
- Steps with no inbound source references from other steps are root tasks

## Examples

### Input
```yaml
cwlVersion: v1.0
class: Workflow
inputs:
  input_file: File
steps:
  step_align:
    run: align.cwl
    in:
      reads: input_file
    out: [aligned]
  step_sort:
    run: sort.cwl
    in:
      input: step_align/aligned
    out: [sorted]
outputs:
  final:
    type: File
    outputSource: step_sort/sorted
```

### Expected Output
```json
{
  "name": "workflow.cwl",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "step_align",
          "id": "step_align",
          "parents": [],
          "children": ["step_sort"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "step_sort",
          "id": "step_sort",
          "parents": ["step_align"],
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
          "id": "step_align",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "align", "arguments": []},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "step_sort",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "sort", "arguments": []},
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
