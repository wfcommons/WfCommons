# Skill: Nextflow

## Triggers
process, channel, nextflow, .nf, workflow {, emit:, input:, output:, script:, params.

## Description
Nextflow workflow system: process-based, channel-driven dataflow programming.
Workflow files use the .nf extension. Processes define tasks; channels connect them;
the workflow block defines execution order.

## Domain Knowledge
### Nextflow Concepts
- **Processes**: Independent tasks. Each `process NAME { ... }` block becomes a WfFormat task.
- **Channels**: Data flow pipes connecting processes. Created via `channel.of(...)`, `channel.fromPath(...)`, `channel.value(...)`.
- **Input/Output declarations**: Inside processes, `input:` and `output:` blocks define data ports.
  - `val` inputs/outputs: simple values
  - `path` or `file` inputs/outputs: file-based data (map to inputFiles/outputFiles)
  - `tuple` inputs/outputs: grouped data
- **Piping syntax**: `channel | PROCESS` or `PROCESS1 | PROCESS2` defines execution flow.
- **Emit**: `PROCESS.out` or named outputs like `PROCESS.out.results` reference process outputs.

### Mapping to WfFormat
- Each `process` block -> one task in specification (name = process name, id = process name)
- Channel connections and piping -> parent/children relationships
- `input:` declarations with `path`/`file` type -> inputFiles
- `output:` declarations with `path`/`file` type -> outputFiles
- `val` type inputs/outputs -> no files, just dependency edges
- `script:` block content -> command.program (typically the first command) and command.arguments
- If no explicit script, use "echo" as the default program

### Dependency Inference
- `A | B` means A is parent of B
- `B(A.out)` means A is parent of B
- `channel.of(...) | A` means A is a root task (no parents)
- If process B takes input from process A's output channel, A is parent of B
- Linear chains: `A | B | C` means A -> B -> C
- Fork: one output consumed by multiple processes = one parent, multiple children
- Join: multiple outputs consumed by one process = multiple parents, one child

## Examples

### Input (hello.nf)
```nextflow
process SAYHELLO {
    input:
    val name

    output:
    val "Hello, $name!" into GREETING

    // script: echo "Hello, $name!"
}

process GREETING {
    input:
    val greeting

    script:
    """
    echo printing the greeting: $greeting from SAYHELLO process
    """
}

workflow {
    channel.of('Alice', 'Bob', 'Charlie') | SAYHELLO
}
```

### Expected Output
```json
{
  "name": "hello.nf",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "SAYHELLO",
          "id": "SAYHELLO",
          "parents": [],
          "children": ["GREETING"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "GREETING",
          "id": "GREETING",
          "parents": ["SAYHELLO"],
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
          "id": "SAYHELLO",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "echo", "arguments": []},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "GREETING",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "echo", "arguments": []},
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
