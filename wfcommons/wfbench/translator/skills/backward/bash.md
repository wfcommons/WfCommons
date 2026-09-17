# Skill: Bash

## Triggers
#!/bin/bash, #!/bin/sh, #!/usr/bin/env bash

## Description
Bash/shell script workflows: sequential or piped commands in a shell script.
Each significant command becomes a task; execution order and pipes define dependencies.

## Domain Knowledge
### Bash Workflow Concepts
- **Commands**: Each standalone command or pipeline stage is a potential task.
- **Sequential execution**: Commands on separate lines run in order (implicit dependency).
- **Pipes**: `cmd1 | cmd2` means cmd1's output feeds into cmd2 (cmd1 is parent).
- **Variables**: `VAR=$(cmd)` captures output; later use of `$VAR` implies dependency.
- **Conditionals/loops**: `if`, `for`, `while` blocks represent control flow, not separate tasks unless they contain distinct computation steps.
- **Functions**: Named shell functions can be treated as task definitions.

### Mapping to WfFormat
- Each significant command -> one task in specification
- Command name -> task name and id (use a descriptive name if possible)
- Sequential order -> parent/children (previous command is parent of next)
- Pipe chains: `a | b | c` -> a is parent of b, b is parent of c
- The first command in the script -> root task (no parents)
- Program name (first word of command) -> command.program
- Remaining arguments -> command.arguments
- File arguments -> inputFiles/outputFiles if identifiable

### What Counts as a Task
- Major computation commands (not `cd`, `echo` for logging, `export`, `set`)
- Pipeline stages that process data
- Function calls that perform work
- Script invocations (`python script.py`, `./run.sh`)

## Examples

### Input
```bash
#!/bin/bash
samtools sort input.bam -o sorted.bam
samtools index sorted.bam
bcftools call -m sorted.bam -o variants.vcf
```

### Expected Output
```json
{
  "name": "pipeline.sh",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "samtools_sort",
          "id": "samtools_sort",
          "parents": [],
          "children": ["samtools_index"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "samtools_index",
          "id": "samtools_index",
          "parents": ["samtools_sort"],
          "children": ["bcftools_call"],
          "inputFiles": [],
          "outputFiles": []
        },
        {
          "name": "bcftools_call",
          "id": "bcftools_call",
          "parents": ["samtools_index"],
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
          "id": "samtools_sort",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "samtools", "arguments": ["sort", "input.bam", "-o", "sorted.bam"]},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "samtools_index",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "samtools", "arguments": ["index", "sorted.bam"]},
          "coreCount": 1,
          "avgCPU": 0,
          "readBytes": 0,
          "writtenBytes": 0,
          "memoryInBytes": 0,
          "machines": ["unknown"]
        },
        {
          "id": "bcftools_call",
          "runtimeInSeconds": 0,
          "executedAt": "1970-01-01T00:00:00Z",
          "command": {"program": "bcftools", "arguments": ["call", "-m", "sorted.bam", "-o", "variants.vcf"]},
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
