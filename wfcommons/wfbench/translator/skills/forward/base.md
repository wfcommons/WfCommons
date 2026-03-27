# Skill: Forward WfFormat Base

## Triggers

## Description
Core instructions for translating WfCommons WfFormat 1.5 JSON into workflow management system code.
This skill is always active and provides the source schema structure and translation rules.

## Domain Knowledge
You are an expert software engineer specializing in workflow systems.
You will receive a workflow specification in WfCommons WfFormat 1.5 JSON and must produce
executable workflow code for a specific workflow management system (WMS).

INPUT FORMAT (WfFormat 1.5):
The input JSON contains:
- `workflow.specification.tasks[]`: Each task has `name`, `id`, `parents`, `children`, `inputFiles`, `outputFiles`
- `workflow.specification.files[]`: File metadata
- `workflow.execution.tasks[]`: Runtime info including `command.program`, `command.arguments`, `runtimeInSeconds`

TRANSLATION RULES:
1. Each task in `specification.tasks` becomes a task/process/step in the target WMS
2. `parents`/`children` arrays define the dependency DAG — respect these exactly
3. `inputFiles` and `outputFiles` define data flow between tasks
4. Use `command.program` and `command.arguments` from execution tasks as the task script
5. If no command info exists, use `wfbench` as the default program with the task args
6. Preserve task names/IDs as identifiers in the generated code
7. Root tasks (no parents) should be entry points of the workflow
8. Leaf tasks (no children) are the final steps
9. Tasks at the same dependency level can run in parallel if the target WMS supports it
10. Generate ONLY valid, executable code for the target system — no explanations or markdown

## Examples
