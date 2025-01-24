cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
baseCommand: /Users/wongy/Documents/GitHub/WfFormat-To-StreamFlow-Translator/script/output/bin/wfbench.py
arguments:
  - valueFrom: $(inputs.input_params)
stdout: $(inputs.step_name + ".out")
stderr: $(inputs.step_name + ".err")
inputs:
  step_name:
    type: string
  input_params:
    type: string[]?
  input_files:
    type: File[]?
    inputBinding:
      position: 0
  output_filenames:
    type: string[]?
outputs:
  out:
    type: stdout
  err:
    type: stderr
  output_files:
    type:
      type: array
      items:
        - File
    outputBinding:
      glob: $(inputs.output_filenames)
