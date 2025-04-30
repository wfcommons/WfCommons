cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
stdout: $(inputs.step_name + ".out")
stderr: $(inputs.step_name + ".err")

inputs:
  command:
    type: string
    inputBinding:
      position: 1
      shellQuote: false
  input_files:
    type: File[]?
  output_filenames:
    type: string[]?
  step_name:
    type: string

outputs:
  out:
    type: stdout
  err:
    type: stderr
  output_files:
    type: File[]
    outputBinding:
      glob: $(inputs.output_filenames)
