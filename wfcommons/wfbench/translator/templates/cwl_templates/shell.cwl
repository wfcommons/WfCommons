cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
stdout: $("logs/" + inputs.step_name + ".out")
stderr: $("logs/" + inputs.step_name + ".err")

arguments:
  - position: 1
    valueFrom: >
      ${
        var cmd = inputs.command;
        if (inputs.input_files) {
          for (var i = 0; i < inputs.input_files.length; i++) {
            cmd = cmd.replace(new RegExp(inputs.input_files[i].basename, 'g'), inputs.input_files[i].path);
          }
        }
        return cmd;
      }
    shellQuote: false

inputs:
  command:
    type: string
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
