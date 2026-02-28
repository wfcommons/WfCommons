cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

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
        cmd = cmd + " > " + runtime.outdir + "/" + inputs.step_name + ".out 2> " + runtime.outdir + "/" + inputs.step_name + ".err";
        cmd = cmd + " ; echo '-- end of stdout for " + inputs.step_name + " --' >> " + runtime.outdir + "/" + inputs.step_name + ".out";
        cmd = cmd + " ; echo '-- end of stderr for " + inputs.step_name + " --' >> " + runtime.outdir + "/" + inputs.step_name + ".err";
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
    type: File
    outputBinding:
      glob: $(inputs.step_name + ".out")
  err:
    type: File
    outputBinding:
      glob: $(inputs.step_name + ".err")
  output_files:
    type: File[]
    outputBinding:
      glob: $(inputs.output_filenames)
