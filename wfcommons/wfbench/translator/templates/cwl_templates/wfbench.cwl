# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
baseCommand: wfbench
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
      itemSeparator: " "
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
