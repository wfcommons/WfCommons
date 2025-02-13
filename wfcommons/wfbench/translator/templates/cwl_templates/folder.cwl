# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

cwlVersion: v1.2
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}

inputs:
  item: 
    type:
      - File
      - Directory
      - type: array
        items:
          - File
          - Directory
  name: string

outputs:
    out: Directory

expression: "${
    if (inputs.item.class == 'Directory'){
        return {
            'out': {
                'class': 'Directory',
                'basename': inputs.name,
                'listing': [inputs.item]
            }
        }
    };
    if (inputs.item.class == 'File'){
        var arr = [inputs.item];
        }
    else {
        var arr = inputs.item;
    }
    return {
        'out': {
            'class': 'Directory',
            'basename': inputs.name,
            'listing': arr
        }
    }
}"
