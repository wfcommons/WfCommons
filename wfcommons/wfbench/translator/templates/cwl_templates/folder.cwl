# A collection of RNA-Seq data analysis tools wrapped in CWL scripts
# Copyright (C) 2019 Alessandro Pio Greco, Patrick Hedley-Miller, Filipe Jesus, Zeyu Yang
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env cwl-runner

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
