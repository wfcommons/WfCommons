#!/bin/bash -euo pipefail
OptiTypePipeline.py -i mapped_1.bam mapped_2.bam -e 1 -b 0.009 \
    -p "example_pe" -c config.ini --dna --outdir example_pe
