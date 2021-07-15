#!/bin/bash -euo pipefail
export HOME="${PWD}/HOME"

qiime demux summarize 		--i-data demux.qza 		--o-visualization demux.qzv

qiime tools export --input-path demux.qzv --output-path demux
