#!/bin/bash

set -e

for backend in "dask" "parsl" "nextflow" "airflow" "bash" "taskvine" "cwl"; do
  echo "Building $backend Docker image..."
  docker build --platform linux/amd64 -t wfcommons/wfcommons-testing-$backend -f Dockerfile.$backend .
done
