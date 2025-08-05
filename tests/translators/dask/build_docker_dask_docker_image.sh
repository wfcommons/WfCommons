#!/bin/bash

set -e

docker build --platform linux/amd64 -t wfcommons/wfcommons-testing-dask -f Dockerfile_Dask .
