#!/bin/bash -euo pipefail
awk -vFS="[:-]" '{
  name = sprintf("%s_%d-%d", $1, $2, $3);
  printf("%s\t%d\t%d\n", $1, $2-1, $3) > name ".bed"
}' small.intervals
