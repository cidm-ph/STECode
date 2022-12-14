#!/bin/bash

# Downloads data needed for tests
#
# NB: currently requires access to a specific server

set -e

DATA_FILES=(
  "test1_trimmed.paired_R1.fq.gz"
  "test1_trimmed.paired_R2.fq.gz"
  "test2_trimmed.paired_R1.fq.gz"
  "test2_trimmed.paired_R2.fq.gz"
)

get_file() {
  rsync "hpc.sydney.edu.au:/project/pgraw/stec/fastq/$1" "tests/data/$1"
}

if ! test -d "tests/data"; then
  "ERROR: run from main STECode directory" >&2
  exit 1
fi

for file in "${DATA_FILES[@]}"; do
  if test -e "tests/data/$file"; then
    echo "OK:  $file"
  else
    echo "GET: $file"
    get_file "$file"
  fi
done
