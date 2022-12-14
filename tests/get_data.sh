#!/bin/bash

set -e

DATA_FILES=(
  "tests/data/test1_trimmed.paired_R1.fq.gz"
  "tests/data/test1_trimmed.paired_R2.fq.gz"
  "tests/data/test2_trimmed.paired_R1.fq.gz"
  "tests/data/test2_trimmed.paired_R2.fq.gz"
)

get_file() {
  rsync "hpc.sydney.edu.au:/project/pgraw/stec/fastq/$1" "tests/data/$1"
}

for file in "${DATA_FILES[@]}"; do
  if ! test -e "$file"; then
    name="$(basename "$file")"
    echo "FETCHING: $name"
    get_file "$name"
  fi
done
