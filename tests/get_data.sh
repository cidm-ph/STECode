#!/bin/bash

set -e

DATA_FILES=(
  "tests/data/test1_trimmed.paired_R1.fq.gz"
  "tests/data/test1_trimmed.paired_R2.fq.gz"
  "tests/data/test2_trimmed.paired_R1.fq.gz"
  "tests/data/test2_trimmed.paired_R2.fq.gz"
)

for file in "${DATA_FILES[@]}"; do
  if ! test -e "$file"; then
    echo "MISSING: $file"
  fi
done
