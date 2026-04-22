#!/bin/bash

set -e

# Transform into a fasta
i=1
while IFS= read -r line; do
    echo ">seq${i}"
    echo "$line"
    ((i++))
done
