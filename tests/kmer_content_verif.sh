#!/bin/bash

set -e

kmc=$(which kmc)       || { echo "kmc not found"; exit 1; }
kmc_tools=$(which kmc_tools) || { echo "kmc_tools not found"; exit 1; }

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fasta_file> [k] [m]" >&2
    exit 1
fi

fasta_file="$1"
k="${2:-32}"
m="${3:-31}"

if [ "$m" -ge "$k" ]; then
    echo "Error: m ($m) must be strictly smaller than k ($k)" >&2
    exit 1
fi

prefix="$(basename "${fasta_file%.*}")_${k}_${m}"
tmpdir="tmp_${prefix}"
mkdir -p "$tmpdir"

# Create a sorted skmer list from the input fasta
./bin/sskm construct -f "$fasta_file" -k "$k" -m "$m" --ascii > "${tmpdir}/${prefix}_skmers.txt"

cat "${tmpdir}/${prefix}_skmers.txt" | cut -d' ' -f1 | ./sequences_2_fa.sh > "${tmpdir}/${prefix}_skmers.fa"

# Count uniq kmers using kmc
mkdir -p "${tmpdir}/kmc_tmp" && $kmc -k"$k" -ci1 -fa "${tmpdir}/${prefix}_skmers.fa" "${tmpdir}/kmc_${prefix}" "${tmpdir}/kmc_tmp" && rm -r "${tmpdir}/kmc_tmp"
$kmc_tools transform "${tmpdir}/kmc_${prefix}" dump /dev/stdout | cut -f1 | sort > "${tmpdir}/kmc_${prefix}.txt"
rm "${tmpdir}/kmc_${prefix}.kmc_pre" "${tmpdir}/kmc_${prefix}.kmc_suf"

# Count uniq kmers from input file
mkdir -p "${tmpdir}/kmc_tmp" && $kmc -k"$k" -ci1 -fm "$fasta_file" "${tmpdir}/kmc_raw_${prefix}" "${tmpdir}/kmc_tmp" && rm -r "${tmpdir}/kmc_tmp"
$kmc_tools transform "${tmpdir}/kmc_raw_${prefix}" dump /dev/stdout | cut -f1 | sort > "${tmpdir}/kmc_raw_${prefix}.txt"
rm "${tmpdir}/kmc_raw_${prefix}.kmc_pre" "${tmpdir}/kmc_raw_${prefix}.kmc_suf"

# Compare counts
diff "${tmpdir}/kmc_raw_${prefix}.txt" "${tmpdir}/kmc_${prefix}.txt"
