#!/usr/bin/env bash
set -euo pipefail

# Simple acceptance test for post-stitch chimera gating semantics
# Usage:
#   ./test.sh /path/to/STAR
# If no argument is provided, the script will try to use 'STAR' from PATH.
#
# This script will create a genome index from chim_lib.fa if it doesn't exist.
# Reads: chim_R1.fa, chim_R2.fa in the same directory as this script

STAR_BIN="${1:-STAR}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

GENOME_DIR="${GENOME_DIR:-${SCRIPT_DIR}/star_index_chimera2}"
GENOME_FA="${SCRIPT_DIR}/chim_lib.fa"
R1="${R1:-${SCRIPT_DIR}/chim_R1.fa}"
R2="${R2:-${SCRIPT_DIR}/chim_R2.fa}"

echo "Using STAR executable: ${STAR_BIN}"

# Create genome index if it doesn't exist
if [ ! -d "${GENOME_DIR}" ] || [ ! -f "${GENOME_DIR}/SA" ]; then
    echo "Creating genome index at ${GENOME_DIR} from ${GENOME_FA}..."
    rm -rf "${GENOME_DIR}"
    mkdir -p "${GENOME_DIR}"
    "${STAR_BIN}" --runMode genomeGenerate \
        --genomeDir "${GENOME_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --genomeSAindexNbases 4 \
        --genomeChrBinNbits 8 \
        --runThreadN 2 > "${GENOME_DIR}/genomeGenerate.log" 2>&1
    echo "Genome index created."
else
    echo "Using existing genome index at ${GENOME_DIR}"
fi

echo "Using genomeDir=${GENOME_DIR}, reads=(${R1} ${R2})"
echo "Old behavior checks (without --chimScoreUsePostStitch):"

rm -rf test_old_21 test_old_22
mkdir -p test_old_21 test_old_22

set +e
"${STAR_BIN}" --runThreadN 2 --genomeDir "${GENOME_DIR}" --readFilesIn "${R1}" "${R2}" \
  --outFileNamePrefix test_old_21/ --outSAMtype BAM Unsorted \
  --chimMultimapNmax 1 --chimSegmentMin 1 --chimScoreDropMax 21 --chimOutType Junctions > test_old_21/run.log 2>&1
rc21=$?
set -e
echo "Exit code (old, DropMax=21): ${rc21} (expected non-zero or no chim junctions)"

"${STAR_BIN}" --runThreadN 2 --genomeDir "${GENOME_DIR}" --readFilesIn "${R1}" "${R2}" \
  --outFileNamePrefix test_old_22/ --outSAMtype BAM Unsorted \
  --chimMultimapNmax 1 --chimSegmentMin 1 --chimScoreDropMax 22 --chimOutType Junctions > test_old_22/run.log 2>&1
echo "Completed old behavior with DropMax=22 (expected chimera present)"

echo
echo "New behavior checks (with --chimScoreUsePostStitch 1):"

rm -rf test_new_6
mkdir -p test_new_6

"${STAR_BIN}" --runThreadN 2 --genomeDir "${GENOME_DIR}" --readFilesIn "${R1}" "${R2}" \
  --outFileNamePrefix test_new_6/ --outSAMtype BAM Unsorted \
  --chimMultimapNmax 1 --chimSegmentMin 1 \
  --chimScoreDropMax 6 \
  --chimScoreUsePostStitch 1 \
  --chimScorePreStitchAllowance 32 \
  --chimOutType Junctions > test_new_6/run.log 2>&1
echo "Completed new behavior with DropMax=6 and chimScoreUsePostStitch=1 (expected chimera present)"

echo
echo "Inspect reports:"
for d in test_old_21 test_old_22 test_new_6; do
  if [[ -f "$d/Chimeric.out.junction" ]]; then
    echo "== $d/Chimeric.out.junction =="
    cat "$d/Chimeric.out.junction"
  else
    echo "== $d: no Chimeric.out.junction produced =="
  fi
done
