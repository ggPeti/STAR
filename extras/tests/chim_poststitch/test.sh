#!/usr/bin/env bash
set -euo pipefail

STAR_BIN="${1:-STAR}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

GENOME_DIR="${GENOME_DIR:-${SCRIPT_DIR}/star_index_chimera2}"
GENOME_FA="${SCRIPT_DIR}/chim_lib.fa"
R1="${R1:-${SCRIPT_DIR}/chim_R1.fa}"
R2="${R2:-${SCRIPT_DIR}/chim_R2.fa}"

echo "Using STAR executable: ${STAR_BIN}"

# Create genome index
echo "Creating genome index at ${GENOME_DIR} from ${GENOME_FA}..."
rm -rf "${GENOME_DIR}"
mkdir -p "${GENOME_DIR}"
"${STAR_BIN}" --runMode genomeGenerate \
    --genomeDir "${GENOME_DIR}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --genomeSAindexNbases 10 \
    --genomeChrBinNbits 8 \
    --runThreadN 2 > "${GENOME_DIR}/genomeGenerate.log" 2>&1
echo "Genome index created."

rm -rf test_new_6
mkdir -p test_new_6

"${STAR_BIN}" --runThreadN 2 --genomeDir "${GENOME_DIR}" --readFilesIn "${R1}" "${R2}" \
  --outFileNamePrefix test_new_6/ --outSAMtype BAM Unsorted \
  --chimMultimapNmax 10 --chimSegmentMin 1 \
  --chimMultimapScoreRange 0 \
  --chimNonchimScoreDropMin 30 \
  --outFilterScoreMin 100 \
  --chimScoreMin 100 \
  --chimScoreDropMax 180 \
  --chimScoreUsePostStitch 1 \
  --chimJunctionOverhangMin 0 \
  --chimScorePreStitchAllowance 50 \
  --chimScoreJunctionNonGTAG 0 \
  --chimAllowSubsetTranscripts 1
