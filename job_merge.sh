#!/bin/bash
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -N garfield_merge
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp
set -euo pipefail

cd "$PBS_O_WORKDIR"
TAG=${RUN_TAG:?RUN_TAG required}
OUTDIR="${OUTBASE:-root}/${TAG}"

module purge || true
source /np1a/phanes/opt/root/6.32.06/bin/thisroot.sh

hadd -f "${OUTDIR}/grid_times.merged.root" "${OUTDIR}"/grid_times.shard*.root
echo "[MERGE] done: ${OUTDIR}/grid_times.merged.root"

====
echo "array job id: ${ARRAY_ID}"

# 配列の全要素がOK後にマージ
qsub -W depend=afterokarray:${ARRAY_ID} \
  -N garf_merge_${TAG} \
  -v RUN_TAG=${TAG},OUTBASE=root \
  job_merge.sh