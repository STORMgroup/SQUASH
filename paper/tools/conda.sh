#!/usr/bin/env bash
# paper/tools/conda.sh NAME: the conda environment envs/NAME.yml (exact versions and builds) into TOOLS_DIR/envs/NAME/.
# NAME: flye, hifiasm, miniasm, samtools, minimap2, quast, remora, build (the compilers of the two source builds).
# An environment that already exists is left as it is (FORCE=1 rebuilds it).
# Settings: TOOLS_DIR; CONDA_CMD (default mamba when on PATH, else conda; micromamba works too).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
NAME=${1:?usage: conda.sh NAME (a file envs/NAME.yml)}
YML=$HERE/envs/$NAME.yml
[[ -s $YML ]] || { echo "conda.sh: no $YML" >&2; exit 2; }
if [[ -z ${CONDA_CMD:-} ]]; then
  if command -v mamba > /dev/null; then CONDA_CMD=mamba; elif command -v conda > /dev/null; then CONDA_CMD=conda
  elif command -v micromamba > /dev/null; then CONDA_CMD=micromamba
  else echo "conda.sh: no conda, mamba or micromamba on PATH (or set CONDA_CMD)" >&2; exit 1; fi
fi
P=$TOOLS_DIR/envs/$NAME
if [[ -e $P/conda-meta/history && ${FORCE:-0} != 1 ]]; then echo "conda.sh: $P exists"; exit 0; fi
rm -rf -- "$P"; mkdir -p "$TOOLS_DIR/envs"
echo "conda.sh: $CONDA_CMD env create -p $P -f $YML"
case $(basename "$CONDA_CMD") in
  micromamba) "$CONDA_CMD" create -y -p "$P" -f "$YML" ;;
  *) "$CONDA_CMD" env create -y -p "$P" -f "$YML" ;;
esac
