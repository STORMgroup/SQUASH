#!/usr/bin/env bash
# paper/tools/install.sh [TOOL ...]: every tool of basecall/, assemble/ and the evaluation, pinned to the paper's versions
# (../versions.tsv), into TOOLS_DIR (default paper/tools/opt; tools/env.sh), then check.sh. TOOL (default: all, in this
# order): dorado-0.9.2 dorado-1.4.0 dorado-2.1.2 flye hifiasm miniasm samtools minimap2 quast remora build rawhash2
# metamdbg. Each step is skipped when its result exists. Needs: bash, curl, git, tar, sha256sum and conda / mamba /
# micromamba (CONDA_CMD); about 21 GB in TOOLS_DIR (Dorado's three releases 13 GB, the models 0.8 GB, the environments
# 6.5 GB).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
ALL=(dorado-0.9.2 dorado-1.4.0 dorado-2.1.2 flye hifiasm miniasm samtools minimap2 quast remora build rawhash2 metamdbg)
(( $# )) && STEPS=("$@") || STEPS=("${ALL[@]}")
for s in "${STEPS[@]}"; do
  echo "=== install.sh: $s ($(date '+%F %T'))"
  case $s in
    dorado-*) bash "$HERE/dorado.sh" "${s#dorado-}" ;;
    rawhash2) bash "$HERE/rawhash2.sh" ;;
    metamdbg) bash "$HERE/metamdbg.sh" ;;
    flye|hifiasm|miniasm|samtools|minimap2|quast|remora|build) bash "$HERE/conda.sh" "$s" ;;
    *) echo "install.sh: unknown tool $s (one of: ${ALL[*]})" >&2; exit 2 ;;
  esac
done
bash "$HERE/check.sh"
