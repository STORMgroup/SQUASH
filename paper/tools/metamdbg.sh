#!/usr/bin/env bash
# paper/tools/metamdbg.sh: nanoMDBG = metaMDBG 1.2 (`metaMDBG --version`), the paper's build: the GaetanBenoitDev/metaMDBG
# repository at tag metaMDBG-1.2 (commit 6a9e7ee), CMake Release build, GCC 15.1.0 of envs/build.yml (its libraries
# found through the binary's RPATH; tools/conda.sh build first). -> TOOLS_DIR/metaMDBG/ (binary build/bin/metaMDBG).
# metaMDBG calls minimap2 from PATH: the paper's runs had minimap2 2.30-r1287 there, and assemble/nanomdbg.sh puts
# MINIMAP2's folder first on PATH (bioconda's metamdbg 1.2 package requires minimap2 2.28, hence the source build).
# Settings: TOOLS_DIR, BUILD_THREADS (8).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
REPO=https://github.com/GaetanBenoitDev/metaMDBG
TAG=metaMDBG-1.2 COMMIT=6a9e7ee8c483a0029c5dc2274faf94646ba870e1
BENV=$TOOLS_DIR/envs/build
[[ -x $BENV/bin/x86_64-conda-linux-gnu-g++ ]] || { echo "metamdbg.sh: no build environment in $BENV (tools/conda.sh build)" >&2; exit 1; }
export PATH=$BENV/bin:$PATH
S=$TOOLS_DIR/metaMDBG
[[ -d $S/.git ]] || git clone "$REPO" "$S"
git -C "$S" checkout -q "$TAG"
[[ $(git -C "$S" rev-parse HEAD) == "$COMMIT" ]] || { echo "metamdbg.sh: tag $TAG is not commit $COMMIT" >&2; exit 1; }
if [[ ! -x $S/build/bin/metaMDBG ]]; then
  cmake -S "$S" -B "$S/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$BENV" \
    -DCMAKE_C_COMPILER="$BENV/bin/x86_64-conda-linux-gnu-gcc" -DCMAKE_CXX_COMPILER="$BENV/bin/x86_64-conda-linux-gnu-g++" \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$BENV/lib"
  cmake --build "$S/build" -j "${BUILD_THREADS:-8}"
fi
"$S/build/bin/metaMDBG" --version 2>&1 | grep -i version
