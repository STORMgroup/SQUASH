#!/usr/bin/env bash
# paper/tools/rawhash2.sh: RawHash2 (Rawsamble) at commit 92366af (`rawhash2 --version`: 2.5), the paper's build: the
# STORMgroup/RawHash2 repository at that commit with the submodules a POD5-only build needs (kmer_models, zstd; at the
# commits the superproject pins), CMake Release build with POD5 support only (the default; CMake downloads POD5's C
# library 0.3.36), GCC 15.1.0 of envs/build.yml (its libraries found through the binary's RPATH; tools/conda.sh build
# first). -> TOOLS_DIR/RawHash2/ (binary build/src/rawhash2).
# Also checks the two R10.4.1 9-mer tables taken from it:
#   extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt   Rawsamble's pore model (-p, RAWHASH2_PORE)
#   extern/local_kmer_models/uncalled_r1041_model_only_means.txt   the levels of Remora's refinement (REMORA_LEVELS)
# The build uses -march=native: build on (a node like) the one that runs it. Settings: TOOLS_DIR, BUILD_THREADS (8).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
REPO=https://github.com/STORMgroup/RawHash2
COMMIT=92366afb09df9fa791e7705414e1372f487912ed
PORE_SHA=e5e154511d50c7e288c020a933fa66cc15d87ec06f903baa739e0a5c390cf23c
LEVELS_SHA=0b28b79af1d63443aabdcf937834597d5d170ae91426a4638a7a033d0084ecf2
BENV=$TOOLS_DIR/envs/build
[[ -x $BENV/bin/x86_64-conda-linux-gnu-g++ ]] || { echo "rawhash2.sh: no build environment in $BENV (tools/conda.sh build)" >&2; exit 1; }
export PATH=$BENV/bin:$PATH
S=$TOOLS_DIR/RawHash2
if [[ ! -d $S/.git ]]; then
  git clone "$REPO" "$S"
fi
git -C "$S" checkout -q "$COMMIT"
git -C "$S" submodule update --init extern/kmer_models extern/zstd
[[ $(git -C "$S" rev-parse HEAD) == "$COMMIT" ]] || { echo "rawhash2.sh: $S is not at $COMMIT" >&2; exit 1; }
if [[ ! -x $S/build/src/rawhash2 ]]; then
  cmake -S "$S" -B "$S/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$BENV/bin/x86_64-conda-linux-gnu-gcc" -DCMAKE_CXX_COMPILER="$BENV/bin/x86_64-conda-linux-gnu-g++" \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$BENV/lib"
  cmake --build "$S/build" -j "${BUILD_THREADS:-8}"
fi
"$S/build/src/rawhash2" --version
for f in "extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt $PORE_SHA" \
         "extern/local_kmer_models/uncalled_r1041_model_only_means.txt $LEVELS_SHA"; do
  read -r p want <<< "$f"
  got=$(sha256sum "$S/$p" | cut -d' ' -f1)
  [[ $got == "$want" ]] || { echo "rawhash2.sh: $p: SHA-256 $got, expected $want" >&2; exit 1; }
  echo "rawhash2.sh: $p as the paper's"
done
