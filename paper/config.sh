# paper/config.sh: the settings of paper/'s scripts, read by every one of them. Edit the values here, or set them in
# the environment (the environment wins). Settings of one machine that should not be committed go in
# config.local.sh next to this file (git-ignored), which is read first and wins over both.
PAPER=${PAPER:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
[[ -f $PAPER/config.local.sh ]] && source "$PAPER/config.local.sh"

# The squash binary (`make` at the repository's root builds it here).
: "${SQUASH:=$PAPER/../target/release/squash}"
# Downloaded sources (data/download/), the rebuilt set POD5s (data/extract.sh) and SQUASH's outputs (squash.sh).
: "${SOURCES_DIR:=$PAPER/data/sources}"
: "${POD5_DIR:=$PAPER/data/pod5}"
: "${OUT_DIR:=$PAPER/out}"
# SQUASH's large temporary files (--tmp-dir; a local disk helps on the large sets); empty: OUT_DIR/KEY/work/tmp.
: "${TMP_DIR:=}"
# Threads for SQUASH (its outputs do not depend on them) and the second count of `make threads`; workers of the POD5 cut.
: "${THREADS:=$(nproc)}"
: "${THREADS2:=4}"
: "${CUT_THREADS:=4}"
# A python with requirements.txt's packages (or environment.yml): the POD5 cut and digest, the FAST5 conversion (pod5
# 0.3.39), the tables and figures (matplotlib).
: "${PYTHON:=python3}"
# The download command.
: "${CURL:=curl -sSfL --retry 5 --retry-delay 30 --connect-timeout 60}"
# An existing copy of a source instead of a download: SRC_<SOURCE>=FOLDER, SOURCE as in data/sources.tsv (upper case,
# '-' as '_'), a folder of the source's POD5 files; see data/download/lib.sh. For example:
#   SRC_NVB_CJ=/data/nanovarbench/ATCC_33560__202309
# curation/ (remaking the read lists; not needed to rebuild the sets): its working folder; the tools of the two
# pool lists (dorado 1.4.0 with its hac v5.2.0 model, minimap2 2.30-r1287, samtools 1.23.1); the folder of the six
# D. melanogaster sequencing summaries (dmel_bkim2023.sh SUMMARIES=1 downloads them).
: "${CURATION_DIR:=$OUT_DIR/curation}"
: "${DORADO_MODEL:=dna_r10.4.1_e8.2_400bps_hac@v5.2.0}"
: "${MINIMAP2:=minimap2}"
: "${SAMTOOLS:=samtools}"
: "${DMEL_SUMMARIES:=$SOURCES_DIR/dmel_bkim2023/summaries}"
export PAPER SQUASH SOURCES_DIR POD5_DIR OUT_DIR TMP_DIR THREADS THREADS2 CUT_THREADS PYTHON CURL
export CURATION_DIR DORADO DORADO_MODEL MINIMAP2 SAMTOOLS DMEL_SUMMARIES
# The tools of basecall/, assemble/ and eval/ (one variable per tool, default the copy tools/install.sh made, else PATH)
# and their run settings: tools/env.sh.
source "$PAPER/tools/env.sh"
# The evaluation (eval/): its threads, temporary files and references (made by eval/refs.sh).
: "${EVAL_THREADS:=$THREADS}"
: "${EVAL_TMP:=$OUT_DIR/eval/tmp}"
: "${EVAL_REFS:=$SOURCES_DIR/references/eval}"
export EVAL_THREADS EVAL_TMP EVAL_REFS
: "${DORADO:=$DORADO_1_4_0}"
export DORADO
