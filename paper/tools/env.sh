# paper/tools/env.sh: where the tools of basecall/ and assemble/ (and the evaluation's) are, sourced after config.sh.
# One variable per tool. A value set in the environment, config.sh or config.local.sh wins; otherwise the copy that
# tools/install.sh put in TOOLS_DIR when it is there; otherwise the tool's name, found on PATH.
#   TOOLS_DIR                            where tools/install.sh installs (default paper/tools/opt)
#   DORADO_0_9_2 DORADO_1_4_0 DORADO_2_1_2   the three Dorado releases (bin/dorado of each)
#   DORADO_MODELS                        the folder of the Dorado models (one folder per model, as `dorado download` writes)
#   METAMDBG FLYE HIFIASM RAWHASH2 MINIASM MINIMAP2 SAMTOOLS QUAST   the programs
#   RAWHASH2_PORE                        Rawsamble's R10.4.1 9-mer pore model (RawHash2 -p)
#   REMORA_LEVELS                        the R10.4.1 9-mer levels of Remora's refinement (evaluation)
#   REMORA_PYTHON                        a python with ont-remora 3.3.0 (evaluation)
: "${PAPER:?source paper/config.sh first}"
: "${TOOLS_DIR:=$PAPER/tools/opt}"
export TOOLS_DIR

# tool_default VAR INSTALLED NAME: VAR unset (or set to the bare NAME, config.sh's default) -> INSTALLED if it exists, else NAME
tool_default() {
  local var=$1 path=$2 name=$3
  if [[ -z ${!var:-} || ${!var} == "$name" ]]; then
    if [[ -e $path ]]; then printf -v "$var" '%s' "$path"; else printf -v "$var" '%s' "$name"; fi
  fi
  export "${var?}"
}
tool_default DORADO_0_9_2 "$TOOLS_DIR/dorado-0.9.2-linux-x64/bin/dorado" dorado
tool_default DORADO_1_4_0 "$TOOLS_DIR/dorado-1.4.0-linux-x64/bin/dorado" dorado
tool_default DORADO_2_1_2 "$TOOLS_DIR/dorado-2.1.2-linux-x64/bin/dorado" dorado
: "${DORADO_MODELS:=$TOOLS_DIR/dorado-models}"
tool_default METAMDBG "$TOOLS_DIR/metaMDBG/build/bin/metaMDBG" metaMDBG
tool_default FLYE "$TOOLS_DIR/envs/flye/bin/flye" flye
tool_default HIFIASM "$TOOLS_DIR/envs/hifiasm/bin/hifiasm" hifiasm
tool_default RAWHASH2 "$TOOLS_DIR/RawHash2/build/src/rawhash2" rawhash2
tool_default MINIASM "$TOOLS_DIR/envs/miniasm/bin/miniasm" miniasm
tool_default MINIMAP2 "$TOOLS_DIR/envs/minimap2/bin/minimap2" minimap2
tool_default SAMTOOLS "$TOOLS_DIR/envs/samtools/bin/samtools" samtools
tool_default QUAST "$TOOLS_DIR/envs/quast/bin/quast.py" quast.py
tool_default REMORA_PYTHON "$TOOLS_DIR/envs/remora/bin/python" python3
: "${RAWHASH2_PORE:=$TOOLS_DIR/RawHash2/extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt}"
: "${REMORA_LEVELS:=$TOOLS_DIR/RawHash2/extern/local_kmer_models/uncalled_r1041_model_only_means.txt}"
export DORADO_MODELS RAWHASH2_PORE REMORA_LEVELS

# dorado_bin VERSION: the variable of that Dorado release (0.9.2, 1.4.0, 2.1.2)
dorado_bin() {
  local var=DORADO_${1//./_}
  [[ -n ${!var:-} ]] || { echo "no Dorado $1 (variable $var)" >&2; return 1; }
  echo "${!var}"
}
