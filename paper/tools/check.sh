#!/usr/bin/env bash
# paper/tools/check.sh [TOOL ...]: each tool (tools/env.sh's variables) prints the version the paper used, and the Dorado
# models and the two 9-mer tables equal the paper's files. One line per tool: status (same / DIFFERS / MISSING), tool,
# expected, printed, path. Writes the same lines to TOOLS_DIR/versions.tsv when every tool is checked; exit 1 when any
# differs or is missing. TOOL: dorado-0.9.2 dorado-1.4.0 dorado-2.1.2 metamdbg flye hifiasm rawhash2 miniasm minimap2
# samtools quast remora models (default: all).
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
ALL=(dorado-0.9.2 dorado-1.4.0 dorado-2.1.2 metamdbg flye hifiasm rawhash2 miniasm minimap2 samtools quast remora models)
(( $# )) && TOOLS=("$@") || TOOLS=("${ALL[@]}")
bad=0 out=()
line() {  # line TOOL EXPECTED PRINTED PATH: same when PRINTED starts with EXPECTED
  local st=same
  if [[ -z $3 ]]; then st=MISSING; elif [[ $3 != "$2"* ]]; then st=DIFFERS; fi
  [[ $st == same ]] || bad=$((bad + 1))
  out+=("$(printf '%s\t%s\t%s\t%s\t%s' "$st" "$1" "$2" "$3" "$4")")
  printf '%-8s %-14s expected %-22s printed %-28s %s\n' "$st" "$1" "$2" "$3" "$4"
}
v() {  # the last non-empty line a command prints (nothing when the command is not there)
  command -v "$1" > /dev/null 2>&1 || return 0
  "$@" 2>&1 | tr -d '\r' | grep -v '^\s*$' | tail -1 | xargs
}
for t in "${TOOLS[@]}"; do
  case $t in
    dorado-*) b=$(dorado_bin "${t#dorado-}" 2>/dev/null) || b=""
              line "$t" "${t#dorado-}+" "$( [[ -n $b ]] && v "$b" --version)" "$b" ;;
    metamdbg) line metamdbg "1.2" "$("$METAMDBG" --version 2>&1 | sed -n 's/.*Version: *//p' | head -1)" "$METAMDBG" ;;
    flye)     line flye "2.9.6-b1802" "$(v "$FLYE" --version)" "$FLYE" ;;
    hifiasm)  line hifiasm "0.25.0-r726" "$(v "$HIFIASM" --version)" "$HIFIASM" ;;
    rawhash2) c=$(git -C "$(dirname "$(command -v "$RAWHASH2")")/../.." rev-parse --short=7 HEAD 2>/dev/null || echo "?")
              r=$(v "$RAWHASH2" --version); line rawhash2 "2.5 (commit 92366af)" "${r:+$r (commit $c)}" "$RAWHASH2" ;;
    miniasm)  line miniasm "0.3-r179" "$(v "$MINIASM" -V)" "$MINIASM" ;;
    minimap2) line minimap2 "2.30-r1287" "$(v "$MINIMAP2" --version)" "$MINIMAP2" ;;
    samtools) line samtools "samtools 1.23.1 htslib 1.23.1" \
                "$("$SAMTOOLS" --version 2>/dev/null | sed -n '1p;2p' | sed 's/Using //' | xargs)" "$SAMTOOLS" ;;
    quast)    line quast "QUAST v5.2.0" "$(v "$QUAST" --version)" "$QUAST" ;;
    remora)   line remora "3.3.0" "$("$REMORA_PYTHON" -c 'import remora; print(remora.__version__)' 2>/dev/null)" "$REMORA_PYTHON" ;;
    models)
      n=$(grep -vc '^#' "$HERE/models.sha256")
      if (cd "$DORADO_MODELS" 2>/dev/null && grep -v '^#' "$HERE/models.sha256" | sha256sum -c --quiet > /dev/null 2>&1); then
        line models "$n files" "$n files" "$DORADO_MODELS"; else line models "$n files" "differ" "$DORADO_MODELS"; fi
      for f in "RAWHASH2_PORE e5e154511d50c7e288c020a933fa66cc15d87ec06f903baa739e0a5c390cf23c" \
               "REMORA_LEVELS 0b28b79af1d63443aabdcf937834597d5d170ae91426a4638a7a033d0084ecf2"; do
        read -r var want <<< "$f"
        got=$(sha256sum "${!var}" 2>/dev/null | cut -d' ' -f1)
        line "$var" "sha256 ${want:0:12}" "${got:+sha256 ${got:0:12}}" "${!var}"
      done ;;
    *) echo "check.sh: unknown tool $t (one of: ${ALL[*]})" >&2; exit 2 ;;
  esac
done
if (( $# == 0 )) && [[ -d $TOOLS_DIR ]]; then
  { printf 'status\ttool\texpected\tprinted\tpath\n'; printf '%s\n' "${out[@]}"; } > "$TOOLS_DIR/versions.tsv"
fi
(( bad == 0 )) || { echo "check.sh: $bad tool(s) differ from the paper's or are missing"; exit 1; }
echo "check.sh: every tool as the paper's"
