# paper/basecall/lib.sh: sourced by the basecalling scripts (config.sh, tools/env.sh, a set's row of basecall/sets.tsv).
set -euo pipefail
BC_HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$BC_HERE/../config.sh"
source "$BC_HERE/../tools/env.sh"
# a local disk for staged inputs and temporary files (default TMPDIR, else /tmp)
: "${LOCAL_DIR:=${TMPDIR:-/tmp}}"
# CPU threads of the CPU basecalling route (the paper's: 64) and of the BAM / FASTQ steps
: "${BASECALL_CPU_THREADS:=64}"
: "${BC_THREADS:=8}"
# DRY_RUN=1: print the Dorado command of each step instead of running it
: "${DRY_RUN:=0}"
export LOCAL_DIR BASECALL_CPU_THREADS BC_THREADS DRY_RUN

# bc_col KEY COLUMN: a value of basecall/sets.tsv
bc_col() {
  awk -F'\t' -v k="$1" -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h[c]); f = 1}
    END {exit f ? 0 : 1}' "$BC_HERE/sets.tsv"
}
# set_pod5 KEY: the set's POD5 (data/extract.sh KEY makes it)
set_pod5() {
  local p=$POD5_DIR/$1.pod5
  [[ -s $p || $DRY_RUN == 1 ]] || { echo "no $p (make -C paper data SET=$1)" >&2; return 1; }
  echo "$p"
}
# model_dir NAME: the model's folder in DORADO_MODELS (tools/dorado.sh downloads it)
model_dir() {
  [[ -s $DORADO_MODELS/$1/config.toml || $DRY_RUN == 1 ]] || { echo "no model $1 in $DORADO_MODELS (tools/install.sh)" >&2; return 1; }
  echo "$DORADO_MODELS/$1"
}
# first_cpus N: the first N CPUs this process may use, comma-separated (empty when fewer are allowed)
first_cpus() {
  local all
  all=$(taskset -pc $$ 2>/dev/null | sed 's/.*: //' | tr ',' '\n' | awk -F- 'NF == 2 {for (i = $1; i <= $2; i++) print i; next} {print}')
  (( $(wc -l <<< "$all") >= $1 )) && head -n "$1" <<< "$all" | paste -sd, || true
}
# fastq_lines: samtools fastq output as `@id`, sequence, `+`, qualities (the paper's FASTQ format: the header's first word)
fastq_lines() { awk 'NR % 4 == 1 {split($1, a, /[ \t]/); print a[1]; next} NR % 4 == 3 {print "+"; next} {print}'; }
# gz_out FILE: stdin gzipped into FILE (pigz when on PATH), written as FILE.part and renamed
gz_out() {
  if command -v pigz > /dev/null; then pigz -p "$BC_THREADS" -c > "$1.part"; else gzip -c > "$1.part"; fi
  mv -- "$1.part" "$1"
}
