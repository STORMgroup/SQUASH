# paper/assemble/lib.sh: sourced by the assembler scripts. The paper's benchmark protocol for every tool: the tool pinned
# to the first ASM_THREADS CPUs (taskset; without pinning when fewer are allowed), its input staged to a local disk
# before the clock (untimed, logs/staging.time.txt), its working and output folders on that disk, /usr/bin/time -v of
# each step of the tool alone (logs/STEP.time.txt), the outputs copied to OUT_DIR/assembly/KEY/TOOL/ afterwards (untimed).
#   ASM_THREADS  threads of every assembler (the paper's runs: 64)
#   LOCAL_DIR    the local disk for the working folders (default TMPDIR, else /tmp)
#   STAGE_DIR    where the staged input goes (default LOCAL_DIR; the paper staged the large sets' inputs in RAM: /dev/shm)
#   SUP_FASTQ, HAC_FASTQ   the reads (default OUT_DIR/reads/KEY/sup.fastq.gz, hac.fastq.gz from basecall/; plain or .gz)
#   SET_POD5     Rawsamble's input (default POD5_DIR/KEY.pod5)
#   DRY_RUN=1    print each timed command instead of running it
set -euo pipefail
AS_HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$AS_HERE/../config.sh"
source "$AS_HERE/../tools/env.sh"
: "${ASM_THREADS:=64}"
: "${LOCAL_DIR:=${TMPDIR:-/tmp}}"
: "${STAGE_DIR:=$LOCAL_DIR}"
: "${DRY_RUN:=0}"
export ASM_THREADS LOCAL_DIR STAGE_DIR

# as_col KEY COLUMN: a value of assemble/options.tsv
as_col() {
  awk -F'\t' -v k="$1" -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h[c]); f = 1}
    END {exit f ? 0 : 1}' "$AS_HERE/options.tsv"
}
# as_pin SCRIPT ARGS...: run the calling script again under taskset -c <first ASM_THREADS CPUs> (once)
as_pin() {
  [[ -n ${AS_PINNED:-} || $DRY_RUN == 1 ]] && return 0
  local all cpus
  all=$(taskset -pc $$ 2>/dev/null | sed 's/.*: //' | tr ',' '\n' | awk -F- 'NF == 2 {for (i = $1; i <= $2; i++) print i; next} {print}')
  if (( $(wc -l <<< "$all") >= ASM_THREADS )); then
    cpus=$(head -n "$ASM_THREADS" <<< "$all" | paste -sd,)
    AS_PINNED=$cpus exec taskset -c "$cpus" bash "$@"
  fi
  export AS_PINNED="none ($(wc -l <<< "$all") CPUs allowed, fewer than $ASM_THREADS)"
}
# as_start KEY TOOL: O (the output folder), LOCAL (working folder), IN (staged input folder); skip when O/done exists
as_start() {
  KEY=$1 TOOL=$2
  as_col "$KEY" key > /dev/null || { echo "$TOOL: no set $KEY in $AS_HERE/options.tsv" >&2; exit 2; }
  O=$OUT_DIR/assembly/$KEY/$TOOL
  if [[ $DRY_RUN == 1 ]]; then   # nothing run: the commands go to fd 3 (stdout), the logs to a temporary folder
    O=$(mktemp -d "${TMPDIR:-/tmp}/asm_dry.XXXXXX"); mkdir -p "$O/logs" "$O/asm"; trap 'rm -rf -- "$O"' EXIT
    cd "$O" && mkdir -p LOCAL STAGE; LOCAL=LOCAL IN=STAGE; exec 3>&1; return 0
  fi
  [[ ! -e $O/done ]] || { echo "$TOOL: $O done before"; exit 0; }
  rm -rf -- "$O"; mkdir -p "$O/logs" "$O/asm" "$LOCAL_DIR" "$STAGE_DIR"
  LOCAL=$(mktemp -d "$LOCAL_DIR/asm.$TOOL.$KEY.XXXXXX")
  IN=$LOCAL; [[ $STAGE_DIR == "$LOCAL_DIR" ]] || IN=$(mktemp -d "$STAGE_DIR/asm.$TOOL.$KEY.XXXXXX")
  trap 'rm -rf -- "$LOCAL"; [[ $IN == "$LOCAL" ]] || rm -rf -- "$IN"' EXIT
  echo "=== $TOOL $KEY on $(hostname) ($(nproc) CPUs allowed, pinned: ${AS_PINNED:-no}), $ASM_THREADS threads, local $LOCAL ($(df --output=avail -BG "$LOCAL" | tail -1 | tr -dc 0-9) GB free), input in $IN $(date '+%F %T')"
}
# timed NAME CMD...: /usr/bin/time -v of CMD -> O/logs/NAME.time.txt (DRY_RUN=1: the command printed)
timed() {
  local name=$1; shift
  if [[ $DRY_RUN == 1 ]]; then echo "$name: $*" >&3; return 0; fi
  /usr/bin/time -v -o "$O/logs/$name.time.txt" "$@"
}
# stage_reads SRC: the reads as one plain FASTQ, IN/reads.fastq (LOCAL/reads.fastq links to it), untimed
stage_reads() {
  local src=$1
  [[ $DRY_RUN == 1 ]] && return 0
  [[ -s $src ]] || { echo "$TOOL: no reads $src (basecall/sup.sh, hac.sh $KEY)" >&2; exit 1; }
  if [[ $src == *.gz ]]; then
    /usr/bin/time -v -o "$O/logs/staging.time.txt" bash -c 'gzip -dc "$1" > "$2"' _ "$src" "$IN/reads.fastq"
  else
    /usr/bin/time -v -o "$O/logs/staging.time.txt" cp -L "$src" "$IN/reads.fastq"
  fi
  [[ $IN == "$LOCAL" ]] || ln -s "$IN/reads.fastq" "$LOCAL/reads.fastq"
}
# as_done: the run's marker (O/done)
as_done() { date '+%F %T' > "$O/done"; echo "=== $TOOL $KEY done $(date '+%F %T')"; }
