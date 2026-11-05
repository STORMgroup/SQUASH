#!/usr/bin/env bash
# paper/basecall/timed.sh KEY sup|hac gpu|cpu: the runtime tables' basecalling of a paper set, timed as the paper timed it:
# FASTQ only (dorado basecaller --emit-fastq, Dorado's default read splitting; what the assemblers read), Dorado 2.1.2
# sup@v5.2.0 / hac@v5.2.0 (D. melanogaster, 4 kHz: Dorado 0.9.2 sup@v4.1.0 / hac@v4.1.0; basecall/sets.tsv). Before the
# clock (untimed): Dorado (bin/dorado and lib/), the model and the POD5 copied to LOCAL_DIR (the POD5 read in place when
# LOCAL_DIR lacks room for it). Timed: /usr/bin/time -v of Dorado alone.
#   gpu   one GPU (the paper: one NVIDIA RTX A5000), the whole set.
#   cpu   --device cpu, pinned to the first BASECALL_CPU_THREADS CPUs (the paper: 64 cores of an AMD EPYC 7502) by
#         taskset, under `timeout CPU_CAP` (3h). The large sets on a subset of about one hour of signal
#         (cpu_subset.py: reads in MD5 order of their id until 3.67e9 samples at 5 kHz, 2.33e9 at 4 kHz; the bacteria
#         whole: cpu_subset_samples `all`), cut by data/pod5_cut.py before the clock; rate = the subset's samples /
#         Dorado's "Finished in" time, the set's time = the set's samples / rate (set_est_s).
# Writes OUT_DIR/timing/KEY/: basecall_KIND_DEV.time.txt (/usr/bin/time -v), basecall_KIND_DEV.dorado.log,
# basecall_KIND_DEV.tsv (key, host, device, dorado, model, stage, reads_out, bases_out, wall_s; cpu also subset_reads,
# subset_samples, set_reads, set_samples, dorado_ms, samples_per_s, cpu_percent, max_rss_kb, set_est_s); cpu also
# KEY.cpu_subset.ids.txt / .tsv. The FASTQ is not kept (hac.sh writes the HAC reads).
# Settings: POD5_DIR, OUT_DIR, LOCAL_DIR, BASECALL_CPU_THREADS (64), CPU_CAP (3h), CUT_THREADS, PYTHON (pod5 0.3.39);
# DRY_RUN=1 prints the timed command only.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
KEY=${1:-} KIND=${2:-} DEV=${3:-}
[[ -n $KEY && $KIND =~ ^(sup|hac)$ && $DEV =~ ^(gpu|cpu)$ ]] || { echo "usage: timed.sh KEY sup|hac gpu|cpu" >&2; exit 2; }
V=$(bc_col "$KEY" timed_dorado) || { echo "timed.sh: no set $KEY in $BC_HERE/sets.tsv" >&2; exit 2; }
MNAME=$(bc_col "$KEY" "timed_${KIND}_model")
MODEL=$(model_dir "$MNAME")
DROOT=$(cd "$(dirname "$(command -v "$(dorado_bin "$V")")")/.." && pwd)
POD5=$(set_pod5 "$KEY")
T=$OUT_DIR/timing/$KEY NAME=basecall_${KIND}_${DEV}
mkdir -p "$T" "$LOCAL_DIR"
S=$(mktemp -d "$LOCAL_DIR/bc_${KEY}_XXXXXX"); trap 'rm -rf -- "$S"' EXIT
IN=$POD5
if [[ $DEV == cpu && $DRY_RUN == 1 ]]; then
  [[ $(bc_col "$KEY" cpu_subset_samples) == all ]] || IN=$S/subset.pod5
elif [[ $DEV == cpu ]]; then
  target=$(bc_col "$KEY" cpu_subset_samples)
  "$PYTHON" "$BC_HERE/cpu_subset.py" "$KEY" "$target" "$POD5" "$T"
  if [[ $target != all ]]; then
    "$PYTHON" "$PAPER/data/pod5_cut.py" cut "$T/$KEY.cpu_subset.ids.txt" "$S/subset.pod5" "$POD5" --procs "$CUT_THREADS"
    IN=$S/subset.pod5
  fi
fi
# staging (untimed): Dorado, the model, the POD5
STAGE=local
if [[ $DRY_RUN == 1 ]]; then [[ $IN == "$S/subset.pod5" ]] || IN=$S/reads.pod5
else
mkdir -p "$S/dorado/bin"
cp "$DROOT/bin/dorado" "$S/dorado/bin/"; cp -r "$DROOT/lib" "$S/dorado/"; cp -r "$MODEL" "$S/dorado/bin/"
need=$(( $(stat -L -c %s "$IN") / 1000000000 * 3 / 2 + 10 ))   # the POD5 + the FASTQ (< 0.5 x) + 10 GB
free=$(df --output=avail -BG "$S" | tail -1 | tr -dc 0-9)
if [[ $IN == "$S/subset.pod5" ]]; then :
elif (( free >= need )); then cp "$IN" "$S/reads.pod5"; IN=$S/reads.pod5
else STAGE=in_place; fi
fi
cmd=("$S/dorado/bin/dorado" basecaller)
if [[ $DEV == cpu ]]; then
  cpus=$(first_cpus "$BASECALL_CPU_THREADS")
  [[ -n $cpus || $DRY_RUN == 1 ]] || { echo "timed.sh: fewer than $BASECALL_CPU_THREADS CPUs allowed" >&2; exit 1; }
  [[ -n $cpus ]] || cpus="0-$((BASECALL_CPU_THREADS - 1))"
  cmd=(timeout "${CPU_CAP:-3h}" taskset -c "$cpus" "${cmd[@]}" --device cpu)
fi
cmd+=(--emit-fastq "$S/dorado/bin/$MNAME" "$IN")
echo "timed.sh: $KEY $KIND $DEV on $(hostname): ${cmd[*]} (stage $STAGE)"
[[ $DRY_RUN != 1 ]] || exit 0
( cd "$S" && /usr/bin/time -v -o "$T/$NAME.time.txt" "${cmd[@]}" > "$S/reads.fastq" 2> "$T/$NAME.dorado.log" )
read -r nreads nbases < <(awk 'NR % 4 == 2 {n++; b += length($0)} END {print n + 0, b + 0}' "$S/reads.fastq")
host=$(hostname)
if [[ $DEV == gpu ]]; then
  dev=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || true)
else
  dev=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)" x $BASECALL_CPU_THREADS"
fi
python3 - "$KEY" "$T/$NAME" "$host" "${dev:-unknown}" "$V" "$MNAME" "$STAGE" "$nreads" "$nbases" "$DEV" "$T/$KEY.cpu_subset.tsv" <<'PY'
import re, sys
key, base, host, dev, dorado, model, stage, nreads, nbases, kind, subset = sys.argv[1:12]
tv = open(base + ".time.txt").read()
def tf(label):
    m = re.search(rf"{label}: (.+)", tv)
    return m.group(1).strip() if m else ""
wall = 0.0
for x in tf(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\)").split(":"):
    wall = wall * 60 + float(x)
hdr = ["key", "host", "device", "dorado", "model", "stage", "reads_out", "bases_out", "wall_s"]
row = [key, host, dev, dorado, model, stage, nreads, nbases, f"{wall:.0f}"]
if kind == "cpu":
    sub = dict(zip(*[l.rstrip("\n").split("\t") for l in open(subset)]))
    ms = re.findall(r"Finished in \(ms\): (\d+)", open(base + ".dorado.log").read())
    if not ms:
        sys.exit("timed.sh: no 'Finished in' in the Dorado log")
    rate = int(sub["subset_samples"]) / (int(ms[-1]) / 1000)
    hdr += ["subset_reads", "subset_samples", "set_reads", "set_samples", "dorado_ms", "samples_per_s", "cpu_percent",
            "max_rss_kb", "set_est_s"]
    row += [sub["subset_reads"], sub["subset_samples"], sub["set_reads"], sub["set_samples"], ms[-1], f"{rate:.0f}",
            tf("Percent of CPU this job got").rstrip("%"), tf(r"Maximum resident set size \(kbytes\)"),
            f"{int(sub['set_samples']) / rate:.0f}"]
open(base + ".tsv", "w").write("\t".join(hdr) + "\n" + "\t".join(row) + "\n")
print("\t".join(f"{a}={b}" for a, b in zip(hdr, row)))
PY
