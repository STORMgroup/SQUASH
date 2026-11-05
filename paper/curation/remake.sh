#!/usr/bin/env bash
# paper/curation/remake.sh KEY [KEY ...]: remake a paper set's read list from its public source with select_reads.py and
# check it against the committed one (../data/sets/KEY/reads.txt.gz, compared as a set: that list is in POD5 order,
# select_reads.py's in md5 order) and the counts of the set's README (reads seen, dropped by the noise rule, kept,
# signal-estimated bases).
#   1. the source (../data/download/SOURCE.sh, as ../data/extract.sh runs it; or the local copy SRC_<SOURCE>), with
#      every read select_reads.py needs: D. melanogaster every read and the runs' sequencing summaries (SUMMARIES=1; a source
#      folder cut to a set's reads is refused), human the reads of the pool list (hg002_giab2025.sh LIST);
#   2. select_reads.py with the set's pool and rules (below), G and coverage from ../sets.tsv -> CURATION_DIR/KEY/; a set
#      whose pool is in another set's table (CURATION_DIR/OTHER/selection.tsv: the same pool, or a longer prefix of the
#      same order) is selected from that table, without reading the signal again;
#   3. the check (exit 1 when the list or a count differs).
# Pools and rules (every set: reads >= 10 kb signal-estimated at 400 bases/s, noise cut 1.15 x the pool's median, md5
# order to COVERAGE x G):
#   eck60, eck100     the deposit; pool = md5-order prefix of its >= 10-kb reads to 150 x G
#   cjk60, vpk60      the deposit; pool prefix to 100 x G        cjk100, vpk100: to 160 x G
#   u9k60, u9k100     the run, reads of ../data/pools/sdys/ut9728.txt.gz; pool prefix to 100 x G (k60) / 160 x G (k100)
#   ark60, ark100     the seven strain deposits; every >= 10-kb read
#   awk60             the run; pool prefix to 100 x G
#   dmk60, dmk100     the six runs, reads with a summary duration x 400 >= 9 kb; lengths from durations at 6 digits
#   h2k60             the reads of ../data/pools/hs13q/hg002_chr13q.txt.gz; lengths from durations at 6 digits
# Settings (paper/config.sh, config.local.sh or the environment): SOURCES_DIR, OUT_DIR, PYTHON (pod5 0.3.39), THREADS;
# CURATION_DIR (default OUT_DIR/curation); DMEL_SUMMARIES (a folder with the six D. melanogaster sequencing summaries,
# default SOURCES_DIR/dmel_bkim2023/summaries).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PAPER=$(cd "$HERE/.." && pwd)
[[ -f $PAPER/config.sh ]] && source "$PAPER/config.sh"
SETS_TSV=${SETS_TSV:-$PAPER/sets.tsv} SOURCES_DIR=${SOURCES_DIR:-$PAPER/data/sources} OUT_DIR=${OUT_DIR:-$PAPER/out}
CURATION_DIR=${CURATION_DIR:-$OUT_DIR/curation} PYTHON=${PYTHON:-python3} THREADS=${THREADS:-$(nproc)}
POOLS=$PAPER/data/pools DL=$PAPER/data/download
export SOURCES_DIR PYTHON
(( $# >= 1 )) || { echo "usage: $0 KEY [KEY ...] (keys: column 1 of $SETS_TSV)" >&2; exit 2; }

col() {  # col KEY COLUMN: a value of sets.tsv
  awk -F'\t' -v k="$1" -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h[c]); f = 1}
    END {exit f ? 0 : 1}' "$SETS_TSV"
}

for KEY in "$@"; do
  src=$(col "$KEY" sources) || { echo "remake.sh: no set $KEY in $SETS_TSV" >&2; exit 2; }
  G=$(col "$KEY" G) COV=$(col "$KEY" coverage)
  OUT=$CURATION_DIR/$KEY
  # the set's pool and rules; FROM = sets whose table holds this set's pool
  case $KEY in
    eck60_both60) OPTS=(--pool-coverage 150) FROM=eck100_both100 ;;
    eck100_both100) OPTS=(--pool-coverage 150) FROM=eck60_both60 ;;
    cjk60_both60) OPTS=(--pool-coverage 100) FROM=cjk100_both100 ;;
    vpk60_both60) OPTS=(--pool-coverage 100) FROM=vpk100_both100 ;;
    cjk100_both100|vpk100_both100) OPTS=(--pool-coverage 160) FROM= ;;
    u9k60_both60) OPTS=(--pool-ids "$POOLS/sdys/ut9728.txt.gz" --pool-coverage 100) FROM=u9k100_both100 ;;
    u9k100_both100) OPTS=(--pool-ids "$POOLS/sdys/ut9728.txt.gz" --pool-coverage 160) FROM= ;;
    ark60_both60) OPTS=() FROM=ark100_both100 ;;
    ark100_both100) OPTS=() FROM=ark60_both60 ;;
    awk60_both60) OPTS=(--pool-coverage 100) FROM= ;;
    dmk60_both60) OPTS=(--summary-min-length 9000 --duration-digits 6) FROM=dmk100_both100 ;;
    dmk100_both100) OPTS=(--summary-min-length 9000 --duration-digits 6) FROM=dmk60_both60 ;;
    h2k60_both60) OPTS=(--pool-ids "$POOLS/hs13q/hg002_chr13q.txt.gz" --duration-digits 6) FROM= ;;
    *) echo "remake.sh: no pool rule for $KEY" >&2; exit 2 ;;
  esac
  echo "=== $KEY: G $G, ${COV}x, $src ($(date '+%F %T'))"
  table=
  if [[ -n $FROM && -s $CURATION_DIR/$FROM/selection.tsv ]]; then table=$CURATION_DIR/$FROM/selection.tsv; fi
  if [[ -n $table ]]; then
    "$PYTHON" "$HERE/select_reads.py" --table "$table" -G "$G" --coverage "$COV" "${OPTS[@]}" --out "$OUT"
  else
    case $src in
      hg002_giab2025)
        d=$SOURCES_DIR/$src
        if [[ -e $d/.complete && -s $d/.keep_ids.txt ]] && \
           [[ -n $(comm -23 <(gzip -dc "$POOLS/hs13q/hg002_chr13q.txt.gz") "$d/.keep_ids.txt" | head -1) ]]; then
          echo "remake.sh: $d was cut for a list without every pool read: use another SOURCES_DIR" >&2; exit 1
        fi
        bash "$DL/$src.sh" "$POOLS/hs13q/hg002_chr13q.txt.gz" ;;
      dmel_bkim2023)
        d=$SOURCES_DIR/$src SUMM=${DMEL_SUMMARIES:-$SOURCES_DIR/dmel_bkim2023/summaries}
        if [[ -s $d/.keep_ids.txt ]]; then
          echo "remake.sh: $d holds only some sets' reads: the selection needs every read" >&2; exit 1
        fi
        if ! ls "$SUMM"/*sequencing_summary* > /dev/null 2>&1 && [[ ! -L $d && -e $d/.complete ]]; then
          rm -f -- "$d/.complete"  # a finished download without the summaries: fetch them (not the FAST5s again)
        fi
        SUMMARIES=1 bash "$DL/$src.sh"
        OPTS+=(--summary "$SUMM"/*sequencing_summary*) ;;
      *) bash "$DL/$src.sh" ;;
    esac
    "$PYTHON" "$HERE/select_reads.py" --pod5 "$SOURCES_DIR/$src" -G "$G" --coverage "$COV" "${OPTS[@]}" --out "$OUT" \
      --threads "$THREADS"
  fi
  "$PYTHON" - "$PAPER/data/sets/$KEY" "$OUT" <<'PY'
import gzip, re, sys
from pathlib import Path
sets, out = Path(sys.argv[1]), Path(sys.argv[2])
want = {l.strip() for l in gzip.open(sets / "reads.txt.gz", "rt") if l.strip()}
got = {l.strip() for l in gzip.open(out / "reads.txt.gz", "rt") if l.strip()}
p = dict(l.rstrip("\n").split("\t", 1) for l in open(out / "params.tsv"))
m = re.search(r"Counts: ([\d,]+) pool reads seen in md5 order, ([\d,]+) dropped by the noise rule, ([\d,]+) kept: "
              r"([\d,]+) signal-estimated", (sets / "README.md").read_text())
readme = tuple(int(x.replace(",", "")) for x in m.groups())
mine = tuple(int(p[k]) for k in ("seen", "dropped_noise", "kept", "kept_bases"))
ok = got == want and mine == readme
print(f"{sets.name}: {len(got):,} reads, committed {len(want):,}, equal as a set: {got == want} "
      f"({len(got - want):,} only here, {len(want - got):,} only committed); seen / noise / kept / bases {mine} "
      f"vs README {readme}: {'EQUAL' if mine == readme else 'DIFFERENT'}")
sys.exit(0 if ok else 1)
PY
done
