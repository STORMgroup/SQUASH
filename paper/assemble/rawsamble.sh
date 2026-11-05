#!/usr/bin/env bash
# paper/assemble/rawsamble.sh KEY: Rawsamble + miniasm on a paper set's raw signal as the paper ran it (assemble/lib.sh's
# protocol): the set's POD5 staged as STAGE/in.pod5, then three timed steps,
#   rawsamble_index    rawhash2 -x r10ava -t ASM_THREADS -p RAWHASH2_PORE -d LOCAL/ava.ind STAGE/in.pod5
#   rawsamble_overlap  rawhash2 -x r10ava -t ASM_THREADS LOCAL/ava.ind STAGE/in.pod5 > STAGE/ava.paf
#   miniasm            miniasm STAGE/ava.paf > LOCAL/ava.gfa
# (-x r10ava = Rawsamble's `-x ava --r10`; RAWHASH2_PORE = the R10.4.1 9-mer pore model, tools/env.sh). The paper ran it
# on the bacterial sets, A. rabiei and human, not on D. melanogaster and A. thaliana (assemble/options.tsv), where it would
# take days. Input SET_POD5 (default POD5_DIR/KEY.pod5). Output OUT_DIR/assembly/KEY/rawsamble/: asm/ava.gfa (the
# assembly graph), asm/ava.paf.gz (the overlaps); logs/rawsamble_index.time.txt, rawsamble_overlap.time.txt,
# miniasm.time.txt, their .out / .err, staging.time.txt; stats.tsv (overlaps, segments, links, SHA-256 of the overlaps);
# done. STAGE_DIR=/dev/shm puts the POD5 and the overlaps in RAM (the paper's large sets).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
as_pin "$0" "$@"
as_start "${1:?usage: rawsamble.sh KEY}" rawsamble
[[ $(as_col "$KEY" paper_tools) == *rawsamble* ]] || echo "rawsamble.sh: note: the paper did not run Rawsamble on $KEY" >&2
POD5=${SET_POD5:-$POD5_DIR/$KEY.pod5}
if [[ $DRY_RUN != 1 ]]; then
  [[ -s $POD5 ]] || { echo "rawsamble.sh: no $POD5 (make -C paper data SET=$KEY)" >&2; exit 1; }
  /usr/bin/time -v -o "$O/logs/staging.time.txt" cp -L "$POD5" "$IN/in.pod5"
fi
timed rawsamble_index "$RAWHASH2" -x r10ava -t "$ASM_THREADS" -p "$RAWHASH2_PORE" -d "$LOCAL/ava.ind" "$IN/in.pod5" \
  > "${O:?}/logs/rawsamble_index.out" 2> "$O/logs/rawsamble_index.err"
timed rawsamble_overlap "$RAWHASH2" -x r10ava -t "$ASM_THREADS" "$LOCAL/ava.ind" "$IN/in.pod5" > "$IN/ava.paf" \
  2> "$O/logs/rawsamble_overlap.err"
timed miniasm "$MINIASM" "$IN/ava.paf" > "$LOCAL/ava.gfa" 2> "$O/logs/miniasm.err"
[[ $DRY_RUN == 1 ]] && exit 0
printf 'overlaps\t%s\nsegments\t%s\nlinks\t%s\n' "$(wc -l < "$IN/ava.paf")" "$(grep -c '^S' "$LOCAL/ava.gfa" || true)" \
  "$(grep -c '^L' "$LOCAL/ava.gfa" || true)" > "$O/stats.tsv"
printf 'sha256_overlaps\t%s\n' "$(sha256sum < "$IN/ava.paf" | cut -c1-64)" >> "$O/stats.tsv"
if command -v pigz > /dev/null; then pigz -c -p "$ASM_THREADS" "$IN/ava.paf" > "$O/asm/ava.paf.gz"; else gzip -c "$IN/ava.paf" > "$O/asm/ava.paf.gz"; fi
cp "$LOCAL/ava.gfa" "$O/asm/ava.gfa"
as_done
