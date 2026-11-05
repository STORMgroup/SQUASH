#!/usr/bin/env bash
# paper/curation/sdys_pool.sh: the pool list of the S. dysgalactiae UT9728 sets (../data/pools/sdys/ut9728.txt.gz), made
# with basecalls: run FBA61885 (SQK-RBK114-24) holds two strains on 24 barcodes, and only basecalls tell the barcodes.
#   1. the run (../data/download/sdys_zenodo17428855.sh, 199.9 GB);
#   2. barcodes: dorado 1.4.0, dna_r10.4.1_e8.2_400bps_hac@v5.2.0, `basecaller --kit-name SQK-RBK114-24
#      --disable-read-splitting` (a GPU; PARTS=N runs it on every N-th file in N parts, one after the other; the
#      barcode calls do not depend on the split) -> WORK/demux/part_P.bam; read id, BC tag ('unclassified' when none)
#      and basecalled length -> WORK/read_bc.tsv;
#   3. barcode -> strain: per barcode the first 3,000 reads of >= 5,000 basecalled bases (pools.py sample) mapped with
#      minimap2 2.30 `-x map-ont --secondary=no` to the depositors' two hybrid assemblies (Zenodo 15558489:
#      UT9728_IPa_Hyb.fasta, UT10237_JJa_Hyb.fasta; their one contig each renamed UT9728 / UT10237); a barcode is a
#      strain's when >= 50% of its sampled reads map and >= 80% of those to that strain (pools.py barcodes)
#      -> WORK/barcodes.tsv (the paper's: barcodes 01-06, 09-14, 18 and 20 = UT9728; 07, 08, 15, 16, 17, 24 = UT10237);
#   4. the UT9728 reads (pools.py strain), then those of >= 10 kb signal-estimated length (select_reads.py --pool-out: rules
#      1-2 of the selection, which can drop no other read) -> WORK/ut9728.txt.gz (gzip -9 -n), compared with the
#      committed list.
# Settings (paper/config.sh, config.local.sh or the environment): SOURCES_DIR, OUT_DIR, PYTHON (pod5 0.3.39), THREADS;
# CURATION_DIR (default OUT_DIR/curation; WORK = CURATION_DIR/sdys), DORADO (the dorado 1.4.0 binary), DORADO_MODEL
# (the model folder or its name), MINIMAP2, SAMTOOLS (the paper's: minimap2 2.30-r1287, samtools 1.23.1), PARTS (1).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../data/download/lib.sh"
CUR=$(cd "$DL_HERE/../../curation" && pwd)
W=${CURATION_DIR:-$OUT_DIR/curation}/sdys
DORADO=${DORADO:-dorado} DORADO_MODEL=${DORADO_MODEL:-dna_r10.4.1_e8.2_400bps_hac@v5.2.0}
MINIMAP2=${MINIMAP2:-minimap2} SAMTOOLS=${SAMTOOLS:-samtools} PARTS=${PARTS:-1} THREADS=${THREADS:-$(nproc)}
mkdir -p "$W/demux" "$W/ref"

# 1. the run
bash "$DL_HERE/sdys_zenodo17428855.sh"
SRC=$SOURCES_DIR/sdys_zenodo17428855
mapfile -t FILES < <(find -L "$SRC" -name '*.pod5' | sort)
(( ${#FILES[@]} )) || die "no POD5 in $SRC"

# 2. barcodes
for ((p = 0; p < PARTS; p++)); do
  out=$W/demux/part_$p.bam
  [[ -s $out ]] && continue
  d=$W/demux/files_$p; rm -rf -- "$d"; mkdir -p "$d"
  for ((i = p; i < ${#FILES[@]}; i += PARTS)); do ln -s "${FILES[i]}" "$d/$(basename "${FILES[i]}")"; done
  log "dorado part $p of $PARTS: $(ls "$d" | wc -l) files"
  "$DORADO" basecaller --kit-name SQK-RBK114-24 --disable-read-splitting "$DORADO_MODEL" "$d" > "$out.part"
  mv "$out.part" "$out"; rm -rf -- "$d"
done
if [[ ! -s $W/read_bc.tsv ]]; then
  for ((p = 0; p < PARTS; p++)); do "$SAMTOOLS" view "$W/demux/part_$p.bam"; done | awk -F'\t' '{bc = "unclassified"
    for (i = 12; i <= NF; i++) if ($i ~ /^BC:Z:/) bc = substr($i, 6)
    print $1 "\t" bc "\t" length($10)}' > "$W/read_bc.tsv.part"
  mv "$W/read_bc.tsv.part" "$W/read_bc.tsv"
fi

# 3. barcode -> strain
Z=https://zenodo.org/api/records/15558489/files
fetch "$Z/UT9728_IPa_Hyb.fasta/content" "$W/ref/UT9728_IPa_Hyb.fasta" 2093271 1680c782c933f4b9227eb610229e918a
fetch "$Z/UT10237_JJa_Hyb.fasta/content" "$W/ref/UT10237_JJa_Hyb.fasta" 2140157 00dc82d814d4a937581d0fdabd591afa
{ awk 'NR == 1 {print ">UT9728"; next} 1' "$W/ref/UT9728_IPa_Hyb.fasta"
  awk 'NR == 1 {print ">UT10237"; next} 1' "$W/ref/UT10237_JJa_Hyb.fasta"; } > "$W/ref/both.fa"
if [[ ! -s $W/sample.paf ]]; then
  "$PYTHON" "$CUR/pools.py" sample "$W/read_bc.tsv" 3000 5000 "$W/sample_ids.txt"
  for ((p = 0; p < PARTS; p++)); do
    "$SAMTOOLS" view -b -N "$W/sample_ids.txt" "$W/demux/part_$p.bam" | "$SAMTOOLS" fastq - 2>/dev/null
  done > "$W/sample.fq"
  "$MINIMAP2" -x map-ont --secondary=no -t "$THREADS" "$W/ref/both.fa" "$W/sample.fq" 2> "$W/sample.minimap2.log" \
    > "$W/sample.paf.part"
  mv "$W/sample.paf.part" "$W/sample.paf"
fi
"$PYTHON" "$CUR/pools.py" barcodes "$W/read_bc.tsv" "$W/sample.paf" 3000 5000 "$W/barcodes.tsv"

# 4. the pool list
"$PYTHON" "$CUR/pools.py" strain "$W/read_bc.tsv" "$W/barcodes.tsv" UT9728 "$W/ut9728_all.txt"
"$PYTHON" "$CUR/select_reads.py" --pod5 "$SRC" --pool-ids "$W/ut9728_all.txt" --min-length 10000 --bases-per-second 400 \
  --pool-out "$W/ut9728.txt" --threads "$THREADS"
gzip -9 -n -f "$W/ut9728.txt"
want=$DL_HERE/../pools/sdys/ut9728.txt.gz
if cmp -s <(gzip -dc "$W/ut9728.txt.gz") <(gzip -dc "$want"); then
  log "sdys: $W/ut9728.txt.gz equals the committed pool list"
else
  n=$(comm -3 <(gzip -dc "$W/ut9728.txt.gz") <(gzip -dc "$want") | wc -l)
  log "sdys: $W/ut9728.txt.gz differs from $want: $n reads in one only"
fi
