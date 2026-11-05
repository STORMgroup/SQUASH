#!/usr/bin/env bash
# paper/curation/hs13q_pool.sh: the pool list of the human chr13 q-arm set (../data/pools/hs13q/hg002_chr13q.txt.gz),
# made with ONT's basecalls and a reference: the reads of HG002's maternal chr13 q-arm in ONT's two flow cells.
#   1. ONT's sup v5.0.0 BAMs of PAW70337 and PAW71238, aligned to GRCh38 (../data/download/manifests/
#      hg002_giab2025_bams.tsv, 328.7 GB; each deleted once its chr13 records are out, unless KEEP_BAMS=1): the primary
#      records on chr13 (`samtools view -F 0x904 -e 'rname=="chr13"'`) -> WORK/chr13.FC.bam;
#   2. candidates: those records' reads not split by the basecaller (no pi tag) with duration x 400 >= 9,000 (du tag)
#      -> WORK/cand_ids.txt;
#   3. the mapping reference: HG002 Q100 v1.2 (hg002v1.2.fasta.gz, the release references.sh also reads), one haplotype:
#      every *_MATERNAL chromosome, chrY_PATERNAL and chrM (a diploid reference gives MAPQ 0 where the haplotypes
#      agree), indexed with minimap2 `-x map-ont`;
#   4. the candidates' basecalls (samtools fastq of their records) mapped with minimap2 2.30 `-x map-ont
#      --secondary=no`; the cohort: reads whose primary record with the highest chaining score s1 has MAPQ >= 10 on
#      chr13_MATERNAL and starts at >= 17,240,000 (the q-arm after the centromere; pools.py cohort)
#      -> WORK/hg002_chr13q.txt.gz (gzip -9 -n), compared with the committed list.
# The selection then needs these reads' POD5s: ../data/download/hg002_giab2025.sh with this list (the whole run,
# 3.28 TB, streamed; ~100 GB kept).
# Settings (paper/config.sh, config.local.sh or the environment): SOURCES_DIR, OUT_DIR, THREADS; CURATION_DIR (default
# OUT_DIR/curation; WORK = CURATION_DIR/hs13q), MINIMAP2, SAMTOOLS (the paper's: minimap2 2.30-r1287, samtools 1.23.1),
# PYTHON, KEEP_BAMS.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../data/download/lib.sh"
CUR=$(cd "$DL_HERE/../../curation" && pwd)
W=${CURATION_DIR:-$OUT_DIR/curation}/hs13q
MINIMAP2=${MINIMAP2:-minimap2} SAMTOOLS=${SAMTOOLS:-samtools} THREADS=${THREADS:-$(nproc)}
mkdir -p "$W/bam" "$W/ref"

# 1. chr13 primary records of ONT's BAMs
while IFS=$'\t' read -r key bytes fc; do
  [[ -s $W/chr13.$fc.bam ]] && continue
  fetch "https://ont-open-data.s3.amazonaws.com/$key" "$W/bam/$fc.bam" "$bytes"
  log "chr13 primary records of $fc"
  "$SAMTOOLS" view -@ "$THREADS" -b -F 0x904 -e 'rname=="chr13"' -o "$W/.chr13.$fc.bam" "$W/bam/$fc.bam"
  mv "$W/.chr13.$fc.bam" "$W/chr13.$fc.bam"
  [[ ${KEEP_BAMS:-0} == 1 ]] || rm -f -- "$W/bam/$fc.bam"
done < <(tail -n +2 "$MANIFESTS/hg002_giab2025_bams.tsv")
mapfile -t CHR13 < <(tail -n +2 "$MANIFESTS/hg002_giab2025_bams.tsv" | cut -f3 | sed "s|.*|$W/chr13.&.bam|")

# 2. candidates
if [[ ! -s $W/cand_ids.txt ]]; then
  for b in "${CHR13[@]}"; do "$SAMTOOLS" view "$b"; done | awk -F'\t' '{pi = 0; du = -1
    for (i = 12; i <= NF; i++) {
      if (substr($i, 1, 5) == "pi:Z:") pi = 1
      if (substr($i, 1, 5) == "du:f:") du = substr($i, 6) + 0
    }
    if (!pi && du * 400 >= 9000) print $1}' | LC_ALL=C sort -u > "$W/cand_ids.txt.part"
  mv "$W/cand_ids.txt.part" "$W/cand_ids.txt"
fi
log "hs13q: $(wc -l < "$W/cand_ids.txt") candidates"

# 3. the one-haplotype reference
H=$SOURCES_DIR/references/.dl/hg002v1.2.fasta
if [[ ! -s $H ]]; then
  mkdir -p "$(dirname "$H")"
  $CURL -o "$H.gz" https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.2.fasta.gz
  gzip -dc "$H.gz" > "$H.part"; mv "$H.part" "$H"; rm -f -- "$H.gz"
fi
if [[ ! -s $W/ref/hg002_hap.mmi ]]; then
  [[ -s $H.fai ]] || "$SAMTOOLS" faidx "$H"
  mapfile -t ctg < <(cut -f1 "$H.fai" | grep -E '_MATERNAL$|^chrY_PATERNAL$|^chrM$')
  "$SAMTOOLS" faidx "$H" "${ctg[@]}" > "$W/ref/hg002_hap.fa"
  "$MINIMAP2" -x map-ont -t "$THREADS" -d "$W/ref/.hg002_hap.mmi" "$W/ref/hg002_hap.fa" 2> "$W/ref/index.log"
  mv "$W/ref/.hg002_hap.mmi" "$W/ref/hg002_hap.mmi"
fi

# 4. the cohort
if [[ ! -s $W/cand.paf ]]; then
  for b in "${CHR13[@]}"; do "$SAMTOOLS" view -b -N "$W/cand_ids.txt" "$b" | "$SAMTOOLS" fastq - 2>/dev/null; done \
    > "$W/cand.fq"
  "$MINIMAP2" -x map-ont --secondary=no -t "$THREADS" "$W/ref/hg002_hap.mmi" "$W/cand.fq" 2> "$W/cand.minimap2.log" \
    > "$W/cand.paf.part"
  mv "$W/cand.paf.part" "$W/cand.paf"
fi
"$PYTHON" "$CUR/pools.py" cohort "$W/cand.paf" chr13_MATERNAL 17240000 10 "$W/hg002_chr13q.txt"
gzip -9 -n -f "$W/hg002_chr13q.txt"
want=$DL_HERE/../pools/hs13q/hg002_chr13q.txt.gz
if cmp -s <(gzip -dc "$W/hg002_chr13q.txt.gz") <(gzip -dc "$want"); then
  log "hs13q: $W/hg002_chr13q.txt.gz equals the committed pool list"
else
  n=$(comm -3 <(gzip -dc "$W/hg002_chr13q.txt.gz") <(gzip -dc "$want") | wc -l)
  log "hs13q: $W/hg002_chr13q.txt.gz differs from $want: $n reads in one only"
fi
