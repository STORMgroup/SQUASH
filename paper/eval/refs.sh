#!/usr/bin/env bash
# paper/eval/refs.sh [KEY ...]: the evaluation references of the sets (default: all of eval/sets.tsv) -> EVAL_REFS
# (default SOURCES_DIR/references/eval/), each made from the public references (data/download/references.sh, run here
# for the ones missing) and checked against the paper's file (eval/references.tsv: every record's name, length and
# sequence SHA-256; the BED and UCSC files by their bytes):
#   ecoli_chr.fa               the E. coli chromosome CP149500.1 (the first record of ecoli_CP149500.fa)
#   u9.fa                      UT9728_IPa_Hyb.fasta with its record named u9_UT9728
#   cj.fa vp.fa me14.fa chr13q_MAT.fa chr13q_hap.fa dmel_r6_nuclear.fa colcc.fa   references.sh's files as they are
#   organelles.fa              A. thaliana chloroplast NC_000932.1 and mitochondrion NC_037304.1 (NCBI efetch)
#   colcc_org.fa               Col-CC's five chromosomes then organelles.fa (D8's reference)
#   colcc_org.masked.fa        colcc_org.fa with the chromosomes' tandem arrays as N (mask_arrays.py; colcc_arrays.bed)
#   dmel_r6_nuclear.masked.fa  release 6 with the rDNA locus, chrY and the pericentromeric heterochromatin as N
#                              (mask_dmel.py on UCSC's dm6 cytoBandIdeo and rmsk tables; dmel_mask.bed)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
R=$SOURCES_DIR/references D=$EVAL_REFS
mkdir -p "$D"
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
COLCC_CHR=CP116280.1,CP116281.2,CP116282.1,CP116283.2,CP116284.1

public() {  # public FILE REF: references.sh's FILE (REF = its row of data/sources.tsv), downloaded if missing
  [[ -s $R/$1 ]] || bash "$EVAL/../data/download/references.sh" "$2"
  need "$R/$1"
}
check() {  # check FILE: the paper's (eval/references.tsv: a FASTA's records in order, each its name, length and the
  # SHA-256 of its upper-case sequence, so line width and header text do not matter; another file by its bytes)
  "$REMORA_PYTHON" - "$EVAL/references.tsv" "$D/$1" "$1" <<'PY' || die "$D/$1 is not the paper's"
import hashlib, sys
man, path, name = sys.argv[1:4]
want = [tuple(l.rstrip("\n").split("\t")[1:]) for l in open(man) if l.split("\t")[0] == name]
if want and want[0][0] == "-":
    b = open(path, "rb").read()
    got = [("-", str(len(b)), hashlib.sha256(b).hexdigest())]
else:
    got, cur, seq = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                got.append((cur, "".join(seq)))
            cur, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if cur is not None:
        got.append((cur, "".join(seq)))
    got = [(r, str(len(s)), hashlib.sha256(s.encode()).hexdigest()) for r, s in got]
if got != want:
    sys.exit(f"{name}: {len(got)} records / content differ from the paper's ({len(want)} records)")
print(f"refs.sh: {name} equals the paper's " + ("(bytes)" if got[0][0] == "-" else f"({len(got)} records)"), file=sys.stderr)
PY
}
declare -A CHECKED
put() {  # put FILE: FILE made by make_FILE if missing, then checked (once)
  [[ -n ${CHECKED[$1]:-} ]] && return 0
  [[ -s $D/$1 ]] || { log "$1"; "make_${1//[.-]/_}"; }
  check "$1"; CHECKED[$1]=1
}
link() { public "$1" "$2"; ln -sfn "$R/$1" "$D/$1"; }

make_ecoli_chr_fa() { public ecoli_CP149500.fa ref_ec; awk '/^>/ {n++} n == 1' "$R/ecoli_CP149500.fa" > "$D/.ecoli_chr.fa"; mv "$D/.ecoli_chr.fa" "$D/ecoli_chr.fa"; }
make_u9_fa() { public UT9728_IPa_Hyb.fasta ref_sdys; sed '1s/^>.*/>u9_UT9728/' "$R/UT9728_IPa_Hyb.fasta" > "$D/.u9.fa"; mv "$D/.u9.fa" "$D/u9.fa"; }
make_cj_fa() { link cj.fa ref_cj; }
make_vp_fa() { link vp.fa ref_vp; }
make_me14_fa() { link me14.fa ref_arabiei; }
make_chr13q_MAT_fa() { link chr13q_MAT.fa ref_hg002; }
make_chr13q_hap_fa() { link chr13q_hap.fa ref_hg002; }
make_dmel_r6_nuclear_fa() { link dmel_r6_nuclear.fa ref_dmel; }
make_colcc_fa() { link colcc.fa ref_athal; }
make_organelles_fa() { $CURL -o "$D/.organelles.fa" "${EFETCH}NC_000932.1,NC_037304.1"; mv "$D/.organelles.fa" "$D/organelles.fa"; }
make_colcc_org_fa() { put colcc.fa; put organelles.fa; cat "$D/colcc.fa" "$D/organelles.fa" > "$D/.colcc_org.fa"; mv "$D/.colcc_org.fa" "$D/colcc_org.fa"; }
make_colcc_org_masked_fa() {
  put colcc_org.fa
  "$REMORA_PYTHON" "$EVAL/mask_arrays.py" "$D/colcc_org.fa" "$D/colcc_org.masked.fa" "$D/colcc_arrays.bed" --mask "$COLCC_CHR"
  check colcc_arrays.bed
}
make_cytoBandIdeo_txt_gz() { $CURL -o "$D/.cytoBandIdeo.txt.gz" "$UCSC/cytoBandIdeo.txt.gz"; mv "$D/.cytoBandIdeo.txt.gz" "$D/cytoBandIdeo.txt.gz"; }
make_rmsk_txt_gz() { $CURL -o "$D/.rmsk.txt.gz" "$UCSC/rmsk.txt.gz"; mv "$D/.rmsk.txt.gz" "$D/rmsk.txt.gz"; }
make_dmel_r6_nuclear_masked_fa() {
  put dmel_r6_nuclear.fa; put cytoBandIdeo.txt.gz; put rmsk.txt.gz
  "$REMORA_PYTHON" "$EVAL/mask_dmel.py" "$D/dmel_r6_nuclear.fa" "$D/cytoBandIdeo.txt.gz" "$D/rmsk.txt.gz" \
    "$D/dmel_r6_nuclear.masked.fa" "$D/dmel_mask.bed"
  check dmel_mask.bed
}

keys=("$@"); (( ${#keys[@]} )) || mapfile -t keys < <(tail -n +2 "$EVAL/sets.tsv" | cut -f1)
for k in "${keys[@]}"; do
  set_paths "$k"
  for f in "$REF" "$PLACE_REF" "$GENOME_REF" ${HAP/#-/} ${MASKED/#-/}; do put "$(basename "$f")"; done
done
