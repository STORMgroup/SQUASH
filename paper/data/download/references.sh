#!/usr/bin/env bash
# The evaluation references (QUAST; never read by SQUASH): references.sh [REF ...] (default: all), REF = a `ref_*` row
# of ../sources.tsv. Each is fetched from its public record (NCBI nuccore efetch, NCBI Datasets, Zenodo, the HG002 Q100
# release), reduced to the records the paper scored (renamed / sliced as noted) and checked record by record against
# manifests/references.tsv (length and SHA-256 of the upper-case sequence, so line width and header text do not
# matter). -> SOURCES_DIR/references/FILE (the file names the evaluation uses).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
R=$SOURCES_DIR/references
mkdir -p "$R/.dl"

efetch() {  # efetch OUT ACC[,ACC...]
  [[ -s $1 ]] && return 0
  $CURL -o "$1.part" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$2&rettype=fasta&retmode=text"
  mv "$1.part" "$1"
}
datasets() {  # datasets ACC -> $R/.dl/ACC.fna (the assembly's genomic FASTA from NCBI Datasets)
  local acc=$1
  [[ -s $R/.dl/$acc.fna ]] && return 0
  $CURL -o "$R/.dl/$acc.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
  (cd "$R/.dl" && rm -rf "$acc.unzip" && mkdir "$acc.unzip" && cd "$acc.unzip" && unzip -q "../$acc.zip")
  find "$R/.dl/$acc.unzip" -name '*.fna' -exec cat {} + > "$R/.dl/$acc.fna.part"
  mv "$R/.dl/$acc.fna.part" "$R/.dl/$acc.fna"; rm -rf -- "$R/.dl/$acc.unzip" "$R/.dl/$acc.zip"
}
pick() {  # pick IN OUT SPEC...: SPEC = NAME[=NEWNAME][@START] (0-based start of a suffix slice); records in SPEC order
  "$PYTHON" - "$@" <<'PY'
import sys
src, out, specs = sys.argv[1], sys.argv[2], sys.argv[3:]
want = {}
for s in specs:
    name, _, start = s.partition("@")
    old, _, new = name.partition("=")
    want[old] = (new or old, int(start or 0))
seqs, cur = {}, None
for line in open(src):
    if line.startswith(">"):
        cur = line[1:].split()[0]
        cur = cur if cur in want else None
        if cur:
            seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
missing = [w for w in want if w not in seqs]
if missing:
    sys.exit(f"{src}: no record {missing}")
with open(out + ".part", "w") as f:
    for old, (new, start) in want.items():
        s = "".join(seqs[old])[start:]
        f.write(f">{new}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")
import os
os.replace(out + ".part", out)
PY
}
check() {  # check REF FILE: every record of FILE as manifests/references.tsv lists it, and no other
  "$PYTHON" - "$MANIFESTS/references.tsv" "$1" "$R/$2" "$2" <<'PY'
import hashlib, sys
man, ref, path, name = sys.argv[1:5]
want = [l.rstrip("\n").split("\t") for l in open(man)][1:]
want = [(w[2], int(w[3]), w[4]) for w in want if w[0] == ref and w[1] == name]
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
got = [(r, len(s), hashlib.sha256(s.encode()).hexdigest()) for r, s in got]
if got != want:
    for g, w in zip(got, want):
        if g != w:
            print(f"{name}: record {g[0]} ({g[1]} bp) differs from the paper's {w[0]} ({w[1]} bp)", file=sys.stderr)
    sys.exit(f"{name}: {len(got)} records, the paper's file has {len(want)}; not the paper's reference")
print(f"{name}: {len(got)} records equal the paper's", file=sys.stderr)
PY
}

ref_ec()  { efetch "$R/ecoli_CP149500.fa" CP149500.1,CP149501.1,CP149502.1,CP149503.1,CP149504.1; check ref_ec ecoli_CP149500.fa; }
ref_cj()  { efetch "$R/cj.fa" CP149499.1; check ref_cj cj.fa; }
ref_vp()  { efetch "$R/vp.fa" CP149505.1,CP149506.1; check ref_vp vp.fa; }
ref_sdys() {
  fetch https://zenodo.org/api/records/15558489/files/UT9728_IPa_Hyb.fasta/content "$R/UT9728_IPa_Hyb.fasta" 2093271 1680c782c933f4b9227eb610229e918a
  check ref_sdys UT9728_IPa_Hyb.fasta
}
ref_athal() { datasets GCA_028009825.2
  pick "$R/.dl/GCA_028009825.2.fna" "$R/colcc.fa" CP116280.1 CP116281.2 CP116282.1 CP116283.2 CP116284.1; check ref_athal colcc.fa; }
ref_dmel() { datasets GCF_000001215.4
  pick "$R/.dl/GCF_000001215.4.fna" "$R/dmel_r6_nuclear.fa" NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NC_004353.4 NC_024512.1
  check ref_dmel dmel_r6_nuclear.fa; }
ref_arabiei() { datasets GCA_004011695.2
  pick "$R/.dl/GCA_004011695.2.fna" "$R/me14.fa" $(seq -f 'CP0%g.1' 95288 95310); check ref_arabiei me14.fa; }
ref_hg002() {  # HG002 Q100 v1.2; the q-arms start after the centromeric alpha satellite
  [[ -s $R/.dl/hg002v1.2.fasta ]] || { $CURL -o "$R/.dl/hg002v1.2.fasta.gz" https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.2.fasta.gz
    gzip -dc "$R/.dl/hg002v1.2.fasta.gz" > "$R/.dl/hg002v1.2.fasta.part"; mv "$R/.dl/hg002v1.2.fasta.part" "$R/.dl/hg002v1.2.fasta"; rm -f -- "$R/.dl/hg002v1.2.fasta.gz"; }
  pick "$R/.dl/hg002v1.2.fasta" "$R/chr13q_MAT.fa" chr13_MATERNAL=chr13q_MATERNAL@17240000
  pick "$R/.dl/hg002v1.2.fasta" "$R/chr13q_hap.fa" chr13_MATERNAL=chr13q_MATERNAL@17240000 chr13_PATERNAL=chr13q_PATERNAL@13120000
  check ref_hg002 chr13q_MAT.fa; check ref_hg002 chr13q_hap.fa; }

refs=("$@"); (( ${#refs[@]} )) || refs=(ref_ec ref_cj ref_vp ref_sdys ref_athal ref_dmel ref_arabiei ref_hg002)
for r in "${refs[@]}"; do
  declare -F "$r" > /dev/null || die "no reference $r"
  log "$r"; "$r"
done
