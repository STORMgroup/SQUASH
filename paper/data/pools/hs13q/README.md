# pools/hs13q: the human HG002 chromosome 13 q-arm pool (set h2k60_both60)

ONT open data giab_2025.01, HG002 flow cells PAW70337 and PAW71238 (R10.4.1, 5 kHz, pass and fail reads). Which reads
come from the maternal chr13 q-arm takes basecalls and a reference, so this pool is committed as a list.

- `hg002_chr13q.txt.gz`: the cohort's read ids, sorted, gzip -9 -n. 323,862 read ids; SHA-256 of the text
  e0e8c2daa096f04f839337300e766c876bcb7f3f2f8bde5b2ce4aa7cb26648fe. Every one of them is in the two flow cells' POD5s.
- Made by `../../../curation/hs13q_pool.sh` (ONT's BAMs, 328.7 GB, and the HG002 Q100 v1.2 assembly):
  1. ONT's sup v5.0.0 basecalls aligned to GRCh38 (`basecalling/sup/HG002/FC/calls.sorted.bam`): the primary records on
     chr13 (`samtools view -F 0x904 -e 'rname=="chr13"'`, samtools 1.23.1);
  2. candidates: those reads not split by the basecaller (no `pi` tag) with duration x 400 >= 9,000 (`du` tag);
  3. their basecalls mapped with minimap2 2.30 `-x map-ont --secondary=no` to one haplotype of HG002 Q100 v1.2 (every
     `*_MATERNAL` chromosome, chrY_PATERNAL and chrM: a diploid reference gives MAPQ 0 wherever the haplotypes agree,
     and paternal chr13 reads map to chr13_MATERNAL);
  4. the cohort: reads whose primary record with the highest chaining score (s1) has MAPQ >= 10 on chr13_MATERNAL and
     starts at >= 17,240,000 (the q-arm after the centromeric satellites): 323,862 of the 349,349 candidates (349,344
     mapped).
- The set: `select_reads.py` on these reads' POD5s with lengths from durations at 6 significant digits, every >= 10-kb read
  in the pool (314,536), 60 x G, G = 96,096,965 = chr13_MATERNAL from 17,240,000 (`../../../curation/README.md`).
  The POD5s: `../../download/hg002_giab2025.sh hg002_chr13q.txt.gz` streams the whole run (3.28 TB) and keeps these
  reads (~100 GB). Every read of h2k60_both60 is in this list.
