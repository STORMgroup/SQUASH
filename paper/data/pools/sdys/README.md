# pools/sdys: the S. dysgalactiae UT9728 pool (sets u9k60_both60, u9k100_both100)

Run FBA61885 (Schababerle et al. 2025; Zenodo 17428855, SQK-RBK114-24, every read) holds two strains, UT9728 and
UT10237, on 24 rapid barcodes. Which read is whose takes basecalls, so this pool is committed as a list.

- `ut9728.txt.gz`: the run's UT9728 reads with a signal-estimated length of at least 10 kb (samples x 400 / sample
  rate, to the nearest base), sorted by id, gzip -9 -n. 237,941 read ids (of the strain's 2,786,491; the shorter reads
  can never pass the selection's length rule, so they are left out to keep the list small). SHA-256 of the text
  f63defdd658d3da411a30551d599ab50018f9d4cf8ffb5577424931236d586e6.
- Made by `../../../curation/sdys_pool.sh` (a GPU and the depositors' two hybrid assemblies):
  1. dorado 1.4.0, dna_r10.4.1_e8.2_400bps_hac@v5.2.0, `basecaller --kit-name SQK-RBK114-24 --disable-read-splitting`
     on every POD5 of the run: each read's barcode;
  2. per barcode the first 3,000 reads of >= 5,000 basecalled bases, mapped with minimap2 2.30 `-x map-ont
     --secondary=no` to UT9728_IPa_Hyb and UT10237_JJa_Hyb (Zenodo 15558489); a barcode is a strain's when >= 50% of its
     sampled reads map and >= 80% of those to that strain. UT9728: barcodes 01-06 and 09-14 (12 libraries, 98.6-99.9%
     of the mapped sampled reads) and 18 and 20 (25 reads in all, 4 of them in this list, none in the sets); UT10237:
     07, 08, 15, 16 (98.5-99.9%), 17 and 24; 19, 21-23 and unclassified reads: neither;
  3. every read of a UT9728 barcode, then the reads of >= 10 kb (`select_reads.py --pool-out`).
- The sets: `select_reads.py` on the run with this list, the md5-order prefix of the pool to 100 x G (u9k60) or 160 x G
  (u9k100), G = 2,093,215 (`../../../curation/README.md`). Every read of both sets is in this list.
- Dorado writes reads in no fixed order, so a new run may sample other reads for step 2; the barcode calls and the
  strain of each barcode are what the list depends on.
