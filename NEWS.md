# SQUASH release notes

## 1.0.0 (2026-11)

The assembler of the paper "SQUASH: Assembling raw nanopore signals with minimizer-space de Bruijn graphs": the
pipeline of the paper's runs, step for step (every output identical to the paper's runs; `paper/PARITY.md`), in one
binary with two presets (`-x small`, `-x large`) and every value an option.
`paper/` reproduces the paper from public data: its 14 read sets, SQUASH's runs, the other tools' runs (basecalling,
nanoMDBG, Flye, hifiasm, Rawsamble + miniasm), the evaluation and the tables.
