# Bioinformatics 1 — Assignment 1 (2026 Spring)

서울대 생물정보학 및 실습 1 — Term Project Week 1.

## Environment
- WSL2 Ubuntu + Miniconda
- conda env `binfo1`: subread, samtools, pandas, numpy, matplotlib, scipy, jupyterlab

## Data (not in repo)
- `binfo1-datapack1.tar` from https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar
- `gencode.vM27.annotation.gtf.gz` from EBI

## Notebooks
- `week1_setup.ipynb` — environment setup, featureCounts, Figure 4D/5B reproduction

## Results (Week 1)
- Cell 33: Figure 4D-like scatter, n=9145, Pearson r=+0.45
- Cell 37: Figure 5B-like, colored by localization
  - integral membrane: median log2(rden_change) = +0.401 (n=1211)
  - nucleus: -0.458 (n=2682), cytoplasm: -0.354 (n=1186)
- Integral membrane proteins show strongest translation upregulation upon LIN28A Knockdown.
