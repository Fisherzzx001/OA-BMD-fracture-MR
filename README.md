# OA-BMD-fracture-MR

**Genetically higher bone mineral density raises knee osteoarthritis but protects against fracture: Mendelian randomisation of 13 skeletal phenotypes**

This repository holds the analysis code and reproduction material for the manuscript submitted to *Osteoarthritis and Cartilage* (`docs/08_投稿稿_OA_BMD_fracture.md`).

## Repository contents

- `code/` — all Python (3.13) and R (4.6.1) scripts that produce the tables and figures.
  - LDSC genetic correlation (`r_g` across 78 phenotype pairs)
  - Two-sample MR (IVW fixed/random, MR-Egger, weighted median, MR-PRESSO, leave-one-out)
  - Sample-overlap correction (MRlap v0.0.3.3)
  - Effect-scale reconstruction (placeholder N=400,000 → true ≈31,873; ×0.2823)
  - Bayesian colocalisation (coloc.abf; WNT16, SOX6, GDF5 + all 31 BMD instrument loci)
  - Multivariable MR conditioning on adult height (Yengo 2022) — grid n_h ∈ {0,100,300,1,000}, EUR and EUR_noUKB
- `code/output/` — intermediate outputs (clumping logs, MRlap outputs, the MVMR height-grid CSV).
- `docs/08_投稿稿_OA_BMD_fracture.md` — the OAC-compliant full-length manuscript (IMRaD body 2,359 words; abstract 207 words; 5 Highlights).
- `docs/Supplementary/` — supplementary tables S8–S11, the STROBE-MR checklist, and the SOS-theory predictions/falsifying tests.
- `Figures/` — high-resolution PNGs for figures S1–S8 and the Graphical Abstract (copied from the figure package at submission time).

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy pandas scipy matplotlib
# MRlap requires R 4.x with the mrlap package:
#   install.packages("remotes"); remotes::install_github("nerfaxl/mrlap")
python code/49_oac_wordcount.py docs/08_投稿稿_OA_BMD_fracture.md   # OAC compliance audit
```

Summary GWAS are publicly available (see Data availability in the manuscript) and were not redistributed.

## License

Code: MIT (see `LICENSE`). Manuscript text is governed by the publisher. Summary statistics are governed by their original consortia.

## Citation

Zheng Z, Chen S, Huang S, Qian C, Yuan H, Bao H, Jiang Q. Genetically higher bone mineral density raises knee osteoarthritis but protects against fracture: Mendelian randomisation of 13 skeletal phenotypes. *Osteoarthritis and Cartilage*. (submitted).
