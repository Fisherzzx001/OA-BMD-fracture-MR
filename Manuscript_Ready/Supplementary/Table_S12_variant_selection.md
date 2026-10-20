# Supplementary Table S12. Exposure (BMD-site) selection and instrument counts

Genome-wide-significant (GWS, p<5×10⁻⁸) locus counts per BMD site were counted from the full candidate-SNP clump-input files and validated against `strict_clump_summary.csv` (514/246/170/277/301 match exactly). `n_strict_instruments` is the count after strict PLINK-equivalent clumping (r²<0.001/10 Mb, 1000G EUR); 'NA' = strict clumping not completed at submission.

The cranial head site (UKB DXA field 23226) was excluded a priori, before any clumping: it is not a load-bearing skeletal region, and head BMD is genetically correlated with intracranial aneurysm in the source GWAS (Qian et al. 2025), a non-skeletal pleiotropy risk. Because the head site was never clumped, its signal is reported as genome-wide-significant SNPs from the deposited top-hits file (20 at p<5×10⁻⁸), documenting that the exclusion is principled rather than instrument-driven.

| BMD_site | n_gws_loci_p5e-8 | n_strict_instruments | entered_MR | reason |
| --- | --- | --- | --- | --- |
| Total body | 514 | 35 | Yes | Strict-clumped; >=15 independent instruments; entered MR |
| Femoral neck | 170 | 28 | Yes | Strict-clumped; >=15 independent instruments; entered MR |
| Femoral total | 262 | NA | No | Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set |
| Lumbar spine | 301 | 23 | Yes | Strict-clumped; >=15 independent instruments; entered MR |
| Spine (total) | 291 | NA | No | Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set |
| Arms | 246 | 20 | Yes | Strict-clumped; >=15 independent instruments; entered MR |
| Legs | 277 | 26 | Yes | Strict-clumped; >=15 independent instruments; entered MR |
| Pelvis | 345 | NA | No | Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set |
| Rib | 229 | NA | No | Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set |
| Trunk | 344 | NA | No | Strong GWS signal present but strict PLINK-equivalent clumping not completed at submission (1000G EUR panel unavailable in compute environment); deferred to avoid non-prespecified ad hoc instrument set |
| Head (excluded a priori) | 20 | NA | No | Excluded a priori: not a load-bearing skeletal region; head BMD is genetically correlated with intracranial aneurysm in the source GWAS (Qian et al. 2025), a non-skeletal pleiotropy risk. Genome-wide signal was present (20 SNPs at p<5e-8), so the exclusion is principled rather than instrument-driven. |
