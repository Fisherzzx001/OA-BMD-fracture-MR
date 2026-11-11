<!-- Standalone supplementary table extracted verbatim from
     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.
     Source slice sha256: 603f08ff569543735f834acf931c26361ede6a22133ec59073c0d000c02c8136 . No statistic was rewritten. -->

## Supplementary Table S25. Bone-versus-cartilage tissue enrichment — competitive gene-set test (F3, this revision)

Tissue specificity of the 13-phenotype genetic architecture was tested by a competitive gene-set enrichment (MAGMA-style): each GWAS SNP was assigned to the nearest gene body within ±50 kb; gene-level association = the minimum p across assigned SNPs; each GO gene set was tested against the remaining genes by a one-sided Mann–Whitney U on −log10 p. Positive enrichment Δ means the set's genes are more associated than background. Bone-side sets: bone development, ossification, skeletal-system development. Cartilage-side sets: cartilage development, chondrocyte differentiation, extracellular-matrix organisation, skeletal-joint morphogenesis. Source GWAS: OA knee, total-body BMD, any-fracture (sumstats in OA_FLAGSHIP/07_reference/ldsc; code/171).

| Trait | Gene set | Tissue | n (set) | mean −log10p set | mean −log10p bg | enrichment Δ | MWW p (1-sided) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| OA_knee | bone_development | bone | 83 | 3.415 | 2.610 | 0.805 | 2.21e-04 |
| OA_knee | cartilage_development | cartilage | 72 | 5.036 | 2.606 | 2.430 | 1.45e-08 |
| OA_knee | ossification | bone | 71 | 3.741 | 2.609 | 1.131 | 2.78e-03 |
| OA_knee | chondrocyte_differentiation | cartilage | 81 | 4.932 | 2.605 | 2.326 | 1.58e-08 |
| OA_knee | bone_morphogenesis | other | 9 | 2.084 | 2.613 | -0.528 | 6.07e-01 |
| OA_knee | extracellular_matrix_organization | cartilage | 37 | 3.084 | 2.612 | 0.472 | 3.27e-02 |
| OA_knee | skeletal_system_development | bone | 42 | 5.354 | 2.608 | 2.746 | 9.96e-07 |
| OA_knee | skeletal_joint_morphogenesis | cartilage | 82 | 4.057 | 2.608 | 1.449 | 4.58e-10 |
| OA_knee | bone_remodeling | other | 1 | 4.853 | 2.612 | 2.241 | NA |
| BMD_total | bone_development | bone | 83 | 2.121 | 1.671 | 0.450 | 3.07e-03 |
| BMD_total | cartilage_development | cartilage | 72 | 2.198 | 1.671 | 0.526 | 1.35e-02 |
| BMD_total | ossification | bone | 71 | 2.285 | 1.671 | 0.614 | 1.00e-01 |
| BMD_total | chondrocyte_differentiation | cartilage | 81 | 2.277 | 1.671 | 0.606 | 2.11e-02 |
| BMD_total | bone_morphogenesis | other | 9 | 2.124 | 1.673 | 0.452 | 1.74e-02 |
| BMD_total | extracellular_matrix_organization | cartilage | 37 | 1.903 | 1.673 | 0.230 | 1.52e-02 |
| BMD_total | skeletal_system_development | bone | 42 | 3.403 | 1.670 | 1.733 | 1.26e-04 |
| BMD_total | skeletal_joint_morphogenesis | cartilage | 82 | 2.530 | 1.670 | 0.860 | 3.16e-07 |
| BMD_total | bone_remodeling | other | 1 | 3.645 | 1.673 | 1.972 | NA |
| fracture | bone_development | bone | 83 | 1.599 | 1.455 | 0.144 | 3.14e-02 |
| fracture | cartilage_development | cartilage | 72 | 1.895 | 1.454 | 0.441 | 2.16e-04 |
| fracture | ossification | bone | 71 | 1.538 | 1.455 | 0.083 | 3.52e-01 |
| fracture | chondrocyte_differentiation | cartilage | 81 | 1.812 | 1.454 | 0.357 | 8.38e-04 |
| fracture | bone_morphogenesis | other | 9 | 1.590 | 1.455 | 0.135 | 1.81e-01 |
| fracture | extracellular_matrix_organization | cartilage | 37 | 1.696 | 1.455 | 0.241 | 5.33e-04 |
| fracture | skeletal_system_development | bone | 42 | 1.825 | 1.455 | 0.370 | 5.42e-04 |
| fracture | skeletal_joint_morphogenesis | cartilage | 82 | 1.902 | 1.454 | 0.448 | 1.36e-06 |
| fracture | bone_remodeling | other | 1 | 1.445 | 1.455 | -0.010 | NA |
