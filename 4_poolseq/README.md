### Data availability

Raw data is available at NCBI SRA under the BioProject number `PRJNA849360`.

### Data processing

Steps below are computationally intensive and intermediate files are several GB in size, so in the folder `0_data` here we provide all the inputs to reproduce graphs and analyses.

**1.** [organize data](/4_poolseq/1_organize)

**2.** [mapping data](/4_poolseq/2_mapping)

**3.** [generate sync files](/4_poolseq/3_sync)

**4.** [generate vcf files](/4_poolseq/4_vcf)

**5.** [calculate Fst](/4_poolseq/5_fst)

**6.** [run CLEAR](/4_poolseq/6_clear)

**7.** [run GLMM](/4_poolseq/7_glmm)

**8.** [run CMH](/4_poolseq/8_cmh)

### Visualization and analyses

**9.** [run CMH](/4_poolseq/9_manhattan_plots)

**10.** [run CMH](/4_poolseq/10_snp_frequency)
