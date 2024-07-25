### Data availability

Raw data is available at NCBI SRA under the BioProject number `PRJNA849360`.

### Data processing

Steps below are computationally intensive and intermediate files are several GB in size, so in the folder `00_data` here we provide all the inputs to reproduce graphs and analyses.

**1.** [organize data](/4_poolseq/01_organize)

**2.** [mapping data](/4_poolseq/02_mapping)

**3.** [generate sync files](/4_poolseq/03_sync)

**4.** [generate vcf files](/4_poolseq/04_vcf)

**5.** [calculate Fst](/4_poolseq/05_fst)

**6.** [run CLEAR](/4_poolseq/06_clear)

**7.** [run GLMM](/4_poolseq/07_glmm)

**8.** [run CMH](/4_poolseq/08_cmh)

### Visualization and analyses

**9.** [Manhattan plots](/4_poolseq/09_manhattan_plots)

**10.** [SNP frequency plots](/4_poolseq/10_snp_frequency)
