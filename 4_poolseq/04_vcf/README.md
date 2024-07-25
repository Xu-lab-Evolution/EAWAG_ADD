### Generate VCF files

For this first part, create the following folders:

```bash
mkdir $PRJDIR/6_gatk
mkdir $PRJDIR/6_gatk/scripts
```
Inside the folder `6_gatk` create the file `chr.list` and paste the list below:

```bash
D_magna_CH1_L
D_magna_CH1_R
D_magna_CH2_L
D_magna_CH2_R
D_magna_CH3_L
D_magna_CH3_R
D_magna_CH4_L
D_magna_CH4_R
D_magna_CH5_L
D_magna_CH5_R
D_magna_CH6
D_magna_CH7_L
D_magna_CH7_R
D_magna_CH8_L
D_magna_CH8_R
D_magna_CH10_L
D_magna_CH10_R
D_magna_CH9_R
D_magna_CH9_L
```

Create gatk dictionary:

```bash
gatk CreateSequenceDictionary -R genome.fasta
```

Run gatk for each sample:

```bash
SAMPLE=############### ADD SAMPLE NAME, e.g. SRR19749256 ###############

NTHREADS=32
DATADIR=$PRJDIR/3_bwa/mapped
REFDIR=$PRJDIR/reference_genome
OUTDIR=$PRJDIR/6_gatk

gatk HaplotypeCaller \
-I $DATADIR/${SAMPLE}.bam \
-R $REFDIR/genome.fasta \
-O $OUTDIR/${SAMPLE}.vcf \
-ERC GVCF -G AS_StandardAnnotation -G StandardAnnotation
```

Combine VCFs and generate input for the GLMM:

```bash
NTHREADS=32
DATADIR=$PRJDIR/3_bwa/mapped
REFDIR=$PRJDIR/reference_genome
OUTDIR=$PRJDIR/6_gatk

cd $OUTDIR

for file in *.vcf; do SAMPLE=$(basename $file ".vcf"); echo $SAMPLE; echo $file; done | paste - - > mapping.txt

gatk GenomicsDBImport \
-R $REFDIR/genome.fasta \
-L $OUTDIR/chr.list \
--sample-name-map $OUTDIR/mapping.txt \
--genomicsdb-workspace-path $OUTDIR/daph_database

gatk GenotypeGVCFs \
-R $REFDIR/genome.fasta \
-V gendb://$OUTDIR/daph_database \
--output $OUTDIR/combined.gvcf.gz

gatk VariantFiltration \
-V $OUTDIR/combined.gvcf.gz \
-R $REFDIR/genome.fasta \
-O $OUTDIR/combinedFILT.gvcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60"

bcftools view -m2 -M2 -v snps -f 'PASS' -O z -o $OUTDIR/final.vcf.gz $OUTDIR/combinedFILT.gvcf.gz

vcftools --gzvcf final.vcf.gz \
--min-alleles 2 --max-alleles 2 --min-meanDP 20 --max-meanDP 400 --maf 0.01 --max-maf 0.99 \
--recode --recode-INFO-all --out final.dp20_400

gatk VariantsToTable -V final.dp20_400.recode.vcf \
-F CHROM -F POS -F ID -F REF -F ALT -GF AD -O final.dp20_400.AD.txt

sed -i 's/\.AD//g' final.dp20_400.AD.txt
```

The file `final.vcf.gz` is the finalized VCF file, while `final.dp20_400.AD.txt` is the input for the GLMM.
