### Map reads towards the reference genome

For this first part, create the following folders:

```bash
mkdir $PRJDIR/3_bwa
mkdir $PRJDIR/reference_genome
mkdir $PRJDIR/3_bwa/mapped
mkdir $PRJDIR/3_bwa/host_removed
```

Place the reference genome `fasta` file inside the folder `$PRJDIR/reference_genome`. This file is available inside the folder `data`.

Inside the folder `3_bwa` create the file `SRR_Acc_List.txt` and paste the list below:

```bash
SRR19749256
SRR19749257
SRR19749258
SRR19749259
SRR19749260
SRR19749261
SRR19749262
SRR19749263
SRR19749264
SRR19749265
SRR19749266
SRR19749268
SRR19749269
SRR19749270
SRR19749271
SRR19749272
SRR19749273
SRR19749274
SRR19749275
SRR25742225
SRR25742226
SRR25742227
SRR25742228
SRR25742229
SRR25742230
SRR25742232
SRR25742233
SRR25742234
SRR25742235
SRR25742238
SRR25742239
SRR25742240
```

Run `bwa`:

```bash
DATADIR=$PRJDIR/2_trimgalore
REFDIR=$PRJDIR/reference_genome
OUTDIR=$PRJDIR/3_bwa
OUTPUT_mapped=$OUTDIR/mapped
OUTPUT_unmapped=$OUTDIR/host_removed

NTHREADS=32
listacc="SRR_Acc_List.txt"

cd $REFDIR
bwa index -p daphnia_bwa Dmagna_Dieter_genome.fasta

cd $OUTDIR

for SAMPLE in $(cat acc_lists/$listacc); do
echo "$SAMPLE bwa"
bwa mem -R '@RG\tID:'"${SAMPLE}"'\tSM:'"${SAMPLE}"'\tPL:ILLUMINA\tPI:330' -M -t $NTHREADS $REFDIR/daphnia_bwa $DATADIR/${SAMPLE}_1_val_1.fq.gz $DATADIR/${SAMPLE}_2_val_2.fq.gz | \
samtools view -@ $NTHREADS -hbS - > $OUTDIR/${SAMPLE}.bam

echo "$SAMPLE process unmapped reads"
samtools view -@ $NTHREADS -b -f 12 -F 256 $OUTDIR/${SAMPLE}.bam | \
samtools sort -@ $NTHREADS -n - | \
samtools fastq -@ $NTHREADS - -1 $OUTPUT_unmapped/${SAMPLE}_host_removed_R1.fastq.gz -2 $OUTPUT_unmapped/${SAMPLE}_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n

echo "$SAMPLE process mapped reads"
samtools view -@ $NTHREADS -b -f 3 -q 30 $OUTDIR/${SAMPLE}.bam | \
samtools sort -@ $NTHREADS -n - | \
samtools fixmate -@ $NTHREADS -m -O bam - - | \
samtools sort -@ $NTHREADS - | \
samtools markdup -@ $NTHREADS -r - $OUTPUT_mapped/${SAMPLE}.bam

echo "$SAMPLE index reads"
samtools index -@ $NTHREADS $OUTPUT_mapped/${SAMPLE}.bam

done

rm *.bam
```
