### Generate SYNC files

For this first part, create the following folders:

```bash
mkdir $PRJDIR/4_mpileup_sync
mkdir $PRJDIR/4_mpileup_sync/mpileup
mkdir $PRJDIR/4_mpileup_sync/sync

cd $PRJDIR/4_mpileup_sync
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/popoolation2/popoolation2_1201.zip
unzip popoolation2_1201.zip
```

Then, run PoPoolation2:

```bash
DATADIR=$PRJDIR/3_bwa/mapped
REFDIR=$PRJDIR/reference_genome
OUTDIR=$PRJDIR/4_mpileup_sync
POPOOLATION2=$OUTDIR/popoolation2_1201

cd $OUTDIR

find $DATADIR -type f -name "*.bam" -printf "%p\n" | sort > $OUTDIR/filelist.txt

function processSYNC {
chr=$1
chrNEW=$(echo ${chr} | sed s/"D_magna_"//)
PRJDIR=######### ADD PATH TO PROJECT DIRECTORY HERE =#########
REFDIR=$PRJDIR/reference_genome
OUTDIR1=$PRJDIR/4_mpileup_sync/mpileup
OUTDIR2=$PRJDIR/4_mpileup_sync/sync
POPOOLATION2=$PRJDIR/4_mpileup_sync/popoolation2_1201
echo "Processing chromosome $chr"
samtools mpileup -b $PRJDIR/4_mpileup_sync/filelist.txt -f $REFDIR/genome.fasta -o $OUTDIR1/daphnia.${chr}.mpileup -r $chr
java -jar $POPOOLATION2/mpileup2sync.jar --input $OUTDIR1/daphnia.${chr}.mpileup --output $OUTDIR2/daphnia.${chrNEW}.sync --fastq-type sanger --min-qual 20 --threads 2
rm $OUTDIR1/daphnia.${chr}.mpileup
}

export -f processSYNC

chrArray=("D_magna_CH1_L" "D_magna_CH1_R" "D_magna_CH2_L" "D_magna_CH2_R" "D_magna_CH3_L" "D_magna_CH3_R" "D_magna_CH4_L" "D_magna_CH4_R" "D_magna_CH5_L" "D_magna_CH5_R" "D_magna_CH6" "D_magna_CH7_L" "D_magna_CH7_R" "D_magna_CH8_L" "D_magna_CH8_R" "D_magna_CH10_L" "D_magna_CH10_R" "D_magna_CH9_R" "D_magna_CH9_L")
parallel -j 12 processSYNC ::: "${chrArray[@]}"

cd $OUTDIR/sync

cat *.sync > daphnia.sync
```
