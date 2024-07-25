### Download data from NCBI

For this first part, create the following folders:

```bash
mkdir $PRJDIR/0_data
mkdir $PRJDIR/2_trimgalore
```

Inside the folder `0_data` create the file `SRR_Acc_List.txt` and paste the list below:

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

Then dowload the data using the following code:

```bash
listacc="SRR_Acc_List.txt"
vdb-config --prefetch-to-cwd

for i in $(cat $listacc); do
echo $i
prefetch -X 9999999999999 ${i} -O ${i} -f yes
done

find . -mindepth 1 -maxdepth 1 -type d | parallel -j 16 fasterq-dump --split-3 {}/{}/{}\.sra

find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} \;
```

And run some data cleanup using TrimGalore:

```bash
INDIR=$PRJDIR/0_data
OUTDIR=$PRJDIR/2_trimgalore

cd $INDIR

find -name "*_1.fastq" | rev | cut -c 9- | rev | parallel -j 32 trim_galore --illumina --paired --fastqc --gzip -o $OUTDIR/ {}\_1.fastq {}\_2.fastq
```
