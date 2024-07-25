### Run CMH

For this first part, create the following folders:

```bash
mkdir $PRJDIR/4_mpileup_sync/cmh
```

The file `snplist.tx` is in the folder `0_data`.

Run PoPoolation2:

```bash
POPOOLATION2=$PRJDIR/4_mpileup_sync/popoolation2_1201
DATADIR=$PRJDIR/4_mpileup_sync/sync
OUTDIR=$PRJDIR/4_mpileup_sync/cmh

cd $OUTDIR

perl $POPOOLATION2/cmh-test.pl --input $DATADIR/daphnia2.sync --output cmh_test2021.txt \
 --min-count 2 --min-coverage 20 --max-coverage 400 \
 --population 4-3,5-6,2-1,16-17,13-12,14-15,9-10,9-11

perl $POPOOLATION2/cmh-test.pl --input $DATADIR/daphnia2.sync --output cmh_test2022.txt \
 --min-count 2 --min-coverage 20 --max-coverage 400 \
 --population 26-32,31-26,22-23,25-24,20-21,30-21,27-28,27-29

awk -F'\t' -v OFS='\t' '{ $(NF+1)=$1"."$2 ; print}' cmh_test2021.sim.bed |
awk '{$3 = ""; print}' |
sed 's/D_magna_//g' > tmp.2021.txt

awk -F'\t' -v OFS='\t' '{ $(NF+1)=$1"."$2 ; print}' cmh_test2022.sim.bed |
awk '{$3 = ""; print}' |
sed 's/D_magna_//g' >  tmp.2022.txt

awk -F' ' 'NR==FNR{ids[$1]; next} $4 in ids' snplist.txt tmp.2021.txt > filt.cmh2021.txt

awk -F' ' 'NR==FNR{ids[$1]; next} $4 in ids' snplist.txt tmp.2022.txt > filt.cmh2022.txt
```

Files `filt.cmh2021.txt` and `filt.cmh2022.txt` are available inside the folder `00_data`.
