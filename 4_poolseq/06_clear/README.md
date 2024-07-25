### Run CLEAR

For this first part, create the following folders:

```bash
mkdir $PRJDIR/4_mpileup_sync/sync_chr
mkdir $PRJDIR/4_mpileup_sync/sync_chr/clear_in
mkdir $PRJDIR/4_mpileup_sync/sync_chr/clear_out
```

Dowload CLEAR and edit the code to output -log10(p):

```bash
git clone --recursive https://github.com/airanmehr/CLEAR.git
# NEED TO EDIT CLEAR.PY!!!!!!!!!
# Open CLEAR.py
# locate a.to_pickle(options.out) towards the end of the script
# change to a.to_csv(options.out)
# also add these two lines after L65 and 70 (the ifelse function before the output)
# f=lambda x: x.alt-x.null
# a=f(a[0.5])
```

Create file `daphnia.sync.ctrl.pops`:
```
0,1     0,2     0,3     0,4     0,5     0,6     0,7     0,8     6,1     6,2     6,3     6,4     6,5     6,6     6,8     16,1    16,2    16,3    16,4    16,5    16,6    16,8
```

Create file `daphnia.sync.herb.pops`:
```
0,1     0,2     0,3     0,4     0,5     0,6     0,7     0,8     6,1     6,2     6,3     6,4     6,5     6,6     6,7     6,8     16,1    16,5    16,3    16,7    16,4    16,8
```

Organize data:

```bash
function processCLEAR {
chr=$1
chrNEW=$(echo ${chr} | sed s/"D_magna_"//)
PRJDIR=####### ADD PATH TO PROJECT DIRECTORY #######
INDIR=$PRJDIR/4_mpileup_sync/sync
OUTDIR=$PRJDIR/4_mpileup_sync/sync_chr/clear_in

echo "Processing chromosome $chr"
awk -v var="$chr" '{ if ($1 == var) { print } } ' $INDIR/daphnia2.sync > $OUTDIR/daphnia.${chrNEW}.sync
cd $OUTDIR
awk '{print $1,$2,$3,$22,$21,$11,$10,$22,$21,$11,$10,$9,$6,$4,$20,$18,$15,$14,$13,$35,$24,$27,$32,$26,$31}' daphnia.${chrNEW}.sync > daphnia.${chrNEW}.temp.sync
sed 's/ /\t/g' daphnia.${chrNEW}.temp.sync > ${chrNEW}.herb.sync
rm daphnia.${chrNEW}.temp.sync
awk '{print $1,$2,$3,$22,$21,$11,$10,$22,$21,$11,$10,$8,$7,$5,$19,$17,$16,$12,$28,$34,$23,$29,$33,$25,$30}' daphnia.${chrNEW}.sync > daphnia.${chrNEW}.temp.sync
sed 's/ /\t/g' daphnia.${chrNEW}.temp.sync > ${chrNEW}.ctrl.sync
rm daphnia.${chrNEW}.temp.sync
rm daphnia.${chrNEW}.sync
}

export -f processCLEAR
chrArray=("D_magna_CH1_L" "D_magna_CH1_R" "D_magna_CH2_L" "D_magna_CH2_R" "D_magna_CH3_L" "D_magna_CH3_R" "D_magna_CH4_L" "D_magna_CH4_R" "D_magna_CH5_L" "D_magna_CH5_R" "D_magna_CH6" "D_magna_CH7_L" "D_magna_CH7_R" "D_magna_CH8_L" "D_magna_CH8_R" "D_magna_CH10_L" "D_magna_CH10_R" "D_magna_CH9_R" "D_magna_CH9_L")
parallel -j $NTHREADS processCLEAR ::: "${chrArray[@]}"
```

Run CLEAR for each chromosome and treatment:

```bash
PRJDIR=####### ADD PATH TO PROJECT DIRECTORY #######
DATADIR=$PRJDIR/4_mpileup_sync/sync_chr/clear_in
OUTDIR=$PRJDIR/4_mpileup_sync/sync_chr/clear_out
CLEARDIR=$PRJDIR/4_mpileup_sync/CLEAR

chr="XXX" # REPLACE XXX WITH THE CHROMOSOME NAME
type="YYY" # REPLACE YYY WITH EITHER herb OR ctrl ACCORDING TO THE TREATMENT
chrNEW=$(echo ${chr} | sed s/"D_magna_"//)

cd $DATADIR
echo "Processing chromosome $chr -- $type ..."
cp $PRJDIR/clear_in/daphnia.sync.${type}.pops $DATADIR/${chrNEW}.${type}.sync.pops
python $CLEARDIR/CLEAR.py --sync ${chrNEW}.${type}.sync --out=$OUTDIR/${chrNEW}.${type}.clear
```

Organize output:

```bash
cat *.ctrl.clear > tmp.ctrl.clear
awk -F',' -v OFS=',' '{ $(NF+1)=$1"."$2 ; print}' tmp.ctrl.clear |
sed 's/D_magna_//g' |
sed '1s/^/CHR,POS,P,SNPID\n/' > final.ctrl.clear

cat *.herb.clear > tmp.herb.clear
awk -F',' -v OFS=',' '{ $(NF+1)=$1"."$2 ; print}' tmp.herb.clear |
sed 's/D_magna_//g' |
sed '1s/^/CHR,POS,P,SNPID\n/' > final.herb.clear

rm tmp.ctrl.clear tmp.herb.clear

awk -F',' 'NR==FNR{ids[$1]; next} $4 in ids' ../sync_fst/snplist.txt final.ctrl.clear > filt.ctrl.clear
rm final.ctrl.clear

awk -F',' 'NR==FNR{ids[$1]; next} $4 in ids' ../sync_fst/snplist.txt final.herb.clear > filt.herb.clear
rm final.herb.clear
```

Files `filt.ctrl.clear` and `filt.herb.clear`  are available inside the folder `0_data`.
