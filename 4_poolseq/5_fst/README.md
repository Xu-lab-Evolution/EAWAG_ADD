### Calculate Fst

For this first part, create the following folders:

```bash
mkdir $PRJDIR/4_mpileup_sync/sync_fst
mkdir $PRJDIR/4_mpileup_sync/fst_res
```

Organize data:

```bash
input=$PRJDIR/4_mpileup_sync/sync/daphnia2.sync
OUTDIR=$PRJDIR/4_mpileup_sync/sync_fst

awk '{print $1,$2,$3,$35,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/A1.sync
awk '{print $1,$2,$3,$24,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/A3.sync
awk '{print $1,$2,$3,$27,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/B2.sync
awk '{print $1,$2,$3,$32,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/B6.sync
awk '{print $1,$2,$3,$26,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/C2.sync
awk '{print $1,$2,$3,$31,$28,$34,$23,$29,$33,$25,$30}' $input > $OUTDIR/C6.sync

awk '{print $1,$2,$3,$9,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T1A.sync
awk '{print $1,$2,$3,$6,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T1D.sync
awk '{print $1,$2,$3,$4,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T2B.sync
awk '{print $1,$2,$3,$20,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T2C.sync
awk '{print $1,$2,$3,$18,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T3A.sync
awk '{print $1,$2,$3,$15,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T3D.sync
awk '{print $1,$2,$3,$14,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T6B.sync
awk '{print $1,$2,$3,$13,$8,$7,$5,$19,$17,$16,$12}' $input > $OUTDIR/T6C.sync
```

Run the following R code:

```R
library("poolfstat")
library("tidyverse")
library("ggtext")
library("ggpubr")

calc.fst <- function(file, timepoint, sample){
  pnames <- as.character(c('A1', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'))
  psizes <- as.numeric(c('500', '500', '500', '500', '500', '500', '500', '500'))
  SG.pooldata <- popsync2pooldata(sync.file = file,
                                  poolsizes = psizes, poolnames = pnames,
                                  min.rc = 4, min.cov.per.pool = 20, max.cov.per.pool = 400,
                                  min.maf = 0.01, noindel = TRUE, nlines.per.readblock = 1e+06)
  
  PairwiseFST <- compute.pairwiseFST(SG.pooldata, output.snp.values = T)
  
  df <- data.frame(chr = SG.pooldata@snp.info$Chromosome,
                   pos = SG.pooldata@snp.info$Position,
                   fst = rowMeans((PairwiseFST@PairwiseSnpFST)[,c(1:7)]),
                   timepoint = timepoint,
                   sample = sample)
  return(df)
}

a1 <- calc.fst(file = "sync_fst/A1.sync", timepoint = "T2", sample = "A1")
write.csv(a1, "fst_res/a1.txt")

a3 <- calc.fst(file = "sync_fst/A3.sync", timepoint = "T2", sample = "A3")
write.csv(a3, "fst_res/a3.txt")

b2 <- calc.fst(file = "sync_fst/B2.sync", timepoint = "T2", sample = "B2")
write.csv(b2, "fst_res/b2.txt")

b6 <- calc.fst(file = "sync_fst/B6.sync", timepoint = "T2", sample = "B6")
write.csv(b6, "fst_res/b6.txt")

c2 <- calc.fst(file = "sync_fst/C2.sync", timepoint = "T2", sample = "C2")
write.csv(c2, "fst_res/c2.txt")

c6 <- calc.fst(file = "sync_fst/C6.sync", timepoint = "T2", sample = "C6")
write.csv(c6, "fst_res/c6.txt")


a1 <- calc.fst(file = "sync_fst/T1A.sync", timepoint = "T1", sample = "A1")
write.csv(a1, "fst_res/ta1.txt")

a3 <- calc.fst(file = "sync_fst/T3A.sync", timepoint = "T1", sample = "A3")
write.csv(a3, "fst_res/ta3.txt")

b2 <- calc.fst(file = "sync_fst/T2B.sync", timepoint = "T1", sample = "B2")
write.csv(b2, "fst_res/tb2.txt")

c2 <- calc.fst(file = "sync_fst/T2C.sync", timepoint = "T1", sample = "C2")
write.csv(c2, "fst_res/tc2.txt")

d1 <- calc.fst(file = "sync_fst/T1D.sync", timepoint = "T1", sample = "D1")
write.csv(d1, "fst_res/td1.txt")

d3 <- calc.fst(file = "sync_fst/T3D.sync", timepoint = "T1", sample = "D3")
write.csv(d3, "fst_res/td3.txt")

b6 <- calc.fst(file = "sync_fst/T6B.sync", timepoint = "T1", sample = "B6")
write.csv(b6, "fst_res/tb6.txt")

c6 <- calc.fst(file = "sync_fst/T6C.sync", timepoint = "T1", sample = "C6")
write.csv(c6, "fst_res/tc6.txt")

#############################################################################################################

a1 <- read.csv("fst_res/a1.txt", row.names = 1) %>% mutate(snp = paste0(chr, ".", pos))
a3 <- read.csv("fst_res/a3.txt") %>% mutate(snp = paste0(chr, ".", pos))
b2 <- read.csv("fst_res/b2.txt") %>% mutate(snp = paste0(chr, ".", pos))
b6 <- read.csv("fst_res/b6.txt") %>% mutate(snp = paste0(chr, ".", pos))
c2 <- read.csv("fst_res/c2.txt") %>% mutate(snp = paste0(chr, ".", pos))
c6 <- read.csv("fst_res/c6.txt") %>% mutate(snp = paste0(chr, ".", pos))

snplist1 <- Reduce(intersect, list(a1$snp,
                                   a3$snp,
                                   b2$snp,
                                   b6$snp,
                                   c2$snp,
                                   c6$snp))  

a1 <- a1 %>% filter(snp %in% snplist1)
a3 <- a3 %>% filter(snp %in% snplist1)
b2 <- b2 %>% filter(snp %in% snplist1)
b6 <- b6 %>% filter(snp %in% snplist1)
c2 <- c2 %>% filter(snp %in% snplist1)
c6 <- c6 %>% filter(snp %in% snplist1)

fst.df.t2 <- data.frame(snp = snplist1,
                     chr = a1$chr,
                     pos = a1$pos,
                     fst = rowMeans(cbind(a1$fst,
                                          a3$fst,
                                          b2$fst,
                                          b6$fst,
                                          c2$fst,
                                          c6$fst)))

fst.df.t2$chr <- gsub("D_magna_", "", fst.df.t2$chr)
fst.df.t2$snp <- gsub("D_magna_", "", fst.df.t2$snp)
fst.df.t2 <- fst.df.t2 %>% na.omit()
write.csv(fst.df.t2, "fst_res/fst_df_t2.txt")


a1 <- read.csv("fst_res/ta1.txt", row.names = 1) %>% mutate(snp = paste0(chr, ".", pos))
a3 <- read.csv("fst_res/ta3.txt") %>% mutate(snp = paste0(chr, ".", pos))
b2 <- read.csv("fst_res/tb2.txt") %>% mutate(snp = paste0(chr, ".", pos))
c2 <- read.csv("fst_res/tc2.txt") %>% mutate(snp = paste0(chr, ".", pos))
d1 <- read.csv("fst_res/td1.txt") %>% mutate(snp = paste0(chr, ".", pos))
d3 <- read.csv("fst_res/td3.txt") %>% mutate(snp = paste0(chr, ".", pos))
b6 <- read.csv("fst_res/tb6.txt") %>% mutate(snp = paste0(chr, ".", pos))
c6 <- read.csv("fst_res/tc6.txt") %>% mutate(snp = paste0(chr, ".", pos))

snplist2 <- Reduce(intersect, list(a1$snp,
                                   a3$snp,
                                   b2$snp,
                                   c2$snp,
                                   d1$snp,
                                   d3$snp,
                                   b6$snp,
                                   c6$snp))  

a1 <- a1 %>% filter(snp %in% snplist2)
a3 <- a3 %>% filter(snp %in% snplist2)
b2 <- b2 %>% filter(snp %in% snplist2)
c2 <- c2 %>% filter(snp %in% snplist2)
d1 <- d1 %>% filter(snp %in% snplist2)
d3 <- d3 %>% filter(snp %in% snplist2)
b6 <- b6 %>% filter(snp %in% snplist2)
c6 <- c6 %>% filter(snp %in% snplist2)

fst.df.t1 <- data.frame(snp = snplist2,
                     chr = a1$chr,
                     pos = a1$pos,
                     fst = rowMeans(cbind(a1$fst,
                                          a3$fst,
                                          b2$fst,
                                          c2$fst,
                                          d1$fst,
                                          d3$fst,
                                          b6$fst,
                                          c6$fst)))

fst.df.t1$chr <- gsub("D_magna_", "", fst.df.t1$chr)
fst.df.t1$snp <- gsub("D_magna_", "", fst.df.t1$snp)
fst.df.t1 <- fst.df.t1 %>% na.omit()
write.csv(fst.df.t1, "fst_res/fst_df_t1.txt")

snplist <- intersect(fst.df.t2$snp, fst.df.t1$snp)
write.table(snplist, "fst_res/snplist.txt", row.names = F, col.names = F, quote=FALSE)
```

Files `fst_df_t1.txt`, `fst_df_t2.txt`, and `snplist.txt` are available inside the folder `0_data`.
