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

Generate 3D plot

```R
EAWG.poolseq <- read.pcadapt("./data/final.dp20_400.recode.vcf.gz", type = "vcf");
x <- pcadapt(input = EAWG.poolseq, K = 20);
(x$singular.values**2/sum(x$singular.values**2))[1:5]

poplist.int<-read.table(file="./data/SRRInfo.txt",sep="\t",header=T)
poplist.int<-poplist.int[order(as.numeric(poplist.int$ID)),]

data4plot<-data.frame(unlist(x$scores))
colnames(data4plot)<-paste(rep("PC",20),seq(1:20),sep="");
#dim(data4plot)
data4plot$SampleID<-poplist.int$Run
data4plot$Treatment<-poplist.int$Treatment
data4plot$Year<-poplist.int$Time
data4plot$Year<-as.factor(data4plot$Year)
data4plot$color<-as.numeric(as.factor(data4plot$Treatment))
data4plot$color[data4plot$color==1]<-"darkred";
data4plot$color[data4plot$color==2]<-"#56B4E9";
data4plot$color[data4plot$color==3]<-"grey";
data4plot$shape<-as.numeric(as.factor(data4plot$Year))
data4plot$shape[data4plot$shape==1]<-19
data4plot$shape[data4plot$shape==2]<-11
data4plot$shape[data4plot$shape==3]<-15

pdf(file = "./Results/PCA_plot.pdf", width = 8, height = 6)
EAWAG.3d<-scatterplot3d( data4plot[1:3],
               color = data4plot$color,
              type="h",
  pch=data4plot$shape,
  xlab="PC1 (19.5%)", ylab="PC2 (13.0%)", zlab="PC3 (7.9%)")

legend(EAWAG.3d$xyz.convert(0.6, 0.2, 0.2), legend = c("Control","Aphid","Initial"),
       col =  c('#56B4E9', 'darkred','grey'), pch = 16)
legend(EAWAG.3d$xyz.convert(0.8, 0.2, 0.2), legend = c("Initial","2021","2022"),
        pch = c(19,11,15))
dev.off()
```

Generate Fst plots

```R
poplist.int<-read.table(file="./data/SRRInfo.txt",sep="\t",header=T);
poplist.int<-poplist.int[order(as.numeric(poplist.int$ID)),];

sampleID2021<-poplist.int$Run[poplist.int$year==2021]

sampleID2022<-poplist.int$Run[poplist.int$year==2022]

sampleIDInital<-poplist.int$Run[poplist.int$year==0]

EAWAG.poolseqAll<-vcf2pooldata(vcf.file = "./data/final.dp20_400.recode.vcf.gz", 
                                poolsizes =rep(250,nrow(poplist.int)) ,
                                poolnames = poplist.int$Run, min.cov.per.pool =20,
                                nlines.per.readblock = 1e+06, 
                                verbose = TRUE )

EAWAG.poolseq.allvsall<-compute.fstats(EAWAG.poolseqAll, computeQmat = F,
                                       output.pairwise.div=F,
                                   nsnp.per.bjack.block = 100,computeF4 = F,
                                   verbose = TRUE)

colnames(EAWAG.poolseq.allvsall@pairwise.fst)<-poplist.int$SampleID
row.names(EAWAG.poolseq.allvsall@pairwise.fst)<-poplist.int$SampleID

col = list(year = c("0" = "lightgrey", "2021" = "darkgrey", "2022" = "black"),
          Treatment = c("Aphid" = "darkred", "Control" = "#56B4E9","Initial"="lightgrey") )

Treatment.ann <- HeatmapAnnotation(
  year = poplist.int$Time, Treatment = poplist.int$Treatment,
  col = col
)

fst.hm <- Heatmap(EAWAG.poolseq.allvsall@pairwise.fst,
cluster_rows =TRUE,cluster_columns=TRUE,name="Fst",
column_title = "Fst=(Q1-Q2)/(1 - Q2)",row_split=3,column_split = 3,
top_annotation=Treatment.ann)

save_pdf(
  fst.hm,"./Results/Fst_heatmap.pdf"
)

EAWAG.poolseq.pwfst<-melt(EAWAG.poolseq.allvsall@pairwise.fst) 

colnames(EAWAG.poolseq.pwfst)<-c("Ind1","Ind2","Fst");

row.names(poplist.int)<-poplist.int$SampleID
EAWAG.poolseq.pwfst$Ind1_Treatment<-poplist.int[EAWAG.poolseq.pwfst$Ind1,"Treatment"]
EAWAG.poolseq.pwfst$Ind2_Treatment<-poplist.int[EAWAG.poolseq.pwfst$Ind2,"Treatment"]
EAWAG.poolseq.pwfst$Time<-paste(poplist.int[EAWAG.poolseq.pwfst$Ind1,"Time"],poplist.int[EAWAG.poolseq.pwfst$Ind2,"Time"],sep="_");

EAWAG.poolseq.pwfst.combined<-EAWAG.poolseq.pwfst %>%
  mutate(Fst_type=ifelse(Ind1_Treatment==Ind2_Treatment,"WithinTreatment","BetweenTreatment")) %>%
  filter(!Ind1==Ind2)

EAWAG.poolseq.pwfst.combined$Fst_type[EAWAG.poolseq.pwfst.combined$Ind1_Treatment=="Control" & EAWAG.poolseq.pwfst.combined$Ind2_Treatment=="Control"]<-"WithinControl";

EAWAG.poolseq.pwfst.combined$Fst_type[EAWAG.poolseq.pwfst.combined$Ind1_Treatment=="Aphid" & EAWAG.poolseq.pwfst.combined$Ind2_Treatment=="Aphid"]<-"WithinAphid";
head(EAWAG.poolseq.pwfst.combined)

head(EAWAG.poolseq.pwfst.combined)

EAWAG.poolseq.pwfst.combined.between<-EAWAG.poolseq.pwfst.combined[EAWAG.poolseq.pwfst.combined$Fst_type=="BetweenTreatment" ,]
summary(lm(Fst~Time,data=EAWAG.poolseq.pwfst.combined.between[EAWAG.poolseq.pwfst.combined.between$Time=="2021_2021"| EAWAG.poolseq.pwfst.combined.between$Time=="2022_2022",]))

EAWAG.poolseq.pwfst.combined.within<-EAWAG.poolseq.pwfst.combined[grep("Within",EAWAG.poolseq.pwfst.combined$Fst_type),] %>%
    filter(Time=="Initial_Initial" | Time=="2021_2021"| Time=="2022_2022") %>%
  mutate(Group=as.factor(paste(Fst_type,Time,sep="_")))

model<-lm(Fst~Group,data=EAWAG.poolseq.pwfst.combined.within)
emmeans(model, pairwise ~ Group)
summary(glht(model, linfct = mcp(Group = 'Tukey')))

EAWAG.poolseq.pwfst.combined.betweenYears<-EAWAG.poolseq.pwfst.combined %>%
  filter(Time=="0_2021"|Time=="0_2022"|Time=="2021_2022") 
EAWAG.poolseq.pwfst.combined.betweenYears$Time<-factor(EAWAG.poolseq.pwfst.combined.betweenYears$Time)
model<-lm(Fst~Time,data=EAWAG.poolseq.pwfst.combined.betweenYears)
summary(glht(model, linfct = mcp(Time = 'Tukey')))

EAWAG.poolseq.pwfst.combined.sum<-summarySE(data=EAWAG.poolseq.pwfst.combined[!EAWAG.poolseq.pwfst.combined$Ind1=="C3"|EAWAG.poolseq.pwfst.combined$Ind2=="C3", ],
                                                    measurevar = "Fst",
                                                  groupvars = c("Fst_type","Time"));

EAWAG.poolseq.pwfst.combined.within.sum<-EAWAG.poolseq.pwfst.combined.sum %>%
    filter(Time=="2021_2021"|Time=="2022_2022"|Time=="0_0)") %>%
    filter(!Fst_type=="BetweenTreatment")

EAWAG.poolseq.pwfst.combined.sum<-summarySE(data=EAWAG.poolseq.pwfst.combined,
                                            measurevar = "Fst",
                                            groupvars = c("Fst_type","Time"));
EAWAG.poolseq.pwfst.combined.within.sum<-EAWAG.poolseq.pwfst.combined.sum %>%
  filter(Time=="2021_2021"|Time=="2022_2022"|Time=="0_0") %>%
  filter(!Fst_type=="BetweenTreatment")


EAWAG.poolseq.pwfst.combined.betweenTreatment.sum<-EAWAG.poolseq.pwfst.combined.sum %>%
  filter(Time=="2021_2021"|Time=="2022_2022") %>%
  filter(Fst_type=="BetweenTreatment")

EAWAG.poolseq.pwfst.combined.betweenYear.sum<-summarySE(data=EAWAG.poolseq.pwfst.combined,
                                                        measurevar = "Fst",
                                                        groupvars = c("Time")) %>%
  filter(Time=="0_2021"|Time=="0_2022"|Time=="2021_2022") %>%
  mutate(Fst_type="BetweenYear");

EAWAG.poolseq.pwfst4plot<-rbind(EAWAG.poolseq.pwfst.combined.within.sum,
                                EAWAG.poolseq.pwfst.combined.betweenTreatment.sum,EAWAG.poolseq.pwfst.combined.betweenYear.sum)

EAWAG.poolseq.pwfst4plot$Fst_group<-paste(EAWAG.poolseq.pwfst4plot$Fst_type,EAWAG.poolseq.pwfst4plot$Time,sep="_");
EAWAG.poolseq.pwfst4plot$Fst_group<-factor(EAWAG.poolseq.pwfst4plot$Fst_group,levels=c("WithinTreatment_0_0","WithinControl_2021_2021","WithinControl_2022_2022","WithinAphid_2021_2021",
  "WithinAphid_2022_2022","BetweenTreatment_2021_2021","BetweenTreatment_2022_2022","BetweenYear_0_2021","BetweenYear_0_2022","BetweenYear_2021_2022"),ordered=T)

data.poolstat.summary.plot<- ggplot(EAWAG.poolseq.pwfst4plot, aes(x=Fst_group, y=mean)) + 
  geom_bar( aes(x=Fst_group, y=mean), stat="identity", fill="grey", alpha=0.5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.0,
                position=position_dodge(0.01)) +
  labs(title="Fst", x="Year", y = "Fst")+
  theme_classic() 

ggsave(data.poolstat.summary.plot,file="./Results/Fst_summaryplot.pdf")

```
