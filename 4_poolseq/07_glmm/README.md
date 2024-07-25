### Run GLMM

Input files `final.dp20_400.AD.txt`, `SequenceID_info.txt`, and `snplist.txt` are inside the `00_data` folder. This is computationally intensive, so the output `glmmTMB.rds` is inside the `00_data` folder as well.

This is the R code to run the model using the inputs above:

```R
library("data.table")
library("tidyr")
library("lme4")
library("dplyr")
library("broom.mixed")
library("aod")
library("glmmTMB")
library("parallel")
library("pbmcapply")

nthreads <- ########## ADD NUMBER OF THREADS ##########

EAWAG.count<-fread(file="final.dp20_400.AD.txt",sep="\t")

# remove inital populations
drop.cols <- c("SRR19749275","SRR19749274","SRR19749263","SRR19749262")
EAWAG.count.clean<-EAWAG.count %>% select(-all_of(drop.cols))
EAWAG.count.clean$CHROM <- gsub("D_magna_", "", EAWAG.count.clean$CHROM)

sampleInfo<-read.table(file="SequenceID_info.txt",sep="\t",row.names=2,header = T)
data<-EAWAG.count.clean

## before lunch this, check the file size. it migth take long time to finish the run.
data.SNP<-paste(data$CHROM,data$POS,sep=".");
data.sub<-data[,-c(1:5)];
data.sub$SNP<-data.SNP
data.sub.long<-data.sub %>%
  pivot_longer(!SNP, names_to = "SampleID", values_to = "Count");
AD.info<-unlist(strsplit(data.sub.long$Count,","));
data.sub.long$REF<-as.numeric(AD.info[(1:length(AD.info)) %%2 ==0]);
data.sub.long$ALT<-as.numeric(AD.info[(1:length(AD.info)) %%2 ==1]);
head(data.sub.long)
data.sub.long$SampleID <- gsub(".AD", "", data.sub.long$SampleID, fixed = TRUE)
data.sub.long$Treatment<-sampleInfo[data.sub.long$SampleID,"Treatment"];
data.sub.long$Year<-sampleInfo[data.sub.long$SampleID,"Year"];
data.sub.long$Block<-sampleInfo[data.sub.long$SampleID,"Block"];
# rescale the count using NEF approach
data.sub.long.rescale <- data.sub.long %>%
  mutate(REF_scal = round( ((500*(REF+ALT)-1)/(500+(REF+ALT))) * (REF/(REF+ALT)),digits = 0),
         ALT_scal =  round(((500*(REF+ALT)-1)/(500+(REF+ALT))) * (ALT/(REF+ALT)),digits = 0))
head(data.sub.long.rescale)
## run glmTMP model for each SNP.

snpfilter <- read.table("snplist.txt", header = F)
freq.df <- data.sub.long.rescale
#freq.df$SNP <- gsub("_", ".", freq.df$SNP, fixed = TRUE)
freq.df <- freq.df %>% filter(SNP %in% snpfilter$V1)
snp.list <- unique(freq.df$SNP)
#snp.list <- sample(snp.list, 100)
freq.df <- freq.df %>% filter(SNP %in% snp.list)

test.fun <- function(data, snp){
  df <- data[data$SNP==snp,]
  set.seed(100)
  glmmTMB.out<- tidy(glmmTMB( (cbind(REF_scal,ALT_scal)) ~ Treatment + (1|Year) + (1|Block), family=betabinomial(link = "logit"), data=df), effects = "fixed")
  glmmTMB.out.clean<-glmmTMB.out %>%
    filter(term=="TreatmentControl")
  output<- data.frame(snp = snp,
                      estimate = glmmTMB.out.clean$estimate,
                      statistic = glmmTMB.out.clean$statistic,
                      p.value = glmmTMB.out.clean$p.value)
  return(output)
}

results <- pbmclapply(snp.list, test.fun, data = freq.df, mc.cores = nthreads, mc.preschedule = FALSE)
results2 <- do.call(rbind, results)

saveRDS(results2, file = paste0("glmmTMB.rds"))
```
