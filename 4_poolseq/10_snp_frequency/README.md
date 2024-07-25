## Generate frequency plots

Input files are located inside the `0_data` folder.

```R
library("tidyverse")
library("lme4")
library("car")
library("parallel")
library("ggpubr")

EAWAG.count<-data.table::fread(file="data/final.dp20_400.AD.txt",sep="\t")
EAWAG.count.clean<-EAWAG.count
EAWAG.count.clean$CHROM <- gsub("D_magna_", "", EAWAG.count.clean$CHROM)

sampleInfo<-read.table(file="data/SequenceID_info.txt",sep="\t",row.names=2,header = T)
data<-EAWAG.count.clean

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
data.sub.long.rescale <- data.sub.long %>%
  mutate(REF_scal = round( ((500*(REF+ALT)-1)/(500+(REF+ALT))) * (REF/(REF+ALT)),digits = 0),
         ALT_scal =  round(((500*(REF+ALT)-1)/(500+(REF+ALT))) * (ALT/(REF+ALT)),digits = 0))
head(data.sub.long.rescale)

snpfilter <- read.table("data/snplist.txt", header = F)
freq.df <- data.sub.long.rescale
freq.df <- freq.df %>% filter(SNP %in% snpfilter$V1)
snp.list <- unique(freq.df$SNP)
freq.df <- freq.df %>% filter(SNP %in% snp.list)

freq.df$Year[freq.df$Year == "Initial"] <- 0
freq.df$Year[freq.df$Year == "2021"] <- 6
freq.df$Year[freq.df$Year == "2022"] <- 16

colnames(freq.df) <- c("snp", "sampleID", "count", "REF", "ALT", "treatment", "generation", "block", "REF_scal", "ALT_scal")

plot.fun <- function(chromosome, snip, pvaluemodel, newchr){
  snp_id <- paste0(chromosome, ".", snip)
  df <- freq.df %>% filter(snp %in% snp_id) %>%
    mutate(freq = ALT_scal/(REF_scal+ALT_scal))
  freq.tab.summ <- df %>% group_by(treatment, generation) %>% summarise(mean = mean(freq), sd = sd(freq), n = n()) %>% mutate(se=sd/sqrt(n))
  freq.tab.summ$generation <- factor(freq.tab.summ$generation, levels=c("0", "6", "16"))
  
  plot <- ggplot(freq.tab.summ, aes(x=generation, y=mean, group=treatment, color=treatment)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                  position=position_dodge(0.05)) +
    theme_classic() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_text(),
          axis.title.y = element_text(),
          axis.text.x = element_text()) +
    geom_segment(aes(x=1,xend=2,
                     y= (freq.tab.summ[which(freq.tab.summ$treatment == "Initial" & freq.tab.summ$generation == "0"),]$mean),
                     yend= (freq.tab.summ[which(freq.tab.summ$treatment == "Control" & freq.tab.summ$generation == "6"),]$mean)),
                 color = "#56b4e8") +
    geom_segment(aes(x=1,xend=2,
                     y= mean(freq.tab.summ[which(freq.tab.summ$treatment == "Initial" & freq.tab.summ$generation == "0"),]$mean),
                     yend=mean(freq.tab.summ[which(freq.tab.summ$treatment == "Aphid" & freq.tab.summ$generation == "6"),]$mean)),
                 color = "#8c0000") +
    scale_fill_manual(name = "Legend", values=c("#252525", "#56b4e8", "#8c0000"), labels = c("Initial", "Control", "Aphid"), breaks = c("Initial", "Control", "Aphid")) +
    scale_color_manual(name = "Legend", values=c("#252525", "#56b4e8", "#8c0000"), labels = c("Initial", "Control", "Aphid"), breaks = c("Initial", "Control", "Aphid")) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_discrete(breaks = c("0", "6", "16"), labels = c("initial", "2021", "2022")) +
    #ggtitle(paste0(snp_id)) +
    annotate("text", x = 0.5, y= c(0.25, 0.15), hjust = 0, label = c(paste0(newchr, ":", snip),
                                                      paste0("P = ", pvaluemodel))) +
    xlab("time") +
    ylab("rel. frequency")
  return(plot)
}

plot01 <- plot.fun("CH4_R", "6748952", "0.005", "CHR04")
plot02 <- plot.fun("CH4_R", "6753940", "0.001", "CHR04")
plot03 <- plot.fun("CH5_L", "6945069", "0.0002", "CHR05")
plot04 <- plot.fun("CH5_L", "6938831", "0.005", "CHR05")


px <- ggarrange(plot01, plot02, plot03, plot04, 
                nrow = 1, ncol = 4, align = "hv", common.legend = T)
```
