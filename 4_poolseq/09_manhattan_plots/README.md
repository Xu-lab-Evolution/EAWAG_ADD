## Generate Manhattan plots

Input files are located inside the `0_data` folder.

### General functions

```R
library("tidyverse")
library("ggtext")
library("ggpubr")

prep.data.plot <- function(df){
  data_cum <- df %>%
    mutate(chr=fct_relevel(chr,c("CH1_L", "CH1_R", "CH2_L", "CH2_R", "CH3_L", "CH3_R", "CH4_L", "CH4_R", "CH5_L", "CH5_R", "CH6", "CH7_L", "CH7_R", "CH8_L", "CH8_R",  "CH9_R", "CH9_L", "CH10_L", "CH10_R"))) %>%
    group_by(chr) %>%
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(chr, bp_add)
  gwas_data <- df %>%
    inner_join(data_cum, by = "chr") %>%
    mutate(bp_cum = pos + bp_add)
  gwas_data$chrOLD <- gwas_data$chr
  gwas_data$chr <- recode(gwas_data$chr,
    CH1_L = "CHR01", CH1_R = "CHR01",
    CH2_L = "CHR02", CH2_R = "CHR02",
    CH3_L = "CHR03", CH3_R = "CHR03",
    CH4_L = "CHR04", CH4_R = "CHR04",
    CH5_L = "CHR05", CH5_R = "CHR05",
    CH6 = "CHR06",
    CH7_L = "CHR07", CH7_R = "CHR07",
    CH8_L = "CHR08", CH8_R = "CHR08",
    CH9_L = "CHR09", CH9_R = "CHR09",
    CH10_L = "CHR10", CH10_R = "CHR10"
  )
  return(gwas_data)
}
```

### Fst

```R
plot.man <- function(df, timepoint){
  axis_set <- df %>% group_by(chr) %>% dplyr::summarize(center = mean(bp_cum))
  chr_ref <- df %>% group_by(chrOLD) %>% dplyr::summarize(bp_add = mean(bp_add))
  chr_ref_fun <- function(x){chr_ref %>% filter(chrOLD == x) %>% select(bp_add) %>% as.numeric()}
  
  manhplot <- ggplot(df, aes(x = bp_cum, y = fst, color = as_factor(chr))) +
    geom_hline(yintercept = min(filter(df, fst > quantile(fst, .999, na.rm = T))$fst), color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = min(filter(df, fst > quantile(fst, .99, na.rm = T))$fst), color = "grey40", linetype = "dotted") +
    geom_vline(xintercept = mean((6592517 + chr_ref_fun("CH4_R")), (6735201 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6732660 + chr_ref_fun("CH4_R")), (6761441 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((788984 + chr_ref_fun("CH7_R")), (785984 + chr_ref_fun("CH7_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6392699 + chr_ref_fun("CH5_R")), (8350017 + chr_ref_fun("CH5_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6813376 + chr_ref_fun("CH5_L")), (7243143 + chr_ref_fun("CH5_L"))), linetype="dotted", color = "orange", size=1) +
    geom_point(alpha = 0.75) +
    geom_point(data = filter(df, fst > quantile(fst, .99, na.rm = T)), aes(x = bp_cum, y = fst), color = "#fc9272") +
    geom_point(data = filter(df, fst > quantile(fst, .999, na.rm = T)), aes(x = bp_cum, y = fst), color = "#ef3b2c") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
    scale_color_manual(values = rep(c("#cccccc", "#636363"), unique(length(axis_set$chr)))) +
    labs(x = NULL, y = paste0("F<sub>ST</sub>", timepoint)) +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_markdown(),
          axis.text.x = element_text()
    )
  return(manhplot)
}

fst.df <- read.csv("data/fst_df_t1.txt", row.names = 1)
t1.plot <- prep.data.plot(fst.df)
snplist1 <- t1.plot$snp
t1.plot <- plot.man(t1.plot, " [2021]")

fst.df <- read.csv("data/fst_df_t2.txt", row.names = 1)
t2.plot <- prep.data.plot(fst.df)
snplist2 <- t2.plot$snp
t2.plot <- plot.man(t2.plot, " [2022]")

px <- ggarrange(t1.plot,
                t2.plot,
                nrow = 2,  align = "hv", common.legend = F)
```
![03_fst](https://github.com/user-attachments/assets/aa843444-060c-42d0-a832-9e1803769d2e)

### GLMM

```R
a1 <- readRDS("data/glmmTMB.rds")
a1 <- a1 %>% separate_wider_delim(cols = snp, delim = ".", names = c("chr", "pos"), too_few = "debug")
a1 <- a1 %>% drop_na() %>% select(chr, pos, p.value)
colnames(a1) <- c("chr", "pos", "p")
a1$pos <- as.numeric(a1$pos)
a1$p <- as.numeric(a1$p)

lm.plot <- prep.data.plot(a1) %>% filter(-log10(p) < 30)

plot.man <- function(df){
  axis_set <- df %>% group_by(chr) %>% dplyr::summarize(center = mean(bp_cum))
  chr_ref <- df %>% group_by(chrOLD) %>% dplyr::summarize(bp_add = mean(bp_add))
  chr_ref_fun <- function(x){chr_ref %>% filter(chrOLD == x) %>% select(bp_add) %>% as.numeric()}

  manhplot <- ggplot(df, aes(x = bp_cum, y = -log10(p), color = as_factor(chr))) +
    geom_hline(yintercept = -log10(0.05/nrow(df)), color = "grey40", linetype = "dashed") +
    geom_vline(xintercept = mean((6592517 + chr_ref_fun("CH4_R")), (6735201 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6732660 + chr_ref_fun("CH4_R")), (6761441 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((788984 + chr_ref_fun("CH7_R")), (785984 + chr_ref_fun("CH7_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6392699 + chr_ref_fun("CH5_R")), (8350017 + chr_ref_fun("CH5_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6813376 + chr_ref_fun("CH5_L")), (7243143 + chr_ref_fun("CH5_L"))), linetype="dotted", color = "orange", size=1) +
    geom_point(alpha = 0.75) +
    geom_point(data = filter(df, p < 0.05/nrow(df)), aes(x = bp_cum, y = -log10(p)), color = "#ef3b2c") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
    scale_color_manual(values = rep(c("#cccccc", "#636363"), unique(length(axis_set$chr)))) +
    labs(x = NULL, y = "-log<sub>10</sub>(p) [GLMM]") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_markdown(),
          axis.text.x = element_text()
    )
  return(manhplot)
}

lm.plot.glmm <- plot.man(lm.plot)
```
![01_lm](https://github.com/user-attachments/assets/ac92de70-eab4-42ae-9e75-b76308ac7d55)

### CMH

```R
a1 <- read.table("data/filt.cmh2021.txt", stringsAsFactors = FALSE, sep = " ")[,c(1,2,4)]
names(a1) <- c("chr", "pos", "p") 
cmh.plot.t1 <- prep.data.plot(a1) %>% mutate(log = -log10(p))

a2 <- read.table("data/filt.cmh2022.txt", stringsAsFactors = FALSE, sep = " ")[,c(1,2,4)]
names(a2) <- c("chr", "pos", "p") 
cmh.plot.t2 <- prep.data.plot(a2) %>% mutate(log = -log10(p))

plot.man <- function(df, timepoint){
  axis_set <- df %>% group_by(chr) %>% dplyr::summarize(center = mean(bp_cum))
  chr_ref <- df %>% group_by(chrOLD) %>% dplyr::summarize(bp_add = mean(bp_add))
  chr_ref_fun <- function(x){chr_ref %>% filter(chrOLD == x) %>% select(bp_add) %>% as.numeric()}
  
  manhplot <- ggplot(df, aes(x = bp_cum, y = -log10(p), color = as_factor(chr))) +
    geom_hline(yintercept = -log10(0.01/nrow(df)), color = "grey40", linetype = "dashed") +
    geom_vline(xintercept = mean((6592517 + chr_ref_fun("CH4_R")), (6735201 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6732660 + chr_ref_fun("CH4_R")), (6761441 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((788984 + chr_ref_fun("CH7_R")), (785984 + chr_ref_fun("CH7_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6392699 + chr_ref_fun("CH5_R")), (8350017 + chr_ref_fun("CH5_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6813376 + chr_ref_fun("CH5_L")), (7243143 + chr_ref_fun("CH5_L"))), linetype="dotted", color = "orange", size=1) +
    geom_point(alpha = 0.75) +
    geom_point(data = filter(df, p < 0.01/nrow(df)), aes(x = bp_cum, y = -log10(p)), color = "#ef3b2c") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values = rep(c("#cccccc", "#636363"), unique(length(axis_set$chr)))) +
    labs(x = NULL, y = paste0("-log<sub>10</sub>(p) [CMH", timepoint, "]")) +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_markdown(),
          axis.text.x = element_text()
    )
  return(manhplot)
}

cmh.plot.t1 <- plot.man(cmh.plot.t1, " 2021")
cmh.plot.t2 <- plot.man(cmh.plot.t2, " 2022")

px <- ggarrange(cmh.plot.t1,
                cmh.plot.t2,
                nrow = 2,  align = "hv", common.legend = F)
```
![02_cmh](https://github.com/user-attachments/assets/d61fa549-6258-49f7-b601-83d14bd4411a)

### CLEAR

```R
a1 <- read.table("data/filt.ctrl.clear", skip = -1, stringsAsFactors = FALSE, sep = ",")
names(a1) <- c("chr", "pos", "p", "snp") 

a2 <- read.table("data/filt.herb.clear", skip = -1, stringsAsFactors = FALSE, sep = ",")
names(a2) <- c("chr", "pos", "p", "snp")

plot.man <- function(df, label){
  axis_set <- df %>% group_by(chr) %>% dplyr::summarize(center = mean(bp_cum))
  chr_ref <- df %>% group_by(chrOLD) %>% dplyr::summarize(bp_add = mean(bp_add))
  chr_ref_fun <- function(x){chr_ref %>% filter(chrOLD == x) %>% select(bp_add) %>% as.numeric()}
  
  manhplot <- ggplot(df, aes(x = bp_cum, y = p, color = as_factor(chr))) +
    geom_hline(yintercept = -log10(0.05/nrow(df)), color = "grey40", linetype = "dashed") +
    geom_vline(xintercept = mean((6592517 + chr_ref_fun("CH4_R")), (6735201 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6732660 + chr_ref_fun("CH4_R")), (6761441 + chr_ref_fun("CH4_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((788984 + chr_ref_fun("CH7_R")), (785984 + chr_ref_fun("CH7_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6392699 + chr_ref_fun("CH5_R")), (8350017 + chr_ref_fun("CH5_R"))), linetype="dotted", color = "orange", size=1) +
    geom_vline(xintercept = mean((6813376 + chr_ref_fun("CH5_L")), (7243143 + chr_ref_fun("CH5_L"))), linetype="dotted", color = "orange", size=1) +
    geom_point(alpha = 0.75) +
    geom_point(data = filter(df, p > -log10(0.05/nrow(df))), aes(x = bp_cum, y = p), color = "#ef3b2c") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(1, 20)) +
    scale_color_manual(values = rep(c("#cccccc", "#636363"), unique(length(axis_set$chr)))) +
    labs(x = NULL, y = label) +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_markdown(),
          axis.text.x = element_text()
    )
  return(manhplot)
}

clear.plot.ctrl <- prep.data.plot(a1)
clear.plot.herb <- prep.data.plot(a2)

clear.plot.ctrl <- plot.man(clear.plot.ctrl, "-log<sub>10</sub>(p) [CLEAR ctrl]")
clear.plot.herb <- plot.man(clear.plot.herb, "-log<sub>10</sub>(p) [CLEAR herb]")

px <- ggarrange(clear.plot.ctrl,
                clear.plot.herb,
                nrow = 2,  align = "hv", common.legend = F)
```
![04_clear](https://github.com/user-attachments/assets/460bdf8b-d93c-4189-a757-5f3808b9f921)
