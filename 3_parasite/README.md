
Input files are in the `data` folder.

```R
dat1a <-read.csv("./data/resistotypes_binary.csv")
treat <-read.csv("./data/treatments.csv")

names(dat1a)
names(treat)

dat1a$host_group <- substring(dat1a$host,1,10)
dat1a$host_group2 <- substring(dat1a$host,1,11)

dat1a$group <- "unknown"
dat1a$group <-  ifelse(dat1a$host_group == "CH-H-2014-", "Swisspondpanel",
                 ifelse(dat1a$host_group == "CH-H-2016-", "Swisspondpanel",      
                  ifelse(dat1a$host_group == "CH-H-2019-", "Swisspondpanel", 
                   ifelse(dat1a$host_group == "CH-H-EEE0-", "EEE_start",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-1A", "EEE_1A",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-1B", "EEE_1B",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-1C", "EEE_1C",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-1D", "EEE_1D",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-2A", "EEE_2A",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-2B", "EEE_2B",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-2C", "EEE_2C",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-2D", "EEE_2D",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-3A", "EEE_3A",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-3B", "EEE_3B",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-3C", "EEE_3C",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-3D", "EEE_3D",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-6A", "EEE_6A",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-6B", "EEE_6B",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-6C", "EEE_6C",
                          ifelse(dat1a$host_group2 == "CH-H-EEE-6D", "EEE_6D", 
                                 "unknow2"))))))))))))))))))))
     
table(dat1a$group, exclude=NULL)                                                        
table(treat$group, exclude=NULL)                                                        

dat1 <- merge(dat1a, treat, by="group")


dat1$host_group <- NULL
dat1$host_group2 <- NULL

# produce binomial variables (replace R and S with 0 and 1).

dat1$C1_b <- ifelse(dat1$C1 == "R", 0,
                    ifelse(dat1$C1 == "S", 1, 
                           ifelse(is.na(dat1$C1), NA, NA)))
dat1$C19_b <- ifelse(dat1$C19 == "R", 0,
                    ifelse(dat1$C19 == "S", 1, 
                           ifelse(is.na(dat1$C19), NA, NA)))
dat1$P15_b <- ifelse(dat1$P15 == "R", 0,
                    ifelse(dat1$P15 == "S", 1, 
                           ifelse(is.na(dat1$P15), NA, NA)))
dat1$P20_b <- ifelse(dat1$P20 == "R", 0,
                    ifelse(dat1$P20 == "S", 1, 
                           ifelse(is.na(dat1$P20), NA, NA)))
dat1$P21_b <- ifelse(dat1$P21 == "R", 0,
                    ifelse(dat1$P21 == "S", 1, 
                           ifelse(is.na(dat1$P21), NA, NA)))

dat1$P38_b <- ifelse(dat1$P38 == "R", 0,
                     ifelse(dat1$P38 == "S", 1, 
                            ifelse(is.na(dat1$P38), NA, NA)))

dat1$P50_b <- ifelse(dat1$P50 == "R", 0,
                     ifelse(dat1$P50 == "S", 1, 
                            ifelse(is.na(dat1$P50), NA, NA)))
dat1$P54_b <- ifelse(dat1$P54 == "R", 0,
                     ifelse(dat1$P54 == "S", 1, 
                            ifelse(is.na(dat1$P54), NA, NA)))


# Summarize data and make graphics

names(dat1)
table(dat1$aphid, dat1$EEE_sample)
dat2 = subset(dat1, select = c(EEE_sample, group, aphid, C1_b, C19_b, P15_b, P20_b, P21_b, 
                               P38_b, P50_b, P54_b))

aggdata <-aggregate(dat2, by=list(dat2$aphid, dat2$EEE_sample, dat2$group), FUN=mean, na.rm=TRUE)
names(aggdata)
aggdata$aphid <- aggdata$Group.1
aggdata$EEE_sample <- aggdata$Group.2
aggdata$group <- aggdata$Group.3
aggdata$Group.1 <- NULL
aggdata$Group.2 <- NULL
aggdata$Group.3 <- NULL
aggdata$block <- NULL

aggdata$time_treat <- paste(aggdata$EEE_sample, aggdata$aphid)
print(aggdata) 

table(aggdata$time_treat)
aggdata$time_treat<-factor(aggdata$time_treat,levels=c("2021/6 zero", "2021/10 no", "2021/10 yes",
                                                      "2022/7 no", "2022/7 yes" ), ordered=T)

p1 <- ggplot(aggdata, aes(x=time_treat, y=C1_b)) + 
  geom_bar(position = 'dodge', fun="mean", stat = 'summary', fill = c('yellow','#56B4E9','darkred','#56B4E9','darkred')) +
  ylim(0,1) +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0) +
  geom_point(shape=16,col="grey") +
  theme_minimal() +
  ggtitle("P. ramosa C1") + 
  theme(text = element_text(size = 8)) +
  labs(x = "Treatment", y ="Parasite attachment")


p2 <- ggplot(aggdata, aes(x=time_treat, y=C19_b)) + 
  geom_bar(position = 'dodge', fun="mean", stat = 'summary', fill=c('yellow','#56B4E9','darkred','#56B4E9','darkred')) +
  ylim(0,1) +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0) +
  geom_point(shape=16,col="grey") +
  theme_minimal() +
  ggtitle("P. ramosa C19") + 
  theme(text = element_text(size = 8)) +
   labs(x = "Treatment", y ="Parasite attachment")

p3 <- ggplot(aggdata, aes(x=time_treat, y=P15_b)) + 
  geom_bar(position = 'dodge', fun="mean", stat = 'summary', fill = c('yellow','#56B4E9','darkred','#56B4E9','darkred')) +
  ylim(0,1) +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0) +
  geom_point(shape=16,col="grey") +
  theme_minimal() +
  ggtitle("P. ramosa P15") + 
  theme(text = element_text(size = 8)) +
   labs(x = "Treatment", y ="Parasite attachment")


p4 <- ggplot(aggdata, aes(x=time_treat, y=P20_b)) + 
  geom_bar(position = 'dodge', fun="mean", stat = 'summary', fill = c('yellow','#56B4E9','darkred','#56B4E9','darkred')) +
  ylim(0,1) +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0) +
  geom_point(shape=16,col="grey") +
  theme_minimal() +
  ggtitle("P. ramosa P20") + 
  theme(text = element_text(size = 8)) +
  labs(x = "Treatment", y ="Parasite attachment")


p5 <- ggplot(aggdata, aes(x=time_treat, y=P21_b)) + 
  geom_bar(position = 'dodge', fun="mean", stat = 'summary', fill=c('yellow','#56B4E9','darkred','#56B4E9','darkred')) +
  ylim(0,1) +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0) +
  geom_point(shape=16,col="grey") +
  theme_minimal() +
  ggtitle("P. ramosa P21") + 
  theme(text = element_text(size = 8)) +
   labs(x = "Treatment", y ="Parasite attachment")


p.all<-grid.arrange(p1, p2, p3,p4, p5, ncol = 5, nrow =1)
ggsave(file="./Results/Figure2.pdf",p.all)



##################################################################
# statistical analysis using GLM;
# binomial regression

dat_time1 <- dat1[(!dat1$EEE_sample=="2021/6"),]
head(dat_time1)

mod_C1  <- glmer(C1_b ~ aphid + (1+EEE_sample|block), 
                data=dat_time1, family=binomial(link = "logit"))
summary(mod_C1)



mod_C19  <- glmer(C19_b ~ aphid + (1+EEE_sample|pond_id), 
                 data=dat_time1, family=binomial(link = "logit"))

summary(mod_C19)

mod_P15  <- glmer(P15_b ~ aphid  + (1+EEE_sample|pond_id), 
                  data=dat_time1, family=binomial(link = "logit"))
summary(mod_P15)

mod_P20  <- glmer(P20_b ~ aphid  + (1+EEE_sample|pond_id), 
                  data=dat_time1, family=binomial(link = "logit"))

summary(mod_P20)

mod_P21  <- glmer(P21_b ~ aphid  + (1+EEE_sample|pond_id), 
                  data=dat_time1, family=binomial(link = "logit"))


summary(mod_P21)
```
