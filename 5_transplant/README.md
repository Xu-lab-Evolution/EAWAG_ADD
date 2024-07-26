### Transplant experiment

```R
## analyze the transplant experiment for duckweed

## read data from KW25 and KW30
data.feedback2022.duckweed.kw30<-readxl::read_excel("./data/FeedbackExperiment2022.xlsx",sheet = "KW30_Duckweed")
data.feedback2022.duckweed.kw25<-readxl::read_excel("./data/FeedbackExperiment2022.xlsx",sheet = "KW25_Duckweed")
PondInfo<-read.table(file="./data/PondInf.txt",sep="\t", header=T,row.names=1);
data<-data.feedback2022.duckweed.kw25
data<-data.feedback2022.duckweed.kw30

Preprocessing<-function(data=data,PondInfo=PondInfo){
  CalResistance<-function(data=NULL, Var=NULL){
    data.out<-matrix(ncol=5,nrow=0);
    Pond<-levels(as.factor(data$Pond));
    Genotype<-levels(as.factor(data$Genotype));
    for(pd in Pond){
      for(gt in Genotype){
        data.sub<-data[data$Pond==pd & data$Genotype==gt,];
        if(nrow(data.sub)==2){
          Resistance<-data.sub[data.sub$Treatment=="+Aphid",Var]-data.sub[data.sub$Treatment=="Control",Var];
          PondTreatment<-as.vector(data.sub[1,"PondTreatment"]);
          Block<-as.vector(data.sub[1,"Block"]);
          data.out<-rbind(data.out,c(pd,Block,gt,PondTreatment, Resistance=Resistance))
        }else{
          #  print ("Error");
          # print (data.sub);
        }
      }
    }
    return(data.out); 
  }
  dat=data.frame(data);
  if( ncol(dat) == 9){
    dat<-dat[is.na(dat$Note),];
    dat<-dat[,-9];
  }else{
  dat<-dat[!is.na(dat$End_FrondN),];
  }
  dat<-dat[!dat$Pond=="6A",];
  dat<-dat[!dat$Pond=="1D",];
  dat<-dat[!dat$Pond=="3D",];
  dat$GR_N<-(log(1+round(dat$End_FrondN,0)) - log(1+round(dat$Start_number,0))) /round(dat$Start_number,0);
  dat$GR_W<-(log (1+round(dat$End_weight,0)) - log(1+round(dat$Start_weight,0)))/round(dat$Start_weight,0);
  dat$PondTreatment<-PondInfo[dat$Pond,"Treatment"];
  dat$Genotype<-as.factor(dat$Genotype);
  dat.ResistN<-CalResistance(dat,Var="GR_N");
  dat.ResistW<-CalResistance(dat,Var="GR_W");
  dat.ResistN<-data.frame(dat.ResistN)
  colnames(dat.ResistN)<-c("Pond","Block","Genotype","PondTreatment","Resistance");
  dat.ResistN$Resistance<-as.numeric(dat.ResistN$Resistance)
  dat.ResistW<-data.frame(dat.ResistW)
  colnames(dat.ResistW)<-c("Pond","Block","Genotype","PondTreatment","Resistance");
  dat.ResistW$Resistance<-as.numeric(dat.ResistW$Resistance);
  dat.resist<-merge(dat.ResistN,dat.ResistW, by.x=c("Pond","Block","Genotype","PondTreatment"),by.y=c("Pond","Block","Genotype","PondTreatment"),all.x = TRUE, all.y = TRUE);
  colnames(dat.resist)[5:6]<-c("ResistN","ResistW");
  return(list(dat,dat.resist));
}

data.feedback2022.duckweed.kw30.GR<-Preprocessing(data=data.feedback2022.duckweed.kw30,PondInfo=PondInfo)[[1]]
data.feedback2022.duckweed.kw30.RS<-Preprocessing(data=data.feedback2022.duckweed.kw30,PondInfo=PondInfo)[[2]]

data.feedback2022.duckweed.kw25.GR<-Preprocessing(data=data.feedback2022.duckweed.kw25,PondInfo=PondInfo)[[1]]
data.feedback2022.duckweed.kw25.RS<-Preprocessing(data=data.feedback2022.duckweed.kw25,PondInfo=PondInfo)[[2]]


## combine the two experiments together
data.feedback2022.duckweed.combined.GR<-rbind(data.feedback2022.duckweed.kw25.GR,data.feedback2022.duckweed.kw30.GR);
data.feedback2022.duckweed.combined.GR$KW<-c(rep('KW25',nrow(data.feedback2022.duckweed.kw25.GR)),rep('KW30',nrow(data.feedback2022.duckweed.kw30.GR)))
head(data.feedback2022.duckweed.combined.GR)

head(data.feedback2022.duckweed.combined.GR)
data.feedback2022.duckweed.combined.GR$Treatment<-factor(data.feedback2022.duckweed.combined.GR$Treatment,levels=c("Control","+Aphid"), ordered = T)
data.feedback2022.duckweed.combined.GR$PondTreatment<-factor(data.feedback2022.duckweed.combined.GR$PondTreatment,levels=c("Control","Aphid"), ordered = T)


data.feedback2022.duckweed.combined.RS<-rbind(data.feedback2022.duckweed.kw25.RS,data.feedback2022.duckweed.kw30.RS);
data.feedback2022.duckweed.combined.RS$KW<-c(rep('KW25',nrow(data.feedback2022.duckweed.kw25.RS)),rep('KW30',nrow(data.feedback2022.duckweed.kw30.RS)))
head(data.feedback2022.duckweed.combined.RS)

## perform statistical analysis for growth rate data without aphid;
head(data.feedback2022.duckweed.combined.GR)
data.feedback2022.duckweed.GR.model1<-lmer(GR_N~PondTreatment*Genotype + (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="Control",])

data.feedback2022.duckweed.GR.model2<-lmer(GR_N~Genotype+PondTreatment+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="Control",])

data.feedback2022.duckweed.GR.model3<-lmer(GR_N~Genotype+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="Control",])


(anova(data.feedback2022.duckweed.GR.model1,data.feedback2022.duckweed.GR.model2))

(anova(data.feedback2022.duckweed.GR.model2,data.feedback2022.duckweed.GR.model3))

Anova(data.feedback2022.duckweed.GR.model2, test.statistic = "F");


data.feedback2022.duckweed.GR.model4<-lmer(GR_N~PondTreatment+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="Control" & 
                                                                                                                            data.feedback2022.duckweed.combined.GR$Genotype=="102",])

Anova(data.feedback2022.duckweed.GR.model4, test.statistic = "F");

## perform statistical analysis for growth rate data with aphid;
data.feedback2022.duckweed.GR.aphid.model1<-lmer(GR_N~PondTreatment*Genotype + (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="+Aphid",])

data.feedback2022.duckweed.GR.aphid.model2<-lmer(GR_N~Genotype+PondTreatment+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="+Aphid",])

data.feedback2022.duckweed.GR.aphid.model3<-lmer(GR_N~Genotype+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="+Aphid",])


(anova(data.feedback2022.duckweed.GR.aphid.model1,data.feedback2022.duckweed.GR.aphid.model2))

(anova(data.feedback2022.duckweed.GR.aphid.model2,data.feedback2022.duckweed.GR.aphid.model3))

Anova(data.feedback2022.duckweed.GR.aphid.model2, test.statistic = "F");


data.feedback2022.duckweed.GR.aphid.model4<-lmer(GR_N~PondTreatment+ (1|Block) + (1|KW), data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="+Aphid" & 
                                                                                                                            data.feedback2022.duckweed.combined.GR$Genotype=="102",])

Anova(data.feedback2022.duckweed.GR.aphid.model4, test.statistic = "F");



# plot GR data, seperate two time points;
PlotCommunityEffects<-function(data=NULL, Var=NULL, ylab=NULL){
  summarySE <- function(data=NULL, measurevar=NULL, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    require(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                     c( N    = length2(xx[,col], na.rm=na.rm),
                        mean = mean   (xx[,col], na.rm=na.rm),
                        sd   = sd     (xx[,col], na.rm=na.rm),
                        sum  = sum     (xx[,col], na.rm=na.rm)
                     )
                   },
                   measurevar,
                   na.rm
    )
    
    # Rename the "mean" column
    # datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
    dat<-data;
    id.num<-grep(Var,colnames(dat));
    colnames(dat)[id.num]<-"TestVar";
  dat.sum<-summarySE(data=dat,measurevar = "TestVar",
                     groupvars = c("PondTreatment","Genotype","KW"));
  dat.sum$PondTreatment<-factor(dat.sum$PondTreatment,levels=c("Control","Aphid"),ordered=T);
  dat.sum$Genotype<-factor(dat.sum$Genotype,levels=c("102","56","58","65"),ordered = T);
  dat.sum.kw25<-dat.sum[dat.sum$KW=="KW25",];
  dat.sum.kw30<-dat.sum[dat.sum$KW=="KW30",];
    
    p.kw25<-ggplot(dat.sum.kw25, aes(x = PondTreatment, y = mean, group=Genotype)) + 
      geom_point(aes(col=PondTreatment),size=5)  +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0,
                    show.legend = FALSE) +
      labs(x="Pond Treatment", y=ylab) +
      geom_line(size=1) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values=c('#00AFBB','darkred'))+
      facet_grid( ~ Genotype)
    p.kw30<-ggplot(dat.sum.kw30, aes(x = PondTreatment, y = mean, group=Genotype)) + 
      geom_point(aes(col=PondTreatment),size=5)  +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0,
                    show.legend = FALSE) +
      labs(x="Pond Treatment", y=ylab) +
      geom_line(size=1) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values=c('#00AFBB','darkred'))+
      facet_grid( ~ Genotype)    
    
    return(list(p.kw25,p.kw30));
      
}

PlotCommunityEffects_MeanKW<-function(data=NULL, Var=NULL, ylab=NULL){
  summarySE <- function(data=NULL, measurevar=NULL, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    require(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                     c( N    = length2(xx[,col], na.rm=na.rm),
                        mean = mean   (xx[,col], na.rm=na.rm),
                        sd   = sd     (xx[,col], na.rm=na.rm),
                        sum  = sum     (xx[,col], na.rm=na.rm)
                     )
                   },
                   measurevar,
                   na.rm
    )
    
    # Rename the "mean" column
    # datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  dat<-data;
  id.num<-grep(Var,colnames(dat));
  colnames(dat)[id.num]<-"TestVar";
  dat.sum<-summarySE(data=dat,measurevar = "TestVar",
                     groupvars = c("PondTreatment","Genotype"));
  dat.sum$PondTreatment<-factor(dat.sum$PondTreatment,levels=c("Control","Aphid"),ordered=T);
  dat.sum$Genotype<-factor(dat.sum$Genotype,levels=c("102","56","58","65"),ordered = T);
#  dat.sum.kw25<-dat.sum[dat.sum$KW=="KW25",];
#  dat.sum.kw30<-dat.sum[dat.sum$KW=="KW30",];
  
 
  p.kw_mean<-ggplot(dat.sum, aes(x = PondTreatment, y = mean, group=Genotype)) + 
    geom_point(aes(col=PondTreatment),size=5)  +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0,
                  show.legend = FALSE) +
    labs(x="Pond Treatment", y=ylab) +
    geom_line(size=1) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values=c('#56B4E9','darkred'))+
    facet_grid( ~ Genotype)    
  
  return(p.kw_mean);
  
}

GR.control_KWmean<-PlotCommunityEffects_MeanKW(data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="Control",],Var="GR_N",ylab = "Growth rate (frond number)");
GR.aphid_KWmean<-PlotCommunityEffects_MeanKW(data=data.feedback2022.duckweed.combined.GR[data.feedback2022.duckweed.combined.GR$Treatment=="+Aphid",],Var="GR_N",ylab = "Growth rate (frond number)");

ggsave(file="./Results/FeedbackGrowthRate_control.pdf",GR.control_KWmean);
ggsave(file="./Results/FeedbackGrowthRate_aphid.pdf",GR.aphid_KWmean);


## analyze the feedback data of the aphid growth rate

data.feedback2022.Aphid.kw25<-readxl::read_excel("FeedbackExperiment2022.xlsx",sheet = "KW25_Aphid")
data.feedback2022.Aphid.kw30<-readxl::read_excel("FeedbackExperiment2022.xlsx",sheet = "KW30_Aphid")
data.feedback2022.Aphid.kw25<-data.frame(data.feedback2022.Aphid.kw25)
data.feedback2022.Aphid.kw30<-data.frame(data.feedback2022.Aphid.kw30)


PreprocessingAphid<-function(data=data,PondInfo=PondInfo){
  dat=data.frame(data);
  dat<-dat[!dat$Pond=="6A",];
  dat<-dat[!dat$Pond=="1D",];
  dat<-dat[!dat$Pond=="3D",];
  dat$GR_N<-log( 1+ ((round(dat$AphidN,0) -log(5))/5));
  dat$PondTreatment<-PondInfo[dat$Pond,"Treatment"];
  dat$PondTreatment<-factor(dat$PondTreatment, levels=c("Control","Aphid"),ordered=T)
  dat$Genotype<-factor(dat$Genotype,levels=c("102","56","58","65"),ordered=T)
  return(dat);
}

data.feedback2022.Aphid.kw25<-PreprocessingAphid(data=data.feedback2022.Aphid.kw25,PondInfo=PondInfo)
data.feedback2022.Aphid.kw30<-PreprocessingAphid(data=data.feedback2022.Aphid.kw30,PondInfo=PondInfo)
data.feedback2022.Aphid.combined<-rbind(data.feedback2022.Aphid.kw25,data.feedback2022.Aphid.kw30);
data.feedback2022.Aphid.combined$KW<-c(rep('KW25',nrow(data.feedback2022.Aphid.kw25)),rep('KW30',nrow(data.feedback2022.Aphid.kw30)))
## perform statistical analysis
aphid.model1<-lmer(GR_N~PondTreatment*Genotype + (1|Block) + (1|KW), data=data.feedback2022.Aphid.combined);
aphid.model2<-lmer(GR_N~PondTreatment+Genotype + (1|Block) + (1|KW), data=data.feedback2022.Aphid.combined);
anova(aphid.model1,aphid.model2)
Anova(aphid.model2, test.statistic = "F")


data.feedback2022.Aphid<-data.feedback2022.Aphid.combined
head(data.feedback2022.Aphid);
# remove containminations and 6A
data.feedback2022.Aphid<-data.feedback2022.Aphid[!data.feedback2022.Aphid$Pond=="6A",];
data.feedback2022.Aphid<-data.feedback2022.Aphid[!data.feedback2022.Aphid$Pond=="1D",];
data.feedback2022.Aphid<-data.feedback2022.Aphid[!data.feedback2022.Aphid$Pond=="3D",];


data.feedback2022.Aphid$PondTreatment<-PondInfo[data.feedback2022.Aphid$Pond,"Treatment"];
data.feedback2022.Aphid$Genotype=as.factor(data.feedback2022.Aphid$Genotype)

data.feedback2022.Aphid$AphidN<-(log(round(data.feedback2022.Aphid$AphidN,0))-log(5))/35;
data.feedback2022.Aphid$PondTreatment<-factor(data.feedback2022.Aphid$PondTreatment, levels=c("Control","Aphid"),ordered=T)
data.feedback2022.Aphid.sum<-summarySE(data=data.feedback2022.Aphid,measurevar = "AphidN", groupvars = c("Genotype","PondTreatment"))
data.feedback2022.Aphid.sum$Genotype<-factor(data.feedback2022.Aphid.sum$Genotype,levels=c("102","56","58","65"),ordered=T)

p<-ggplot(data.feedback2022.Aphid.sum, aes(x = PondTreatment, y = mean, group=Genotype)) + 
  geom_point(aes(col=PondTreatment),size=5)  +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0,
                show.legend = FALSE) +
  labs(x="Pond Treatment", y="AphidGrowth(Frond number)") +
  geom_line(size=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c('#56B4E9','darkred'))+
  facet_grid( ~ Genotype)
p

ggsave(p, file="./Results/AphidGrowth_KWmean.pdf", height  =5.7,width  = 12.1/2);

```
