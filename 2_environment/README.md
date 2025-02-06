Input files are in the folder `data`.

Create general functions:
```R
# List of packages for session
.packages = c("ggplot2", "plyr", "data.table","EnvStats","tibble","tidyr","qqman","stringr",
              "vcfR","poolfstat","aod","glmmTMB","scatterplot3d","lme4",
             "tidyverse","CMplot","tidyr","broom.mixed","parallel","hrbrthemes","reshape2",
             "multcomp",   "emmeans","devtools","clusterProfiler","ontologyIndex","vegan","car","gridExtra","pcadapt","ComplexHeatmap","pheatmap","dendextend","tidyHeatmap",  "lmerTest")


# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) BiocManager::install(.packages[!.inst],update = T,ask = F)

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

# load required funcitons.

AphidEffects<-function(data=data, Treatment=NULL, Pond=NULL, Block=NULL,Time=NULL, Year=NULL, Week=NULL, TimeRemove=NULL, logTransform=TRUE, Parameters="Unknown"){
  library(readxl)
  library(vegan)
  library(ggplot2)
  library(lme4)
  library("car")
  library("stringr")
  TimeRemove=TimeRemove;
  summarySE <- function(data=NULL, measurevar=NULL, groupvars=NULL, na.rm=T,
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

  makePlot<-function(data.tmp=NULL, TimeRemove=NULL){
       data.tmp$Block<-as.factor(data.tmp$Block)
       data.tmp$Year<-as.factor(data.tmp$Year) 
    if (!is.null(TimeRemove)){
      data.sub<-data.tmp[data.tmp$Time>TimeRemove,];
      ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
      Aphid.model<-lmer(Response~Treatment+Year + poly(Week, 2) + (Week|Block), data=data.sub,REML = T,control = ctrl);
      Aphid.model.p<- format(Anova(Aphid.model,test="F")$`Pr(>F)`,scientific=T, digits=2)[1];
      Aphid.model.DF<-Anova(Aphid.model,test="F")$Df[1];
      Aphid.model.Fvalue<-sprintf("%.2f", Anova(Aphid.model,test="F")$F[1])
      
    }else{
      ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
      Aphid.model<-lmer(Response~Treatment+Year + poly(Week, 2) + (Week|Block), data=data.tmp,REML = T,control = ctrl);
      Aphid.model<-lmer(Response~Treatment+Time + (1 + Pond|Time), data=data.tmp);
      Aphid.model.p<- format(Anova(Aphid.model,test="F")$`Pr(>F)`,scientific=T, digits=2)[1];
      Aphid.model.DF<-Anova(Aphid.model,test="F")$Df[1];
      Aphid.model.Fvalue<-sprintf("%.2f", Anova(Aphid.model,test="F")$F[1])
    }
    
    # make the plot;
    # put the information into the plot
    ID<-(str_split(Resp.name, "/",n=3 )[[1]][3])
    title<-paste(ID,"P =",Aphid.model.p,"F-value=",Aphid.model.Fvalue,sep=" ")
    data.tmp.sum<-summarySE(data=data.tmp,groupvars = c("Time","Treatment"), measurevar = "Response");
    p<-ggplot(data=data.tmp.sum, aes(x=Time, y = mean, colour = Treatment, shape=Treatment,
                                     group = Treatment)) + 
      geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.0, 
                    position=position_dodge(0.05)) +
      geom_line(linewidth=2) +
      geom_point(size=6)+
      ggtitle(title) +
      labs(x="Calendar Week", size=22)+
      xlim(-1, 75)+
      theme_classic()
    p<-p + scale_color_manual(values=c('darkred','#56B4E9'))
    p
    if(is.na(match(Resp.name,"\\."))){
      filename=paste(Resp.name,"pdf",sep="."); 
    }else{
      filename<-paste(str_split(Resp.name, "\\.",n=2 )[[1]][1],"pdf",sep=".");
    }
    ggsave(p,file=filename,width = 8, height = 6)
  }
  
  if(is.null(ncol(data))){
    ## if only a vector is given;
    Resp.name<-Parameters;
    data.tmp<-data.frame(Response=as.numeric(as.character(data)),Treatment=Treatment, Block=Block, Week=Week, Year=Year, Pond=Pond, Time=Time);
    data.tmp<-na.omit(data.tmp);
    makePlot(data.tmp,TimeRemove);
    
  }else{
    for(i in 1:ncol(data)){
      # calculate P-value;
      Resp.name<-colnames(data)[i];
      check.na<-is.na(data[,Resp.name]);
      if(is.na(max(data[,Resp.name]))){ # if there is NA, remove them;
        data.sub<-data[!check.na,];
      }else{
        data.sub<-data;
      }
      if(max(data.sub[,Resp.name])!=min(data.sub[,Resp.name])){
        if(logTransform){ # log transform the data if needed;
          data.tmp<-data.frame(Response=log10(as.numeric(as.character(data.sub[,Resp.name]+1))),Treatment=Treatment[!check.na], Pond=Pond[!check.na], Time=Time[!check.na]);
        }else{
          data.tmp<-data.frame(Response=as.numeric(as.character(data.sub[,Resp.name])),Treatment=Treatment[!check.na], Pond=Pond[!check.na], Time=Time[!check.na]);
        }
        data.tmp<-na.omit(data.tmp);
        makePlot(data.tmp,TimeRemove);
      }
    }
  }
  
}
  
summarySE <- function(data=NULL, measurevar=NULL, groupvars=NULL, na.rm=T,
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

```

Read and prepare data:

```R
## read all data into one excell sheet;

Data.all<- sapply(readxl::excel_sheets("./data/ADD-Monitoring_v2.xlsx"), simplify = F, USE.NAMES = F,
                  function(X) readxl::read_excel("./data/ADD-Monitoring_v2.xlsx", sheet = X))
PondInfo<-read.table(file="./data/PondInf.txt",sep="\t", header=T,row.names=1);
## 26 data sheet
head(Data.all)
  ADD.Monitor<-data.frame();
# read all data, except the phytoplantkon data (in a different format).
for(i in 1 :(length(Data.all)-1)){
  data.tmp<-data.frame(Data.all[i]);
  colnames(data.tmp)<-data.tmp[1,];
  data.tmp<-data.tmp[-1,];
  colnames(data.tmp)[1]<-"Pond";
  data.tmp<-melt(data.tmp, id=c("Pond"));  
  colnames(data.tmp)<-c("Pond","Time",colnames(data.frame(Data.all[i]))[1]);
  data.tmp$Time<-as.numeric(as.character( data.tmp$Time))
  if(nrow(ADD.Monitor)==0){
    ADD.Monitor<-data.tmp;
  }else{
    ADD.Monitor<-merge(ADD.Monitor,data.tmp,by.x=c("Pond","Time"),by.y=c("Pond","Time"),all.x = TRUE, all.y = TRUE);
  }
}

dim(ADD.Monitor)
head(ADD.Monitor)

# merge the data with phytoplankton
phytoplanton<-data.frame(Data.all[length(Data.all)]);
phytoplanton$Time<-as.numeric(phytoplanton$Week);
head(phytoplanton)
phytoplanton<-phytoplanton[,-c(4,5)];
phytoplanton$SimpsonDiv<-diversity(log(phytoplanton[,c(5:16)]+1),index = "simpson");


ADD.Monitor.all<-merge(ADD.Monitor,phytoplanton,by.x=c("Pond","Time"),by.y=c("Pond","Time"),all.x = TRUE, all.y = TRUE);
ADD.Monitor.all$Treatment<-PondInfo[ADD.Monitor.all$Pond,"Treatment"];


# removed the data from pond 6A, because snails were not established in this pond.
ADD.Monitor.all<-ADD.Monitor.all[!ADD.Monitor.all$Pond=="6A",]

# for statistical analysis, we removed the data from the first two weeks (as aphid population remained very low) and the data after 73th week (when the last sample from Daphnia was taken). Duckweed coverage was monitored for a longer time.

ADD.Monitor.all.sub<-ADD.Monitor.all[ADD.Monitor.all$Time<=73,];
Block<-substr(ADD.Monitor.all.sub$Pond,start = 1,stop = 1)

ADD.Monitor.all.sub$Year=ifelse((ADD.Monitor.all.sub$Time>=27),2022,2021);
ADD.Monitor.all.sub$Week=ifelse((ADD.Monitor.all.sub$Year==2021),ADD.Monitor.all.sub$Time+25,ADD.Monitor.all.sub$Time+25-52);
```

Plot:

```R

AphidEffects(data=ADD.Monitor.all.sub$DuckweedCoverage,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2, Parameters  = "./results/DuckweedCoverage");

AphidEffects(data=ADD.Monitor.all.sub$ChlA..µg.L.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/ChlA");

AphidEffects(data=ADD.Monitor.all.sub$AphidDensity,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/AphidDensity")

AphidEffects(data=ADD.Monitor.all.sub$Phosphat..µg.L.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Phosphat")

AphidEffects(data=ADD.Monitor.all.sub$Total.Phosphor..µg.L.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/TotalP")

AphidEffects(data=ADD.Monitor.all.sub$Ammonium..µg.L.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Ammonium")

AphidEffects(data=ADD.Monitor.all.sub$Temp...C.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Temperature")

AphidEffects(data=ADD.Monitor.all.sub$O2..mg.L.,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/O2")

AphidEffects(data=log(ADD.Monitor.all.sub$Chlorophyta..cells.l.+1),Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Chlorophyta")

AphidEffects(data=log(ADD.Monitor.all.sub$Bacillariophyceae..cells.l.+1),Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Bacillariophyceae")

AphidEffects(data=log(ADD.Monitor.all.sub$Streptophyta..cells.l.+1),Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Streptophyta")

AphidEffects(data=log(ADD.Monitor.all.sub$TotalDaphina_magna+1),Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/TotalDaphina_magna")

AphidEffects(data=ADD.Monitor.all.sub$pH,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/PH")

AphidEffects(data=ADD.Monitor.all.sub$SimpsonDiv,Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/Diversity")

AphidEffects(data=(as.numeric(ADD.Monitor.all.sub$ChlA..µg.L.)),Block = Block,Pond=ADD.Monitor.all.sub$Pond,Time =  ADD.Monitor.all.sub$Time,Week=ADD.Monitor.all.sub$Week, Year=ADD.Monitor.all.sub$Year, Treatment = ADD.Monitor.all.sub$Treatment,TimeRemove = 2,Parameters  = "./results/ChlA")



```
