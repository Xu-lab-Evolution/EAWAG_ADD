```{r load function and packages}

setwd("./")

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

AphidEffects <- function(data, 
                         Treatment = NULL, 
                         Pond = NULL, 
                         Block = NULL, 
                         Time = NULL, 
                         Year = NULL, 
                         Week = NULL, 
                         TimeRemove = NULL, 
                         logTransform = TRUE, 
                         Parameters = "Unknown") {
  # Create results folder if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  ##########################################
  # Helper: summarySE function
  ##########################################
  summarySE <- function(data, measurevar, groupvars, na.rm = TRUE,
                        conf.interval = 0.95, .drop = TRUE) {
    # Function to calculate count, mean, sd, se, and confidence interval.
    length2 <- function(x, na.rm = FALSE) {
      if (na.rm) sum(!is.na(x)) else length(x)
    }
    datac <- ddply(data, groupvars, .drop = .drop,
                   .fun = function(xx, col, na.rm) {
                     c(N = length2(xx[, col], na.rm = na.rm),
                       mean = mean(xx[, col], na.rm = na.rm),
                       sd = sd(xx[, col], na.rm = na.rm),
                       sum = sum(xx[, col], na.rm = na.rm))
                   },
                   measurevar, na.rm)
    datac$se <- datac$sd / sqrt(datac$N)  # Standard error
    ciMult <- qt(conf.interval / 2 + 0.5, datac$N - 1)
    datac$ci <- datac$se * ciMult
    return(datac)
  }
  
  ##########################################
  # Helper: run_model to fit model and report stats
  ##########################################
  run_model <- function(data_tmp, TimeRemove, response_label) {
    # Set lmer control options
    ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    
    # Subset data based on TimeRemove if provided
    if (!is.null(TimeRemove)) {
      data_sub <- data_tmp[data_tmp$Time > TimeRemove, ]
    } else {
      data_sub <- data_tmp
    }
    data_sub$Week <- as.numeric(data_sub$Week)
    data_sub$Year <- factor(data_sub$Year, ordered = TRUE)
    
    # Fit the model
    model <- lmer(Response ~ Treatment * Year + poly(Week, 2) + (Time | Pond),
                  data = data_sub, REML = TRUE, control = ctrl)
    
    # Extract ANOVA table with denominator degrees of freedom
    anova_table <- anova(model)
    
    # Initialize storage for effects and calculate partial eta squared
    effects <- rownames(anova_table)
    effect_size <- rep(NA, nrow(anova_table))
    
    for (i in 1:nrow(anova_table)) {
      Fval <- anova_table$`F value`[i]
      df_effect <- anova_table$NumDF[i]
      df_error <- anova_table$DenDF[i]
      if (!is.na(Fval) && Fval > 0 && !is.na(df_error) && (Fval * df_effect + df_error) != 0) {
        effect_size[i] <- (Fval * df_effect) / (Fval * df_effect + df_error)
      } else {
        effect_size[i] <- NA
      }
    }
    
    # Combine into summary table
    result_table <- data.frame(
      Effect = effects,
      NumDF = anova_table$NumDF,
      DenDF = anova_table$DenDF,
      F_value = anova_table$`F value`,
      p_value = anova_table$`Pr(>F)`,
      Partial_Eta_Squared = effect_size,
      row.names = NULL
    )
    
    # Round numeric columns to two digits, except for p_value
    num_cols <- sapply(result_table, is.numeric)
    for (col in names(result_table)[num_cols]) {
      if (col != "p_value") {
        result_table[[col]] <- round(result_table[[col]], 2)
      }
    }
    # Format p_value column to preserve scientific notation
    result_table$p_value <- sapply(result_table$p_value, function(x) {
      if (is.na(x)) return(NA)
      format(x, scientific = TRUE, digits = 2)
    })
    
    return(list(model = model, stats = result_table))
  }
  
  ##########################################
  # Helper: makePlot to generate and save plot
  ##########################################
  makePlot <- function(data_tmp, stats_table, Resp_name) {
    # Convert Block and Year to factors
    data_tmp$Block <- as.factor(data_tmp$Block)
    data_tmp$Year <- as.factor(data_tmp$Year)
    
    # Build a title that reports the Treatment effect statistics (if available)
    treat_row <- stats_table[grep("Treatment", stats_table$Effect), ]
    if (nrow(treat_row) > 0) {
       title <- sprintf("%s: Treatment effect: F(%d, %d) = %.2f, p = %s", 
                 Resp_name, 
                 as.integer(treat_row$NumDF[1]), 
                 as.integer(treat_row$DenDF[1]), 
                 treat_row$F_value[1], 
                 formatC(as.numeric(treat_row$p_value[1]), format = "e", digits = 2))
    } else {
      title <- Resp_name
    }
    
    # Summarize data for plotting
    data_sum <- summarySE(data = data_tmp, groupvars = c("Time", "Treatment"), measurevar = "Response")
    
    # Create plot
    p <- ggplot(data = data_sum, aes(x = Time, y = mean, colour = Treatment, shape = Treatment, group = Treatment)) +
      geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0, position = position_dodge(0.05)) +
      geom_line(size = 1.2) +
      geom_point(size = 4) +
      ggtitle(title) +
      labs(x = "Calendar Week") +
      xlim(min(data_tmp$Time, na.rm = TRUE) - 1, max(data_tmp$Time, na.rm = TRUE) + 1) +
      theme_classic() +
      scale_color_manual(values = c('darkred', '#56B4E9'))
    
    # Determine file name based on Resp_name (remove extension if necessary) and save in the results folder
    file_name <- ifelse(grepl("\\.", Resp_name),
                        paste0("results/", str_split(Resp_name, "\\.", n = 2)[[1]][1], ".pdf"),
                        paste0("results/", Resp_name, ".pdf"))
    
    ggsave(filename = file_name, plot = p, width = 8, height = 6)
    return(p)
  }
  
  ##########################################
  # Main analysis
  ##########################################
  if (is.null(ncol(data))) { 
    # When data is a vector (single response)
    Resp_name <- Parameters
    data_tmp <- data.frame(Response = as.numeric(as.character(data)), 
                           Treatment = Treatment, 
                           Block = Block, 
                           Week = Week, 
                           Year = Year, 
                           Pond = Pond, 
                           Time = Time)
    data_tmp <- na.omit(data_tmp)
    
    analysis <- run_model(data_tmp, TimeRemove, response_label = Resp_name)
    p <- makePlot(data_tmp, analysis$stats, Resp_name)
    
    # Print the statistical table
    print(analysis$stats)
    
  } else {
    # When data is a data.frame with multiple response columns
    for (i in 1:ncol(data)) {
      Resp_name <- colnames(data)[i]
      # Remove NA values for the response variable
      check_na <- is.na(data[, Resp_name])
      data_sub <- if (any(check_na)) data[!check_na, ] else data
      
      # Only proceed if there is variation in the response
      if (max(data_sub[, Resp_name], na.rm = TRUE) != min(data_sub[, Resp_name], na.rm = TRUE)) {
        if (logTransform) {
          data_tmp <- data.frame(
            Response = log10(as.numeric(as.character(data_sub[, Resp_name])) + 1),
            Treatment = Treatment[!check_na],
            Pond = Pond[!check_na],
            Time = Time[!check_na],
            Block = Block[!check_na],
            Year = Year[!check_na],
            Week = Week[!check_na]
          )
        } else {
          data_tmp <- data.frame(
            Response = as.numeric(as.character(data_sub[, Resp_name])),
            Treatment = Treatment[!check_na],
            Pond = Pond[!check_na],
            Time = Time[!check_na],
            Block = Block[!check_na],
            Year = Year[!check_na],
            Week = Week[!check_na]
          )
        }
        data_tmp <- na.omit(data_tmp)
        
        analysis <- run_model(data_tmp, TimeRemove, response_label = Resp_name)
        p <- makePlot(data_tmp, analysis$stats, Resp_name)
        
        # Print the statistical table for this response
        print(analysis$stats)
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

## set working directory and load the data


```{r Load monitoring data, echo=T}

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

MW_NH4 <- 18.04   # ammonium
MW_PO4 <- 94.97   # phosphate

#ADD.Monitor.all.sub$Time
#head(ADD.Monitor.all)

#is.numeric(ADD.Monitor.all.sub$Time)
```


```{r run statistics , echo=T}

# Clear the results table each time you run if you want a fresh table

results_table <- data.frame(
  Response     = character(),
  Explanatory  = character(),
  NumDF           = numeric(),
  DenDF           = numeric(),
  F_value      = numeric(),
  P_value      = numeric(),
  Effect_size  = numeric(),
  stringsAsFactors = FALSE
)

ADD.Monitor.all.sub<-ADD.Monitor.all.sub %>%
  mutate(LightIntensity=log(as.numeric(Light_weekly_average))) %>%
  mutate(temperature=as.numeric(Temp...C._weekly.average_HOBO.Logger)) %>% 
  mutate(ChlA=as.numeric(ChlA..µg.L.))%>%
  mutate(O2=as.numeric(O2..mg.L.)) %>%
  mutate(Chlorophyta=log(as.numeric(Chlorophyta..cells.l.)+1)) %>%
  mutate(Bacillariophyceae=log(as.numeric(Bacillariophyceae..cells.l.)+1))%>%
  mutate(Streptophyta=log(as.numeric(Streptophyta..cells.l.)+1))%>%
  mutate(dmagna=log(as.numeric(TotalDaphina_magna)+1))%>%
  mutate(ammonium=as.numeric(Ammonium..µg.L.)) %>%
  mutate(phosphate=as.numeric(Phosphat..µg.L.))%>%
  mutate(TotalCarbon=as.numeric(TotalerC..mg.L.)) %>%
  mutate(Cyanobacteria=log(as.numeric(Cyanobacteria..cells.l.)+1)) %>%
  filter(!is.na(ammonium)) %>%
  mutate(NPratio=ammonium * MW_PO4 /(phosphate * MW_NH4)) 
    

colnames(ADD.Monitor.all.sub)
col.export<-c("Pond", "Time",	"Year",	"Week",	"Treatment","AphidDensity",	"DuckweedCoverage",
              "ammonium","phosphate","LightIntensity","temperature",
              "chlA","rotifer","dmagna","cyanobacteria")

write.table(ADD.Monitor.all.sub[,col.export], file="./data/environment.txt",sep="\t", quote = F,row.names = F)
  
  ## remove the data point when temperature is below 15 degrees, as no duckweed and aphids.

# List of variables to be analyzed

Phytoplankton <- c(  "Chlorophyta",
                     "Bacillariophyceae",
                     "Streptophyta",
                     "Cyanobacteria","SimpsonDiv","ChlA")

ADD<-c("AphidDensity","DuckweedCoverage" )

AbioticFactor<-c("ammonium",
                     "phosphate",
                     "TotalCarbon",
                     "temperature",
                     "O2..mg.L.",
                     "pH","LightIntensity","NPratio")

# Loop through each variable and store the returned row in results_table
# for abiotic responses and D. magna, remove the first 7 data points from statistical analysis,when duckweed population started showing differences.
colnames(ADD.Monitor.all.sub)
for (var_name in c(ADD,AbioticFactor,"dmagna")) {
  # Extract the numeric data
  response_vec <- as.numeric(ADD.Monitor.all.sub[[var_name]])
  # Call the function
  res_row <- AphidEffects(data = response_vec,
                          Treatment = ADD.Monitor.all.sub$Treatment,
                          Pond = ADD.Monitor.all.sub$Pond,
                          Block = substr(ADD.Monitor.all.sub$Pond,1,1),
                          Time = ADD.Monitor.all.sub$Time,
                          Week = ADD.Monitor.all.sub$Week,
                          Year = ADD.Monitor.all.sub$Year,
                          TimeRemove = 7, #exclude the data from the firs 7 weeks from statistical analysis, when duckweed showed different population size.
                          Parameters = var_name)
  res_row$Response<-var_name;
  # Append to our results table
  results_table <- rbind(results_table, res_row)
}

# for phytoplankton, remove the first 2 data points from statistics, as the response is quick.
for (var_name in c(Phytoplankton)) {
  # Extract the numeric data
  response_vec <- as.numeric(ADD.Monitor.all.sub[[var_name]])
  # Call the function
  res_row <- AphidEffects(data = response_vec,
                          Treatment = ADD.Monitor.all.sub$Treatment,
                          Pond = ADD.Monitor.all.sub$Pond,
                          Block = substr(ADD.Monitor.all.sub$Pond,1,1),
                          Time = ADD.Monitor.all.sub$Time,
                          Week = ADD.Monitor.all.sub$Week,
                          Year = ADD.Monitor.all.sub$Year,
                          TimeRemove = 2, #exclude the data from the firs 7 weeks from statistical analysis, when duckweed showed different population size.
                          Parameters = var_name)
  res_row$Response<-var_name;
  # Append to our results table
  results_table <- rbind(results_table, res_row)
}


write.table(results_table,file="AphidEffect.txt",sep="\t",quote = F,col.names =T)



# Inspect the final table
print(results_table)
#ADD.Monitor.all.sub$Time<-factor(ADD.Monitor.all.sub$Time, ordered = T)

## analyse the N:P ratio
ADD.Monitor.all.sub.clean<-ADD.Monitor.all.sub %>%
  mutate(Phosphate=as.numeric(Phosphat..µg.L.)) %>%
  filter(!is.na(Phosphate)) %>%
  mutate(ammonium=as.numeric(Ammonium..µg.L.)) %>%
  filter(!is.na(ammonium)) 

ADD.Monitor.all.sub.clean$NP_ratio=as.numeric(ADD.Monitor.all.sub.clean$ammonium) * MW_PO4 /(as.numeric(ADD.Monitor.all.sub.clean$Phosphate) * MW_NH4)

ADD.Monitor.all.sub.clean$N_MC<-14*as.numeric(ADD.Monitor.all.sub.clean$ammonium)/MW_NH4
ADD.Monitor.all.sub.clean$P_MC<-31*as.numeric(ADD.Monitor.all.sub.clean$Phosphate)/MW_PO4
Block<-substr(ADD.Monitor.all.sub.clean$Pond,start = 1,stop = 1)

ADD.Monitor.all.sub.clean[ADD.Monitor.all.sub.clean$Time=="57",c("Treatment","NP_ratio")]

AphidEffects(data=(log(ADD.Monitor.all.sub.clean$NP_ratio)),Block,Pond=ADD.Monitor.all.sub.clean$Pond,Time =  ADD.Monitor.all.sub.clean$Time,Week=ADD.Monitor.all.sub.clean$Week, Year=ADD.Monitor.all.sub.clean$Year, Treatment = ADD.Monitor.all.sub.clean$Treatment,TimeRemove = 2,Parameters  = "NP_ratio")

AphidEffects(data=((ADD.Monitor.all.sub.clean$N_MC)),Block,Pond=ADD.Monitor.all.sub.clean$Pond,Time =  ADD.Monitor.all.sub.clean$Time,Week=ADD.Monitor.all.sub.clean$Week, Year=ADD.Monitor.all.sub.clean$Year, Treatment = ADD.Monitor.all.sub.clean$Treatment,TimeRemove = 2,Parameters  = "N_MC")


ADD.Monitor.all.sub.clean[ADD.Monitor.all.sub.clean$Time>=40,c("Week","NP_ratio","Pond","Treatment")]


## analyze the changes of phytoplankoton between week 30 and 34 in 2021.
ADD.Monitor.all.sub.kw30<-ADD.Monitor.all.sub %>%
  filter(Year=="2021") %>%
  filter(Week>=30 &Week <=34)
 colnames(ADD.Monitor.all.sub.kw30)
 var_name<-"temperature"
 response_vec <- as.numeric(ADD.Monitor.all.sub.kw30[[var_name]])
  # Call the function
  data.tmp<-data.frame(Response=response_vec,
                       Treatment = ADD.Monitor.all.sub.kw30$Treatment,
                          Pond = ADD.Monitor.all.sub.kw30$Pond,
                          Block = substr(ADD.Monitor.all.sub.kw30$Pond,1,1),
                          Time = ADD.Monitor.all.sub.kw30$Time,
                          Week = ADD.Monitor.all.sub.kw30$Week,
                          Year = ADD.Monitor.all.sub.kw30$Year)
      ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

    model <- lmer(Response ~ Treatment  +  (Time | Pond),
                  data = data.tmp, REML = TRUE, control = ctrl)
anova(model)

## calculate the average changes of abiotic factors in 2022
#ADD.Monitor.all.sub$Week
ADD.Monitor.all.sub.clean.2022<-ADD.Monitor.all.sub %>%
  filter(Year==2022) %>%
  filter(Week>30) %>%
  mutate(TotalCarbon=as.numeric(TotalerC..mg.L.)) %>%
  mutate(Temperature=as.numeric(Temp...C._weekly.average_HOBO.Logger)) %>%
  mutate(pH=as.numeric(pH)) %>%
  mutate(Phosphate=as.numeric(Phosphat..µg.L.)) %>%
  mutate(ammonium=as.numeric(Ammonium..µg.L.)) 

ADD.Monitor.all.sub.clean.2022$LightIntensity
summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "ammonium",groupvars = "Treatment")[1,3]/summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "ammonium",groupvars = "Treatment")[2,3] -1

summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "Phosphate",groupvars = "Treatment")[1,3]/summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "Phosphate",groupvars = "Treatment")[2,3] -1

summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "TotalCarbon",groupvars = "Treatment")[1,3]/summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "TotalCarbon",groupvars = "Treatment")[2,3] -1

summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "LightIntensity",groupvars = "Treatment")[1,3]/summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "LightIntensity",groupvars = "Treatment")[2,3] -1
summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "pH",groupvars = "Treatment")[1,3] - summarySE(ADD.Monitor.all.sub.clean.2022,measurevar = "pH",groupvars = "Treatment")[2,3] 


```

