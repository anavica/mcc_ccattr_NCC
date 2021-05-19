##############################################################################
###### The burden of heat-related mortality attributable to recent 
###### human-induced climate - Vicedo-Cabrera et al. 2021 Nature Clim Change
##############################################################################


## This code has been implemented to reproduce the main analysis applied in the 
## study referenced above. It uses temperature-mortality data of 10 hypothetical 
## cities and the corresponding simulated temperature for 10 models defining the 
## two scenarios. 
## 
## The analysis is divided in three stages: (1) 1st stage analysis (time-series 
## analysis), (2) Meta-regression (to extract the BLUPs and MMTs), 
## (3) Quantification of impacts under the two scenarios (factual / counterfactual)

# LOAD LIBRARIES
library(dlnm) ; library(mixmeta) ; library(splines) 
library(lubridate) ; library(MASS) ;library(ggplot2) ;library(patchwork)
library(RColorBrewer) ; library(ggtext)

# LOAD OBSERVED DATA
load("10cities_data.RData")

##############################################################################
####### - 1st STAGE ANALYSIS
##############################################################################

### 1) SET PARAMETERS OF THE MODEL
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 10
lagnk <- 2

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 4

# DEGREE OF FREEDOM FOR TREND
dftrend <- 1

# MODEL FORMULA (SUMMER FORMULA)
formula <- deaths ~ cb + dow + ns(yday,df=dfseas):factor(year) + ns(date,df=round(length(unique(year))/dftrend/10))

### 2) CREATE EMPTY OBJECTS TO STORE THE COEFF-VCOV & OTEHR CITY-SPECIFIC INFO
ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)
coefall <- matrix(NA,nrow(metatable),ncoef,dimnames=list(metatable$city))
vcovall <- vector("list",nrow(metatable))
names(vcovall) <- metatable$city

# LOOP ACROSS CITIES
for(i in seq(length(listdata))) {

  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  datacity <- listdata[[i]]
  
  # SUBSET TO THE SUMMER PERIOD
  data <- subset(datacity, month %in% listmonthmat[i,])

  # CREATE YDAY
  data$yday <- yday(data$date)

  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun,knots=quantile(data$tmean,varper/100,na.rm=T), 
                 Bound=range(data$tmean,na.rm=T))
  arglag <- list(knots=logknots(lag,lagnk))

  cb <- crossbasis(data$tmean,lag=lag,argvar=argvar,arglag=arglag, group=data$indsummer)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")

  # REDUCTION TO OVERALL CUMULATIVE
  redall <- crossreduce(cb,model,cen=mean(data$tmean,na.rm=T))
  coefall[i,] <- coef(redall)
  vcovall[[i]] <- vcov(redall)
}


##############################################################################
####### - 2nd STAGE ANALYSIS
##############################################################################

# DEFINE AVERAGE ESTIMATES OF THE META-PREDICTORS BY CITY
avgtmean <- sapply(listdata,function(x) mean(x$tmean,na.rm=T))
rangetmean <- sapply(listdata,function(x) IQR(x$tmean,na.rm=T))

# RUN META-REGRESSION MODEL (with mixmeta package)
mvmlall <- mixmeta(coefall ~ avgtmean + rangetmean,
  vcovall,metatable,control=list(showiter=T, igls.inititer=10),method="reml")

# OBTAIN BLUPS
blupall <- blup(mvmlall,vcov=T)

# ESTIMATION MMT FROM BLUPS
minperccity <- mintempcity <- rep(NA,length(listdata))
names(mintempcity) <- names(minperccity) <- metatable$city

# DEFINE MINIMUM MORTALITY VALUES
for(i in seq(length(listdata))) {
  # EXTRACT THE DATA
  datacity <- listdata[[i]]
  # SUBSET TO THE SUMMER PERIOD
  data <- subset(datacity, month %in% listmonthmat[i,])
  predvar <- quantile(data$tmean,25:98/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (25:98)[which.min((bvar%*%blupall[[i]]$blup))]
  mintempcity[i] <- quantile(data$tmean,minperccity[i]/100,na.rm=T)
}


##############################################################################
####### - QUANTIFICATION OF THE IMPACTS
##############################################################################

# LOAD SIMULATED TEMPERATURE FOR THE TWO SCENARIOS
load("TmeanModelled.RData")

# DEFINE SCENARIOS
scenario <- c("fct", "cnfact")

# DEFINE MODELS
gcm <- substr(unique(tmean_hist$model),1,3)

# NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION OF ATTRIBUTABLE RISK
nsim <- 1000

# CREATE AN ARRAY TO STORE THE ATTRIBUTABLE DEATHS (WITH CI)
# STRATIFIED BY EST/CI, RANGE, PERIOD, RCP
ancity <- afcity <- array(0,dim=c(nrow(metatable),length(gcm)+1,3,3)
    ,dimnames=list(metatable$city,c(gcm, "ensemble"),c("est","ci.l","ci.u"),
    c(scenario,"dif")))
antotal <- aftotal <- array(0,dim=c(length(gcm)+1,3,3),
     dimnames=list(c(gcm, "ensemble"),c("est","ci.l","ci.u"),c(scenario,"dif")))

# CREATE ARRAYS TO STORE INFO ON PROJECTED AVERAGE TEMPERATURE
tmeanyear <- array(NA,c(nrow(metatable),length(gcm),length(1991:2018),length(scenario)),
  dimnames=list(metatable$city,gcm,1991:2018,scenario))

# FUNCTIONS USED FOR CALIBRATION
source("fhempel_main.R") ## to derive the corr parameters
source("fhempel_corr.R")

# CREATE A VECTOR TO STORE THE PROJECTED MORTALITY IN EACH PERIOD
deathbase <- rep(0,nrow(metatable))
names(deathbase) <- metatable$city

# CREATE TEMPORARY OBJECT TO STORE SIMULATED RESULTS USED FOR UNCERTAINTY
ancitysim <- array(0,dim=c(nrow(metatable),length(gcm),3,nsim+1),
  dimnames=list(metatable$city,gcm,c(scenario,"dif"),
  c("est",paste0("sim",seq(nsim)))))

# LOOP ACROSS GCM
for (g in seq(length(gcm))){
    
    # PRINT
    cat("\n\n",gcm[g],"\n")

    # LOOP ACROSS CITIES
    for(i in seq(nrow(metatable))) {
    
          # PRINT
          cat(i,"")
    
          # SELECT OBSERVED DATA
          datacity <- listdata[[i]]
    
          # SUBSET TO THE SUMMER PERIOD
          data <- subset(datacity, month %in% listmonthmat[i,])    
  
          # DEFINE ARGVAR (AS IN ESTIMATION), CENTERING AND COEF-VCOV
          argvar <- list(fun=varfun,knots=quantile(data$tmean,varper/100,
            na.rm=T),Bound=range(data$tmean,na.rm=T))
          cen <- mintempcity[i]
  
          # SELECT PROJECTED SERIES FOR CITY i
          # NB: transform into data.frame, otherwise the fhempel function does not work
          hist.sel <- data.frame(tmean_hist[tmean_hist$city==names(listdata)[i] & 
                                   tmean_hist$model==unique(tmean_hist$model)[grep(gcm[g],unique(tmean_hist$model))],])
          nat.sel <- data.frame(tmean_nat[tmean_nat$city==names(listdata)[i] & 
                                 tmean_nat$model==unique(tmean_nat$model)[grep(gcm[g],unique(tmean_nat$model))],])
          
          # CHECK NO REPETITION IN OVERLAP HIST-SSP
          if (nrow(hist.sel)!=length(unique(tmean_hist$date))){
                      hist.sel <- aggregate(hist.sel$ta, 
                      by = list(hist.sel$date), 
                      FUN = mean) 
          colnames(hist.sel) <- c("date", "ta")
          }

          # ESTIMATE CORR FACTORS BETWEEN HISTORICAL AND OBSERVED
          corr <- fhempel_main(listdata[[i]][c("date","tmean")],hist.sel[c("date","ta")],output="correction")

          # APPLY CALIBRATION & SELECT THE SERIES
          ind1 <- substr(hist.sel$date,6,10)!="02-29"
          hist.sel <- hist.sel[ind1,]
          hist.proj <- fhempel_corr(hist.sel[c("date","ta")], corr) 
          hist.proj <- subset(hist.proj, month(hist.proj$date) %in% listmonthmat[i,])   
          
          ind1 <- substr(nat.sel$date,6,10)!="02-29" 
          nat.sel <- nat.sel[ind1,]
          nat.proj <- fhempel_corr(nat.sel[c("date","ta")], corr) 
          nat.proj <- subset(nat.proj, month(nat.proj$date) %in% listmonthmat[i,]) 

          # INFO ON PROJECTED TEMPERATURE
          tmeanyear[i,,,"fct"] <- tapply(hist.proj[,"ta"],
            substr(hist.proj[,"date"],1,4),mean)
          
          tmeanyear[i,,,"cnfact"] <- tapply(nat.proj[,"ta"],
            substr(nat.proj[,"date"],1,4),mean)
        
          # DEFINE PROJECTED MORTALITY SERIES
          ind1 <- substr(data$date,6,10)!="02-29"
          deathobs <- data$deaths[ind1]
          deathdoy <- tapply(deathobs,yday(data$date)[ind1],mean,na.rm=T)
          while(any(isna <- is.na(deathdoy)))
            deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
          deathproj <- rep(deathdoy,length=nrow(hist.proj))

          # STORE THE BASELINE TOTAL MORTALITY
          deathbase[i] <- sum(deathdoy)*28
  
          # EXTRACT PARAMETERS
          coef <- blupall[[i]]$blup
          vcov <- blupall[[i]]$vcov
  
          # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
          set.seed(13041975+g)
          coefsim <- mvrnorm(nsim,coef,vcov)
  
          # FACTUAL
          
            # DERIVE THE CENTERED BASIS
            bvar <- do.call(onebasis,c(list(x=hist.proj[,2]),argvar))
            cenvec <- do.call(onebasis,c(list(x=cen),argvar))
            bvarcen <- scale(bvar,center=cenvec,scale=F)
      
            # INDICATOR FOR COLD/HEAT DAYS
            indheat <- hist.proj[,2]>cen
      
            # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
            an <- (1-exp(-bvarcen%*%coef))*deathproj
      
            # SUM BY RANGE AND PERIOD, STORE BEFORE THE ITERATIONS
            # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
            ancitysim[i,g,"fct",1] <- sum(an[indheat])
        
            # LOOP ACROSS ITERATIONS
            for(s in seq(nsim)) {
          
                # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
                an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
                
                # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
                # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
                ancitysim[i,g,"fct",s+1] <- sum(an[indheat])
              }
          
          # COUNTERFACTUAL
          
            # DERIVE THE CENTERED BASIS
            bvar <- do.call(onebasis,c(list(x=nat.proj[,2]),argvar))
            cenvec <- do.call(onebasis,c(list(x=cen),argvar))
            bvarcen <- scale(bvar,center=cenvec,scale=F)
      
            # INDICATOR FOR COLD/HEAT DAYS
            indheat <- nat.proj[,2]>cen
      
            # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
            an <- (1-exp(-bvarcen%*%coef))*deathproj
      
            # SUM BY RANGE AND PERIOD, STORE BEFORE THE ITERATIONS
            # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
            ancitysim[i,g,"cnfact",1] <- sum(an[indheat])
        
            # LOOP ACROSS ITERATIONS
            for(s in seq(nsim)) {
          
                # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
                an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
                
                # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
                # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
                ancitysim[i,g,"cnfact",s+1] <- sum(an[indheat])
  
            }
    }
}


################################################################################
# DIFFERENCE BETWEEN SCENARIOS

  ancitysim[,,3,] <- ancitysim[,,1,] - ancitysim[,,2,]

################################################################################
# SUM ACCROSS SIMULATIONS AT A HIGHER LEVEL
  
  # TOTAL
  antotalsim <- apply(ancitysim[,,,],2:4,sum)

################################################################################
# ATTRIBUTABLE DEATHS - BY GCM 
  
  # BY CITY, WITH CI
  ancity[,-11,1,] <- apply(ancitysim[,,,-1],1:3,mean)
  ancity[,-11,2,] <- apply(ancitysim[,,,-1],1:3,quantile,0.025)
  ancity[,-11,3,] <- apply(ancitysim[,,,-1],1:3,quantile,0.975)

  # TOTAL, WITH CI
  antotal[-11,1,] <- antotalsim[,,1]
  antotal[-11,2,] <- apply(antotalsim[,,-1],1:2,quantile,0.025,na.rm=T)
  antotal[-11,3,] <- apply(antotalsim[,,-1],1:2,quantile,0.975,na.rm=T)

# ATTRIBUTABLE DEATHS - ENSEMBLE 
  
  # BY CITY, WITH CI
  ancity[,11,1,] <- apply(ancitysim[,,,1],c(1,3),mean)
  ancity[,11,2,] <- apply(ancitysim[,,,-1],c(1,3),quantile,0.025)
  ancity[,11,3,] <- apply(ancitysim[,,,-1],c(1,3),quantile,0.975)

  # TOTAL, WITH CI
  antotal[11,1,] <- apply(antotalsim[,,1],2,mean)
  antotal[11,2,] <- apply(antotalsim[,,-1],2,quantile,0.025,na.rm=T)
  antotal[11,3,] <- apply(antotalsim[,,-1],2,quantile,0.975,na.rm=T)
  
  
################################################################################
# ATTRIBUTABLE FRACTIONS 

  afcity[,,,] <- ancity[,,,]/deathbase*100
  aftotal[,,] <- antotal[,,]/
    as.numeric(sum(deathbase))*100
################################################################################
# PROPORTION OF MORTALITY ATTRIBUTED TO ANTHROPOGENIC CLIMATE CHANGE

apcity <- (afcity[,11,1,3] / afcity[,11,1,1])*100
aptotal <- (aftotal[11,1,3] /aftotal[11,1,1] )*100
  
rm(ancitysim, antotalsim, deathdoy, deathbase, deathobs, deathproj, 
   nat.sel, nat.proj, data, an, datacity, hist.sel, hist.proj)
  
##############################################################################
####### - SUMMARIES & PLOTS
##############################################################################


# 1. TS PLOT FACTUAL - COUNTERFACTUAL
# (soon)
# 2. EXPOSURE-RESPONSE FUNCTIONS
listerplot <- list()

for(i in seq(length(names(listdata)))) {
  datacity <- listdata[[i]]
  data <- subset(datacity, month %in% listmonthmat[i,])
  argvar <- list(x=data$tmean,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  pred2 <- crosspred(bvar,coef=blupall[[i]]$blup,
                     vcov=blupall[[i]]$vcov,
  model.link="log",cen=mintempcity[i],by=0.1, from=range(data$tmean, na.rm=T)[1],
    to=quantile(data$tmean, 0.997,na.rm=T))
  res <- data.frame(cbind(pred2$predvar, pred2$allRRfit, pred2$allRRlow, pred2$allRRhigh))
  names(res) <- c("temp", "rr", "rrlow", "rrhigh")

 listerplot[[i]] <- ggplot(res, aes(x=temp, y=rr)) +
      geom_line(color="red", size=1.5) + 
      geom_ribbon(aes(ymin=rrlow, ymax=rrhigh), linetype=2, alpha=0.2)+
      ylim(c(0.6, 6))+
      labs(y="Relative Risk", x="Temperature (ºC)", title=names(listdata[i]))+
      geom_vline(xintercept=quantile(data$tmean,0.99, na.rm=T),
                 linetype="dashed", color = "red") +  
      geom_hline(yintercept=1,
                 linetype="solid", color = "grey40") + 
      theme_minimal()+
      theme(axis.line = element_line(size =0.5),
            axis.text = element_text(size=12), 
            axis.title = element_text(size=13),
            plot.title = element_text(size=18, face="bold"),
            plot.subtitle = element_text(size=15, face="italic"),
            axis.ticks = element_line(size=0.5 ),
            axis.ticks.length = unit(4, "pt"))
}

erplot <- wrap_plots(listerplot, ncol=4)

erplot




# 3. BAR PLOT + DOTS

res_af <- data.frame(afcity[,11,1,c(1,2)])
res_af$id <- c(1:nrow(res_af))
label_data <- rownames(res_af)
col_reg <- brewer.pal(10, "Paired")

# BAR PLOT WITH THE TWO SCENARIOS
barplot <- ggplot(data=res_af) +  
                geom_crossbar(aes(x = id, y = 0, ymin = 0, 
                                  ymax = fct, fill = factor(id)),
                              fatten = 0, colour = "white") + 
                geom_crossbar(aes(x = id, y = 0, ymin = 0, 
                                  ymax = cnfact),fatten = 0, 
                              width = .8, 
                              fill = "white",
                              alpha = .5,
                              colour = "transparent",
                              show.legend = FALSE) +
                geom_hline(yintercept = 0) +
                coord_flip(ylim=c(-2,9)) +
                scale_y_continuous("", breaks = seq(0,10, 1.5))+
                labs(y = "Excess mortality (%)", x = "",
                             title = "Heat-related mortality by scenario")+
                scale_x_reverse() + 
                geom_text( aes(x=id, y=-0.1, label=rownames(res_af), hjust=1),
                            alpha=0.8, size=5, inherit.aes = FALSE) +
                scale_fill_manual(values = col_reg) +
                theme_minimal() +
                theme(legend.position="none",
                              plot.margin=unit(c(2,-1,0,0), "lines"),
                              plot.title = element_textbox(size = 17, hjust = 0.5, face ="bold"),
                              axis.text.y = element_blank(),
                              axis.text.x = element_text(size=12),
                              axis.title.x = element_text(size=12),
                              axis.ticks = element_blank(),
                              panel.grid.major.x = element_line(size=1),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor = element_blank())



# DOT/WISKER PLOT WITH DIFFERENCE
res_afdif <- data.frame(diff=afcity[,11,1,3],diff_low=afcity[,11,2,3], diff_high=afcity[,11,3,3])
res_afdif$id <- c(1:nrow(res_afdif))
label_data <- rownames(res_afdif)

dot_plot <- ggplot(res_afdif, aes(x = id, y = diff)) +
                    geom_hline(yintercept=00, linetype="solid", color = "black") +
                    geom_errorbar(aes(ymin = diff_low, ymax=diff_high, colour = factor(id)),
                                  width = 0.2, size=1, show.legend = FALSE) +
                    geom_point(aes(y=diff,colour=factor(id)), size = 3, shape = 21, 
                               fill = "white",
                               show.legend = FALSE) +
                    labs(y = "", x = "", title = "Human-induced CC heat-related mortality")+
                    coord_flip() +
                    scale_y_continuous("", breaks = seq(-1,5,1))+
                    scale_color_manual(values = col_reg) +
                    scale_x_reverse()+ 
                    theme_minimal() +
                    theme(legend.position="none",
                          plot.margin = unit(c(2, 2, 0, -2), units = "lines"),
                          plot.title = element_textbox(size = 17, hjust = 0.5, face ="bold"),
                          axis.text.x = element_text(size=12), 
                          axis.text.y = element_blank(),
                          panel.grid.major.x = element_line(size=1),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank())



# COMBINE PLOTS 
barplot + dot_plot + plot_layout(widths = c(2.5,2.5), ncol=2) +  
  theme(plot.margin = unit(c(1, 1, 0, 1.5), units = "lines"))

### TO BE INCLUDED: LEGEND, X AXIS LABEL

save.image("Results.RData")

