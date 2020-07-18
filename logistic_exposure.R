##--------------------------SETUP------------------------------------------

##packages
packages <- c("MuMIn", "ggplot2", "plyr", "readxl",
              "AICcmodavg", "viridis", "viridisLite", "Hmisc",
              "rms", "ResourceSelection", "arm", "lme4",
              "regclass")
lapply(packages, require, character.only=TRUE)

##functions
expand.nest.data <- function( nestdata, 
                              FirstFound="FirstFound", 
                              LastChecked="LastChecked",
                              LastPresent="LastPresent",
                              Fate="Fate",
                              AgeDay1="AgeDay1"){
  require(plyr)
  nestdata2 <- plyr::adply(nestdata, 1, function(x){
    # expand for days when nest is known to be alive
    # create data frame because tibbles can cause problems
    x <- as.data.frame(x)
    #browser()
    orig.x <- x
    x <- x[rep(1,x[,LastPresent]-x[,FirstFound]),]  
    #browser()
    if(nrow(x)>0){
      x$Day <- x[,FirstFound][1]:(x[,LastPresent][1]-1)
      x$Exposure<- 1
      x$Survive <- 1
      #browser()
      if(AgeDay1 %in% names(orig.x)) x$NestAge <- orig.x[,AgeDay1] + x$Day -1 
    }
    
    # now for the final record where the nest fails in interval
    # We define the date of the last interval as the midpoint of the interval
    # We define the nestage in the last interval as the age in the midpoint of the interval as well.
    #browser()
    if(orig.x[,LastChecked] > orig.x[,LastPresent]){
      x2 <- orig.x
      x2$Exposure <- orig.x[,LastChecked] - orig.x[,LastPresent]
      x2$Survive  <- 1-orig.x[,Fate]
      x2$Day      <- x2[,LastPresent] + x2$Exposure/2
      if(AgeDay1 %in% names(orig.x)) x2$NestAge  <- orig.x[,AgeDay1] + x2$Day -1
      x <- rbind(x, x2)
      #browser()
    }
    # remove any records with 0 exposure
    x <- x[ !x$Exposure==0,]
    # return the expanded data set
    x
  })
  nestdata2}
# define the modification for the logit() link to account for exposure
logexp <- function(exposure = 1)
{
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  ## FIXME: is there some trick we can play here to allow
  ##   evaluation in the context of the 'data' argument?
  linkinv <- function(eta)  plogis(eta)^exposure
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps,
           exp(eta)/(1+exp(eta))^2)
    ## OR .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
  }
  mu.eta <- function(eta) {       
    exposure * plogis(eta)^(exposure-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

logit <- function(x){log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}
# expand the data to generate the effective sample size
# See the GIM, Chapter 17, page 17-8
# We assume that the variable names on the nest record match those
# required by RMark, but these can be changed
# Returns an expanded dataset with variables Exposure, Survive, Day, NestAge


##--------------------------IMPORT-----------------------------------------

##working directory
setwd("R:/Common Grackles (Noah Horsley)/Analysis/Nest Survival (Logistic Exposure)/")

##import data
nests <- read.csv("Nests.csv", header = T)
View(nests)
nests <- nests[-c(197, 198, 199, 200), ]
nests$X <- NULL


##--------------------------ANALYSIS---------------------------------------

##EVERY RUN
##remove CA nests
nests <- subset(nests, site != "CA")
##remove low confidence nests
nests <- subset(nests, Confidence %in% c("Low", "Moderate", "High"))
##expand data to effective sample size
nests <- expand.nest.data(nests)
##create Age 
nests$Age <- nests$Day - nests$Age_Adj
##store models
candidate.mods <- list()
##drop unused levels
nests <- droplevels(nests)

##SUBSET LEVELS
##create height subset
nests <- subset(nests, !is.na(nests$height))
##create nestlings subset
nests <- subset(nests, nestlings > 0)
##create 2018 subset
nests2018 <- subset(nests, year == 2018)
##create 2019 subset
nests2019 <- subset(nests, year == 2019)


##--------------------------SINGLE VARIABLE MODELS-------------------------

##NULL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR <- glm(Survive~1, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR)
-2*logLik(modDSR) ##DSR on LOGIT scale
##store
candidate.mods[[1]] <- modDSR
##logit(DSR) to DSR
DSR <- expit(coef(modDSR))
DSR.se <- arm::se.coef(modDSR)*DSR*(1-DSR)
cat("DSR", DSR, "(SE", DSR.se, ")\n") ##DSR + SE
expit(confint(modDSR)) ##95% confidence interval
##nest survival
days <- 26
NS <- DSR^days
NS.se <- DSR.se * days * DSR^(days-1)
cat("NS ", days," days ", NS, "(SE ", NS.se, ")\n") ##nest survival rate over nesting period

##SURVIVAL ~ YEAR + SITE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_YearSite <- glm(Survive~year+site, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_YearSite)
-2*logLik(modDSR_YearSite) ##DSR on LOGIT scale

##SURVIVAL ~ SEASONALITY (for PLOT) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_SeasonPlot <- glm(Survive~Day, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_SeasonPlot)
-2*logLik(modDSR_SeasonPlot)

##SURVIVAL ~ SEASONALITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_Season <- glm(Survive~Day+year+site, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_Season)
-2*logLik(modDSR_Season)
##store
candidate.mods[[2]] <- modDSR_Season

##SURVIVAL ~ NEST AGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_Age <- glm(Survive~Age+year+site, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_Age)
-2*logLik(modDSR_Age)
##store
candidate.mods[[3]] <- modDSR_Age

##SURVIVAL ~ CLUTCH SIZE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_Eggs <- glm(Survive~eggs+year+site, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_Eggs)
-2*logLik(modDSR_Eggs)
##store
candidate.mods[[4]] <- modDSR_Eggs


##--------------------------SINGLE VARIABLE AICc---------------------------

SV_AICc <- aictab(candidate.mods, modnames = c("NULL", "CONTROLS", "SEASONALITY", "NEST AGE", "CLUTCH SIZE"), sort = T)
SV_AICc

##export
write.csv(SV_AICc, "AIC Table 10.30.19.csv", row.names = F)


##--------------------------ADDITIVE MODELS--------------------------------

##SURVIVAL ~ SEASONALITY + NEST AGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_SeasonAge <- glm(Survive~Day+Age+year+site, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_SeasonAge)
-2*logLik(modDSR_SeasonAge)
##check for collinearity
cor(nests$Age, nests$Day) ##0.6843-0.7200
##VIF Test
VIF(modDSR_SeasonAge) ##<2 
##store
candidate.mods[[5]] <- modDSR_SeasonAge

##SURVIVAL ~ FULL MODEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modDSR_Full <- glm(Survive~Day+year+site+Age+eggs, family = binomial(link = logexp(nests$Exposure)), data = nests)
summary(modDSR_Full)
-2*logLik(modDSR_Full)
##store
candidate.mods[[6]] <- modDSR_Full

##--------------------------FINAL AICc-------------------------------------

MV_AICc <- aictab(candidate.mods, modnames = c("NULL", "SEASONALITY", "NEST AGE", 
                                               "CLUTCH SIZE", "SEASONALITY + NEST AGE", "SEASONALITY + NEST AGE + CLUTCH SIZE"), sort = T)
MV_AICc

##export
write.csv(MV_AICc, "AIC Table 03.13.2020.csv", row.names = F)

##--------------------------PREDICT----------------------------------------

##seasonality
pred.data <- data.frame(expand.grid(Day = min(nests$FirstFound):max(nests$LastChecked)))
logit.dsr.pred.linear <- predict(modDSR_SeasonPlot, newdata = pred.data, se.fit = TRUE)
dsr.linear <- cbind(pred.data, logit.dsr = logit.dsr.pred.linear$fit, logit.dsr.se = logit.dsr.pred.linear$se.fit)
dsr.linear$dsr <- expit(dsr.linear$logit.dsr)
dsr.linear$dsr.se <- dsr.linear$logit.dsr.se * dsr.linear$dsr * (1 - dsr.linear$dsr)
dsr.linear

##old topmod
pred.data <- data.frame(expand.grid(Day = min(nests2019$FirstFound):max(nests2019$LastChecked), year = 2018:2019, eggs = 1:6))
logit.dsr.pred.linear <- predict(modDSR_Eggs.Season.Year, newdata = pred.data, se.fit=TRUE) 
dsr.linear <- cbind(pred.data, logit.dsr=logit.dsr.pred.linear$fit, logit.dsr.se=logit.dsr.pred.linear$se.fit)
dsr.linear$dsr <- expit(dsr.linear$logit.dsr)
dsr.linear$dsr.se <- dsr.linear$logit.dsr.se* dsr.linear$dsr * (1-dsr.linear$dsr)
dsr.linear
dsr.linear2018 <- subset(dsr.linear, year == 2018)
dsr.linear2019 <- subset(dsr.linear, year == 2019)
dsr.linear1egg <- subset(dsr.linear, eggs == 1)
dsr.linear2egg <- subset(dsr.linear, eggs == 2)
dsr.linear3egg <- subset(dsr.linear, eggs == 3)
dsr.linear4egg <- subset(dsr.linear, eggs == 4)
dsr.linear5egg <- subset(dsr.linear, eggs == 5)
dsr.linear6egg <- subset(dsr.linear, eggs == 6)

##--------------------------PLOTS------------------------------------------

##SURVIVAL ~ SEASONALITY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

COGR.seasonmod <- ggplot(data = dsr.linear, aes(x = Day))+
  geom_line(aes(y = dsr))+
  geom_line(aes(y = expit(logit.dsr-1.96*logit.dsr.se)), linetype = "dashed")+
  geom_line(aes(y = expit(logit.dsr+1.96*logit.dsr.se)), linetype = "dashed")+
  xlim(100, 171)+
  ylim(0.7, 1)+
  ggtitle("Influence of Seasonality on Nest Success", subtitle="Logistic Exposure  |  p < 0.001  |  n = 188")+ 
  xlab("Julian Date")+
  ylab("Daily Survival Rate")+
  theme(axis.title.x=element_text(color="Black", size=12),
        axis.title.y=element_text(color="Black", size=12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(color="Black", size=15, hjust=0.5),
        plot.subtitle=element_text(color="Black", size=10, hjust=0.5),
        strip.text=element_text(color="Black", size=12),
        #plot
        panel.grid.minor.y=element_line(colour="grey94"),
        panel.grid.minor.x=element_line(colour="grey94"),
        strip.background=element_rect(color="black", fill="white"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(linetype="solid", fill=NA))
COGR.seasonmod

##SURVIVAL ~ SEASON + YEAR + EGGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
COGR.topmod <- ggplot()+
  ##1egg
  geom_line(data = dsr.linear1egg, aes(x = Day, y = dsr, color = "1"), size = 1)+
  geom_ribbon(data = dsr.linear1egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "1"), alpha=.2)+
  ##2eggs
  geom_line(data = dsr.linear2egg, aes(x = Day, y = dsr, color = "2"), size = 1)+
  geom_ribbon(data = dsr.linear2egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "2"), alpha=.2)+
  ##3eggs
  geom_line(data = dsr.linear3egg, aes(x = Day, y = dsr, color = "3"), size = 1)+
  geom_ribbon(data = dsr.linear3egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "3"), alpha=.2)+
  ##4eggs
  geom_line(data = dsr.linear4egg, aes(x = Day, y = dsr, color = "4"), size = 1)+
  geom_ribbon(data = dsr.linear4egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "4"), alpha=.2)+
  ##5eggs
  geom_line(data = dsr.linear5egg, aes(x = Day, y = dsr, color = "5"), size = 1)+
  geom_ribbon(data = dsr.linear5egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "5"), alpha=.2)+
  ##6eggs
  geom_line(data = dsr.linear6egg, aes(x = Day, y = dsr, color = "6"), size = 1)+
  geom_ribbon(data = dsr.linear6egg, aes(x = Day, ymin=expit(logit.dsr-1.96*logit.dsr.se),
                                         ymax=expit(logit.dsr+1.96*logit.dsr.se), fill = "6"), alpha=.2)+
  xlim(100, 160)+
  ylim(0.5, 1)+
  ggtitle("Common Grackle Nest Survival in Central Illinois", subtitle="Logistic Exposure  |  n = 188")+ 
  xlab("Ordinal Date")+
  ylab("Daily Survival Rate")+
  facet_grid(.~year)+
  #legend
  #scale_color_manual(name="Eggs", values=c("1"="grey85", "2"="grey70", "3"="grey55", 
  #                                                   "4"="grey40", "5"="grey25", "6"="black"))+
  #scale_fill_manual(name="Eggs", values=c("1"="grey85", "2"="grey70", "3"="grey55", 
  #                                         "4"="grey40", "5"="grey25", "6"="black"))+
  scale_color_viridis(name="Clutch Size", option = "viridis", discrete=T, direction=-1)+
  scale_fill_viridis(name="SE", labels=c("0.0273", "0.0174", "0.0109", "0.0072", "0.0057", "0.0051"), option="viridis", discrete=T, direction=-1)+
  guides(color=guide_legend(ncol=1, order=1), fill=guide_legend(ncol=1, order=2))+
  #text
  theme(axis.title.x=element_text(color="Black", size=15),
        axis.title.y=element_text(color="Black", size=15),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(color="Black", size=20, hjust=0.5),
        plot.subtitle=element_text(color="Black", size=12, hjust=0.5),
        strip.text=element_text(color="Black", size=15),
        #plot
        panel.grid.minor.y=element_line(colour="grey94"),
        panel.grid.minor.x=element_line(colour="grey94"),
        strip.background=element_rect(color="black", fill="white"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(linetype="solid", fill=NA),
        #legend
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.position="right",
        legend.direction="vertical",
        legend.key=element_rect(fill="white"),
        legend.box.background=element_rect(),
        legend.box.margin=margin(6,6,6,6))
COGR.topmod