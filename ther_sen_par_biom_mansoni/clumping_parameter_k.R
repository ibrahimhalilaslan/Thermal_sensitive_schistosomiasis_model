
#################################################################################
# BEGINNING of Giulio's latest integration on Feb 20 2023
 
# compute clumping parameters
 
#create a dataframe to store the results
cp_df <- worm_data[, c("School", "year")]
cp_df$Sh_Amean <- NA
cp_df$Sm_Amean <- NA
cp_df$Sh_cp    <- NA
cp_df$Sm_cp    <- NA
 
#create a fucntion for the MLE estimator
fitNegBin <- function(SPar){
 
  mu <- mean(SPar, na.rm=T)
  vari <- var(SPar, na.rm=T)
 
  betfitSh <-  tryCatch(fitdistr(round(na.exclude(SPar)), densfun = dnbinom,
                           start = list(mu = mu, size = mu^2 /(vari-mu)/10),
                           method = "L-BFGS-B",
                           lower = c(1, 0.00001),
                           upper = c(2*mean(SPar, na.rm=T), 10)),
                           error=function(err) NA)
 
  if (length(betfitSh) == 1){
     betfitSh <-  tryCatch(fitdistr(round(na.exclude(SPar)), densfun = dnbinom,
                                   start = list(mu = mu, size = mu^2 /(vari-mu)/5),
                                   method = "L-BFGS-B",
                                   lower = c(1, 0.00001),
                                   upper = c(2*mean(SPar, na.rm=T), 10)),
                          error=function(err) NA)
  }
 
  
    return(ifelse(is.na(betfitSh), NA, betfitSh$estimate[2]))
}
 
for(i in 1:nrow(cp_df)){
    ndf <- filter(human_data, School == cp_df[i,]$School, year ==  cp_df[i,]$year)
 
#   ndf <- filter(human_data, School == si, year ==  yi)
   
    SPar <- ndf$SmW
   
    cp_df[i, 'Sh_Amean'] <- mean(ndf$ShW, na.rm=T)
    cp_df[i, 'Sm_Amean'] <- mean(ndf$SmW, na.rm=T)
    cp_df[i, 'Sh_cp'] <- ifelse(round(mean(ndf$ShW, na.rm=T))>0, fitNegBin(ndf$ShW), NA) 
    cp_df[i, 'Sm_cp'] <- ifelse(round(mean(ndf$SmW, na.rm=T))>0, fitNegBin(ndf$SmW), NA) 
    
  }
 
plot(cp_df$Sh_Amean, cp_df$Sh_cp, log = 'xy')
 
plot(cp_df$Sm_Amean, cp_df$Sm_cp, log = 'xy')
 
hist(cp_df$Sh_Amean)
hist(cp_df$Sm_Amean)
hist(cp_df$Sh_cp)
hist(cp_df$Sm_cp)
 
ggplot(cp_df, aes(x = Sh_Amean)) +
  geom_histogram(binwidth=15, color="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sh_Amean)),color="blue", linetype="dashed", size=1) +
  xlab("mean parasite burden") + ylab("count") +
  ggtitle("S. haematobium")
 
ggplot(na.omit(cp_df, col = "Sm_cp"), aes(x = Sm_Amean)) +
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sm_Amean)),color="blue", linetype="dashed", size=1) +
  xlab("mean parasite burden") + ylab("count") +
  ggtitle("S. mansoni") +
  scale_x_continuous(trans='log10')
 
  ggplot(cp_df, aes(x = Sh_cp)) +
    geom_histogram(color="black", fill="white") +
    geom_vline(aes(xintercept=mean(Sh_cp)),color="blue", linetype="dashed", size=1) +
    xlab("clumping paramater") + ylab("count") +
    ggtitle("S. haematobium")
 
  ggplot(na.omit(cp_df, col = "Sm_cp"), aes(x = Sm_cp)) +
    geom_histogram(color="black", fill="white") +
    geom_vline(aes(xintercept=mean(Sm_cp)),color="blue", linetype="dashed", size=1) +
    xlab("clumping parameter") + ylab("count") +
    ggtitle("S. mansoni")
 
  
  summary(cp_df$Sh_Amean)
  summary(cp_df$Sm_Amean)
 
# linear regression for the two parasites independently, not ested in years
 
regSh <- lm(Sh_cp ~ Sh_Amean, data = cp_df)
summary(regSh)
 
regSm <- lm(Sm_cp ~ Sm_Amean, data = cp_df)
summary(regSm)
 
# log-log regression
 
regSh <- lm(log(Sh_cp) ~ log(Sh_Amean), data = cp_df) #filter(cp_df, year == 2018))
summary(regSh)
 
regSm <- lm(log(Sm_cp) ~ log(Sm_Amean), data = cp_df)
summary(regSm)
 
 
 
ggplot(data = cp_df,aes(x = log(Sh_cp), y = log(Sh_Amean), color = year, group = 1)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab("log(parasites)") + ylab("log(clumping parameter 'k')") +
  ggtitle("S. haematobium")
 
ggplot(data = cp_df,aes(x = log(Sm_cp), y = log(Sm_Amean), color = year, group = 1)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab("log(parasites)") + ylab("log(clumping parameter 'k')") +
  ggtitle("S. mansoni")
 
 
 
 
# note that the results are almost identical
# lets us lmer4
 
library('lme4')
 
cp_df$year <- as.factor(cp_df$year)
 
regSh <- lmer(log(Sh_cp) ~ log(Sh_Amean) + (1|year), data = cp_df )
summary(regSh)
 
regSm <- lmer(log(Sm_cp) ~ log(Sm_Amean) + (1|year), data = cp_df )
summary(regSm)
##############################################################################
# END of Giulio latest integration on Feb 20 2023
##############################################################################