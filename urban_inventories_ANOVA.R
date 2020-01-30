# -------------------------------------------------------#
#  Analysis of the uncertainties in emission inventories #
# ------------ with 1-way ANOVA with blocking -----------#
# -------------------------------------------------------#

# This is the code related to the publication 
# "A statistical model for the uncertainties in emission inventories"

# Authors: Degraeuwe, B. and Peduzzi, E.

# contact: bart.degraeuwe@ec.europa.eu

# The code is available on github 
# Create a local git folder folder and execute in the git Shell:
# $ git clone https://github.com/bd77/InventoryUncertainty.git

# This script reads data for 11 cities, 3 sectors, 6 inventories and 4 pollutants.
# A one-way ANOVA is performed on each city-sector combination. For each inventory
# pair the log-activity and log-emission factor differences are calculated with
# their confidence interval.

# input file: 'inventories.csv'
# 5 columns:
#   - city: Barcelona, Bucarest, Budapest, Katowice, London, Madrid, Milan, Paris, Sofia, Utrecht, Warsaw
#   - sector: ms34 (industry), ms2 (residential), ms7 (road transport)
#   - inventory: ctm4iam, edgar, emep, jrc, macc2, macc3
#   - pollutant: nox, pm25, voc, sox
#   - emission (kton)

# output to determine the sources of uncertainty:
# - Figures
# <city>_<sector>_fig1_pol_logE.png: emissions versus pollutants
# <city>_<sector>_fig2_inv_logE.png: emissions versus inventories
# <city>_<sector>_fig3_logA_ratios.png: log-activity differences between inventory pairs 
# <city>_<sector>_fig4_A_pctdiff.png: percentage activity differences between inventory pairs
# <city>_<sector>_fig5_logEF_<pollutant>_ratios.png: log-EF differences between inventory pairs 
# <city>_<sector>_fig6_logEF_<pollutant>_pctdiff.png: percentage EF differences between inventory pairs
# <city>_<sector>_fig7_residuals_pollutant: residual plot to check for heteroskedasticity
# <city>_<sector>_fig8_residuals_inventory: residual plot to check for heteroskedasticity

# - results table: "inventories_ANOVA_results.csv"

# load libraries
library(plyr)
library(ggplot2)

# clean up
rm(list = ls())
# set the working directory
wd <- "D:/Diamondplot/github/InventoryUncertainty"
setwd(wd)

# read the input data
inv.df <- read.table('inventories.csv', header = T, sep = ",")
city.list <- as.vector(unique(inv.df$city))
sector.list <- as.vector(unique(inv.df$sector))

# plot settings
plot.font.size <- 20

# data frame with all the results
res.df <- data.frame()
ANOVA.results.path <- "ANOVA_results"
if (!dir.exists(ANOVA.results.path)) {dir.create(ANOVA.results.path)}

# loop over all city - sector combinations
for (city in city.list) {
  for (sector in sector.list) {
    print(paste("ANOVA on emissions of", sector, "in", city))
    
    # create an empty data.frame to store results of a city-sector combination
    sc.res.df <- data.frame()
    
    # select just the data of one city-sector combination
    inv.city.sector.df <- inv.df[inv.df$city==city & inv.df$sector==sector,]
    inventory.list <- as.vector(unique(inv.city.sector.df$inventory))
    pollutant.list <- as.vector(unique(inv.city.sector.df$pollutant))
    
    # add a column with log emissions
    inv.city.sector.df$logE <- log(inv.city.sector.df$emission)

    # Line plot of log emission as a function of pollutant. This plot allows a visual
    # diagnostic of the main source of uncertainty (parrallel lines => activity,
    # big spread for one pollutant => emission factor)
    p <- ggplot(inv.city.sector.df, aes(x=pollutant, y=logE, group=inventory, col=inventory)) 
    p <- p + geom_point() + geom_line() 
    p <- p + theme(text = element_text(size=plot.font.size))
    p <- p + labs(title = paste0(city, " ", sector), x="Pollutant", y="log(emission)")
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig1_pol_logE.png")))
    print(p)
    dev.off()
    
    # plot the log-emisssions per inventory
    p <- ggplot(inv.city.sector.df, aes(x=inventory, y=logE, group=pollutant, col=pollutant)) 
    p <- p + geom_line() + geom_point()
    p <- p + labs(title = paste0(city, " ", sector), x="Inventory", y="log(emission)") 
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig2_inv_logE.png")))
    print(p)
    dev.off()
    
    
    # Implementation of ANOVA formulas
    # --------------------------------
    
    # sums of squares
    # Treatment sum of squares, here inventory sum of squares
    NP <- length(unique(inv.city.sector.df$pollutant)) # number of pollutants
    NI <- length(unique(inv.city.sector.df$inventory)) # # number of inventories
    N <- NI*NP
    # total sum of squares
    SS.T <- sum(inv.city.sector.df$logE^2) - (sum(inv.city.sector.df$logE)^2)/N
    # Treatment/Inventory sum of squares
    sum.logE.inv <- ddply(inv.city.sector.df, c("inventory"), summarise, sum.logE.inv = sum(logE))
    SS.I <- 1/NP * sum(sum.logE.inv$sum.logE.inv^2) - (sum(inv.city.sector.df$logE)^2)/N
    MSS.I <- SS.I / (NI-1)
    # block/pollutant sum of squares
    sum.logE.pol <- ddply(inv.city.sector.df, c("pollutant"), summarise, sum.logE.pol = sum(logE))
    SS.P <- 1/NI * sum(sum.logE.pol$sum.logE.pol^2) - (sum(inv.city.sector.df$logE)^2)/N
    # error sum of squares
    SS.E <- SS.T - SS.I - SS.P
    MSS.E <- SS.E / ((NI-1)*(NP-1))
    
    # least significant difference or half the width of the CI of a difference
    # with 0.05 >> 90% CI (one sided 95% or 2-sided 90%)
    # with 0.025 >> 95% CI two sided
    LSD <- qt(0.025, (NI-1)*(NP-1), lower.tail = FALSE) * sqrt(2*MSS.E/NP)
    
    # Calculate the expected value and confidence intervales for the log-activity difference
    # of each inventory pair.
    # The expected value of the activity (+ mean emission factor) of an inventory is the mean of all 
    # emissions (see paper paragraph 2.2).
    mean.logE.inv <- ddply(inv.city.sector.df, c("inventory"), summarise, mean.logE.inv = mean(logE))
    # a data.frame with the results for all inventory pairs
    logA.res.df <- data.frame()
    # loop over all inventory pairs
    for (i1 in 1:(NI-1)) {
      for (i2 in (i1+1):NI) {
        inv1 <- inventory.list[i1]
        inv2 <- inventory.list[i2]
        lr.label <- paste0("Activity_ratio_", inv1, "_", inv2)
        logA1 <- mean.logE.inv$mean.logE.inv[mean.logE.inv$inventory == inv1] 
        logA2 <- mean.logE.inv$mean.logE.inv[mean.logE.inv$inventory == inv2]
        lrA1A2 <- logA1 - logA2
        if (lrA1A2 < 0) {
          inv1 <- inventory.list[i2]
          inv2 <- inventory.list[i1]
          lrA1A2 <- -lrA1A2
        }
        p.value.lrA1A2 <- pt(q = lrA1A2/sqrt(2*MSS.E/NP), df = (NI-1)*(NP-1), lower.tail = FALSE)
        # table with the results for activities
        logA.res.df <- rbind(logA.res.df,
                             data.frame(city = city,
                                        sector = sector,
                                        pollutant = "",
                                        outcome = "logA.ratio",
                                        inventory1 = inv1,
                                        inventory2 = inv2,
                                        EV.log.ratio = lrA1A2,
                                        CI.low = lrA1A2 - LSD,
                                        CI.high = lrA1A2 + LSD,
                                        EV.pct.diff = 100*(exp(lrA1A2) - 1),
                                        CI.low.pct.diff = 100*(exp(lrA1A2 - LSD) - 1),
                                        CI.high.pct.diff = 100*(exp(lrA1A2 + LSD) - 1),
                                        prob.below.zero = p.value.lrA1A2,
                                        significant = if(p.value.lrA1A2 < 0.05) {"1"} else {"0"}))
      }
    }
    # add the log-activity results to the city-sector results data.frame
    sc.res.df <- rbind(sc.res.df,logA.res.df)
    
    # error bar plot for log-activity ratios of each inventory pair
    logA.plot.df <- logA.res.df
    logA.plot.df$inventories <- paste0(logA.plot.df$inventory1, " - ", logA.res.df$inventory2)
    logA.plot.df <- logA.plot.df[order(logA.plot.df$EV.log.ratio),]
    logA.plot.df$inventories <- factor(logA.plot.df$inventories, ordered = T, levels = logA.plot.df$inventories)
    p <- ggplot(data = logA.plot.df, aes(x=EV.log.ratio, y = inventories)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.low, xmax = CI.high))
    p <- p + geom_vline(xintercept = 0, col = "red")
    p <- p + labs(title = paste(city, sector, "\nlog(activity ratio) with 2-sided 95% CI"), 
                  x = "log(A1/A2)", y="Inventory pair")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig3_logA_ratios.png")))
    print(p)
    dev.off()
    
    # percent difference plot for activity ratios of each inventory pair
    p <- ggplot(data = logA.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point() 
    p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
    p <- p + geom_vline(xintercept = 0, col = "red")
    p <- p + labs(title = paste(city, sector, "\nRelative activity diff. with 2-sided 95% CI"), 
                  x = "A1/A2 - 1 (%)", y="Inventory pair")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig4_A_pctdiff.png")))
    print(p)
    dev.off()
    
    # Calculate the expected value and confidence intervals for the log-emission factor difference
    # and relative emission factor differences of each inventory pair.
    inv.city.sector.df <- merge(inv.city.sector.df[,c('city', 'inventory', 'sector', 'emission', 'pollutant', 'logE')], 
                                mean.logE.inv, c("inventory"))
    inv.city.sector.df$logEF <- inv.city.sector.df$logE - inv.city.sector.df$mean.logE.inv
    
    # Log emission factor differences
    logEF.res.df <- data.frame()
    for (p1 in 1:NP) {
      pol.name <- pollutant.list[p1]
      for (i1 in 1:(NI-1)) {
        for (i2 in (i1+1):NI) {
          inv1 <- inventory.list[i1]
          inv2 <- inventory.list[i2]
          # log EF ratio for pollutant p1 between inv1 and inv2
          logEF.p1.inv1 <- inv.city.sector.df$logEF[inv.city.sector.df$inventory == inv1 & 
                                                      inv.city.sector.df$pollutant == pol.name]
          logEF.p1.inv2 <- inv.city.sector.df$logEF[inv.city.sector.df$inventory == inv2 & 
                                                      inv.city.sector.df$pollutant == pol.name]
          lrEF1EF2 <- logEF.p1.inv1 - logEF.p1.inv2
          if (lrEF1EF2 < 0) {
            inv1 <- inventory.list[i2]
            inv2 <- inventory.list[i1]
            lrEF1EF2 <- -lrEF1EF2
          }
          # percent difference of all samples
          pct.diff.EF1EF2 <- 100*(exp(lrEF1EF2)-1)
          lr.label <- paste0(pol.name, "_EF_ratio_", inv1, "_", inv2)
          # p-value, probablity that log ratio < 0
          p.value.EF1EF2 <- pt(q = lrEF1EF2/sqrt(2*MSS.E/NP), df = (NI-1)*(NP-1), lower.tail = FALSE)
          # table with the results for activities
          logEF.res.df <- rbind(logEF.res.df,
                                data.frame(city = city,
                                     sector = sector,
                                     outcome = "logEF.ratio",
                                     pollutant = pol.name,
                                     inventory1 = inv1,
                                     inventory2 = inv2,
                                     EV.log.ratio = lrEF1EF2,
                                     CI.low = lrEF1EF2 - LSD,
                                     CI.high = lrEF1EF2 + LSD,
                                     EV.pct.diff = pct.diff.EF1EF2,
                                     CI.low.pct.diff = 100*(exp(lrEF1EF2 - LSD)-1),
                                     CI.high.pct.diff = 100*(exp(lrEF1EF2 + LSD)-1),
                                     prob.below.zero = p.value.EF1EF2,
                                     significant = if(p.value.EF1EF2 < 0.05) {"1"} else {"0"}))
        }
      }
    }
    sc.res.df <- rbind(sc.res.df,logEF.res.df)
    
    # write a file with the results for the city sector
    write.table(sc.res.df, row.names = F, sep = ",", 
                file = file.path(ANOVA.results.path,
                                 paste0("ratio_results_ANOVA_", city, "_", sector, ".csv")))
    # add the results for the city-sector to all results
    res.df <- rbind(res.df, sc.res.df)
    
    # error bar plot for log-emission factor ratios and relative emission factor differences
    for (i.p in 1:NP) { # loop over all pollutants
      # select data for the pollutant
      logEF.plot.df <- logEF.res.df[logEF.res.df$pollutant == pollutant.list[i.p],]
      # add labels for the inventory pairs
      logEF.plot.df$inventories <- paste0(logEF.plot.df$inventory1, " - ", logEF.plot.df$inventory2)
      # sort from large to small differences
      logEF.plot.df <- logEF.plot.df[order(logEF.plot.df$EV.log.ratio),]
      logEF.plot.df$inventories <- factor(logEF.plot.df$inventories, ordered = T, levels = logEF.plot.df$inventories)
      # plot for log EF ratio
      p <- ggplot(data = logEF.plot.df, aes(x=EV.log.ratio, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low, xmax = CI.high))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutant.list[i.p]), 
                                  "\nlog(EF ratios) with 2-sided 95% CI"), 
                    x = "log(EF1/EF2)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(ANOVA.results.path, 
                    paste0(city, "_", sector, "_fig5", letters[i.p], "_logEF_", pollutant.list[i.p], "_ratios.png")))
      print(p)
      dev.off()
      
      # plot for percentage EF difference
      p <- ggplot(data = logEF.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutant.list[i.p]), 
                                  "\nRelative EF diff. with 2-sided 95% CI"), 
                    x = "EF1/EF2 - 1 (%)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(ANOVA.results.path, 
                    paste0(city, "_", sector, "_fig6", letters[i.p], "_EF_", pollutant.list[i.p], "_pctdiff.png")))
      print(p)
      dev.off()
    }
    
    # residual plots per pollutant and inventory
    # calculate the residuals
    mean.logEF.pol <- ddply(inv.city.sector.df, c("pollutant"), summarise, mean.logEF.pol = mean(logE))
    mean.logE <- mean(inv.city.sector.df$logE)
    residuals.df <- inv.city.sector.df[,c('inventory', 'pollutant', 'logE')]
    residuals.df <- merge(residuals.df, mean.logE.inv)
    residuals.df <- merge(residuals.df, mean.logEF.pol)
    residuals.df$residuals <- residuals.df$logE - residuals.df$mean.logE.inv - residuals.df$mean.logEF.pol + mean.logE
    
    # residuals per pollutant
    p <- ggplot(data = residuals.df, aes(x=pollutant, y=residuals, group=inventory, col=inventory)) 
    p <- p + geom_point() + geom_line()
    p <- p + labs(title = paste(city, sector, "Residuals"))
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig7_residuals_pollutant.png")))
    print(p)
    dev.off()

    # residuals per inventory
    p <- ggplot(data = residuals.df, aes(x=inventory, y=residuals, group=pollutant, col=pollutant)) 
    p <- p + geom_point() + geom_line()
    p <- p + labs(title = paste(city, sector, "Residuals"))
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(ANOVA.results.path, paste0(city, "_", sector, "_fig8_residuals_inventory.png")))
    print(p)
    dev.off()
    
    # Control with the R lm function (gives the same results)
    # ------------------------------
    
    # ANOVA with fixed effect inventory and blocking factor pollutant
    # The order in the formula is important. The confidence intervals given by
    # confint give the CI on the difference with the reference inventory (the first one, ctm4iam)
    # lm.inv <- lm(data = inv.city.sector.df, 
    #              formula = logE ~ inventory + pollutant )
    # summary(lm.inv)  
    # 
    # # compare the residuals
    # resid.lm.df <- data.frame(inv.city.sector.df[,c('inventory', 'pollutant')],
    #                           residuals.lm = residuals.lm <- as.data.frame(proj(lm.inv))[,4])
    # resid.lm.df <- merge(residuals.df, resid.lm.df)
    
    
  } # close the loop over all sectors of a city
} # close the loop over all cities

# write all results to a file
write.table(res.df, row.names = F, sep = ",", 
            file = file.path("inventories_ANOVA_results.csv"))


