# -------------------------------------------------------#
#  Analysis of the uncertainties in emission inventories #
# ------------ with Markov Chain Monte Carlo ------------#
# -------------------------------------------------------#

# This is the code related to the publication 
# "A statistical model for the uncertainties in emission inventories"

# Authors: Degraeuwe, B. and Peduzzi, E.

# contact: bart.degraeuwe@ec.europa.eu

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
# <city>_<sector>_fig2a_StdDev.png: standard deviations of activity and emission factors with CI
# <city>_<sector>_fig2b_StdDev.png: coefficient of variation of activity and emission factors with CI
# <city>_<sector>_fig3_logA_distributions.png: PDF of the activities of the inventories
# <city>_<sector>_fig4_logA_ratios.png: log-activity differences between inventory pairs
# <city>_<sector>_fig5_A_pctdiff.png: percentage activity differences between inventory pairs
# <city>_<sector>_fig6a_logEF_<pollutant>_distributions.png: PDF of the EFs of the inventories
# <city>_<sector>_fig7a_logEF_<pollutant>_ratios.png: log-EF differences between inventory pairs
# <city>_<sector>_fig8a_EF_<pollutant>_pctdiff.png: percentage EF differences between inventory pairs

# - results tables: 
# write the results to a file:
# ratio_results_MCMC.csv: The expected values and confidence intervals of activity and emission factor differences.
# stddev_results_MCMC.csv: The expected values of standard deviations and coefficients of variation an 
#   their confidence intervals
# stddev_compare_results_MCMC.csv: pairwise comparison of activity/emission factor SDs and credibility 
#   that one is bigger than the other based on MCMC samples.



# load libraries
library(nlme)
library(plyr)
library(R2jags)
library(ggplot2)

# clean up
rm(list = ls())
# set the working directory
wd <- "C:/Documenten/InventoryUncertainty"
setwd(wd)

# create a folder for the results
MCMC.results.path <- "MCMC_results"
if (!dir.exists(MCMC.results.path)) {
  dir.create(MCMC.results.path)
}

# credibility interval for log-activity and emission factor differences
cred.level <- 0.95
# plot font size
plot.font.size <- 15
# make trace plots
trace.plots <- TRUE # takes a lot of time
# set the seed for repeatable results
set.seed(666)

# read the input data
inv.df <- read.table('inventories.csv', header = T, sep = ",")
city.list <- as.vector(unique(inv.df$city))
sector.list <- as.vector(unique(inv.df$sector))

# Data exploration
# ----------------

# An overview of the coefficients of variation (relative standard deviations) of the emissions 
# and the standard deviations of the emissions
# These should be related by sd/mean = sqrt(exp(var(log(E)))-1) if E is log-nomrally distributed.
sd.city.sector.pol.overview <- ddply(inv.df, c('city', 'sector', 'pollutant'), summarise, 
                                     rel.sd.E = sd(emission)/mean(emission), 
                                     sd.log.E = sqrt(exp(var(log(emission)))-1))
# some outliers here but the assumption of a log-normal distribution looks ok
p <- ggplot(data = sd.city.sector.pol.overview, aes(x=rel.sd.E, y=sd.log.E)) + geom_point()
p <- p + geom_abline(intercept = 0, slope = 1, col="red")
p <- p + labs(x="Relative sd of the emissions", y="sd of log-emissions")
png(file.path(MCMC.results.path, "1_relative_sd_VS_sd_of_logE.png"))
print(p)
dev.off()
# The sd of the log VOC emissions in Katowice for ms34 are an extreme outlier (10.0)
# Normal values are around 0.5

# boxplots of the standard deviation of log-emissions
p <- ggplot(sd.city.sector.pol.overview, aes(x=pollutant, y=sd.log.E)) + geom_boxplot()
p <- p + lims(y=c(0,2.5)) # There is one big outlier for VOC in Katowice
png(file.path(MCMC.results.path, "2_Boxplot_sd_logE_per_pollutant.png"))
print(p)
dev.off()

# quantile plot: log-emissions verus normal distribution
png(file.path(MCMC.results.path, "3_Normal_quantile_plot_of_logE.png"))
qqplot(rnorm(10*NROW(inv.df)), log(inv.df$emission),
       main = "QQ-plot of log-emissions",
       xlab="Quantiles of a normal distribution",
       ylab="Quantiles of log(emissions)")
dev.off()

# quantile plots per city-sector combination
qqplot.folder <- file.path(MCMC.results.path, "quantile_plots")
if (!dir.exists(qqplot.folder)) {
  dir.create(qqplot.folder)
}
for (city in city.list) {
  for (sector in sector.list) {
    png(file.path(qqplot.folder, 
                  paste0("Normal_quantile_plot_of_logE_for", city, "_", sector, ".png")))
    qqplot(rnorm(10*NROW(inv.df)), log(inv.df$emission[inv.df$city==city & inv.df$sector==sector]),
           xlab="Quantiles of a normal distribution",
           ylab="Quantiles of log(emissions)",
           main = paste("QQ-plot for the log-emissions of", sector, "in", city))
    dev.off()
  }
}

# Markov Chain Monte Carlo
# ------------------------

# The stochastic model used in JAGS
# The response Y are the log-emissions
# Each pollutant has its own emission factor standard deviation (sd.logEF[j])
# The log-activities (logA[i]) of the inventories are normally distributed with
# standard deviation sd.logA.
mcmc_model_string <- "model{
  # Likelihood
  for(i in 1:NI){
    for (j in 1:NP) {
      Y[i,j] ~ dnorm(logA[i] + mu.logEF[j], tau.logEF[j])
    }
  }
  # calculate logEF[i,j], i.e. the residuals
  for (i in 1:NI) {
    for (j in 1:NP) {
      logEF[i,j] <- Y[i,j] - logA[i] 
    }
  }
  # random effects for logA
  for (i in 1:NI) {
    logA[i] ~ dnorm(0, tau.logA)
  }
  
  # Priors
  sd.logA <- (tau.logA)^(-0.5)
  tau.logA ~ dgamma(0.0001, 0.0001)
  for (j in 1:NP) {
    mu.logEF[j] ~ dnorm(0, 0.00001)
    tau.logEF[j] ~ dgamma(0.0001, 0.0001)
    sd.logEF[j] <- (tau.logEF[j])^(-0.5)
  }
}"

# data frame with all the results for activities and emission factors
res.df <- data.frame()
# data frame with all the results for standard deviations
sd.res.df <- data.frame()
# comparison of standard deviations
sd.compare.df <- data.frame()

# loop over cities and sectors
for (city in city.list[1]) {
  for (sector in sector.list[3]) {
    print(paste("MCMC on emissions of", sector, "in", city))
    
    # data.frame for the results of a city-sector combination.
    sc.res.df <- data.frame()
    
    # select just the data of one city-sector combination
    inv.city.sector.df <- inv.df[inv.df$city==city & inv.df$sector==sector,]
    # add a column with log emissions
    inv.city.sector.df$logE <- log(inv.city.sector.df$emission)

    # output folder for the city-sector
    cs.output.folder <- file.path(MCMC.results.path, paste(city, sector, sep = "_"))
    if (!(dir.exists(cs.output.folder))) {
      dir.create(cs.output.folder)
    }
    
    pollutant.list <- as.vector(unique(inv.city.sector.df$pollutant))
    NP <- length(pollutant.list)
    inventory.list <- as.vector(unique(inv.city.sector.df$inventory))
    NI <- length(inventory.list)
      
    # log emissions as a matrix
    Y <- matrix(0, nrow = NI, ncol = NP, dimnames = list(inventory.list, pollutant.list))
    for (i in 1:NI) {
      for (j in 1:NP) {
        Y[i,j] <- inv.city.sector.df$logE[inv.city.sector.df$inventory==inventory.list[i] &
                                            inv.city.sector.df$pollutant == pollutant.list[j]]
      }
    }
    
    # Line plot of log emission as a function of pollutant. This plot allows a visual
    # diagnostic of the main source of uncertainty (parrallel lines => activity,
    # big spread for one pollutant => emission factor)
    p <- ggplot(data = inv.city.sector.df, 
                aes(x = pollutant, y = logE, group = inventory, col = inventory)) 
    p <- p + geom_line() + geom_point()
    p <- p + labs(title = paste("log-emissions of", sector, "in", city))
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig1_pol_logE.png")))
    print(p)
    dev.off()
    
    # -------- JAGS ----------
    # JAGS simulation
    n.chains <- 3
    jags.input.list <- list(Y = Y, NI = NI, NP = NP)
    jagsfit <- jags(data = jags.input.list, n.iter=10000, n.chains=n.chains,
                    parameters.to.save = c('logA', 'logEF', 'mu.logEF', 'sd.logEF', 'sd.logA'),
                    model.file = textConnection(mcmc_model_string))
    
    # Reorganize the JAGS output
    # --------------------------
    jagsfit.mcmc <- as.mcmc(jagsfit)
    output.df <- data.frame()
    for (i in 1:n.chains) {
      chain.df <- cbind(chain = toString(i), iteration = 1:NROW(jagsfit.mcmc[[i]]), 
                        as.data.frame(jagsfit.mcmc[[i]]))
      output.df <- rbind(output.df, chain.df)
    }
    
    # rename columns of output.df and remove brackets
    # "B[2,2]" becomes B2_2
    output.col.names <- names(output.df)
    n.col.output <- NCOL(output.df)
    for (i.col in 1:n.col.output) {
      col.name <- output.col.names[i.col]
      new.col.name <- gsub(',', '_', col.name)
      new.col.name <- gsub('\\[', '', new.col.name)
      new.col.name <- gsub('\\]', '', new.col.name)
      names(output.df)[i.col] <- new.col.name
    }
    
    # Trace plots and PDF per chain
    # -----------------------------
    
    if (trace.plots == TRUE) {
      for (varname in names(output.df)[3:NCOL(output.df)]) {
        # from varname to real names of inventory and pollutant
        if (grepl(pattern = "^(logA)", varname)) {
          inv.id <- as.numeric(gsub(pattern = "logA", "", varname))
          var.label <- paste0('logA_', inventory.list[inv.id])
        } else if (grepl(pattern = "^(logEF)", varname)) {
          inv.pol.id <- gsub(pattern = "logEF", "", varname)
          inv.id <- as.numeric(strsplit(inv.pol.id, "_")[[1]][1])
          pol.id <- as.numeric(strsplit(inv.pol.id, "_")[[1]][2])
          var.label <- paste0('logEF_', inventory.list[inv.id], '_', pollutant.list[pol.id])
        } else if (grepl(pattern = "^(sd.logEF)", varname)) {
          pol.id <- as.numeric(gsub(pattern = "sd.logEF", "", varname))
          var.label <- paste0('sd.logEF_', pollutant.list[pol.id])
        } else {
          var.label <- varname
        }
        
        p <- ggplot(data = output.df,
                    aes_string(x='iteration', y=varname, col='chain')) + geom_line()
        p <- p + labs(title = paste("Traceplot", var.label))
        p <- p + theme(text = element_text(size=plot.font.size))
        png(file.path(cs.output.folder, paste0(city, "_", sector, "_trace_", var.label, ".png")),
            width = 3*480, height = 0.5*480)
        print(p)
        dev.off()
        p <- ggplot(data = output.df,
                    aes_string(x=varname, col='chain')) + geom_density()
        p <- p + labs(title = paste("Probability density", var.label))
        p <- p + theme(text = element_text(size=plot.font.size))
        png(file.path(cs.output.folder, paste0(city, "_", sector, "_pdf_", var.label, ".png")),
            width = 480, height = 480)
        print(p)
        dev.off()
      }
    }
    
    # Results for standard deviations (SD) and coefficient of variation (CV)
    # ----------------------------------------------------------------------
    
    # SD and CV of the activity
    # -------------------------
    sc.sd.res.df<- data.frame()
    sd.A.output <- output.df$sd.logA
    # coefficient of variation (%) or relative standard deviation (%)
    CV.A.output <- 100*sqrt(exp(sd.A.output^2)-1) 
    sc.sd.res.df <- data.frame(city = city, sector = sector,
                               outcome = "sigma.A", pollutant = "", label = "Activity", 
                               output.name = "sd.logA",
                               EV.sd = mean(sd.A.output),
                               CI.sd.low = as.numeric(quantile(sd.A.output, probs = 0.5-cred.level/2)),
                               CI.sd.high = as.numeric(quantile(sd.A.output, probs = 0.5+cred.level/2)),
                               EV.CV = as.numeric(quantile(CV.A.output, probs = 0.5)),
                               CI.CV.low = as.numeric(quantile(CV.A.output, probs = 0.5-cred.level/2)),
                               CI.CV.high = as.numeric(quantile(CV.A.output, probs = 0.5+cred.level/2)))
    
    # SD and CV of the emission factors
    # ---------------------------------
    for (i.p in 1:NP) {
      # standard deviation output
      sd.EF.output <- output.df[, paste0("sd.logEF", i.p)]
      # Coefficient of variation or relative standard deviation
      CV.EF.output <- 100*sqrt(exp(sd.EF.output^2)-1) # coefficient of variation (%)
      sc.sd.res.df <- rbind(sc.sd.res.df,
                            data.frame(city = city, sector = sector,
                                       outcome = "sigma.EF", pollutant = pollutant.list[i.p], 
                                       output.name = paste0("sd.logEF", i.p),
                                       label = paste(toupper(pollutant.list[i.p]), "EF"),
                                       EV.sd = mean(sd.EF.output),
                                       CI.sd.low = as.numeric(quantile(sd.EF.output, probs = 0.5-cred.level/2)),
                                       CI.sd.high = as.numeric(quantile(sd.EF.output, probs = 0.5+cred.level/2)),
                                       EV.CV = as.numeric(quantile(CV.EF.output, probs = 0.5)),
                                       CI.CV.low = as.numeric(quantile(CV.EF.output, probs = 0.5-cred.level/2)),
                                       CI.CV.high = as.numeric(quantile(CV.EF.output, probs = 0.5+cred.level/2))))
      
    }
    
    # ranking of the sigmas
    sc.sd.res.df <- sc.sd.res.df[order(-sc.sd.res.df$EV.sd),]
    sc.sd.res.df$prob.sd.gt.next <- NA
    i.row <- 1
    n.sample <- NROW(output.df)
    for (i.row in 1:(NROW(sc.sd.res.df)-1)) {
      # get the index of each pollutant in the 'pollutant.list' vector
      this.output <- toString(sc.sd.res.df$output.name[i.row])
      next.output <- toString(sc.sd.res.df$output.name[i.row + 1])
      sc.sd.res.df$prob.sd.gt.next[i.row] <- sum(output.df[, this.output] > output.df[, next.output]) / n.sample
    }
    
    # compare every sigma with every other
    sc.sd.compare.df <- data.frame()
    sd.output.vec <- as.vector(unique(sc.sd.res.df$output.name))
    for (sd.output.1 in sd.output.vec) {
      label.1 <- sc.sd.res.df$label[which(sc.sd.res.df$output.name == sd.output.1)]
      for (sd.output.2 in sd.output.vec[sd.output.vec != sd.output.1]) {
        label.2 <- sc.sd.res.df$label[which(sc.sd.res.df$output.name == sd.output.2)]
        prob.sd.1.gt.2 <- sum(output.df[, sd.output.1] > output.df[, sd.output.2]) / n.sample
        sc.sd.compare.df <- rbind(sc.sd.compare.df,
                                  data.frame(city = city, sector = sector,
                                             label.1 = label.1, label.2 = label.2,
                                             prob.sd.1.gt.2 = prob.sd.1.gt.2))
      }
    }
    sd.compare.df <- rbind(sd.compare.df, sc.sd.compare.df)
    
    sc.sd.res.df$label <- factor(sc.sd.res.df$label, levels = rev(sc.sd.res.df$label), ordered = T)
    
    # plot coefficient of variation with error bars
    p <- ggplot(data = sc.sd.res.df, aes(x=EV.sd, y = label)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.sd.low, xmax = CI.sd.high))
    p <- p + labs(title = paste(city, sector, "\nCV with 2-sided", round(cred.level*100), "% CI"), 
                  x = "sd of log(A or EF)", y="")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig2a_StdDev", ".png")))
    print(p)
    dev.off()
    
    
    # plot coefficient of variation with error bars
    p <- ggplot(data = sc.sd.res.df, aes(x=EV.CV, y = label)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.CV.low, xmax = CI.CV.high))
    p <- p + labs(title = paste(city, sector, "\nCV with 2-sided", round(cred.level*100), "% CI"), 
                  x = "Coefficient of variation (%)", y="")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig2b_CV", ".png")))
    print(p)
    dev.off()

    # Distribution of log Activity ratios between pairs of inventories
    # ----------------------------------------------------------------
    # data frame with results of activity ratios (log-activity differences)
    logA.res.df <- data.frame()
    # loop over all inventory pairs
    for (i1 in 1:(NI-1)) {
      for (i2 in (i1+1):NI) {
        inv1 <- inventory.list[i1]
        inv2 <- inventory.list[i2]
        lr.label <- paste0("Activity_ratio_", inv1, "_", inv2)
        # vector with all the samples of the log of the activity ratio
        lrA1A2 <- output.df[,paste0("logA",i1)] - output.df[,paste0("logA",i2)]
        mean.lrA1A2 <- mean(lrA1A2)
        if (mean.lrA1A2 < 0) {
          inv1 <- inventory.list[i2]
          inv2 <- inventory.list[i1]
          lrA1A2 <- -lrA1A2
          mean.lrA1A2 <- -mean.lrA1A2
        }
        
        # # calculate the emission ratios
        # logE1E2 <- Y[inv1,] - Y[inv2,]
        # asymtotes <- sort(c(0, as.numeric(logE1E2)))
        # 
        # # function giving the zeros in the 2 inventory case
        # zero.function <- function(ar, logE1E2) {
        #   1/ar + sum(1/(ar - logE1E2))
        # }
        # 
        # ar.zero.vec <- rep(NA, length(logE1E2))
        # for (i.interval in 1:length(logE1E2)) {
        #   width.interval <- asymtotes[i.interval+1] - asymtotes[i.interval]
        #   search.interval <- c(asymtotes[i.interval], asymtotes[i.interval+1]) + width.interval/10000 * c(1,-1)
        #   uniroot.output <- uniroot(f = zero.function, interval = search.interval, logE1E2) 
        #   ar.zero.vec[i.interval] <- uniroot.output$root
        # }

        # percent difference of all samples
        pct.diff.A1A2 <- 100*(exp(lrA1A2)-1)
        mean.pct.diff.A1A2 <- mean(pct.diff.A1A2)
        prob.below.zero <- sum(lrA1A2 <= 0) / length(lrA1A2)
        
        # PDF of log of the activity ratio
        p <- ggplot(data.frame(lrA1A2 = lrA1A2), aes(x=lrA1A2)) + geom_density()
        p <- p + labs(title = paste0(lr.label, " (Pr<=0 = ", round(prob.below.zero*100,1), ")"))
        p <- p + theme(text = element_text(size=plot.font.size))
        p <- p + geom_vline(xintercept = mean.lrA1A2, col = "green")
        # p <- p + geom_vline(xintercept = asymtotes, col = "red")
        # p <- p + geom_vline(xintercept = c(0, ar.zero.vec), col = "grey")
        # print(p)
        png(file.path(cs.output.folder, paste0(city, "_", sector, "_", lr.label, ".png")))
        print(p)
        dev.off()
        # table with the results for activities
        logA.res.df <- rbind(logA.res.df,
                             data.frame(city = city,
                                  sector = sector,
                                  pollutant = "",
                                  outcome = "logA.ratio",
                                  inventory1 = inv1,
                                  inventory2 = inv2,
                                  EV.log.ratio = mean(lrA1A2),
                                  CI.low = as.numeric(quantile(lrA1A2, probs = 0.5-cred.level/2 )),
                                  CI.high = as.numeric(quantile(lrA1A2, probs = 0.5+cred.level/2 )),
                                  EV.pct.diff = mean.pct.diff.A1A2,
                                  CI.low.pct.diff = as.numeric(quantile(pct.diff.A1A2, probs = 0.5-cred.level/2 )),
                                  CI.high.pct.diff = as.numeric(quantile(pct.diff.A1A2, probs = 0.5+cred.level/2 )),
                                  prob.below.zero = prob.below.zero,
                                  significant = if(prob.below.zero < 0.05) {"1"} else {"0"}))
      }
    }
    # add the results of acitivities to the results data frame
    sc.res.df <- rbind(sc.res.df, logA.res.df)
    
    # one plot with all the logA's
    logA.df <- data.frame()
    for (i1 in 1:NI) {
      logA.df <- rbind(logA.df, data.frame(inventory = inventory.list[i1], 
                                           logA = output.df[, paste0("logA",i1)]))
    }
    p <- ggplot(logA.df, aes(x=logA, col=inventory)) + geom_density()
    p <- p + labs(title = "logA distributions")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig3_logA_distributions.png")))
    print(p)
    dev.off()
    
    # error bar plot for each log(activity ratio), ordered from big to small
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
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig4_logA_ratios.png")))
    print(p)
    dev.off()

    # percent difference plot for activity ratios
    p <- ggplot(data = logA.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
    p <- p + geom_vline(xintercept = 0, col = "red")
    p <- p + labs(title = paste(city, sector, "\nRelative activity diff. with 2-sided 95% CI"), 
                  x = "A1/A2 - 1 (%)", y="Inventory pair")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(MCMC.results.path, paste0(city, "_", sector, "_fig5_A_pctdiff.png")))
    print(p)
    dev.off()
    
    # Distribution of log emission factor ratios
    logEF.res.df <- data.frame()
    for (p1 in 1:NP) {
      pol.name <- pollutant.list[p1]
      for (i1 in 1:(NI-1)) {
        for (i2 in (i1+1):NI) {
          inv1 <- inventory.list[i1]
          inv2 <- inventory.list[i2]
          # log ratio of all the samples
          lrEF1EF2 <- output.df[,paste0("logEF",i1,"_",p1)] - output.df[,paste0("logEF",i2,"_",p1)]
          mean.lrEF1EF2 <- mean(lrEF1EF2)
          if (mean.lrEF1EF2 < 0) {
            inv1 <- inventory.list[i2]
            inv2 <- inventory.list[i1]
            mean.lrEF1EF2 <- -mean.lrEF1EF2
            lrEF1EF2 <- -lrEF1EF2
          }
          # percent difference of all samples
          pct.diff.EF1EF2 <- 100*(exp(lrEF1EF2)-1)
          mean.pct.diff.EF1EF2 <- mean(pct.diff.EF1EF2)
          lr.label <- paste0(pol.name, "_EF_ratio_", inv1, "_", inv2)
          prob.below.zero <- sum(lrEF1EF2 <= 0) / length(lrEF1EF2)
          
          # calculate the emission ratios
          # logE1E2 <- Y[inv1,] - Y[inv2,]
          
          # probability density of the log ratio
          p <- ggplot(data.frame(lrEF1EF2 = lrEF1EF2), aes(x=lrEF1EF2)) + geom_density()
          p <- p + labs(title = paste0(lr.label, "\n(Pr<=0 = ", round(prob.below.zero*100,1), ")"))
          p <- p + theme(text = element_text(size=plot.font.size))
          # p <- p + geom_vline(xintercept = mean.lrEF1EF2, col = "green")
          # p <- p + geom_vline(xintercept = c(0, logE1E2), col = "red")
          # print(p)
          png(file.path(cs.output.folder, paste0(city, "_", sector, "_", lr.label, ".png")))
          print(p)
          dev.off()

          # table with the results for activities
          logEF.res.df <- rbind(logEF.res.df,
                                data.frame(city = city,
                                           sector = sector,
                                           outcome = "logEF.ratio",
                                           pollutant = pol.name,
                                           inventory1 = inv1,
                                           inventory2 = inv2,
                                           EV.log.ratio = mean.lrEF1EF2,
                                           CI.low = as.numeric(quantile(lrEF1EF2, probs = 0.5-cred.level/2 )),
                                           CI.high = as.numeric(quantile(lrEF1EF2, probs = 0.5+cred.level/2 )),
                                           EV.pct.diff = mean.pct.diff.EF1EF2,
                                           CI.low.pct.diff = as.numeric(quantile(pct.diff.EF1EF2, probs = 0.5-cred.level/2 )),
                                           CI.high.pct.diff = as.numeric(quantile(pct.diff.EF1EF2, probs = 0.5+cred.level/2 )),
                                           prob.below.zero = prob.below.zero,
                                           significant = if(prob.below.zero < 0.05) {"1"} else {"0"}))
        }
      }
    }
    # add the results of emission factors to the results data frame
    sc.res.df <- rbind(sc.res.df,logEF.res.df)

    # one plot with all the logEF's
    for (p1 in 1:NP) {
      logEF.df <- data.frame()
      for (i1 in 1:NI) {
        logEF.df <- rbind(logEF.df, data.frame(inventory = inventory.list[i1], 
                                               logEF = output.df[, paste0("logEF", i1, "_", p1)]))
      }
      p <- ggplot(logEF.df, aes(x=logEF, col=inventory)) + geom_density()
      p <- p + labs(title = paste0("logEF ", pollutant.list[p1], " distributions"))
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(MCMC.results.path, 
                    paste0(city, "_", sector, "_fig6", letters[p1], "_logEF_", 
                           toupper(pollutant.list[p1]), "_distributions.png")))
      print(p)
      dev.off()
    }
    write.table(sc.res.df, file = paste0("MCMC_results_", city, "_", sector, ".csv"), row.names = F, sep = ",")
    res.df <- rbind(res.df, sc.res.df)
    
    # error bar plot for logEF ratios
    for (i.p in 1:NP) {
      logEF.plot.df <- sc.res.df[sc.res.df$pollutant == pollutant.list[i.p],]
      logEF.plot.df$inventories <- paste0(logEF.plot.df$inventory1, " - ", logEF.plot.df$inventory2)
      logEF.plot.df <- logEF.plot.df[order(logEF.plot.df$EV.log.ratio),]
      logEF.plot.df$inventories <- factor(logEF.plot.df$inventories, ordered = T, levels = logEF.plot.df$inventories)
      p <- ggplot(data = logEF.plot.df, aes(x=EV.log.ratio, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low, xmax = CI.high))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutant.list[i.p]), 
                                  "\nlog(EF ratios) with 2-sided 95% CI"), 
                    x = "log(EF1/EF2)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(MCMC.results.path, 
                    paste0(city, "_", sector, "_fig7", letters[p1], "_logEF_", 
                           pollutant.list[i.p], "_ratios.png")))
      print(p)
      dev.off()
      
      # error bar plot for percent difference in emission factors
      p <- ggplot(data = logEF.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutant.list[i.p]), 
                                  "\nRelative EF diff. with 2-sided 95% CI"), 
                    x = "EF1/EF2 - 1 (%)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(MCMC.results.path, 
                    paste0(city, "_", sector, "_fig8", letters[p1], "_EF_", pollutant.list[i.p], "_pctdiff.png")))
      print(p)
      dev.off()
    }
    
    # add them to the rest
    sd.res.df <- rbind(sd.res.df, sc.sd.res.df)
  }
}

# write the results to a file:
# - the expected values and confidence intervals of activity and emission factor differences
write.table(res.df, file = "ratio_results_MCMC.csv", row.names = F, sep = ",")
# - the expected values of standard deviations and coefficients of variation an their confidence intervals
write.table(sd.res.df, file = "stddev_results_MCMC.csv", row.names = F, sep = ",")
# - pairwise comparison of activity/emission factor SDs and credibility that one is bigger than the other
# based on MCMC samples.
write.table(sd.compare.df, file = "stddev_compare_results_MCMC.csv", row.names = F, sep = ",")
