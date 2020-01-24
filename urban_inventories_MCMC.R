# -------------------------------------------------------------
# Emission inventories: where does the uncertainty come form?
# Activity or Emission factors?
# -------------------------------------------------------------

# Input: Total emissions for several cities from 6 inventories for
# 3 sectors and 4 pollutants.
# Question: Does the variablility in the total emissions tell us something
# about the ratio of the relative errors of the activity and the emission factor?
# relative error of activity = std. dev. activity / mean activity
# idem for emission factor.

library(nlme)
library(plyr)
library(R2jags)
library(ggplot2)

set.seed(666)

rm(list = ls())
wd <- "D:/DiamondPlot/_code_for_paper/urban_inventories_JAGS_RANDOM_EFFECTS"
setwd(wd)

# 0.95 is chosen because for the log ratios we consider 1-sided credibility intervals
cred.level <- 0.95
# plot font size
plot.font.size <- 20
# make trace plots
trace.plots <- TRUE # takes a lot of time

# The stochastic model used in JAGS
random_effects_string <- "model{
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

# read the input data
inv.table.df <- read.table("cities_emission_inventories.txt", header = T, sep = "\t")
# the Synth inventory is the median of the others, remove it
inv.table.df <- inv.table.df[inv.table.df$inventory != "synth",]
inv.table.df$inventory <- factor(inv.table.df$inventory, levels = as.vector(unique(inv.table.df$inventory)))

# lists of cities and sectors
city.list <- as.vector(unique(inv.table.df$city)) # c("Barcelona", "Katowice") # 
sector.list <- as.vector(unique(inv.table.df$sector)) # c("ms7", "ms34") # 

# reorganize, 4 pollutant columns to 1 with the pollutant name and 1 with the emissions
inv.df <- data.frame()
for (pollutant in c('nox', 'pm25', 'so2', 'voc')) {
  inv.pol.df <- inv.table.df[, c('city', 'inventory', 'sector', pollutant)]
  names(inv.pol.df)[4] <- "emission"
  inv.pol.df$pollutant <- pollutant
  inv.df <- rbind(inv.df, inv.pol.df)
}

# An overview of the coefficients of variation (relative standard deviations) of the emissions 
# and the standard deviations of the emissions
# These should be related by sd/mean = sqrt(exp(var(log(E)))-1) if E is log-nomrally distributed.
sd.city.sector.pol.overview <- ddply(inv.df, c('city', 'sector', 'pollutant'), summarise, 
                                     rel.sd.e = sd(emission)/mean(emission), 
                                     sd.log.e = sqrt(exp(var(log(emission)))-1))
# some outliers here but the assumption of a log-normal distribution looks ok
p <- ggplot(data = sd.city.sector.pol.overview, aes(x=rel.sd.e, y=sd.log.e)) + geom_point()
p <- p + geom_abline(intercept = 0, slope = 1, col="red")
p <- p + labs(x="Relative sd of the emissions", y="sd of log-emissions")
png("relative_sd_VS_sd_of_logE.png")
print(p)
dev.off()
# The sd of the log VOC emissions in Katowice for ms34 are an extreme outlier (10.0)
# Normal values are around 0.5

# boxplots of sd of log-emissions
p <- ggplot(sd.city.sector.pol.overview, aes(x=pollutant, y=sd.log.e)) + geom_boxplot()
p <- p + lims(y=c(0,2.5))
png("Boxplot_sd_logE_per_pollutant.png")
print(p)
dev.off()

# find a distribution that describes the sd of the emissions. It is used for the simulation study.
# log-normal works fine
p <- ggplot() + geom_line(data = data.frame(x=seq(0,10,0.1), y=dlnorm(seq(0,10,0.1), meanlog = -0.59, sdlog = 0.73)),
                          aes(x = x, y = y, col="red")) + geom_density(data=sd.city.sector.pol.overview, aes(x=sd.log.e))
p

# range of standard deviations to be used in the simulation studies
sd.df <- ddply(inv.df, c('city', 'sector', 'pollutant'), summarise, sd.pol = sd(log(emission)))
range(sd.df$sd.pol)
# quantiles (there are a few big outliers)
quantile(sd.df$sd.pol, probs = c(0.05, 0.95))
# 5%       95% 
# 0.1903845 1.2808845 
# as percentage error
sqrt(exp(quantile(sd.df$sd.pol, probs = c(0.05, 0.95))^2)-1)
# 5%       95% 
# 0.1921228 2.0392646

# quantile plot: log-emissions verus normal distribution
png("Normal quantile plot of logE.png")
qqplot(rnorm(10*NROW(inv.df)), log(inv.df$emission),
       xlab="Quantiles of a normal distribution",
       ylab="Quantiles of log(emissions)")
dev.off()
# quantile plots per city-sector combination
if (!dir.exists("quantile_plots")) {dir.create("quantile_plots")}
for (city in city.list) {
  for (sector in sector.list) {
    png(file.path("quantile_plots", 
                  paste0("Normal quantile plot of logE for", city, " ", sector, ".png")))
    qqplot(rnorm(10*NROW(inv.df)), log(inv.df$emission[inv.df$city==city & inv.df$sector==sector]),
           xlab="Quantiles of a normal distribution",
           ylab="Quantiles of log(emissions)",
           main = paste(city, sector))
    dev.off()
  }
}

# data frame with all the results for activities and emission factors
res.df <- data.frame()
# data frame with all the results for standard deviations
sd.res.df <- data.frame()
# comparison of standard deviations
sd.compare.df <- data.frame()


plot.path <- "plots"
if (!dir.exists(plot.path)) {dir.create(plot.path)}

# for testing
city <- "Budapest"
sector <- "ms34"

# loop over cities and sectors
for (city in city.list) {
  for (sector in sector.list) {
    sc.res.df <- data.frame()
    
    print(paste(city, sector))
    # select just the data of one city-sector combination
    inv.city.sector.df <- inv.df[inv.df$city==city & inv.df$sector==sector,]
    # add a column with log emissions
    inv.city.sector.df$logE <- log(inv.city.sector.df$emission)
    # sort the data !!! ORDER IMPORTANT FOR NEXT STEP !!!
    inv.city.sector.df <- inv.city.sector.df[order(inv.city.sector.df$inventory,inv.city.sector.df$pollutant),]
    
    # output folder for the city-sector
    cs.output.folder <- file.path(paste(city, sector, sep = "_"))
    if (!(dir.exists(cs.output.folder))) {
      dir.create(cs.output.folder)
    }
    
    pollutants <- as.vector(unique(inv.city.sector.df$pollutant))
    NP <- length(pollutants)
    inventories <- as.vector(unique(inv.city.sector.df$inventory))
    NI <- length(inventories)
      
    # log emissions as a matrix
    Y <- matrix(0, nrow = NI, ncol = NP, dimnames = list(inventories, pollutants))
    for (i in 1:NI) {
      for (j in 1:NP) {
        Y[i,j] <- inv.city.sector.df$logE[inv.city.sector.df$inventory==inventories[i] &
                                            inv.city.sector.df$pollutant == pollutants[j]]
      }
    }
    
    # Line plot of logE
    p <- ggplot(data = inv.city.sector.df, 
                aes(x = pollutant, y = logE, group = inventory, col = inventory)) 
    p <- p + geom_line() + geom_point()
    p <- p + labs(title = paste(city, sector))
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(plot.path, paste0(city, "_", sector, "_0_logE_LinePlot.png")))
    print(p)
    dev.off()
    
    # -------- JAGS ----------
    # JAGS simulation
    jags.input.list <- list(Y = Y, NI = NI, NP = NP)
    jagsfit <- jags(data = jags.input.list, n.iter=10000,
                    parameters.to.save = c('logA', 'logEF', 'mu.logEF', 'sd.logEF', 'sd.logA'),
                    model.file = textConnection(random_effects_string))
    
    # process the output
    jagsfit.mcmc <- as.mcmc(jagsfit)
    output.df <- data.frame()
    for (i in 1:3) {
      chain.df <- cbind(chain = toString(i), iteration = 1:NROW(jagsfit.mcmc[[i]]), 
                        as.data.frame(jagsfit.mcmc[[i]]))
      output.df <- rbind(output.df, chain.df)
    }
    
    # rename columns of output.df: remove brackets, they're annoying
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
    
    # trace plots and PDF per chain
    if (trace.plots == TRUE) {
      for (varname in names(output.df)[3:NCOL(output.df)]) {
        # from varname to real names of inventory and pollutant
        if (grepl(pattern = "^(logA)", varname)) {
          inv.id <- as.numeric(gsub(pattern = "logA", "", varname))
          var.label <- paste0('logA_', inventories[inv.id])
        } else if (grepl(pattern = "^(logEF)", varname)) {
          inv.pol.id <- gsub(pattern = "logEF", "", varname)
          inv.id <- as.numeric(strsplit(inv.pol.id, "_")[[1]][1])
          pol.id <- as.numeric(strsplit(inv.pol.id, "_")[[1]][2])
          var.label <- paste0('logEF_', inventories[inv.id], '_', pollutants[pol.id])
        } else if (grepl(pattern = "^(sd.logEF)", varname)) {
          pol.id <- as.numeric(gsub(pattern = "sd.logEF", "", varname))
          var.label <- paste0('sd.logEF_', pollutants[pol.id])
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
    
    # Distribution of log Activity ratios between pairs of inventories
    # i1 <- 1 # for testing
    # i2 <- 2 # for testing
    # data frame with results of activity ratios (log-activity differences)
    logA.res.df <- data.frame()
    # loop over all inventory pairs
    for (i1 in 1:(NI-1)) {
      for (i2 in (i1+1):NI) {
        inv1 <- inventories[i1]
        inv2 <- inventories[i2]
        lr.label <- paste0("Activity_ratio_", inv1, "_", inv2)
        # vector with all the samples of the log of the activity ratio
        lrA1A2 <- output.df[,paste0("logA",i1)] - output.df[,paste0("logA",i2)]
        mean.lrA1A2 <- mean(lrA1A2)
        if (mean.lrA1A2 < 0) {
          inv1 <- inventories[i2]
          inv2 <- inventories[i1]
          lrA1A2 <- -lrA1A2
          mean.lrA1A2 <- -mean.lrA1A2
        }
        
        # calculate the emission ratios
        logE1E2 <- Y[inv1,] - Y[inv2,]
        asymtotes <- sort(c(0, as.numeric(logE1E2)))
        
        # function giving the zeros in the 2 inventory case
        zero.function <- function(ar, logE1E2) {
          1/ar + sum(1/(ar - logE1E2))
        }
        
        ar.zero.vec <- rep(NA, length(logE1E2))
        for (i.interval in 1:length(logE1E2)) {
          width.interval <- asymtotes[i.interval+1] - asymtotes[i.interval]
          search.interval <- c(asymtotes[i.interval], asymtotes[i.interval+1]) + width.interval/10000 * c(1,-1)
          uniroot.output <- uniroot(f = zero.function, interval = search.interval, logE1E2) 
          ar.zero.vec[i.interval] <- uniroot.output$root
        }

        # percent difference of all samples
        pct.diff.A1A2 <- 100*(exp(lrA1A2)-1)
        mean.pct.diff.A1A2 <- mean(pct.diff.A1A2)
        prob.below.zero <- sum(lrA1A2 <= 0) / length(lrA1A2)
        
        # PDF of log of the activity ratio
        p <- ggplot(data.frame(lrA1A2 = lrA1A2), aes(x=lrA1A2)) + geom_density()
        p <- p + labs(title = paste0(lr.label, " (Pr<=0 = ", round(prob.below.zero*100,1), ")"))
        p <- p + theme(text = element_text(size=plot.font.size))
        p <- p + geom_vline(xintercept = mean.lrA1A2, col = "green")
        p <- p + geom_vline(xintercept = asymtotes, col = "red")
        p <- p + geom_vline(xintercept = c(0, ar.zero.vec), col = "grey")
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
      logA.df <- rbind(logA.df, data.frame(inventory = inventories[i1], 
                                           logA = output.df[, paste0("logA",i1)]))
    }
    p <- ggplot(logA.df, aes(x=logA, col=inventory)) + geom_density()
    p <- p + labs(title = "logA distributions")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(plot.path, paste0(city, "_", sector, "_logA_distributions.png")))
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
    p <- p + labs(title = paste(city, sector, "Activity", 
                                "\n% diff with 2-sided ", round(cred.level*100), "% CI"), 
                  x = "log(A1/A2)", y="Inventory pair")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(plot.path, paste0(city, "_", sector, "_logA_ratios.png")))
    print(p)
    dev.off()

    # percent difference plot for activity ratios
    p <- ggplot(data = logA.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
    p <- p + geom_vline(xintercept = 0, col = "red")
    p <- p + labs(title = paste(city, sector, "Activity", 
                                "\n% diff with 2-sided ", round(cred.level*100), "% CI"), 
                  x = "A1/A2 - 1 (%)", y="Inventory pair")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(plot.path, paste0(city, "_", sector, "_A_pctdiff.png")))
    print(p)
    dev.off()
    
    # Distribution of log emission factor ratios
    i1<-1
    i2<-2
    logEF.res.df <- data.frame()
    for (p1 in 1:NP) {
      pol.name <- pollutants[p1]
      for (i1 in 1:(NI-1)) {
        for (i2 in (i1+1):NI) {
          inv1 <- inventories[i1]
          inv2 <- inventories[i2]
          # log ratio of all the samples
          lrEF1EF2 <- output.df[,paste0("logEF",i1,"_",p1)] - output.df[,paste0("logEF",i2,"_",p1)]
          mean.lrEF1EF2 <- mean(lrEF1EF2)
          if (mean.lrEF1EF2 < 0) {
            inv1 <- inventories[i2]
            inv2 <- inventories[i1]
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
        logEF.df <- rbind(logEF.df, data.frame(inventory = inventories[i1], 
                                               logEF = output.df[, paste0("logEF", i1, "_", p1)]))
      }
      p <- ggplot(logEF.df, aes(x=logEF, col=inventory)) + geom_density()
      p <- p + labs(title = paste0("logEF ", pollutants[p1], " distributions"))
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(plot.path, paste0(city, "_", sector, "_logEF_", toupper(pollutants[p1]), "_distributions.png")))
      print(p)
      dev.off()
    }
    write.table(sc.res.df, file = paste0("ratio_results_", city, "_", sector, ".csv"), row.names = F, sep = ",")
    res.df <- rbind(res.df, sc.res.df)
    
    # error bar plot for logEF ratios
    for (i.p in 1:NP) {
      logEF.plot.df <- sc.res.df[sc.res.df$pollutant == pollutants[i.p],]
      logEF.plot.df$inventories <- paste0(logEF.plot.df$inventory1, " - ", logEF.plot.df$inventory2)
      logEF.plot.df <- logEF.plot.df[order(logEF.plot.df$EV.log.ratio),]
      logEF.plot.df$inventories <- factor(logEF.plot.df$inventories, ordered = T, levels = logEF.plot.df$inventories)
      p <- ggplot(data = logEF.plot.df, aes(x=EV.log.ratio, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low, xmax = CI.high))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutants[i.p]), 
                                  "\n% diff with 2-sided ", round(cred.level*100), "% CI"), 
                    x = "log(EF1/EF2)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(plot.path, paste0(city, "_", sector, "_logEF_", pollutants[i.p], "_ratios.png")))
      print(p)
      dev.off()
      
      # error bar plot for percent difference in emission factors
      p <- ggplot(data = logEF.plot.df, aes(x=EV.pct.diff, y = inventories)) + geom_point()
      p <- p + geom_errorbarh(aes(xmin = CI.low.pct.diff, xmax = CI.high.pct.diff))
      p <- p + geom_vline(xintercept = 0, col = "red")
      p <- p + labs(title = paste(city, sector, toupper(pollutants[i.p]), 
                                  "\n% diff with 2-sided ", round(cred.level*100), "% CI"), 
                    x = "EF1/EF2 - 1 (%)", y="Inventory pair")
      p <- p + theme(text = element_text(size=plot.font.size))
      png(file.path(plot.path, paste0(city, "_", sector, "_EF_", pollutants[i.p], "_pctdiff.png")))
      print(p)
      dev.off()
    }
    
    # Results for standard deviations and coefficient of variation
    # -------------------------------------------------------------
    
    # Std dev or CV of the activity
    # -----------------------------
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
    
    for (i.p in 1:NP) {
      # standard deviation output
      sd.EF.output <- output.df[, paste0("sd.logEF", i.p)]
      # Coefficient of variation or relative standard deviation
      CV.EF.output <- 100*sqrt(exp(sd.EF.output^2)-1) # coefficient of variation (%)
      sc.sd.res.df <- rbind(sc.sd.res.df,
                            data.frame(city = city, sector = sector,
                                       outcome = "sigma.EF", pollutant = pollutants[i.p], 
                                       output.name = paste0("sd.logEF", i.p),
                                       label = paste(toupper(pollutants[i.p]), "EF"),
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
      # get the index of each pollutant in the 'pollutants' vector
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
    png(file.path(plot.path, paste0(city, "_", sector, "_StdDev", ".png")))
    print(p)
    dev.off()
    
    
    # plot coefficient of variation with error bars
    p <- ggplot(data = sc.sd.res.df, aes(x=EV.CV, y = label)) + geom_point()
    p <- p + geom_errorbarh(aes(xmin = CI.CV.low, xmax = CI.CV.high))
    p <- p + labs(title = paste(city, sector, "\nCV with 2-sided", round(cred.level*100), "% CI"), 
                  x = "Coefficient of variation (%)", y="")
    p <- p + theme(text = element_text(size=plot.font.size))
    png(file.path(plot.path, paste0(city, "_", sector, "_CV", ".png")))
    print(p)
    dev.off()
    
    # add them to the rest
    sd.res.df <- rbind(sd.res.df, sc.sd.res.df)
  }
}

write.table(res.df, file = "ratio_results_MCMC.csv", row.names = F, sep = ",")
write.table(sd.res.df, file = "stddev_results_MCMC.csv", row.names = F, sep = ",")
write.table(sd.compare.df, file = "stddev_compare_results_MCMC.csv", row.names = F, sep = ",")
