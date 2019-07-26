# -------------------------------------------------------------
# Emission inventories: where does the uncertainty come form?
# Activity or Emission factors?
# -------------------------------------------------------------

# Input: Total emissions for several cities from 7 inventories for
# 3 sectors and 4 pollutants.
# Question: Does the variablility in the total emissions tell us something
# about the ratio of the relative errors of the activity and the emission factor.
# relative error of activity = std. dev. activity / mean activity
# idem for emission factor.

rm(list = ls())
wd <- "D:/DiamondPlot"
setwd(wd)

library(nlme)
library(plyr)
library(R2jags)
source("auxiliary functions.R")

# The stochastic model used in JAGS
error_model_string <- "model{
  # Likelihood
  for(i in 1:V){
    for (j in 1:P) {
      Y[(i-1)*P + j] ~ dnorm(mulogE[j], tau.error[j])
    }
  }
  for (j in 1:P) {
    tau.error[j] <- (tau.logAC^(-1) + tau.logEF[j]^(-1))^(-1)
  }
  # Priors
  tau.logAC  ~ dgamma(1.0E-2, 1.0E-2)
  for (j in 1:P) {
    tau.logEF[j] ~ dgamma(1.0E-2, 1.0E-2)
    mulogE[j] ~ dnorm(meanlogE[j], 1.0E-4)
  }
}"


# read the input data
inv.table.df <- read.table("cities_emission_inventories.txt", header = T, sep = "\t")

# reorganize, 4 pollutant columns to 1 with the pollutant name and 1 with the emissions
inv.df <- data.frame()
for (pollutant in c('nox', 'pm25', 'so2', 'voc')) {
  inv.pol.df <- inv.table.df[, c('city', 'inventory', 'sector', pollutant)]
  names(inv.pol.df)[4] <- "emission"
  inv.pol.df$pollutant <- pollutant
  inv.df <- rbind(inv.df, inv.pol.df)
}

CI.results.df <- data.frame()

city <- "London"
sector <- "ms7"
city.list <- as.vector(unique(inv.df$city))
sector.list <- as.vector(unique(inv.df$sector))
for (city in city.list) {
  for (sector in sector.list) {
    print(paste(city, sector))
    inv.city.sector.df <- inv.df[inv.df$city==city & inv.df$sector==sector,]
    inv.city.sector.df$logE <- log(inv.city.sector.df$emission)
    sim.cs.input <- list(n.inv = 7, n.pol = 4, real.AC = NA, real.EF = c(NA, NA, NA, NA), 
                         rel.sigma.AC = NA, rel.sigma.EF = c(NA, 0.1, 0.1, 0.1),
                         output.folder = file.path("urban_inventories", paste(city, sector, sep = "_")))
    # run jags
    jagsfit.sc <- run.jags(sim.cs.input, inv.city.sector.df, error_model_string)
    
    # ranking
    jagsfit.sc.mcmc <- as.mcmc(jagsfit.sc)
    output.df <- data.frame()
    for (i in 1:3) {
      chain.df <- cbind(chain = toString(i), iteration = 1:NROW(jagsfit.sc.mcmc[[i]]), 
                        as.data.frame(jagsfit.sc.mcmc[[i]]))
      output.df <- rbind(output.df, chain.df)
    }
    # rename columns of output.df: remove brackets, they're annoying
    output.col.names <- names(output.df)
    n.col.output <- NCOL(output.df)
    for (i.col in 1:n.col.output) {
      col.name <- output.col.names[i.col]
      new.col.name <- gsub('\\[', '', col.name)
      new.col.name <- gsub('\\]', '', new.col.name)
      names(output.df)[i.col] <- new.col.name
    }
    
    mean.taus <- colMeans(output.df[, c("tau.logAC", "tau.logEF1", "tau.logEF2", "tau.logEF3", "tau.logEF4")])
    ranking.df <- data.frame(varname = names(mean.taus), mean.tau = as.numeric(mean.taus))
    ranking.df <- ranking.df[order(ranking.df$mean.tau),]
    ranking.df$prob.gt.next <- NA
    N.iter <- NROW(output.df)
    ranking.str <- ""
    for (i in 1:(NROW(ranking.df)-1)) {
      this.tau <- toString(ranking.df$varname[i])
      next.tau <- toString(ranking.df$varname[i+1])
      ranking.df$prob.gt.next[i] <- sum(output.df[,this.tau] < output.df[,next.tau]) / N.iter
      ranking.str <- paste0(ranking.str, gsub("tau.log", "se", ranking.df$varname[i]), 
                            ">(", round(100*ranking.df$prob.gt.next[i]), "%)")
    }
    ranking.str <- paste0(ranking.str, gsub("sigma.log", "se", ranking.df$varname[NROW(ranking.df)]))
    # png(file.path(sim.cs.input$output.folder, paste0("Hist_Ranking.png")), width = 480*1.5, height = 480)
    # p <- ggplot(data = sigma.df, aes(x=log(sigma), col=varname)) + geom_density()
    # p <- p + labs(title = paste0("Ranking ", sim.cs.input$output.folder, "\n", ranking.str))
    # p <- p + theme(text = element_text(size=20))
    # print(p)
    # dev.off()
    
    CI.results.df <- rbind(CI.results.df,
                           data.frame(city=city, sector=sector,
                                      ranking.str=ranking.str))
    
  }
}

write.table(CI.results.df, file="cities_uncertainty_ranking.txt", row.names = F, quote = F, sep = "\t")
