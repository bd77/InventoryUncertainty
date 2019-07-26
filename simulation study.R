# --------------------------------------------------------------
# Simulation study: uncertainty analysis of emission inventories
# --------------------------------------------------------------

rm(list = ls())
wd <- "D:/DiamondPlot"
setwd(wd)

library(R2jags)
library(plyr)
library(ggplot2)
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

# 1) One more uncertain EF
# ------------------------
sim.input.1 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 0.1, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.2),
                    output.folder = "1_One_Uncertain_EF")
# create a data set
inv.1.df <- create.dataset(sim.input.1) 
# visualize it
visualize.invenory(sim.input.1, inv.1.df)
# run jags
jagsfit.1 <- run.jags(sim.input.1, inv.1.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.1, inv.1.df, jagsfit.1)

# 2) One very uncertain EF
# ------------------------
sim.input.2 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 0.1, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.5),
                    output.folder = "2_One_Very_Uncertain_EF")
# create a data set
inv.2.df <- create.dataset(sim.input.2) 
# visualize it
visualize.invenory(sim.input.2, inv.2.df)
# run jags
jagsfit.2 <- run.jags(sim.input.2, inv.2.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.2, inv.2.df, jagsfit.2)

# 3) Uncertain AC
# ------------------------
sim.input.3 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 0.2, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.1),
                    output.folder = "3_Uncertain_AC")
# create a data set
inv.3.df <- create.dataset(sim.input.3) 
# visualize it
visualize.invenory(sim.input.3, inv.3.df)
# run jags
jagsfit.3 <- run.jags(sim.input.3, inv.3.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.3, inv.3.df, jagsfit.3)

# 4) Very Uncertain AC
# ------------------------
sim.input.4 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 1, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.1),
                    output.folder = "4_Very_Uncertain_AC")
# create a data set
inv.4.df <- create.dataset(sim.input.4) 
# visualize it
visualize.invenory(sim.input.4, inv.4.df)
# run jags
jagsfit.4 <- run.jags(sim.input.4, inv.4.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.4, inv.4.df, jagsfit.4)

# 5) Extremely Uncertain AC
# ------------------------
sim.input.5 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 2, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.1),
                    output.folder = "5_Extremely_Uncertain_AC")
# create a data set
inv.5.df <- create.dataset(sim.input.5) 
# visualize it
visualize.invenory(sim.input.5, inv.5.df)
# run jags
jagsfit.5 <- run.jags(sim.input.5, inv.5.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.5, inv.5.df, jagsfit.5)

# 6) One more uncertain EF, more inventories
# ------------------------
sim.input.6 <- list(n.inv = 30, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 0.1, rel.sigma.EF = c(0.1, 0.1, 0.1, 0.2),
                    output.folder = "6_Uncertain_EF_Many_Invs")
# create a data set
inv.6.df <- create.dataset(sim.input.6) 
# visualize it
visualize.invenory(sim.input.6, inv.6.df)
# run jags
jagsfit.6 <- run.jags(sim.input.6, inv.6.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.6, inv.6.df, jagsfit.6)

# 7) Very certain activity and one EF
# -----------------------------------
sim.input.7 <- list(n.inv = 7, n.pol = 4, real.AC = 2, real.EF = c(2, 4, 6, 8), 
                    rel.sigma.AC = 0.01, rel.sigma.EF = c(0.01, 0.2, 0.2, 0.2),
                    output.folder = "7_Certain_AC_EF1")
# create a data set
inv.7.df <- create.dataset(sim.input.7) 
# visualize it
visualize.invenory(sim.input.7, inv.7.df)
# run jags
jagsfit.7 <- run.jags(sim.input.7, inv.7.df, error_model_string)
# visualize the output
visualize.jagsfit(sim.input.7, inv.7.df, jagsfit.7)
