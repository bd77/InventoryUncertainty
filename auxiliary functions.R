# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      AUXILIARY FUNCTIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1) create.dataset
# Creates a dataset for 1 sector P pollutants and 

# function to create a synthetic data set
# ---------------------------------------

create.dataset <- function(sim.input) {
  # input: sim.input, a list with
  # - n.inv: number of invenories
  # - n.pol: number of pollutants
  # - real.AC: the real activity of the sector (1 float)
  # - real.EF: the real emission factor for each pollutant (1xP float vector)
  # - rel.sigma.AC: relative standard deviation of the activity sigma(AC)/mean(AC) (1 float)
  # - rel.sigma.AC: relative standard deviation of the emission factors sigma(EF_p)/mean(EF_p)
  # output: a dataset with logE = logAC + logEF_p for n.inv*n.pol combinations. logAC of each
  # inventory is sampled from a normal distribution N(log(real.AC), )
  n.inv <- sim.input$n.inv 
  n.pol <- sim.input$n.pol
  real.AC <- sim.input$real.AC
  real.EF <- sim.input$real.EF
  rel.sigma.AC<- sim.input$rel.sigma.AC
  rel.sigma.EF <- sim.input$rel.sigma.EF
  
  # relative uncertainty > absolute uncertainty
  real.logAC <- log(real.AC)
  sigma.logAC <- sqrt(log(rel.sigma.AC^2+1))
  real.logEF <- log(real.EF)
  sigma.logEF <- sqrt(log(rel.sigma.EF^2+1))
  
  # number of rows
  N <- n.pol*n.inv
  # vectors with pollutant and inventory names
  pollutant.vec <- paste0("pol", 1:n.pol)
  inventory.vec <- paste0("inv", 1:n.inv)
  
  # all combinations of inventories and pollutants
  # !!! ORDER IMPORTANT FOR NEXT STEP !!! 
  # expand.grid sorts first per inventory, than pollutant
  inv.df <- expand.grid(pollutant = pollutant.vec, inventory = inventory.vec)
  inv.df <- inv.df[,c('inventory', 'pollutant')] # looks nicer
  
  # create random activities and emission factors with the given sigma
  # same activity for all pollutants of the same inventory
  inv.df$real.logAC <- real.logAC
  # sample log activity for each inventory
  inv.df$logAC <- real.logAC + rep(rnorm(n.inv, mean = 0, sd = sigma.logAC), each = n.pol)
  
  # the pattern in sigma.ef is repeated n.inv times
  inv.df$real.logEF <- rep(real.logEF, times = n.inv)
  inv.df$logEF <- inv.df$real.logEF + rnorm(N, mean = 0, sd = sigma.logEF)
  
  # calculate log of the emissions
  inv.df$logE <- inv.df$logAC + inv.df$logEF
  
  return(inv.df)
}

# Function to visualize the inventory
# ------------------------------------
visualize.invenory <- function(sim.input, inv.df) {
  # create output folder
  if (!(dir.exists(sim.input$output.folder))) {
    dir.create(sim.input$output.folder)
  }
  # Boxplot of the log activities
  p <- ggplot(data = inv.df, aes(x="", y=logAC)) + geom_boxplot()  
  p <- p + geom_point(aes(col=inventory)) + geom_point(aes(x="", y=real.logAC), shape = 4)
  p <- p + labs(title = paste("Log activities for", sim.input$output.folder), x = "")  
  p <- p + theme(text = element_text(size=20))
  png(file.path(sim.input$output.folder, "Boxplot of logAC.png"))
  print(p)
  dev.off()
  # Boxplot of the log emission factors
  p <- ggplot(data = inv.df, aes(x=pollutant, y=logEF)) + geom_boxplot()  
  p <- p + geom_point(aes(col=inventory)) + geom_point(aes(x=pollutant, y=real.logEF), shape = 4)
  p <- p + labs(title = paste("Log emission factors for", sim.input$output.folder))
  p <- p + theme(text = element_text(size=20))
  png(file.path(sim.input$output.folder, "Boxplot of logEF.png"))
  print(p)
  dev.off()
  # Boxplot of the log emissions
  p <- ggplot(data = inv.df, aes(x=pollutant, y=logE)) + geom_boxplot()  
  p <- p + geom_point(aes(col=inventory)) + geom_point(aes(x=pollutant, y=real.logAC + real.logEF), shape = 4)
  p <- p + labs(title = paste("Log emissions for", sim.input$output.folder))
  p <- p + theme(text = element_text(size=20))
  png(file.path(sim.input$output.folder, "Boxplot of logE.png"))
  print(p)
  dev.off()
}

# run JAGS for the synthetic inventory
# -------------------------------------
run.jags <- function(sim.input, inv.df, model_string) {
  V <- sim.input$n.inv
  P <- sim.input$n.pol
  N <- NROW(inv.df)
  Y <- inv.df$logE
  mean.logE.df <- ddply(inv.df, c('pollutant'), summarise, meanlogE = mean(logE))
  meanlogE <- mean.logE.df$meanlogE
  
  # jags input data
  jags.input.list <- list(Y=Y, V=V, P=P, meanlogE=meanlogE)
  
  # create a function that generates starting values
  sample.var.logE.df <- ddply(inv.df, c('pollutant'), summarise,
                              sample.var.logE = var(logE))
  # the variance of the activity is taken as half of the smallest variance on logE
  # > distribute variance equaly over activity and emission factor
  var.logAC.start <- min(sample.var.logE.df$sample.var.logE)/2
  tau.logAC.start <- 1/var.logAC.start
  
  # the rest of the variance of logE is given to the variance of the EFs
  var.logEF.start <- sample.var.logE.df$sample.var.logE - var.logAC.start
  tau.logEF.start <- 1/var.logEF.start
  
  # initial values for mulogE(pollutant)
  mulogE.start <- meanlogE
  
  create.initial.values <- function(){
    list("tau.logAC" = runif(1, min = tau.logAC.start*0.9, max=tau.logAC.start*1.1),
         "tau.logEF" = runif(4, min = tau.logEF.start*0.9, max=tau.logEF.start*1.1),
         "mulogE" = runif(4, min = mulogE.start - abs(mulogE.start)*0.9, 
                          max = mulogE.start + abs(mulogE.start)*1.1))
  }
  # run jags
  jagsfit <- jags(data = jags.input.list, n.iter=10000,
                  inits = create.initial.values,
                  parameters.to.save = c('mulogE', 'tau.logAC', 'tau.logEF'),
                  model.file = textConnection(model_string))
  
  return(jagsfit)  
}

# visualize jags output
# ---------------------
visualize.jagsfit <- function(sim.input, inv.df, jagsfit) {
  # What happens?
  # - output to dataframe
  # - trace plots
  # - output for mulogE vs meanlogE
  # - output of logAC, logEFs versus real values
  # - uncertainty ratios
  # - uncertainty ranking
  
  # make a data.frame
  # -----------------
  jagsfit.mcmc <- as.mcmc(jagsfit)
  output.df <- data.frame()
  for (i in 1:3) {
    chain.df <- cbind(chain = toString(i), iteration = 1:NROW(jagsfit.mcmc[[i]]), 
                      as.data.frame(jagsfit.mcmc[[i]]))
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
  
  # Make trace plots of all variables
  # ---------------------------------
  for (col.name in names(output.df)[3:n.col.output]) {
    p <- ggplot(data = output.df, aes_string(x="iteration", y=col.name, col="chain")) + geom_line()
    p <- p + theme(text = element_text(size=20))
    png(file.path(sim.input$output.folder, paste0("Trace_", col.name, ".png")), height = 480/2, width = 480*2)
    print(p)
    dev.off()
  }
  
  # Histograms for mulogE(pol) and compare to meanlogE(pol)
  # -------------------------------------------------------
  mean.logE.df <- ddply(inv.df, c('pollutant'), summarise, meanlogE = mean(logE))
  meanlogE <- mean.logE.df$meanlogE
  for (i.pol in  1:sim.input$n.pol) {
    mulogEP <- paste0('mulogE', i.pol)
    p <- ggplot(data = output.df, aes_string(x=mulogEP, col="chain")) + geom_density()
    p <- p + geom_vline(xintercept = meanlogE[i.pol])
    p <- p + theme(text = element_text(size=20))
    png(file.path(sim.input$output.folder, paste0("Hist_mulogE_pol", i.pol, ".png")))
    print(p)
    dev.off()
  }
  
  # Histogram of the sigmas and the real sigma of AC and EFs
  # --------------------------------------------------------
  tau.names <- c("tau.logAC", paste0("tau.logEF", 1:sim.input$n.pol))
  tau <- tau.names[2]
  sigma.df <- data.frame()
  for (tau in tau.names) {
    # get the real data from the simulaton input
    if (tau == "tau.logAC") {
      # Mind the relation between relative sigma X and sigma log(X)!!!
      real.log.sigma <- sqrt(log(sim.input$rel.sigma.AC^2 + 1))
      label <- "AC"
    } else {
      i.pol <- as.integer(gsub("tau.logEF", "", tau))
      real.log.sigma <- sqrt(log(sim.input$rel.sigma.EF[i.pol]^2 + 1))
      label <- paste0("EF", i.pol)
    }
    # credible interval
    sigma.log.var <- sqrt(1/output.df[,tau])
    credible.interval <- as.numeric(quantile(sigma.log.var, probs = c(0.025, 0.5, 0.975)))
    p <- ggplot(data = output.df, aes_string(x=paste0("sqrt(1/", tau, ")"), col="chain")) + geom_density()
    p <- p + geom_vline(xintercept = real.log.sigma)
    p <- p + geom_vline(xintercept = credible.interval, col = 'red')
    p <- p + labs(x = paste("sigma log", label), title = paste(sim.input$output.folder, label))
    p <- p + theme(text = element_text(size=20))
    png(file.path(sim.input$output.folder, paste0("Hist_sigmaLog", label, ".png")))
    print(p)
    dev.off()
    # data.frame for uncertainty ranking
    sigma.df <- rbind(sigma.df,
                      data.frame(iteration = output.df$iteration,
                                 chain = output.df$chain,
                                 varname = gsub("tau", "sigma",tau),
                                 sigma = sigma.log.var))
  }
  
  # Uncertainty ratios
  # ------------------
  var.logAC <- output.df[,"tau.logAC"]^-1
  CI.df <- data.frame()
  i.pol<-4
  output.folder <- sim.input$output.folder
  n.pol <- sim.input$n.pol
  for (i.pol in 1:sim.input$n.pol) {
    real.LogRatio <- log(sim.input$rel.sigma.EF[i.pol] / sim.input$rel.sigma.AC)
    var.logEF <- output.df[,paste0("tau.logEF", i.pol, "")]^-1
    sigma.EF.AC.ratio <- sqrt((exp(var.logEF)-1) / (exp(var.logAC)-1))
    sigma.EF.AC.ratio.df <- data.frame(chain = output.df$chain, sigma.EF.AC.ratio=sigma.EF.AC.ratio)
    credible.interval.Ratio <- quantile(sigma.EF.AC.ratio.df$sigma.EF.AC.ratio, probs = c(0.025, 0.5, 0.975))
    png(file.path(sim.input$output.folder, paste0("LogRatio_seAC_seEF_pol", i.pol, ".png")), width = 480, height = 480)
    p <- ggplot(data=sigma.EF.AC.ratio.df, aes(x=log(sigma.EF.AC.ratio), col=chain)) + geom_density()
    p <- p + geom_vline(xintercept = real.LogRatio)
    p <- p + geom_vline(xintercept = log(credible.interval.Ratio), col = 'red')
    p <- p + labs(title = paste0(sim.input$output.folder, ", Pol ", i.pol))
    p <- p + theme(text = element_text(size=20))
    print(p)
    dev.off()
    CI.pol.df <- as.data.frame(t(credible.interval.Ratio))
    CI.pol.df$pollutant <- paste0("Pol", i.pol)
    CI.df <- rbind(CI.df, CI.pol.df)
    write.table(credible.interval.Ratio, row.names = F, sep = ';',
                file = file.path(output.folder, paste0("credible_interval_Ratio_pol", i.pol, ".txt")))
    
  }
  
  # Uncertainty ranking
  # -------------------
  ranking.df <- ddply(sigma.df, c('varname'), summarise, mean.sigma = mean(sigma))
  ranking.df <- ranking.df[order(-ranking.df$mean.sigma),]
  ranking.df$prob.gt.next <- NA
  N.iter <- NROW(output.df)
  ranking.str <- ""
  for (i in 1:(NROW(ranking.df)-1)) {
    this.tau <- gsub("sigma", "tau", ranking.df$varname[i])
    next.tau <- gsub("sigma", "tau", ranking.df$varname[i+1])
    ranking.df$prob.gt.next[i] <- sum(output.df[,this.tau] < output.df[,next.tau]) / N.iter
    ranking.str <- paste0(ranking.str, gsub("sigma.log", "se", ranking.df$varname[i]), 
                          ">(", round(100*ranking.df$prob.gt.next[i]), "%)")
  }
  ranking.str <- paste0(ranking.str, gsub("sigma.log", "se", ranking.df$varname[NROW(ranking.df)]))
  png(file.path(sim.input$output.folder, paste0("Hist_Ranking.png")), width = 480*1.5, height = 480)
  p <- ggplot(data = sigma.df, aes(x=log(sigma), col=varname)) + geom_density()
  p <- p + labs(title = paste0("Ranking ", sim.input$output.folder, "\n", ranking.str))
  p <- p + theme(text = element_text(size=20))
  print(p)
  dev.off()
}

