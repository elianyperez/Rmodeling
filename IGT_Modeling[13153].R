setwd ("G:/My Drive/EP_Year1/IGT analysis/Rmodeling")


#to change settings for r to run go to https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
source("hBayesDM_model.R")
library(devtools)
library(hBayesDM)
library(rstan)

######## M1 ####################
startTime <- Sys.time()

output<-igt_pvl_delta(data = "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/igt_data2.txt", niter = 4000, nwarmup = 1000,
              nchain = 4, ncore = 1, nthin = 1, inits = "vb",
              indPars = "mean", modelRegressor = FALSE, vb = FALSE,
              inc_postpred = FALSE, adapt_delta = 0.95, stepsize = 1,
              max_treedepth = 10)

deck5<-igt_pvl_delta(data = "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/igt_deck5.txt",niter = 4000, nwarmup = 1000,
                      nchain = 4, ncore = 1, nthin = 1, inits = "vb",
                      indPars = "mean", modelRegressor = FALSE, vb = FALSE,
                      inc_postpred = FALSE, adapt_delta = 0.95, stepsize = 1,
                      max_treedepth = 10)
  # Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
  plot(output, type = "trace")
  
  # Check Rhat values (all Rhat values should be less than or equal to 1.1)
  rhat(output) #checked and good for this model!
  
  # Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
  plot(output)
  
  # Show the WAIC and LOOIC model fit estimates
  printFit(output)
  
  plotInd(output)


  
  ## After model fitting is complete for both groups,
  ## evaluate the group difference (e.g., on the 'pi' parameter) by examining the posterior distribution of group mean differences.
  
  diffDist = output$parVals$mu_lambda - deck5$parVals$mu_lambda  # group1 - group2 
  HDIofMCMC( diffDist )  # Compute the 95% Highest Density Interval (HDI). 
  plotHDI( diffDist )  
  

  library("bayesplot")
  library("rstanarm")
  library("ggplot2")
  posterior <- as.matrix(output$fit)
  
  plot_title <- ggtitle("Posterior distributions",
                        "with medians and 80% intervals")
  mcmc_areas(posterior,
             pars = c("sigma[1]", "A[1]", "alpha[1]"),
             prob = 0.8) + plot_title

  color_scheme_set("darkgray")
  mcmc_scatter(
    as.matrix(output$fit),
    pars = c("cons", "A"),
    np = nuts_params(output$fit),
    np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
  )
  color_scheme_set("darkgray") 
  fit<-as.data.frame(posterior)
  ggplot(fit, aes(mu_cons, mu_A)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)
  
  stan_hist(output$fit)
  (cons_vs_A <- stan_scat(output$fit, pars = c("mu_cons", "mu_lambda"), color = "blue", size = 4))
  cons_vs_A +
    ggplot2::coord_flip() +
    theme(panel.background = ggplot2::element_rect(fill = "black"))
  
  
  setwd("G:/My Drive/EP_Year1/IGT analysis/Rmodeling")
  cons<-as.data.frame(output$allIndPars)
  write.csv(cons, "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/allindpars.csv", row.names=FALSE)
  stan_hist(cons)
  allindpars<-read.csv("G:/My Drive/EP_Year1/IGT analysis/Rmodeling/allindpars.csv")
  Etotal<-read.csv( "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/Etotal.csv")
  merged<-merge(allindpars,Etotal,by.x="ID")
  
  
  library(ggplot2)
  library(ggpubr)
  theme_set(
    theme_minimal() +
      theme(legend.position = "top")
  )
  b <- ggplot(merged, aes(x = E_choice, y = cons))
  # Scatter plot with regression line
  b + geom_point()+
    geom_smooth(method = "lm") 
  
 
  
  
  
  
  
  p<-ggplot(Etotal, aes(x = E_choice, y = Total)) +
    geom_point(aes(color = E_choice))+
    geom_smooth(method=lm,formula=y~x,   # Add linear regression line
                se=TRUE,col="black", size=1)+
    theme(legend.position = "right")+
    labs(title = "Variability in explore deck", 
         x = "Explore deck clicks", y = "Total")
  p+theme(text=element_text(family="Futura"))
######## M2 ####################


output2<-igt_orl(data = "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/igt_data2.txt", niter = 4000, nwarmup = 1000, nchain = 4,
  ncore = 1, nthin = 1, inits = "vb", indPars = "mean",
  modelRegressor = FALSE, vb = FALSE, inc_postpred = FALSE,
  adapt_delta = 0.95, stepsize = 1, max_treedepth = 10)

# Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
plot(output2, type = "trace")

# Check Rhat values (all Rhat values should be less than or equal to 1.1)
rhat(output2)

# Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
plot(output2)

# Show the WAIC and LOOIC model fit estimates
printFit(output2)


######## M3 ###################

output3<-igt_pvl_decay(data = "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/igt_data2.txt", niter = 4000, nwarmup = 1000,
  nchain = 4, ncore = 1, nthin = 1, inits = "vb",
  indPars = "mean", modelRegressor = FALSE, vb = FALSE,
  inc_postpred = FALSE, adapt_delta = 0.95, stepsize = 1,
  max_treedepth = 10)

# Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
plot(output3, type = "trace")

# Check Rhat values (all Rhat values should be less than or equal to 1.1)
rhat(output3)

# Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
plot(output3)

# Show the WAIC and LOOIC model fit estimates
printFit(output3)



######## M4 ####################


output4<-igt_vpp(data = "G:/My Drive/EP_Year1/IGT analysis/Rmodeling/igt_data2.txt", niter = 4000, nwarmup = 1000, nchain = 4,
  ncore = 4, nthin = 1, inits = "vb", indPars = "mean",
  modelRegressor = FALSE, vb = FALSE, inc_postpred = FALSE,
  adapt_delta = 0.95, stepsize = 1, max_treedepth = 10)

# Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
plot(output4, type = "trace")

# Check Rhat values (all Rhat values should be less than or equal to 1.1)
rhat(output4)

# Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
plot(output4)

# Show the WAIC and LOOIC model fit estimates
printFit(output4)


###Compare MODELS and GROUPS
printFit(output, output2, output3, output4)

output$allIndPars
output2$allIndPars
output3$allIndPars
output4$allIndPars

output$fit
output2$fit
output3$fit
output4$fit

output$parVals
output$modelRegressor

plotInd(output, "ep")


## fit example data with the gng_m3 model and run posterior predictive checks

output5<-igt_orl(data = "C:/Users/elian/Desktop/Rmodeling/igt_exampleData.txt", niter = 4000, nwarmup = 1000, nchain = 4,
                 ncore = 1, nthin = 1, inits = "vb", indPars = "mean",
                 modelRegressor = FALSE, vb = TRUE, inc_postpred = TRUE,
                 adapt_delta = 0.95, stepsize = 1, max_treedepth = 10)

## dimension of x$parVals$y_pred
dim(output5$parVals$y_pred)   # y_pred --> 4000 (MCMC samples) x 10 (subjects) x 240 (trials)


y_pred_mean = apply(output5$parVals$y_pred, c(2,3), mean)  # average of 4000 MCMC samples

dim(y_pred_mean)  # y_pred_mean --> 10 (subjects) x 240 (trials)

numSubjs = dim(output5$allIndPars)[1]  # number of subjects

subjList = unique(output5$rawdata$subjID)  # list of subject IDs
maxT = max(table(output5$rawdata$subjID))  # maximum number of trials
true_y = array(NA, c(numSubjs, maxT)) # true data (`true_y`)

## true data for each subject
for (i in 1:numSubjs) {
  tmpID = subjList[i]
  tmpData = subset(output5$rawdata, subjID == tmpID)
  true_y[i, ] = tmpData$choice  # only for data with a 'choice' column
}

## Subject #1
plot(true_y[1, ], type="l", xlab="Trial", ylab="Choice (0 or 1)", yaxt="n")
lines(y_pred_mean[1,], col="red", lty=2)
axis(side=2, at = c(0,1) )
legend("bottomleft", legend=c("True", "PPC"), col=c("black", "red"), lty=1:2)