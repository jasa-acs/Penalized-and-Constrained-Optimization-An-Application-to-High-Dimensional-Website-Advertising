### This file contains the code used to generate the results shown in the ACS author form for the minimal working example
###   in "Penalized and Constrained Optimization: An Application to High-Dimensional Website Advertising".
### Unfortunately, due to the proprietary nature of the comScore data set used in the case study, the authors
###   cannot make the complete data available to reproduce all the results from the paper.
### Instead, we have provided this example of the code used, with two example data sets generated to be similar
###   to the format of the actual comScore data:
###         (1) Page_View_Matrix_Example.csv: a matrix of example page views at 500 Internet websites 
###         (2) 500_Site_Info_Example.csv: the corresponding website information (e.g., CPM cost to advertise) for the 500 websites

### The code has been commented throughout to show the steps used in the analysis. Please see the accompanying
###   README file or the specific functions in PaC_Functions_Revised.R for any specific function details.

## First, read in the source file with the functions for the reproducibility code:
source("PaC_Functions_Revised.R")

## Next, read in the two example data sets, to be used in creating the reproducible working example:
z.500 = read.csv("Page_View_Matrix_Example.csv", header = T)
info.500 = read.csv("500_Site_Info_Example.csv", header = T)

## Isolate the information from the full data to give a vector of CPMs (cost per one thousand impressions) at
##    each website (c.500), a set of tau (total website visits) values (tau.500), a set of clickthrough rates
##    (q.500), and replace any missing values in the page view matrix with 0 values (since this information will
##    not help in predicting viewership/reach/clickthrough)

c.500 = info.500[,3]
tau.500 = info.500[,4]
q.500 = info.500[,5]
z.500[is.na(z.500)] <- 0

# Calculate the gamma values for each site, as detailed in the PaC paper
gamma.500 = 1/(c.500*tau.500)

### First, assume no clickthrough (just reach)
###   The ELMSO function calculates the reach based on the Paulson et al. (2018) paper
###   The compare.to.ELMSO function calculates for PaC to compare to Paulson et al.

trial.ELMSO = ELMSO(z.500[,-1], gamma = gamma.500, step = 0.02, size = 250, C.full = NULL, b = NULL)
trial.PAC = compare.to.ELMSO(z.500[,-1], gamma = gamma.500, step = 0.02, size = 250, C.full = NULL, b = NULL)

### Repeat the function calls above, but for click rate
trial.ELMSO.CTR = ELMSO(z.500[,-1], gamma = gamma.500, step = 0.03, size = 250, C.full = NULL, b = NULL, q = q.500)
trial.PAC.CTR = compare.to.ELMSO(z.500[,-1], gamma = gamma.500, step = 0.03, size = 250, C.full = NULL, b = NULL, q = q.500)

## Extract budget points for each lambda (saved in w.sum) for each of the four models
B.ELMSO = trial.ELMSO$w.sum
B.PAC = trial.PAC$w.sum
B.ELMSO.CTR = trial.ELMSO.CTR$w.sum
B.PAC.CTR = trial.PAC.CTR$w.sum

### Note that all budget vectors here have same length, though different values
###   Create a vector for each model to save the reach and click values
###   In addition, we will also calculate cost-adjusted and naive reach/click rate for comparison

reach.pac = reach.elmso = ctr.pac = ctr.elmso = reach.cost = ctr.cost = reach.naive = ctr.naive = rep(0,250)

## Use the reach.calc function in PaC_Functions_Revised.R to calculate reach (or click rate) for each
##    model, cost-adjusted, and naive setting
## Note to speed up code, we have skipped all runs where the budget is 0, since then reach and click rate
##    is also necessarily 0

for(i in 1:250){
  
  if(trial.PAC$w.sum[i]!=0){reach.pac[i] = reach.calc(trial.PAC$w[,i],z.500[,-1], gamma.500)}
  if(trial.PAC$w.sum[i]!=0){reach.naive[i] = reach.calc(rep(trial.PAC$w.sum[i]/500,500), z.500[,-1], gamma.500)}
  if(trial.PAC$w.sum[i]!=0){reach.cost[i] = reach.calc(trial.PAC$w.sum[i]*(rev(c.500/sum(c.500))), z.500[,-1], gamma.500)}
  if(trial.ELMSO$w.sum[i]!=0){reach.elmso[i] = reach.calc(trial.ELMSO$w[,i], z.500[,-1], gamma.500)}
  
  if(trial.PAC$w.sum[i]!=0){ctr.pac[i] = reach.calc(trial.PAC$w[,i], z.500[,-1], gamma.500, q = q.500)}
  if(trial.PAC$w.sum[i]!=0){ctr.naive[i] = reach.calc(rep(trial.PAC$w.sum[i]/500,500), z.500[,-1], gamma.500, q = q.500)}
  if(trial.PAC$w.sum[i]!=0){ctr.cost[i] = reach.calc(trial.PAC$w.sum[i]*(rev(c.500/sum(c.500))), z.500[,-1], gamma.500, q = q.500)}
  if(trial.ELMSO$w.sum[i]!=0){ctr.elmso[i] = reach.calc(trial.ELMSO$w[,i], z.500[,-1], gamma.500, q = q.500)}
  
}

## Plot results to get the reproducible minimal example from the ACS file. Note that although all models have
##    the same length of run, they do not produce the same budgets due to differences in calculation, so the 
##    indexing is adjusted in the plotting code to plot equivalent budget levels

par(mfrow = c(1,2), xpd = TRUE)

plot(B.ELMSO[1:250]/1000000, reach.elmso[1:250],type = "l", xlab = "Budget (in millions)", ylab = "Reach", col = "blue")
lines(B.PAC[1:203]/1000000, reach.pac[1:203], col = "black")
lines(B.PAC[1:203]/1000000, reach.cost[1:203], col = "purple")
lines(B.PAC[1:203]/1000000, reach.naive[1:203], col = "green")

plot(B.PAC[1:203]/1000000, ctr.pac[1:203], type = "l", xlab = "Budget (in millions)", ylab = "CTR", col = "black")
lines(B.ELMSO[1:250]/1000000, ctr.elmso[1:250], col = "blue")
lines(B.PAC[1:203]/1000000, ctr.cost[1:203], col = "purple")
lines(B.PAC[1:203]/1000000, ctr.naive[1:203], col = "green")

