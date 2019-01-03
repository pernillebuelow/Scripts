

#count of SS neurons for WT and KO
SS_counts <- c(1, 14)

#total number of cells for WT and KO
total_counts <- c(23, 39)

# proportions of SS for WT and KO respectively 
SS_counts/total_counts

#performing 2-sample test to test if proportions are significantly different
prop.test(x=SS_counts, n=total_counts)

## Next step: calculating power for proportions with unequal Ns
library(pwr)

#first calculating effect size 
h = ES.h(0.04, 0.36)

#then calculating the power of the analysis done
pwr.2p2n.test(h = -0.88, n1 = 23, n2 = 39, sig.level = 0.05)

library(MonteCarlo)

ttest<-function(n, loc, scale){
  sample<- rnorm(n, loc, scale)
  stat <- sqrt(n)*mean(sample)/sd(sample)
  desicion<-abs(stat)>1.96
  return(list("decision"=decision))
}

n_grid <-c(1, 22, 14, 25)
loc_grid<-seq(0, 1, 0.2)
scale_grid<-c(1,2)

param_list=list("n"=n_grid, "loc"=loc_grid, "scale"=scale_grid)

MC_result<-MonteCarlo(func=ttest, nrep=1000, param_list=param_list)

library(datasets)
x<-rnorm(10, 20, 2)
x
summary(x)

