
################################
# Monte Carlo and MCMC Methods #
################################
library(TeachingDemos)
library(pscl)




####################################################
##### Sparrow data
####################################################

library(mvtnorm)
# May need to install the mvtnorm package first?
# If so, type at the command line:  install.packages("mvtnorm", dependencies=T)
# while plugged in to the internet.

sparrow.data <- read.table("sparrowdata.txt", header=T)

y <- sparrow.data[,1]
x <- sparrow.data[,2]; xsq <- x^2

X <- cbind(rep(1,times=length(x)), x, xsq)

p <- dim(X)[2]  # number of columns of matrix X

beta.prior.mean <- rep(0, times=p)
beta.prior.sd <- rep(10,times=p)

proposal.cov.matrix <- var(log(y+1/2))*solve(t(X)%*%X)  
# should be reasonable in this case:   should be close to (sigma^2)*(X'X)^{-1}

S <- 10000
beta.current <- rep(0,times=p) # initial value for M-H algorithm
acs <- 0 # will be to track "acceptance rate"

beta.values <- matrix(0,nrow=S,ncol=p)  # will store sampled values of beta vector

for (s in 1:S) {
 beta.proposed <- t(rmvnorm(1, beta.current, proposal.cov.matrix))

 log.accept.ratio <- sum(dpois(y,exp(X %*% beta.proposed), log=T)) - sum(dpois(y,exp(X %*% beta.current), log=T)) +
        sum(dnorm(beta.proposed, beta.prior.mean, beta.prior.sd, log=T)) - 
        sum(dnorm(beta.current, beta.prior.mean, beta.prior.sd, log=T) )

  if (log.accept.ratio > log(runif(1)) ) {
    beta.current <- beta.proposed
    acs <- acs + 1
  }

beta.values[s,] <- beta.current
}

acs/S  # gives the acceptance rate

acf(beta.values[,1])  # plot autocorrelation values for beta_0
acf(beta.values[,2])  # plot autocorrelation values for beta_1
acf(beta.values[,3])  # plot autocorrelation values for beta_2

# Seems to be an issue with serial dependence.

# Thinning out the sampled values by taking every 10th row:

beta.values.thin <- beta.values[10*(1:(S/10) ),]

# Looks much better...

### Posterior summary:

apply(beta.values.thin,2,median)  # Posterior medians for each regression coefficent

cbind(apply(beta.values.thin, 2, quantile, probs=0.025),   apply(beta.values.thin, 2, quantile, probs=0.975) )
# approximate 0.025 and 0.975 quantiles for each regression coefficent
# to get approximate 95% quantile-based intervals

rbind( emp.hpd(beta.values.thin[,1], conf=0.95), 
emp.hpd(beta.values.thin[,2], conf=0.95), 
emp.hpd(beta.values.thin[,3], conf=0.95) ) # approximate 95% HPD intervals


# Plots of posterior median for expected offspring for ages 1,2,3,4,5,6:

my.X <- cbind( rep(1,times=6), (1:6), (1:6)^2 )

y.hats <- exp( my.X %*% as.matrix( apply(beta.values.thin,2,median), ncol=1 ) )

plot( (1:6), y.hats, type='b', xlab='age', ylab='expected offspring')

# Trace plots for the sampled beta_0 and beta_1 values:

plot(beta.values.thin[,1], type='l')

plot(beta.values.thin[,2], type='l')

plot(beta.values.thin[,3], type='l')





