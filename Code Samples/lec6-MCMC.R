####################################################################################################################################
# A Gibbs sampler example
# Assume X_ij ~ N(\theta_i, \sigma^2_i), where X_ij denotes the differential expression of gene i in subject j between two condidtions (say two different tissues)
# Prior: \theta_i \sim N(\mu_\theta, \tau); \sigma^2_i \sim IG(a_\sigma, b_\sigma)
# Data: n genes in k samples. So X = [X_ij]_(n by k)
###################################################################################################################################



############
set.seed(1234)
library(pscl)  # loading pscl package, to use inverse gamma distribution


n <- 50
k <- 50



########## Simulate Gene expression data from a  model a #########
sim.data.sens.model <- function(theta, sd.theta){
        x <- NULL
    
    for(i in 1:n){
        x <- rbind(x, rnorm(k, theta[i], sd.theta[i] ) )
    }
    #x <- t(x)
    return(x)
}
################################################################################################


###### Data matrix x is nxk, each column a gene, each row a sample #####


theta.tr <- sample(c(1,-1,0),n,replace=TRUE) ## the true mean expression is centered at -1, 0, or 1, only three values
##theta.tr
sigma.tr <- c(rep(1^2, n/2), rep(2^2, n/2))  ## the true var expression is 1 for the first half (n/2) genes and 4 for the second half
##sigma.tr


x <- sim.data.sens.model(theta=theta.tr, sd.theta=sqrt(sigma.tr))


#### MCMC Simulation ####
xbar <- apply(x, 1, mean) ## xbar is the row mean of x, which is \sum_j x_{ij} / k

N.sim <- 1000 ## N.sim is the number of iterations in the MCMC



## Prior and hyperparameters a0, b0, and var of delta
a.sigma <- 2; b.sigma <- 5### This gives an inverse gamma prior for sigma^2_i with mean 5 and variance infinity. IG mean b/(a-1); var b^2/(a-1)^2/(a-2) for a>2
mean.theta <- 0
tau <- 10^2


a.tau <- 2; b.tau <- 5 ### This gives an inverse gamma prior for tau with mean 5 and variance infinity. IG mean b/(a-1); var b^2/(a-1)^2/(a-2) for a>2


#### Initialize  #####

sim.sigma <- matrix(rep(1, n*N.sim), ncol=N.sim)
sim.theta <- matrix(rep(0, n*N.sim), ncol=N.sim)




par(mfrow=c(3,3)) ### for plotting


for(NN in 1:(N.sim-1)){

   ##### Assume sim.sigma is first initialized (=1), MCMC proceeds by sampling theta_i, sigma^2_i, mean.theta, tau  ######
    
    ## sample all the theta_i's
    for(i in  1:n){
        w <- ( sim.sigma[i, NN] / k ) / ( sim.sigma[i, NN] / k + tau )
        curr.draw <- rnorm(n=1, mean = w * mean.theta + (1-w) * xbar[i], sd = sqrt(w * tau) )
        sim.theta[i, NN+1] <- curr.draw
    }
    
    ## sample all the sigma^2_i's
    for(i in 1:n){
        ss.x <- sum( (x[i, ] - sim.theta[i, NN+1])^2 )
        curr.draw <- rigamma(n=1, alpha= a.sigma + 0.5 * k , beta= b.sigma + 0.5 * ss.x)
        sim.sigma[i, NN+1] <- curr.draw
    }
    

}

##for(i in 1:n){
plot(1:N.sim, sim.theta[1,], "l")  ## trace plot
abline(theta.tr[1], 0, col=2)

plot(1:N.sim, sim.sigma[1,], "l") ## trace plot
abline(sigma.tr[1], 0, col=2)

hist(apply(sim.theta, 1, mean)-theta.tr)  ## Histogram of the bias between posterior mean and the truth for each gene



###################################
#### Questions
# Q1: which genes are differentially expressed defined as theta_i > 0? See below for code. df.pp contains the pr(theta_i>0 |data)

# df.pp <- NULL
# for(i in 1:n){
#     df.pp <- c(df.pp, mean(sim.theta[i,]>0))
# }
# plot(df.pp)


# Q2: If we want to get more accurate estimates of the gene expression theta_i's, which number should we increase in the experiment?

# k, the sample size

# Q3: If we assume \theta_i \iid~ N(theta_0, \tau), i.e., the prior mean  E(theta_i) = theta_0 and variance var(theta_i) = \tau. Assume theta_0 ~ N(mean.theta0=0, var.theta0=10^2) and \tau ~ Inverse Gamma distribution IG(a.tau, b.tau). Generate MCMC samples now including theta_0 and \tau





###### Data matrix x is nxk, each column a gene, each row a sample #####


theta.tr <- sample(c(1,-1,0),n,replace=TRUE) ## the true mean expression is centered at -1, 0, or 1, only three values
##theta.tr
sigma.tr <- c(rep(1^2, n/2), rep(2^2, n/2))  ## the true var expression is 1 for the first half (n/2) genes and 4 for the second half
##sigma.tr


x <- sim.data.sens.model(theta=theta.tr, sd.theta=sqrt(sigma.tr))


#### MCMC Simulation ####
xbar <- apply(x, 1, mean) ## xbar is the column mean of x, which is \sum_j x_{ij} / k

N.sim <- 1000 ## N.sim is the number of iterations in the MCMC



## Prior and hyperparameters a0, b0, and var of delta
a.sigma <- 2; b.sigma <- 5### This gives an inverse gamma prior for sigma^2_i with mean 5 and variance infinity. IG mean b/(a-1); var b^2/(a-1)^2/(a-2) for a>2
mean.theta <- 0
var.theta <- 10^2


a.tau <- 2; b.tau <- 5 ### This gives an inverse gamma prior for tau with mean 5 and variance infinity. IG mean b/(a-1); var b^2/(a-1)^2/(a-2) for a>2


#### Initialize  #####

sim.sigma <- matrix(rep(1, n*N.sim), ncol=N.sim)
sim.theta <- matrix(rep(0, n*N.sim), ncol=N.sim)


mean.theta0 <- 0; var.theta0 <- 10^2

sim.theta0 <- rep(0, N.sim)
sim.tau <- rep(1, N.sim)



par(mfrow=c(2,3)) ### for plotting
par(mar=c(1,1,1,1))


for(NN in 1:(N.sim-1)){
    
    ##### Assume sim.sigma is first initialized (=1), MCMC proceeds by sampling theta_i, sigma^2_i, mean.theta, tau  ######
    
    ## sample all the theta_i's
    for(i in  1:n){
        w <- ( sim.sigma[i, NN] / k ) / ( sim.sigma[i, NN] / k + sim.tau[NN] )
        curr.draw <- rnorm(n=1, mean = w * sim.theta0[NN] + (1-w) * xbar[i], sd = sqrt(w * sim.tau[NN]) )
        sim.theta[i, NN+1] <- curr.draw
    }
    
    ## sample all the sigma^2_i's
    for(i in 1:n){
        ss.x <- sum( (x[i, ] - sim.theta[i, NN+1])^2 )
        curr.draw <- rigamma(n=1, alpha= a.sigma + 0.5 * k , beta= b.sigma + 0.5 * ss.x)
        sim.sigma[i, NN+1] <- curr.draw
    }
    
    
    ## sample theta0
    theta.bar <- mean(sim.theta[,NN+1])
    w0 <- (sim.tau[NN]/n) / (sim.tau[NN]/n + var.theta0)
    sim.theta0[NN+1] <- rnorm(n=1, mean = w0 * mean.theta0 + (1-w0) * theta.bar, sd = sqrt(w0 * var.theta0) )
    
    ## sample tau
    sim.tau[NN+1] <- rigamma(n=1, alpha = a.tau + 0.5 * n, b.tau + 0.5* sum( (sim.theta[,NN+1] - sim.theta0[NN+1])^2 ) )
    
}

##for(i in 1:n){
plot(1:N.sim, sim.theta[1,], "l")
abline(theta.tr[1], 0, col=2)

plot(1:N.sim, sim.sigma[1,], "l")
abline(sigma.tr[1], 0, col=2)

plot(1:N.sim, sim.theta0, "l")
abline(mean(theta.tr), 0, col=2)

plot(1:N.sim, sim.tau, "l")
abline(var(theta.tr), 0, col=2)

hist(apply(sim.theta, 1, mean)-theta.tr)



