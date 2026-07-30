
#################################
## Bayesian Linear Regression with Reversible Jump ###
## Assume y = mu + beta1* x1 + beta2* x1^2 + epsilon, where epsilon ~ N(0, sig). 
## We assume mu and beta1 follow proper prior
## We assume m=1, 2 is a binary random variable. P(m=1)=0.5
## We assume beta2 | m=1 ~ Normal prior N(0, sig.beta); and beta2 | m=2 ~ point mass at 0. That is, when model indicator m=2, beta2=0.
## We want to use reversible jump MCMC to sample all the parameters, including mu, beta1, beta2, and m. 
## Note that when m=2, we have one fewer parameter beta2 in the model. 
## Yuan Ji, 2020
#################################

set.seed(123)


### Simulate data assuming a true model
x1 <- seq(1,10, by=.05)
x1 <- x1/sd(x1)
x2 <- x1*x1                                  
sig <- sqrt(.5)
y <- rnorm(length(x1), mean = 1 + 2*x1 + 0.5*x2 , sd=sig) ## here the truth is that beta2 = 0.5. So 
#y <- rnorm(length(x1), mean = 1 + 2*x1 + 0*x2 , sd=sig)

##X <- as.matrix( cbind(rep(1, length(x1)), x1, x2) )

par(mfrow=c(3,3))
plot(x1, y)


## Initialize parameter values
n.mc <- 100000
m <- rep(1, n.mc)

h <- matrix(rep(0.5, 4), ncol=2, nrow=2)  ## the probability of transition between model m=1 and model m=2; with probability 0.5 the model can go from 1->1, 1->2, 2->1, and 2->2. Note that it must be true that h(1,1) + h(1,2) = h(2,1) + h(2,2) = 1, since from model 1 (or 2), the algorithm only allows either staying at model 1 (or 2) or jump to model 2 (or 1), respectively. 

beta01 <- matrix(0, ncol=2, nrow=n.mc)
beta2 <- rep(0, n.mc)


### Set up prior parameters: mu, beta1 ~ N(0, sig.beta), beta2 | m=2 ~ N(0, sig.beta), 
sig.u <- sqrt(1)
sig.beta <- sqrt(10) 

## Assume prior P(m=1) = pi.M = 0.5.
pi.M <- 0.5

#### tau is the standard deviation of proposal density for mu, beta1, and beta2
tau <- sqrt(0.05)


m[1] <- 1
beta01[1,] <- c(1,1)
##beta2[1] <- 0.5



### Function samp012 is to sample the full model with mu, beta1, beta2. Here, b01 = c(mu, beta1), and b2 = beta2 ###
samp012 <- function(b01, b2){

    curr <- c(b01, b2)
    epi <- rnorm(3, 0, sd=tau)
    prop <- curr + epi 
    
    like.ratio <- - sum(log(dnorm(y, mean = curr[1]+curr[2]*x1 + curr[3]*x2, sd=sig))) + sum(log(dnorm(y, prop[1]+ prop[2]*x1 + prop[3]*x2, sd=sig)))

    #cat("curr", curr, "prop", prop, "\n")
    #cat("like-ratio", like.ratio, "\n"); readline()
    
    prior.ratio <- - sum(log(dnorm(curr, 0, sig.beta))) + sum(log(dnorm(prop, 0, sig.beta))) 

    acc <- exp(like.ratio + prior.ratio)
    #cat("samp012 acc", acc, "\n"); #readline()

    ind <- (acc > runif(1))

#    res <- c(res, beta012.tmp * ind + c(beta01,beta2) * (1-ind))
    res <- prop * ind + curr * (1-ind)
    #cat("ind", ind, "res", res[i+1, ], "i", i, "\n"); #readline()
    #cat(res); #readline()
    return(res)
}




### Function samp01 is to sample the reduced model with mu and beta1 Here, b01 = c(mu, beta1) ###
samp01 <- function(b01){

    curr <- b01
    epi <- rnorm(2, 0, sd=tau)
    prop <- curr + epi 

    
    like.ratio <- - sum(log(dnorm(y, mean=curr[1]+curr[2]*x1, sd=sig ))) + sum(log(dnorm(y, mean=prop[1]+ prop[2]*x1, sd=sig )))

    prior.ratio <- - sum(log(dnorm(curr, 0, sig.beta))) + sum(log(dnorm(prop, 0, sig.beta)))

    acc <- exp(like.ratio + prior.ratio)

    #cat("samp01 acc", acc, "\n"); #readline()
    
    ind <- (acc > runif(1))
    
    res <- prop * ind + curr * (1-ind)
    ##cat(res)
    
    return(res)
}



#### Set up the RJMCMC

for(sim in 1:(n.mc-1)){

    if(m[sim]==1){

        chg.ind <- (runif(1) < h[1,2])  ## Set chg.ind to TRUE with probability h[1,2] (prob from model 1 to 2). If chg.ind is TRUE, propose a move, which will add beta2 to the model. The move needs to be accepted.
        
        if(chg.ind){

            u <- rnorm(1, 0, sd=sig.u)
            beta2.tmp <- u
            #cat("u", u, "\n"); readline()
            
            like.ratio <- - sum(-2 * beta2.tmp * x2 * (y - beta01[sim,1]- beta01[sim,2] * x1) + beta2.tmp^2 * x2^2) / 2 / sig^2
            prior.ratio <- - log(sqrt(2*pi) * sig.beta) - beta2.tmp^2 / 2 / sig.beta^2 ## + log((1-pi.M)) - log(pi.M)
            proposal <- log(sqrt(2*pi) * sig.u) + u^2 / 2 / sig.u^2

            acc <- exp(like.ratio + prior.ratio + proposal)

            ind <- (acc > runif(1))  ## If ind is TRUE, accept the move from model 1 --> 2. If ind is FALSE, do not move.

            if(ind){
                m[sim+1] <- 2
                beta2[sim+1] <- beta2.tmp
                beta01[sim+1, ] <- beta01[sim, ]
                #beta01[sim+1, ] <- c(1,1)
            }
            if(!ind){
                m[sim+1] <- m[sim]
                beta2[sim+1] <-  0
                beta01[sim+1, ] <- beta01[sim, ]
                #beta01[sim+1, ] <- c(1,1)
            }
        }

        if(!chg.ind){  ## if chg.ind is FALSE, do not propose model 2; resample mu and beta1.

            rest <- samp01(beta01[sim,])
            beta01[sim+1, ] <- rest
            beta2[sim+1] <- 0 #beta2[sim]
            m[sim+1] <- m[sim]
               # beta01[sim+1, ] <- c(1,1)
        }
    }

#    cat(sim, c(beta01[sim+1, ], beta2[sim+1]), "\n")
  
  ### The reversible move from model 2 to model 1.
    
    if(m[sim]==2){

        chg.ind <- (runif(1) < h[2,1]) ## Set chg.ind to TRUE with probability h[2,1] (prob from model 2 to 1). If chg.ind is TRUE, propose a move, which will eliminate beta2 from the model. The move needs to be accepted.
        
        if(chg.ind){

            u <- beta2[sim]
            beta2.tmp <- beta2[sim]
        
            like.ratio <- - sum(-2 * beta2.tmp * x2 * (y - beta01[sim,1]- beta01[sim,2] * x1) + beta2.tmp^2 * x2^2) / 2 / sig^2
            prior.ratio <- - log(sqrt(2*pi) * sig.beta) - beta2.tmp^2 / 2 / sig.beta^2
            proposal <- log(sqrt(2*pi) * sig.u) + u^2 / 2 / sig.u^2

            acc <- exp(-like.ratio - prior.ratio - proposal)

            ind <- (acc > runif(1))
            
            if(ind){
                m[sim+1] <- 1
                beta2[sim+1] <- 0
                beta01[sim+1, ] <- beta01[sim, ]
                #beta01[sim+1, ] <- c(1,1)
            }
            
            if(!ind){
                m[sim+1] <- m[sim]
                beta01[sim+1, ] <-beta01[sim, ]
                beta2[sim+1] <- beta2[sim]
                #beta01[sim+1, ] <- c(1,1)
            }
        }
        
        if(!chg.ind){
           rest <- samp012(beta01[sim,], beta2[sim])
           beta01[sim+1, ] <- rest[1:2]
           beta2[sim+1] <- rest[3]
           m[sim+1] <- m[sim]
                #beta01[sim+1, ] <- c(1,1)
        }
    }
#   cat(sim, c(beta01[sim+1, ], beta2[sim+1]), "\n")
    
    ##readline()   
               
}

plot(1:n.mc, beta01[,1], "l")

for(i in seq(1,n.mc, by=1000)){
    points(i, beta01[i,1], col=2, size=4)
}

plot(1:n.mc, beta01[,2], "l")
for(i in seq(1,n.mc, by=1000)){
    points(i, beta01[i,2], col=3)
}
plot(1:n.mc, beta2, "l")
for(i in seq(1,n.mc, by=1000)){
    points(i, beta2[i], col=4)
}

plot(y, mean(beta01[,1]) + mean(beta01[,2])*x1 + mean(beta2)*x2, "n")
for(i in 1:length(y)){
    points(y[i], mean(beta01[,1]) + mean(beta01[,2])*x1[i] + mean(beta2)*x2[i], col=2)
}

abline(0,1)
plot(y- ( mean(beta01[,1]) + mean(beta01[,2])*x1 + mean(beta2)*x2), ylim=c(-5,5))
abline(0,0)

print(mean(m))
