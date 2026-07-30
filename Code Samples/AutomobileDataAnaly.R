###################################################
### Conjugate Prior/Hierarchical Model (Automobile data)
### Originally written by David Hitchcock 
### Yuan Ji, 2021 for teaching purpose
###################################################

autodata <- read.table("C:/Users/rewin/OneDrive/Documents/STAT_35920/Code Samples/autoregresslarge.txt", header=T)

y <- autodata$mpg  ### response is a car's mile per gallon (mpg)

### three covariates: engine displacement; horsepower; and car weight
x1 <- autodata$displacement 
x2 <- autodata$horsepower 
x3 <- autodata$weight

X <- as.matrix( cbind(rep(1, length(x1)), x1, x2, x3) )

###### Setting up prior specification :

# 3 predictor variables, so we need exactly 3+1=4 hypothetical "prior observations"

# Suppose that based on "expert opinion" we have the following guesses:

xtilde.obs.1 <- 
c(150, 100, 2000)
ytilde.obs.1 <- 25  
# Prior belief:  a car with displacement=150, horsepower=100, weight=2000 should have mpg around 25

xtilde.obs.2 <- 
c(200, 160, 4500)
ytilde.obs.2 <- 10  
# Prior belief:  a car with displacement=200, horsepower=160, weight=4500 should have mpg around 10

xtilde.obs.3 <- 
c(250, 140, 3000)
ytilde.obs.3 <- 20
# Prior belief:  a car with displacement=250, horsepower=140, weight=3000 should have mpg around 20

xtilde.obs.4 <- 
c(100, 80, 1800)
ytilde.obs.4 <- 35
# Prior belief:  a car with displacement=100, horsepower=80, weight=1800 should have mpg around 35

# Making the matrix X~ :
prior.obs.stacked <- rbind(xtilde.obs.1, xtilde.obs.2, xtilde.obs.3, xtilde.obs.4)
xtilde <- cbind(rep(1, times=nrow(prior.obs.stacked) ), prior.obs.stacked)

# Making the vector Y~ :
ytilde <- c(ytilde.obs.1, ytilde.obs.2, ytilde.obs.3, ytilde.obs.4)

# Diagonal matrix D contains weights that indicate how much worth we place on 
# our hypothetical prior observations (note the weights could vary if we are more
# confident about some hypothetical prior observations than about others):

D <- diag(c(4, 4, 4, 4))

# This D yields a D^{-1} that is diag(c(1/4, 1/4, 1/4, 1/4)), which
# assumes that our four hypothetical prior observations TOGETHER are
# "worth" about one actual sample observation.

# The prior mean on the beta vector that we are inducing is:

prior.mean.beta <- as.vector(solve(xtilde) %*% ytilde)
print("Prior mean for beta vector was:"); print(round(prior.mean.beta,3))


### Choosing prior parameters a and b for gamma prior on tau:

# The parameter "a" reflects our confidence in our prior:
# Let's choose "a" to be 0.5, which implies that our prior precision is "low"

a <- 0.5

# Consider prior observation 1:
# Prior belief:  a car with displacement=150, horsepower=100, weight=2000 should have mpg around 25
# Suppose the expert believes the largest realistic mpg for this type of car is 35

# Then set the prior guess for sigma to be (35 - 25)/1.645 = 6.08
# So the prior guess for sigma^2 is (6.08^2) = 36.97
# So the prior guess for tau is 1/36.97 = 0.027

# Since the gamma prior mean is a/b, let
# b = a / tau.guess:

tau.guess <- 0.027
b <- a / tau.guess 

# Here b is around 18.5 ...

##### Posterior information for tau and beta:

# Calculate beta-hat :

betahat <- as.vector(solve(t(X) %*% X + t(xtilde)%*%solve(D)%*%xtilde) %*% (t(X)%*%y + t(xtilde)%*%solve(D)%*%ytilde))


# Calculate s* :

n <- length(y)

sstar <- as.numeric( (t(y-X%*%betahat)%*%(y-X%*%betahat) + t(ytilde-xtilde%*%betahat)%*%solve(D)%*%(ytilde-xtilde%*%betahat) + 2*b)/(n+2*a) )


### Point estimates for tau (and thus for sigma^2):

p.mean.tau <- 1 / sstar
p.mean.sig.sq <- 1 / p.mean.tau

p.median.tau <- qgamma(0.50, shape=(n+2*a)/2, rate=((n+2*a)/2)*sstar)
p.median.sig.sq <- 1 / p.median.tau

print(paste("posterior.mean for sigma^2=", round(p.mean.sig.sq,3), 
      "posterior.median for sigma^2=", round(p.median.sig.sq,3) ))

### Marginal Interval estimate for tau (and sigma^2):

library(TeachingDemos) # loading TeachingDemos package, to use hpd function

hpd.95.tau <- hpd(qgamma, shape=(n+2*a)/2, rate=((n+2*a)/2)*sstar )

hpd.95.sig.sq <- 1 / hpd.95.tau

print("95% posterior interval for sigma^2:")

round(sort(hpd.95.sig.sq), 3)



### A NUMERICAL APPROACH TO GET POINT ESTIMATES FOR BETA:

# Randomly generate many tau values from its posterior distribution:

how.many <-30000

tau.values <- rgamma(n=how.many, shape=(n+2*a)/2, rate=((n+2*a)/2)*sstar )

library(mvtnorm)

beta.values<- matrix(0,nr = how.many, nc=length(betahat))
for (j in 1:how.many){
beta.values[j,] <- rmvnorm(n=1, mean=betahat, sigma= (1/tau.values[j])*solve(t(X)%*%X + t(xtilde)%*%solve(D)%*%xtilde) )
}

# Posterior median for each regression coefficient:

post.medians<-apply(beta.values,2,median)

# 95% posterior interval for each regression coefficient:

post.lower <- apply(beta.values,2,quantile, probs=0.025)
post.upper <- apply(beta.values,2,quantile, probs=0.975)

## Summarizing:

names.predictors <- c("x1", "x2", "x3")
beta.post.summary <- data.frame(cbind(post.lower, post.medians, post.upper), row.names=c("intercept", names.predictors))
names(beta.post.summary) <- c('0.025 Quantile', '0.5 Quantile', '0.975 Quantile')
print(beta.post.summary)


