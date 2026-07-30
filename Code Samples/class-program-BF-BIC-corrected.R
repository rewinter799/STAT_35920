### Lecture 3: Using Bayes Factor for Bayesian hypothesis testing #####


###############################################
# Comparing Two Means:  Bayes Factor Approach #
###############################################

# Data for blood pressure reduction on 21 subjects
# (10 patients took a blood pressure medication (calcium supplement), 11 patients took placebo):

calcium <- c(7,-4,18,17,-3,-5,1,10,11,-2)
placebo <- c(-1,12,-1,-3,3,-5,5,2,-11,-1,-3)

x1 <- calcium
x2 <- placebo

## making a boxplot of the two vectors. Does there appear to be a difference?  ### 
cp.data<-data.frame(cbind(x1, x2))
boxplot(cp.data)

###### Prepare to compute the Bayes factor ######
xbar.1 <- mean(x1); xbar.2 <- mean(x2);
s2.1 <- var(x1); s2.2 <- var(x2);
n1 <- length(x1);  n2 <- length(x2);


s2.p <- ((n1-1)*s2.1 + (n2-1)*s2.2)/(n1+n2-2)

n.alt <-   1 / (1/n1 + 1/n2)

# Usual two-sample t-statistic:
t.star <- (xbar.1 - xbar.2) / sqrt(s2.p/n.alt)

# prior mean and variance for Delta = |mu1 - mu2|/sigma:

mean.Delta <- 0  
var.Delta <- 1/9 ## "objective" prior information

# mean.Delta <- 0.5  
# var.Delta <- 1/9 ## stronger prior belief that the two population means differ

my.pv <- sqrt(1 + n.alt * var.Delta) # scale parameter for noncentral-t
my.ncp <- mean.Delta * sqrt(n.alt) / my.pv # noncentrality parameter for noncentral-t

# Calculating Bayes Factor:

B <- dt(t.star, df = n1+n2-2, ncp = 0) / ( dt(t.star/my.pv, df = n1+n2-2, ncp = my.ncp)/my.pv )

print(B)

# Probability that H_0 is true, given equal prior probabilities for H_0 and H_a:

1 / (1 + (1/B))


##################################
# BICs for the blood pressure Example:  #
##################################

# Data:

calcium <- c(7,-4,18,17,-3,-5,1,10,11,-2)
placebo <- c(-1,12,-1,-3,3,-5,5,2,-11,-1,-3)

my.x <- c(calcium,placebo) # combining both samples into a single vector

# defining a grouping variable:
my.group <- factor(c(rep(1,times=length(calcium)), rep(2,times=length(placebo))) )

# Model 1:  Assuming equal means:

M.1 <- lm(my.x ~ 1);
p1 <- 1 ## there is only one free parameter, \mu_1 = \mu_2 = \mu

# Model 2:  Assuming different means:

M.2 <- lm(my.x ~ my.group);
p2 <- 2 ## there are two free parameters, \mu_1 and \mu_2

# BICs (uses AIC code with k = ln(n) ):

AIC(M.1, M.2, k= log(length(my.x)))
BIC(M.1, M.2)

# Model 1 is SLIGHTLY favored...
# The sample size is MUCH TOO SMALL here to view the difference in BICs as 
# as an approximation to a Bayes Factor.


###############################
# BICs in regression setting: #
###############################

### Oxygen uptake data from measuring people participating different training groups ###


y <- c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)
# y = change in maximal oxygen uptake (positive number = improvement)
x1 <- c(0,0,0,0,0,0,1,1,1,1,1,1)              # 0 = running group, 1 = aerobics group
x2 <- c(23,22,22,25,27,20,31,23,27,28,22,24)  # ages of subjects
x3 <- x1*x2                                   # cross-product term (interaction b/w group and age)

X <- as.matrix( cbind(rep(1, length(x1)), x1, x2, x3) )

# Model 1:  Assuming no useful predictors:
M.1 <- lm(y ~ 1)

# Model 2:  Assuming only group as a useful predictor:
M.2 <- lm(y ~ x1)

# Model 3:  Assuming only age as a useful predictor:
M.3 <- lm(y ~ x2)

# Model 4:  Assuming age and group as useful predictors:
M.4 <- lm(y ~ x1 + x2)

# Model 5:  Assuming age, group, and age-by-group interaction as useful predictors:
M.5 <- lm(y ~ x1 + x2 + x3)

# BICs (uses AIC code with k = ln(n) ):

AIC(M.1, M.2, M.3, M.4, M.5, k= log(nrow(X)) )

# Model 4 is favored, with Model 5 second-best.
#


#################################################
# Comparing Two Means:  Gibbs Sampling Approach #
#################################################

# Data: blood pressure (medicine is just calcium)

calcium <- c(7,-4,18,17,-3,-5,1,10,11,-2)
placebo <- c(-1,12,-1,-3,3,-5,5,2,-11,-1,-3)

x1 <- calcium
x2 <- placebo

boxplot(x1, x2) ## the boxplots of the two vectors

xbar.1 <- mean(x1); xbar.2 <- mean(x2);
s2.1 <- var(x1); s2.2 <- var(x2);
n1 <- length(x1);  n2 <- length(x2);

###

S <- 5000 # number of Gibbs iterations

# Prior parameters for mu: 

mean.mu <- 2  # expect slight overall mean reduction?
var.mu <- 50  # chosen quite wide

# Prior parameters for tau: 

mean.tau <- 0  # not favoring either group a priori
var.tau <- 50  # chosen quite wide again

# Prior parameters for sigma^2: 

nu1 <- 1  
nu2 <- 20

## Vectors to hold values of mu, tau, and sigma^2:
mu.vec <- rep(0, times=S)
tau.vec <- rep(0, times=S)
sigma.sq.vec <- rep(0, times=S)

## Starting values for mu and tau:

mu.vec[1] <- 2
sigma.sq.vec[1] <- 1

## Starting Gibbs algorithm:

for (s in 2:S){

  
# update tau: 

var.tau.cond <- 1/(1/var.tau+((4*n1*n2)/((n1+n2)*sigma.sq.vec[s-1])))
mean.tau.cond <- var.tau.cond*(mean.tau/var.tau+0.5*(mean(x1)-mean(x2))/((n1+n2)*sigma.sq.vec[s-1]/(4*n1*n2)))
tau.vec[s] <- rnorm(1,mean.tau.cond,sqrt(var.tau.cond))

# update sigma.sq:

sigma.sq.vec[s] <- 1/rgamma(1,(nu1+n1+n2)/2,(nu1*nu2+sum((x1-mu.vec[s-1]-tau.vec[s])^2)+
                   sum((x2-mu.vec[s-1]+tau.vec[s])^2))/2)


# update mu:

var.mu.cond <- 1/(1/var.mu+(n1+n2)/sigma.sq.vec[s])
mean.mu.cond <- var.mu.cond*(mean.mu/var.mu+sum(x1-tau.vec[s])/sigma.sq.vec[s] + 
                sum(x2+tau.vec[s])/sigma.sq.vec[s])
mu.vec[s] <- rnorm(1,mean.mu.cond,sqrt(var.mu.cond))

}
## End Gibbs algorithm.

# MCMC diagnostics for mu, tau, sigma^2

par(mfrow=c(2,2))  # set up 3-by-1 array of plots
plot(sigma.sq.vec, type='l')
plot(mu.vec, type='l')
plot(tau.vec, type='l')

#windows()  # open new plotting window
par(mfrow=c(2,2))
acf(sigma.sq.vec)
acf(mu.vec)
acf(tau.vec)

par(mfrow=c(1,1))  # reset window

# Plotting estimated posterior density for mu (after deleting first 100 values as burn-in):

plot(density(mu.vec[-(1:100)]))

#windows()  # open new plotting window

# Plotting estimated posterior density for tau (after deleting first 100 values as burn-in):

plot(density(tau.vec[-(1:100)]))

# Posterior probability that mean for first group is higher than mean for second group:

mean( tau.vec[-(1:100)] > 0)


# Sampling from the posterior predictive distribution:

sample.1.pred <- rnorm(4900, mean = (mu.vec[-(1:100)] + tau.vec[-(1:100)]), sd=sqrt(sigma.sq.vec[-(1:100)]) )

sample.2.pred <- rnorm(4900, mean = (mu.vec[-(1:100)] - tau.vec[-(1:100)]), sd=sqrt(sigma.sq.vec[-(1:100)]) )


# Posterior predictive probability that a random subject from the first group has a greater 
# blood pressure reduction than a random subject from the second group:

mean(sample.1.pred > sample.2.pred)




