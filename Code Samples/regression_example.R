library(rjags)

# set up data
n <- 100
p <- 2
x1 <- rnorm(100)
x2 <- rnorm(100, sd = 0.2)
y <- rnorm(100, 3 + 0.5*x1 - 2*x2, sd = 1)

summary(lm(y ~ x1 + x2))

x <- cbind(rep(1, n), x1, x2)

# compile model
model <- jags.model(file = "regression_example.txt", 
                    data = list(y = y, x = x, n = n, p = p),
                    n.chains = 2)

# do MCMC
posterior_mcmc <- coda.samples(model, 
                               variable.names = c("beta", "sigma", 
                                                  "sigma.sq.inv"), 
                               n.iter = 1000)

# examine the posterior 
str(posterior_mcmc)
summary(posterior_mcmc)
posterior_draws <- do.call(rbind, posterior_mcmc)
head(posterior_draws)

