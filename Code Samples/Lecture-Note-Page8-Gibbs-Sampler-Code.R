
X <- c(1,1,1,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,1)  ## responses from 19 participants
Y <- sum(X)

n.sim <- 10000
n.burn.in <- 1000
th <- rep(0, n.sim)
x20 <- rep(0, n.sim)

## Intialize
x20[1] <- 1
th[1] <- 0.5

## Gibbs
for(j in 2:n.sim){
  th[j] <- rbeta(1, Y+x20[j-1]+1, 20-Y-x20[j-1]+1)
  x20[j] <- rbinom(n=1, size=1, th[j])
}

### Posterior sample ####
x20.final <- x20[(n.burn.in+1):n.sim]
th.final <- th[(n.burn.in+1):n.sim]

plot(1:n.sim, th, "l")
abline(mean(th.final),0,col=2)
acf(th.final)
plot(density(th.final))

mean((th.final>0.5))

####
mean(th.final) ## posteior mean of theta
mean(x20.final) ## pr(X20 = 1 | y)
