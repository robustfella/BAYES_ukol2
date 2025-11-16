library(bayesrules)
library(tidyverse)
# 1 Jenda
# generování dat
set.seed(69)
n <- 30
p_uspech <- 0.35
real_array <- data.frame(mu = rbinom(1,n,p_uspech))
ggplot(data = real_array, aes(x=mu))+
  geom_histogram(aes(y=..density..),color="white", bins=15)
simulated_p <- mean(real_array$mu)/30

one_iteration <- function(a, b, current){
  # STEP 1: Propose the next chain location
  proposal <- rbeta(1, a, b)
  
  # STEP 2: Decide whether or not to go there
  proposal_plaus <- dbeta(proposal, 3, 1) * dbinom(11, 30, proposal)
  proposal_q     <- dbeta(proposal, a, b)
  current_plaus  <- dbeta(current, 3, 1) * dbinom(11, 30, current)
  current_q      <- dbeta(current, a, b)
  alpha <- min(1, proposal_plaus / current_plaus * current_q / proposal_q)
  next_stop <- sample(c(proposal, current), 
                      size = 1, prob = c(alpha, 1-alpha))
  
  return(data.frame(proposal, alpha, next_stop))
}
betabin_tour <- function(N, a, b){
  # 1. Start the chain at location 0.5
  current <- 0.5
  
  # 2. Initialize the simulation
  pi <- rep(0, N)
  
  # 3. Simulate N Markov chain stops
  for(i in 1:N){    
    # Simulate one iteration
    sim <- one_iteration(a = a, b = b, current = current)
    
    # Record next location
    pi[i] <- sim$next_stop
    
    # Reset the current location
    current <- sim$next_stop
  }
  
  # 4. Return the chain locations
  return(data.frame(iteration = c(1:N), pi))
}

a <- 3; b <- 1; y <- real_array$mu
post_a <- 3 + y
post_b <- 1 + n - y
betabin_sim <- betabin_tour(N = 2000, a = a, b = b)
posterior_p <- mean(betabin_sim$pi)
hist(betabin_sim$pi, breaks = 30, main = "Posterior of p", xlab = "p")

ggplot(betabin_sim, aes(x = iteration, y = pi)) + 
  geom_line()
ggplot(betabin_sim, aes(x = pi)) + 
  geom_histogram(aes(y = ..density..), color = "white") + 
  stat_function(fun = dbeta, args = list(post_a, post_b), color = "blue")

# 2 Jenda


# 3 Anička

# 4 Petr a Honza

library(bayesrules)
library(rstan)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(glue)

data(loons)

# 8.21

# 8.21c: Počet datových bodů
n_obs <- nrow(loons)

# průměrný počet loonů na 100 hodin
mean_count_100 <- mean(loons$count_per_100)

print(paste("8.21c - Počet datových bodů (n):", n_obs))
print(paste("8.21c - Průměrný počet loonů na 100 hodin:", mean_count_100))

# 8.21d

alpha0 <- 4
beta0  <- 2

y <- loons$count_per_100
n  <- length(y)
S  <- sum(y)

alpha_post <- alpha0 + S
beta_post  <- beta0 + n

# 95% posteriorní interval
ci_95 <- qgamma(c(0.025, 0.975),
                shape = alpha_post,
                rate  = beta_post)

print("95% posteriorní interval:")
print(ci_95)
