library(bayesrules)
library(tidyverse)
#====================================================================
# 1 Jenda
#====================================================================
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
  current <- 0.3
  
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
#====================================================================
# 2 Jenda
#====================================================================
one_mh_iteration_normal <- function(w, current){
  # STEP 1: Propose the next chain location
  proposal <- runif(1, min = current - w, max = current + w)
  
  # STEP 2: Decide whether or not to go there
  proposal_plaus <- dnorm(proposal, 70, 10) * dnorm(77.88, proposal, 0.77)
  current_plaus  <- dnorm(current, 70, 10) * dnorm(77.88, current, 0.77)
  alpha <- min(1, proposal_plaus / current_plaus)
  next_stop <- sample(c(proposal, current), 
                      size = 1, prob = c(alpha, 1-alpha))
  
  # Return the results
  return(data.frame(proposal, alpha, next_stop))
}
mh_tour_normal <- function(N, w){
  # 1. Start the chain at location 3
  current <- 70
  
  # 2. Initialize the simulation
  mu <- rep(0, N)
  
  # 3. Simulate N Markov chain stops
  for(i in 1:N){    
    # Simulate one iteration
    sim <- one_mh_iteration_normal(w = w, current = current)
    
    # Record next location
    mu[i] <- sim$next_stop
    
    # Reset the current location
    current <- sim$next_stop
  }
  
  # 4. Return the chain locations
  return(data.frame(iteration = c(1:N), mu))
}
data("airquality")
y <- airquality$Temp
y <- na.omit(y)

n <- length(y)
y_bar <- mean(y)
s <- sd(y)
se <- s/sqrt(n)
sigma <- s^2

normalnormal_model = summarize_normal_normal(mean = 70, sd = 10, sigma = s, n=n,y_bar = y_bar)

normalnormal_chain <- mh_tour_normal(N = 2000, w = 1)
posterior_p <- mean(normalnormal_chain$mu)
hist(normalnormal_chain$mu, breaks = 30, main = "Posterior of μ", xlab = "μ")

ggplot(normalnormal_chain, aes(x = iteration, y = mu)) + 
  geom_line()
ggplot(normalnormal_chain, aes(x = mu)) + 
  geom_histogram(aes(y = ..density..), color = "white") + 
  stat_function(fun = dnorm, args = list(mean = normalnormal_model$mean[2], sd = normalnormal_model$sd[2]), color = "blue")
# 3 Anička
# 4 Petr a Honza
