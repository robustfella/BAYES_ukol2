library(bayesrules)
library(tidyverse)
# 1 Jenda
# generování dat
set.seed(69)
n <- 30
p_uspech <- 0.35
real_array <- data.frame(mu = rbinom(5000,n,p_uspech))
ggplot(data = real_array, aes(x=mu))+
  geom_histogram(aes(y=..density..),color="white", bins=15)


one_mh_iteration <- function(w, current, y, n, a, b){
  # STEP 1: Propose the next chain location
  proposal <- runif(1, min = current - w, max = current + w)
  if(proposal <= 0 || proposal >= 1) proposal <- current
  
  # STEP 2: Decide whether or not to go there
  proposal_plaus <- dbeta(proposal, a, b) * dbinom(y, n, proposal)
  current_plaus  <- dbeta(current, a, b) * dbinom(y, n, current)
  alpha <- min(1, proposal_plaus / current_plaus)
  next_stop <- sample(c(proposal, current), 
                      size = 1, prob = c(alpha, 1-alpha))
  
  # Return the results
  return(data.frame(proposal, alpha, next_stop))
}
mh_tour <- function(N, w, y, n, a, b){
  # 1. Start the chain at location 3
  current <- 0.5
  chain <- numeric(N)
  
  # 3. Simulate N Markov chain stops
  for(i in 1:N){    
    # Simulate one iteration
    sim <- one_mh_iteration(w = w, current = current, y, n, a, b)
    
    # Record next location
    chain[i] <- sim$next_stop
    
    # Reset the current location
    current <- sim$next_stop
  }
  
  # 4. Return the chain locations
  return(data.frame(iteration = c(1:N), p=chain))
}

a <- 2; b <- 5; y <- real_array$mu
mh_chain <- mh_tour(N = 5000, w = 0.1, y = y, n = n, a = a, b = b)
mean(mh_chain$p)
hist(mh_chain$p, breaks = 30, main = "Posterior of p", xlab = "p")
# 2 Jenda
# 3 Anička
# 4 Petr a Honza
