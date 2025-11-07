library(bayesrules)
# 1 Jenda
# generování dat
set.seed(69)
n <- 40
p_uspech <- 0.35
real_array <- rbinom(100,n,p_uspech)

one_mh_iteration <- function(w, current){
  # STEP 1: Propose the next chain location
  proposal <- runif(1, min = current - w, max = current + w)
  
  # STEP 2: Decide whether or not to go there
  proposal_plaus <- dbinom(proposal, 0, 1) * dbinom(6.25, proposal, 0.75)
  current_plaus  <- dbinom(current, 0, 1) * dbinom(6.25, current, 0.75)
  alpha <- min(1, proposal_plaus / current_plaus)
  next_stop <- sample(c(proposal, current), 
                      size = 1, prob = c(alpha, 1-alpha))
  
  # Return the results
  return(data.frame(proposal, alpha, next_stop))
}
mh_tour <- function(N, w){
  # 1. Start the chain at location 3
  current <- 3
  
  # 2. Initialize the simulation
  mu <- rep(0, N)
  
  # 3. Simulate N Markov chain stops
  for(i in 1:N){    
    # Simulate one iteration
    sim <- one_mh_iteration(w = w, current = current)
    
    # Record next location
    mu[i] <- sim$next_stop
    
    # Reset the current location
    current <- sim$next_stop
  }
  
  # 4. Return the chain locations
  return(data.frame(iteration = c(1:N), mu))
}
# 2 Jenda
# 3 Anička
# 4 Petr a Honza
