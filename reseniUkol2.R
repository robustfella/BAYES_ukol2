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

#8.19
#b - μ asi 200 mm, možná 140 az 260
# cca 95 % hodnot bude v rozmezi 2*smerodatne odchylky od stredni hodnoty
# 200 - 140 = 60 ... 2*sd ... sd = 30
plot_normal(200,30)

#c
data("penguins_bayes")
adelie<-subset(penguins_bayes,species=="Adelie" & !is.na(flipper_length_mm))
n<-nrow(adelie)
y<-mean(adelie$flipper_length_mm)

#d
sigma<-sd(adelie$flipper_length_mm)
summarize_normal_normal(200,30,sigma,y,n)
plot_normal_normal(200,30,sigma,y,n)
qnorm(c(0.025,0.975),y,0.5320898)

#8.20
#c
pnorm(220, y, sd = 0.5320898) - pnorm(200, y, sd = 0.5320898)

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
ci95 <- qgamma(c(0.025, 0.975),
                shape = alpha_post,
                rate  = beta_post)

print("95% posteriorní interval:")
print(ci95)

# 8.22b
ci95[2] < 1   # je horní mez pod 1?

# 8.22c
post_prob_lambda_lt_1 <- pgamma(1, shape = alpha_post, rate = beta_post)
post_prob_lambda_lt_1

# 8.23a
library(rstan)

stan_code <- "
data {
  int<lower=1> n;
  int<lower=0> y[n];
}
parameters {
  real<lower=0> lambda;
}
model {
  // prior
  lambda ~ gamma(4, 2);      // shape = 4, rate = 2
  
  // likelihood
  y ~ poisson(lambda);
}
"

# připrava dat pro Stan
y <- loons$count_per_100  
stan_data <- list(
  n = length(y),
  y = as.integer(y)
)

# spuštění MCMC
fit_loons <- stan(model_code = stan_code,
                  data  = stan_data,
                  chains = 4,
                  iter   = 10000,
                  warmup = 2000,
                  seed   = 123)

# 8.23b

library(bayesplot)
library(broom.mixed)

# základní souhrn
print(fit_loons, pars = "lambda")

# trace plot
mcmc_trace(as.array(fit_loons), pars = "lambda")

# autocorrelation
mcmc_acf_bar(as.array(fit_loons), pars = "lambda")

# 8.23c

post_draws <- extract(fit_loons, pars = "lambda")$lambda

ci_95_mcmc <- quantile(post_draws, probs = c(0.025, 0.975))
ci_95_mcmc
ci_95

# 8.23d
post_prob_lambda_lt_1_mcmc <- mean(post_draws < 1)
post_prob_lambda_lt_1_mcmc
post_prob_lambda_lt_1

# 8.24a

set.seed(123)
y_rep <- rpois(length(post_draws), lambda = post_draws)

# histogram posteriorně prediktivního rozdělení
library(ggplot2)

ggplot(data.frame(y_rep = y_rep),
       aes(x = y_rep)) +
  geom_histogram(binwidth = 1, boundary = -0.5) +
  labs(x = "Počet loonů v dalším 100-hodinovém období",
       y = "Posteriorní hustota (počty simulací)",
       title = "Posteriorně prediktivní rozdělení Y'")

# 8.24b

pred_summary <- c(
  mean = mean(y_rep),
  sd   = sd(y_rep),
  q025 = quantile(y_rep, 0.025),
  q975 = quantile(y_rep, 0.975)
)
pred_summary

# 8.24c
prob_Yrep_eq_0 <- mean(y_rep == 0)
prob_Yrep_eq_0
