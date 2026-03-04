SIM = function(theta0, model, gamma, n, B = 2e4, llim = 0){
  
  # Progress bar
  pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
  
  # Selection threshold for exponential model
  t = 1/theta0
  
  # Log-prior for exponential model
  log_prior = function(theta) dgamma(theta, shape = 1, rate = 0.01, log = T)
  
  # Log-selection probability
  log_phi = function(theta){
    if (model == "normal"){
      return(pnorm(theta/sqrt(theta0^2 + gamma^2*theta0^4), log.p = T))
    }
    else if (model == "exponential"){
      integrand = Vectorize(function(w) exp(dnorm(w, sd = gamma/sqrt(n), log = T) + pgamma(t - w, shape = n, rate = n*theta, log.p = T)))
      return(log(integrate(integrand, -10*gamma/sqrt(n), 10*gamma/sqrt(n))$value))
    }
  }
  
  # Simulation
  set.seed(123)
  res = numeric(B)
  for (i in 1:B){
    
    # Sample data
    if (model == "normal"){
      cond = F
      while(cond == F){
        z = rnorm(1, sd = theta0)
        w = rnorm(1)
        cond = (z + theta0^2*gamma*w > 0)
      } 
    } else if (model == "exponential"){
      cond = F
      while(cond == F){
        y = mean(rexp(n, rate = theta0))
        w = gamma*rnorm(1)/sqrt(n)
        cond = (mean(y) + w < t)
      } 
    }
    
    # Compute posterior CDF at theta0
    if (model == "normal"){
      post_u = Vectorize(function(theta) exp(dnorm(z, mean = theta, sd = theta0, log = T) - log_phi(theta)))
      K = integrate(post_u, -300, 10)$value
      post = Vectorize(function(theta) post_u(theta)/K)
      res[i] = integrate(post, -300, 0)$value
    } else if (model == "exponential"){
      post_u = Vectorize(function(theta) exp(log_prior(theta) + dgamma(y, shape = n, rate = n*theta, log = T) - log_phi(theta)))
      K = integrate(post_u, llim, 5)$value
      post = Vectorize(function(theta) post_u(theta)/K)
      res[i] = integrate(post, llim, theta0)$value
    }
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  return(res)
  
}

# Results
res_normal_d_50 = SIM(theta0 = 2, gamma = 0.01, n = 50, model = "normal")
res_normal_r_50 = SIM(theta0 = 2, gamma = 1/5, n = 50, model = "normal")
res_exponential_d_50 = SIM(theta0 = 2, gamma = 0.01, n = 50, model = "exponential")
res_exponential_r_50 = SIM(theta0 = 2, gamma = 1/5, n = 50, model = "exponential")
res_normal_d_100 = SIM(theta0 = 2, gamma = 0.01, n = 100, model = "normal")
res_normal_r_100 = SIM(theta0 = 2, gamma = 1/5, n = 100, model = "normal")
res_exponential_d_100 = SIM(theta0 = 2, gamma = 0.01, n = 100, model = "exponential")
res_exponential_r_100 = SIM(theta0 = 2, gamma = 1/5, n = 100, model = "exponential")

# Plot
xs = seq(0, 1, length.out = 1e3)
par(mfrow = c(1, 2))
plot(xs, sapply(xs, function(x) mean(res_normal_d_50 < x)), type = "l", lty = 2, main = "", xlab = "", ylab = "", col = "black", xlim = c(0, 1)) 
lines(xs, sapply(xs, function(x) mean(res_normal_r_50 < x)), lty = 3, col = "black") 
lines(xs, sapply(xs, function(x) mean(res_exponential_d_50 < x)), lty = 2, col = "blue") 
lines(xs, sapply(xs, function(x) mean(res_exponential_r_50 < x)), lty = 3, col = "blue") 
abline(a = 0, b = 1)
plot(xs, sapply(xs, function(x) mean(res_normal_d_100 < x)), type = "l", lty = 2, main = "", xlab = "", ylab = "", col = "black") 
lines(xs, sapply(xs, function(x) mean(res_normal_r_100 < x)), lty = 3, col = "black") 
lines(xs, sapply(xs, function(x) mean(res_exponential_d_100 < x)), lty = 2, col = "blue") 
lines(xs, sapply(xs, function(x) mean(res_exponential_r_100 < x)), lty = 3, col = "blue") 
abline(a = 0, b = 1)
