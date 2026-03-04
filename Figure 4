Posterior_density = function(theta0, y, model, gamma, n, llim = 0){
  
  # z value
  z = theta0^2*sqrt(n)*(1/theta0 - y)
  
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
      integrand = Vectorize(function(w) exp(dnorm(w, sd = gamma/sqrt(n), log = T) + pgamma(t - w, shape = n, rate = max(1e-6, n*theta), log.p = T)))
      return(log(integrate(integrand, -10*gamma/sqrt(n), 10*gamma/sqrt(n))$value))
    }
  }
  
  # Compute posterior PDF
  if (model == "normal"){
    post_u = Vectorize(function(theta) exp(dnorm(z, mean = theta, sd = theta0, log = T) - log_phi(theta)))
    K = integrate(post_u, -300, 10)$value
    post = Vectorize(function(theta) post_u(theta)/K)
  } else if (model == "exponential"){
    post_u = Vectorize(function(theta) exp(log_prior(theta) + dgamma(y, shape = n, rate = max(1e-6, n*theta), log = T) - log_phi(theta)))
    K = integrate(post_u, llim, 5)$value
    post = Vectorize(function(h) post_u(theta0 + h/sqrt(n))/(sqrt(n)*K))
  }
  
  return(sapply(grid, post))
  
}

# Plot
grid = seq(-18, 10, length.out = 1e4)
par(mfrow = c(1, 2))
n = 50
plot(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 0.01, n = n, model = "normal"), type = "l", lty = 1, main = "", xlab = "", ylab = "", col = "black", ylim = c(0, 0.2), lwd = 1.2) 
lines(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 0.01, n = n, model = "exponential"), lty = 1, col = "blue", lwd = 1.2) 
lines(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 1, n = n, model = "exponential"), lty = 1, col = "red", lwd = 1.2) 
plot(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 1/5, n = n, model = "normal"), type = "l", lty = 1, main = "", xlab = "", ylab = "", col = "black", ylim = c(0, 0.2), lwd = 1.2) 
lines(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 1/5, n = n, model = "exponential"), lty = 1, col = "blue", lwd = 1.2) 
lines(grid, Posterior_density(theta0 = 2, y = 0.45, gamma = 1, n = n, model = "exponential"), lty = 1, col = "red", lwd = 1.2) 
