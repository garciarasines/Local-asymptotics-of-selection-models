SIM = function(n, theta, theta0, t, B = 1e4, Nsim = 1e5){
  
  set.seed(123)
  
  simval = numeric(Nsim)
  
  for (i in 1:Nsim){
    
    cond = F
    while(!cond){
      y = rexp(n, rate = theta0)
      w = rnorm(1)
      cond = (sqrt(n)*mean(y) + w < t)
    }
    simval[i] = sqrt(n)*mean(y)
  }
  
  pvals_SLAN = numeric(B)
  pvals_exact = numeric(B)
  
  for (b in 1:B){
    
    cond = F
    while(!cond){
      y = rexp(n, rate = theta)
      w = rnorm(1)
      cond = (sqrt(n)*mean(y) + w < t)
    }
    
    z_obs = (theta0^2)*(1/theta0 - mean(y))*sqrt(n)
    a = 1/theta0
    c = (t - sqrt(n)/theta0)
    z_density = function(z) exp(dnorm(z/theta0, log = T) + pnorm(c + z/theta0^2, log.p = T) - pnorm(c/sqrt(1 + a^2), log.p = T))/theta0
    
    pvals_SLAN[b] = integrate(z_density, z_obs, 15)$value
    pvals_exact[b] = mean(sqrt(n)*mean(y) > simval)
    
  }
  
  if (theta == theta0){
    plot(ecdf(pvals_SLAN), col = 1, lwd = 2, main = paste("n =", n, ", null"), xlab = "p", ylab = "F(p)")
    lines(ecdf(pvals_exact),col = 8, lwd = 2)
    abline(a = 0, b = 1)
  } else {
    plot(ecdf(pvals_SLAN), col = 1, lwd = 2, main = paste("n =", n, ", alt"), xlab = "p", ylab = "F(p)")
    lines(ecdf(pvals_exact),col = 8, lwd = 2)
    abline(a = 0, b = 1)
  }
  
}

pdf("Exponential_pvals.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
SIM(n = 100, theta = 2, theta0 = 2, t = 3.57, B = 1e5)
SIM(n = 100, theta = 2.15, theta0 = 2, t = 3.57, B = 1e5)
dev.off()
