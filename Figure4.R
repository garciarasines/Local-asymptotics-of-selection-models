library(mcmc)

# p_star function 
p_star = function(z, mu, sigma2, t, n1, n2) {
  
  n = n1 + n2
  gamma1 = n1/n
  gamma2 = n2/n
  
  # observed (y, v) from z
  y = mu + z[1]/sqrt(n)
  v = sigma2 + z[2]/sqrt(n) - (z[1]/sqrt(n))^2
  
  # fixed terms (precompute once)
  logfy_y = dnorm(y, sd = 1/sqrt(n), log = TRUE)
  logfv_v = dchisq(n*v, df = n - 1, log = TRUE) + log(n)
  
  # helpers (vectorized in y1)
  y2 = function(y1) (y - gamma1*y1) / gamma2
  
  v2 = function(y1, v1) {
    y2v = y2(y1)
    (v - gamma1*v1 - gamma1*gamma2*(y1 - y2v)^2)/gamma2
  }
  
  logfy1 = function(y1) dnorm(y1, sd = 1/sqrt(n1), log = TRUE)
  logfy2 = function(y2v) dnorm(y2v, sd = 1/sqrt(n2), log = TRUE)
  
  logfv1 = function(v1) dchisq(n1*v1, df = n1 - 1, log = TRUE) + log(n1)
  logfv2 = function(v2v) dchisq(n2*v2v, df = n2 - 1, log = TRUE) + log(n2)
  
  # inner integrand: vector in y1 for a fixed v1
  inner = function(y1, v1) {
    
    y2v = y2(y1)
    v2v = v2(y1, v1)
    
    ok = is.finite(y1) & is.finite(v2v) & (v1 > 0) & (v2v > 0)
    out = numeric(length(y1))
    if (!any(ok)) return(out)
    
    out[ok] = exp(
      logfy1(y1[ok]) +
        logfy2(y2v[ok]) +
        logfv1(v1) +
        logfv2(v2v[ok]) -
        logfy_y -
        logfv_v
    )
    
    out
  }
  
  # g(v1) must accept vector v1 because integrate() passes vectors
  g = function(v1_vec) {
    
    sapply(v1_vec, function(v1) {
      
      lower_y1 = t*sqrt(v1 / n1)
      upper_y1 = y + 5*sqrt(sigma2/n1)
      
      integrate(
        function(y1) inner(y1, v1),
        lower = lower_y1,
        upper = upper_y1,
        rel.tol = 1e-6
      )$value
      
    })
    
  }
  
  integrate(g, lower = 0, upper = sigma2 + 5/sqrt(n1), rel.tol = 1e-6)$value
}

# MLEs
y_bar = function(y) mean(y)
v = function(y) mean((y - mean(y))^2)

# Sample data selectively
test = function(y, t, n1, n2) {
  mean(y[1:n1]) > t*sqrt(v(y[1:n1])/n1)
}
samp = function(mu, sigma2, t, n1, n2) {
  
  n = n1 + n2
  cond = FALSE
  
  while (!cond) {
    y = rnorm(n, mean = mu, sd = sqrt(sigma2))
    cond = test(y, t, n1, n2)
  }
  
  return(y)
}

# p-value
p_val = function(y, mu0, type, t, n1, n2, MCMC_samples) {
  
  n = n1 + n2
  sigma2_est = var(y)
  z1_obs = sqrt(n)*(mean(y) - mu0)
  
  if (n2 == 0) {
    
    if (type == "SLAN") {
      
      a = function(z1) z1^2 / sqrt(n) - sqrt(n) * sigma2_est
      b = function(z1) (sqrt(n) / t^2) * (z1 + sqrt(n) * mu0)^2 + z1^2 / sqrt(n) - sqrt(n) * sigma2_est
      
      integrand = function(z1) {
        dnorm(z1, sd = sqrt(sigma2_est)) * (
          exp(pnorm(b(z1), sd = sqrt(2) * sigma2_est, log.p = TRUE)) -
            exp(pnorm(a(z1), sd = sqrt(2) * sigma2_est, log.p = TRUE))
        )
      }
      
      p = integrate(integrand, z1_obs, 20, rel.tol = 1e-8)$value /
        integrate(integrand, -sqrt(n) * mu0, 20, rel.tol = 1e-8)$value
      
    } else if (type == "plugin") {
      
      den = pnorm(t * sqrt(v(y)) - mu0 * sqrt(n), sd = sqrt(sigma2_est), lower.tail = FALSE)
      num = pnorm(z1_obs, sd = sqrt(sigma2_est), lower.tail = FALSE)
      
      p = num / den
    }
    
  } else if (n2 > 0) {
    
    if (type == "SLAN") {
      
      log_den = function(z) {
        dnorm(z[1], sd = sqrt(sigma2_est), log = T) + dnorm(z[2], sd = 2*sigma2_est, log = T) + log(p_star(z, mu0, sigma2_est, t, n1, n2))
      }
      
      MCMC = metrop(log_den, c(0, 0), nbatch = MCMC_samples, scale = 1)
      
      z1_obs = sqrt(n)*(mean(y) - mu0)
      
      p = mean(MCMC$batch[ , 1] > z1_obs)
      
    }
  }
  
  p
  
}

# Simulation
SIM = function(mu, sigma2, mu0, type, t, n1, n2, B = 1e3, MCMC_samples = 5e3) {
  
  pvals = numeric(B)
  
  pb = txtProgressBar(min = 0, max = B, style = 3)
  on.exit(close(pb), add = TRUE)
  
  for (b in 1:B) {
    
    set.seed(b)
    y = samp(mu, sigma2, t, n1, n2)
    pvals[b] = p_val(y, mu0, type, t, n1, n2, MCMC_samples)

    plot(ecdf(pvals[1:b]))
    
    setTxtProgressBar(pb, b)
  }
   
  pvals
  
}

# Results
ps0 = SIM(mu = 0, sigma2 = 1, mu0 = 0, type = "SLAN", t = 1, n1 = 100, n2 = 50)
ps1 = SIM(mu = 0.15, sigma2 = 1, mu0 = 0, type = "SLAN", t = 1, n1 = 100, n2 = 50)

# Plots
pdf("Gaussian_pvals.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(ecdf(ps0), col = 4, main = "null", xlab = "p", ylab = "F(p)")
abline(a = 0, b = 1)
plot(ecdf(ps1), col = 4, main = "alt", xlab = "p", ylab = "F(p)")
abline(a = 0, b = 1)
dev.off()
