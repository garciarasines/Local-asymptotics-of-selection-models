library(mcmc)
library(mvtnorm)

# Parameters
t = 2
mu0 = 0
sigma20 = 1
n1 = 200
n2 = 100
n = n1 + n2
gamma = n1/n

# Number of samples for MC computation of phi_star
M = 5e4

# Number of evaluations of phi_star for interpolation
phi_samples = 1e4

# Prior
pi = function(mu, sigma2) 1/sqrt(sigma2)

# True selection probability
phi = function(mu, sigma2) {
  log_phi_v = function(v) pnorm(sqrt(n1/sigma2)*(mu - t*sqrt(v/n1)), log.p = T)
  integrand = Vectorize(function(v) exp(dchisq(n1*v/sigma2, df = n1-1, log = T) + log(n1/sigma2) + log_phi_v(v)))
  integrate(integrand, 1e-4, sigma2 + 5/sqrt(n1),stop.on.error = FALSE)$value
}

# True posterior 
log_post = function(h){
  
  mu = mu0 + h[1]/sqrt(n)
  sigma2 = sigma20 + h[2]/sqrt(n)
  
  if (sigma2 < 1e-3) return(-Inf)
  else if (phi(mu, sigma2) < 1e-50) return(-Inf)
  else{
    log(pi(mu, sigma2)) - 0.5*n*log(sigma2) - (0.5/sigma2)*(n1*(v1 + (mu - y1)^2) + n2*(v2 + (mu - y2)^2)) - log(phi(mu, sigma2)) 
  }
  
}

# SLAN approximation
FI_inv = matrix(c(sigma20, 0, 0, 2*sigma20^2), nrow = 2)
y_bar = function(z) mu0 + z[1]/sqrt(n)
sigma_hat = function(z) sqrt(max(sigma20 + z[2]/sqrt(n) - z[1]^2/n, 1e-6))
test = function(y) mean(y) - t*sqrt(n1-1)*sd(y)/n1 > 0

set.seed(123)
# Samples_1: Conditional samples of Y given y_bar and sigma_hat are sigma_hat*Samples_1 + y_bar
Samples_1 = sqrt(n/(n-1))*scale(matrix(rnorm(n*M), nrow = n))

# Samples_2: Samples of Z's for MC estimation of varphi(h) are set as Z = h + FI_inv^(1/2)*Samples_2
Samples_2 = matrix(rnorm(2*M), nrow = 2) 
phi_star = function(h){
  tests = numeric(M)
  for (j in 1:M){
    Sample = y_bar(h + sqrt(FI_inv)%*%Samples_2[ , j]) + sigma_hat(h + sqrt(FI_inv)%*%Samples_2[ , j])*Samples_1[ , j]
    tests[j] = test(Sample[1:n1])
  }
  mean(tests)
}

# Evaluate phi_star in "phi_samples" points and average values at 2 nearest neighbours at other points
h1s = seq(-8, 8, length.out = floor(sqrt(phi_samples)))
h2s = seq(-8, 8, length.out = floor(sqrt(phi_samples)))
hs = expand.grid(h1s, h2s)
phi_app_samples = apply(hs, 1, phi_star)
phi_star_approx = function(h) {
  distances = (hs[ , 1] - h[1])^2 + (hs[ , 2] - h[2])^2
  sorted_indices = order(distances)
  mean(phi_app_samples[sorted_indices[1:2]])
}
log_post_star = function(h){
  zed1=(n1*y1+n2*y2-n*mu0)/sqrt(n)
  zed2=(-n*sigma20+n1*(v1 + (mu0 - y1)^2) + n2*(v2 + (mu0 - y2)^2))/sqrt(n)
  
  if (sigma2 < 1e-3) return(-Inf)
  else if (phi_star_approx(h) < 1e-50) return(-Inf)
  else{
    dnorm(sqrt((1/sigma20)*(h[1]-zed1)^2+(1/(2*sigma20^2))*(h[2]-zed2)^2),log=T) - log(phi_star_approx(h)) 
  }
}

# Simulation
B = 1e3

B = 3

MCMC_samples = 2e4

res_true_mu = numeric(B)
res_app_mu = numeric(B)
res_true_sigma2 = numeric(B)
res_app_sigma2 = numeric(B)
covmu = numeric(B)
covsigma = numeric(B)
set.seed(1123)
i = 1
while (i < B + 1){
  print(i) 
  # Data generation
  cond = F
  while(cond == F){
    v1 = (sigma20/n1)*rchisq(1, df = n1-1)
    v2 = (sigma20/n2)*rchisq(1, df = n2-1)
    y1 = rnorm(1, mean = mu0, sd = sqrt(sigma20/n1))
    y2 = rnorm(1, mean = mu0, sd = sqrt(sigma20/n2))
    cond = (y1 > t*sqrt(v1/n1))
  }
  y = gamma*y1 + (1 + gamma)*y2
  v = (n1*(v1 + (y - y1)^2) + n2*(v2 + (y - y2)^2))/n
  
  # MCMC
  MCMC = tryCatch(metrop(log_post, c(0, 0), nbatch = MCMC_samples, scale = 3.5), error = function(e) "error")
  MCMC_app = tryCatch(metrop(log_post_star_ay, c(0, 0), nbatch = MCMC_samples, scale = 3.0), error = function(e) "error")
  
  # If there is an error, repeat with new sample. Otherwise, store results
  if (length(MCMC) == 1 | length(MCMC_app) == 1){ i = i }
  else {
    res_true_mu[i] = mean(MCMC$batch[ , 1] < 0)
    res_true_sigma2[i] = mean(MCMC$batch[ , 2] < 0)
    res_app_mu[i] = mean(MCMC_app$batch[ , 1] < 0)
    res_app_sigma2[i] = mean(MCMC_app$batch[ , 2] < 0)
    sortedmu  = sort(MCMC$batch[,1])
    sortedsigma = sort(MCMC$batch[,2])
    uppermu=sortedmu[MCMC_samples*0.95]
    lowermu=sortedmu[MCMC_samples*0.05]
    uppersigma=sortedsigma[MCMC_samples*0.95]
    lowersigma=sortedsigma[MCMC_samples*0.05]
    covmu[i]=sum(MCMC_app$batch[ , 1] < uppermu & MCMC_app$batch[ , 1] > lowermu)/MCMC_samples
    covsigma[i]=sum(MCMC_app$batch[ , 2] < uppersigma & MCMC_app$batch[ , 2] > lowersigma)/MCMC_samples
    print(c(covmu[i],covsigma[i]))
    par(mfrow = c(1, 2))
    
    i = i + 1
  }
}

print(c(t, n1, n2, B, mean(covmu), mean(covsigma)))

# FIGURE 5: This plots data from iteration i=1 of the above code.

par(mfrow = c(1, 2))
plot(density(MCMC$batch[ , 1], bw = 0.3), xlab = "h1", main = "Posterior h1", ylim = c(0,0.35), col = "blue")
lines(density(MCMC_app$batch[ , 1], bw = 0.4))
abline(v = -0.3114,lty = 3)
abline(v = 3.8097, lty = 3)
plot(density(MCMC$batch[ , 2], bw = 0.3), xlab = "h2", main = "Posterior h2", ylim = c(0,0.35), col="blue")
lines(density(MCMC_app$batch [ ,2], bw = 0.4))
abline(v = -2.1232, lty = 3)
abline(v = 2.7154, lty = 3)
