# Parameters
t = 1
mu0 = 0
sigma20 = 1
n1 = 20
n2 = 20
n = n1 + n2
gamma1 = n1/n
gamma2 = n2/n
M = 5e4 

# Selection function

# Write (y, v) in terms of (z1, z2) for a true parameter (mu0, sigma20)
z_inv_1 = function(z) mu0 + z[1]/sqrt(n)
z_inv_2 = function(z) sigma20 + z[2]/sqrt(n) - (z_inv_1(z) - mu0)^2  

# Log-densities of y_bar, y_bar^1, y_bar^2
logfy = function(y) dnorm(y, sd = 1/sqrt(n), log = T)
logfy1 = function(y1) dnorm(y1, sd = 1/sqrt(n1), log = T)
logfy2 = function(y2) dnorm(y2, sd = 1/sqrt(n2), log = T)

# Log-densities of v, v1, v2
logfv = function(v) dchisq(n*v, df = n - 1, log = T) + log(n)
logfv1 = function(v1) dchisq(n1*v1, df = n1 - 1, log = T) + log(n1)
logfv2 = function(v2) dchisq(n2*v2, df = n2 - 1, log = T) + log(n2)

y2 = function(y, y1) (y - gamma1*y1)/gamma2
v2 = function(y, v, y1, v1) (v - gamma1*v1 - gamma1*gamma2*(y1 - y2(y, y1))^2)/gamma2

# Conditional density of (y1, v1)|(y, v)
fcond = function(y1, v1, y, v) exp(logfy1(y1) + logfy2(y2(y, y1)) + logfv1(v1) + logfv2(v2(y, v, y1, v1)) - logfy(y) - logfv(v))

# Selection function. Selection condition is y1 > t*sqrt(v1/n1)
p_star = function(z){
  y = z_inv_1(z)
  v = z_inv_2(z)
  conditional_density = function(y1, v1) fcond(y1, v1, y, v)
  g = Vectorize(function(v1) integrate(Vectorize(function(y1) conditional_density(y1, v1)),
                                       lower = t*sqrt(v1/n1), upper = y + 5*sqrt(sigma20/n1))$value)
  integrate(g, 0, sigma20 + 5/sqrt(n1))$value
} 

#  True selection probability 
varphi = function(h) {
  mu = mu0 + h[1]/sqrt(n)
  sigma2 = sigma20 + h[2]/sqrt(n)
  log_phi_v = function(v) pnorm(sqrt(n1/sigma2)*(mu - t*sqrt(v/n1)), log.p = T)
  integrand = Vectorize(function(v) exp(dchisq(n1*v/sigma2, df = n1-1, log = T) + log(n1/sigma2) + log_phi_v(v)))
  integrate(integrand, 1e-4, sigma2 + 5/sqrt(n1), stop.on.error = FALSE)$value
}

varphi_ratio = function(h) varphi(h)/varphi(c(0, 0))

#  Asymptotic selection probability 

# Selection probability
FI_inv = matrix(c(sigma20, 0, 0, 2*sigma20^2), nrow = 2)
y_bar = function(z) mu0 + z[1]/sqrt(n)
sigma_hat = function(z) sqrt(max(sigma20 + z[2]/sqrt(n) - z[1]^2/n, 1e-6))
test = function(y) mean(y) - t*sqrt(n1-1)*sd(y)/n1 > 0
# Approximation of phi*
set.seed(123)
# Samples_1: Conditional samples of Y given y_bar and sigma_hat are sigma_hat*Samples_1 + y_bar
Samples_1 = sqrt(n/(n-1))*scale(matrix(rnorm(n*M), nrow = n))
# Samples_2: Samples of Z's for MC estimation of varphi(h) are set as Z = h + FI_inv^(1/2)*Samples_2
Samples_2 = matrix(rnorm(2*M), nrow = 2) 
varphi_star = function(h){
  tests = numeric(M)
  for (j in 1:M){
    Sample = y_bar(h + sqrt(FI_inv)%*%Samples_2[ , j]) + sigma_hat(h + sqrt(FI_inv)%*%Samples_2[ , j])*Samples_1[ , j]
    tests[j] = test(Sample[1:n1])
  }
  mean(tests)
}

varphi_star0 = varphi_star(c(0, 0))
varphi_star_ratio = function(h) varphi_star(h)/varphi_star0


# Plot
par(mfrow = c(2, 2))

z1s = seq(-1, 5, length.out = 1e3)
z2 = 0
plot(z1s, sapply(z1s, function(z1) p_star(c(z1, z2))), xlab = "z1", ylab = "p*", type = "l", lty = 2, lwd = 1.5)

z2s = seq(-1, 5, length.out = 1e3)
z1 = 1
plot(z2s, sapply(z2s, function(z2) p_star(c(z1, z2))), xlab = "z2", ylab = "p*", type = "l", lty = 2, lwd = 1.5)

h1s = seq(-1, 5, length.out = 50)
h2 = 0
plot(h1s, log(sapply(h1s, function(h1) varphi_ratio(c(h1, h2)))), xlab = "h1", ylab = "r/r*", type = "l", lty = 1, lwd = 1, col = "blue")
lines(h1s, log(sapply(h1s, function(h1) varphi_star_ratio(c(h1, h2)))), lwd = 1)

h2s = seq(-1, 5, length.out = 50)
h1 = 1
plot(h2s, log(sapply(h2s, function(h2) varphi_ratio(c(h1, h2)))), xlab = "h2", ylab = "r/r*", type = "l", lty = 1, lwd = 1, col = "blue")
lines(h2s, log(sapply(h2s, function(h2) varphi_star_ratio(c(h1, h2)))), lwd = 1)
