# Parameters
n = 40
theta = 2
t = 0.5
sigmaW = 1
y_bar = 1/theta 
u = 3.8

# p_star functions
p_star = function(z, type){
  if (type == "det"){
    return(as.numeric(z > theta^2*sqrt(n)*(1/theta - t)))
  } else if (type == "carving"){
    K = 2*(sqrt(n)*theta^2)^(n-1)*(gamma(n)/gamma(n/2)^2)*2^(-n)
    a = (2/theta - 2*z/(sqrt(n)*theta^2))
    f = Vectorize(function(x) (x*(a-x))^(n/2 - 1))
    return(K*(integrate(f, 0, t)$value)/(sqrt(n)*theta - z)^(n-1))
  } else if (type == "rand"){
    return(exp(pnorm((sqrt(n)/sigmaW)*(t - 1/theta) + z/(sigmaW*theta^2), log.p = T)))
  } else if (type == "cond"){
    return(exp(dnorm(u - sqrt(n)/theta + z/theta^2, sd = sigmaW, log = T)))
  }
}

# Selection probability ratios
varphi = function(theta, type){
  if (type == "det"){
    return(pgamma(t, shape = n, rate = theta*n))
  } else if (type == "carving"){
    return(pgamma(t, shape = n/2, rate = theta*n/2))
  } else if (type == "rand"){
    integrand = Vectorize(function(w) pgamma(t - w/sqrt(n), shape = n, rate = n*theta)*dnorm(w, sd = sigmaW))
    return(integrate(integrand, -4*sigmaW, 4*sigmaW)$value)
  } else if (type == "cond"){
    integrand = Vectorize(function(w) exp(dgamma(u - w, shape = n, rate = sqrt(n)*theta, log = T) + dnorm(w, sd = sigmaW, log = T)))
    return(integrate(integrand, -4*sigmaW, 4*sigmaW)$value)
  }
}
varphi_ratio = function(h, type) varphi(theta + h/sqrt(n), type)/varphi(theta, type)
varphi_star = function(h, type){
  integrand = Vectorize(function(z) dnorm(z, mean = h, sd = theta)*p_star(z, type))
  return(integrate(integrand, -3*theta, 3*theta)$value)
}
varphi_star_ratio = function(h, type) varphi_star(h, type)/varphi_star(0, type)

# Plots
par(mfrow = c(2, 2))
hs = seq(-2, 2, length.out = 1e2)
plot(hs, log(sapply(hs, function(h) varphi_ratio(h, "det"))), type = "l", xlab = "h", ylab = "r/r*", col = "blue")
lines(hs, log(sapply(hs, function(h) varphi_star_ratio(h, "det"))))
plot(hs, log(sapply(hs, function(h) varphi_ratio(h, "carving"))), type = "l", xlab = "h", ylab = "r/r*", col = "blue")
lines(hs, log(sapply(hs, function(h) varphi_star_ratio(h, "carving"))))
plot(hs, log(sapply(hs, function(h) varphi_ratio(h, "rand"))), type = "l", xlab = "h", ylab = "r/r*", col = "blue")
lines(hs, log(sapply(hs, function(h) varphi_star_ratio(h, "rand"))))
plot(hs, log(sapply(hs, function(h) varphi_ratio(h, "cond"))), type = "l", xlab = "h", ylab = "r/r*", col = "blue")
lines(hs, log(sapply(hs, function(h) varphi_star_ratio(h, "cond"))))
