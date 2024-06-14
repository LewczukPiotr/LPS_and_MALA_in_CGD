
## This code generates MALA within Gibbs Sampler ###########
## Gamma Shared Frailty Model for the CGD study (1991). ####
## To get the initial values, CGD_LPS needs to be run first.

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(mvtnorm)
library(coda)
source("CGD_LPS.R")  ## necessary only if CGD_LPS has not been run
set.seed(1234)

############# Initial values obtained from LPS #############
lambdaMALA = lambda
gammaMALA = gamma
ZetaMALA = matrix(c(NRA_star$Xi_star, log(gammaMALA)), ncol=1)
XiMALA = matrix(ZetaMALA[1:(K+p)], ncol=1)
deltaMALA = 1
rho = 1

P_tilde = matrix(0, ncol=K+p, nrow=K+p)
P_tilde[1:K, 1:K] = Q[1:K, 1:K]

Sigma_L = diag(ncol=K+p+1, nrow=K+p+1)
Sigma_L[1:(K+p), 1:(K+p)] = -solve(NRA_star$Hessian_star)

B = 50000
accepted = 0
I_clust = length(ni)
resMALA = matrix(0, ncol=nrow(ZetaMALA) + 4, nrow=B)

######################### MALA within Gibbs Loop #########################
for (xx in 1:B) {
  
  deltaMALA  = rgamma(1, shape=nu/2 + a_delta, rate=lambdaMALA*nu/2 + b_delta)
  
  lambdaMALA = rgamma(1, shape=(nu+K)/2, 
               rate=(t(XiMALA) %*% P_tilde %*% XiMALA + nu*deltaMALA)/2)
  
  Q_lambdaMALA = Q
  Q_lambdaMALA[1:K, 1:K] = lambdaMALA * Q_lambdaMALA[1:K, 1:K]
  
  GradientMALA = gradient_for_MALA(ZetaMALA, Q_lambdaMALA)
  l_post_current = log_posterior(ZetaMALA, Q_lambdaMALA)
  
  ZetaMALA_propos = matrix(rmvnorm(1, mean = ZetaMALA + 1/2*rho*Sigma_L %*% GradientMALA, 
                            sigma = rho*Sigma_L), ncol=1)
  
  GradientMALA_propos = gradient_for_MALA(ZetaMALA_propos, Q_lambdaMALA)
  l_post_propos  = log_posterior(ZetaMALA_propos, Q_lambdaMALA)
  
  log_probab = min((l_post_propos - l_post_current
                -1/2 * as.numeric(
                t(GradientMALA_propos + GradientMALA) %*%
                (ZetaMALA_propos - ZetaMALA + rho/4 * Sigma_L %*%
                (GradientMALA_propos - GradientMALA))) ), 0) 
  
  if (log_probab > log(runif(1))) {
    ZetaMALA = ZetaMALA_propos
    accepted = accepted + 1}
  
  zz = sqrt(rho) + 1/xx * (exp(log_probab) - 0.57)
  rho = as.numeric((1e-4*(zz < 1e-4) + zz*(zz > 1e-4 & zz < 1e4) + 1e4*(zz > 1e4))^2)
  
  XiMALA = matrix(ZetaMALA[1:(K+p)], ncol=1)
  gammaMALA = exp(ZetaMALA[K+p+1])
  resMALA[xx, ] = c(t(XiMALA), gammaMALA, lambdaMALA, deltaMALA, rho, accepted/xx)
  if (xx%%1000 == 0) {cat("Iteration: ", xx, "\n")}
}

accepted/B
plot(resMALA[ , ncol(resMALA)], type="l", ylim=c(0,1))
abline(h=0.57)

out = mcmc(resMALA[10001:B, ])
param = (K+1):(K+p+3)
summary(out)$statistics[param, ]
round(HPDinterval(out[ ,param]), 4)
traceplot(out[ ,param])
densplot(out[ ,param])

geweke.plot(out[ ,param])

