
## This code generates LPS Gamma Shared Frailty Model for the CGD study (1991). ##

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(blapsr)

d = read.csv("cgd.csv")
t = d$tstop - d$tstart
delta = d$status
ni = table(d$id)

## Zx = matrix(0, nrow=sum(ni), ncol=1)     ## for a Null model (baseline hazard)
Zx = matrix(c(d$treat, d$sex, d$age, d$height, d$weight), nrow=sum(ni))
p = ncol(Zx)

J = 300
Delta = max(t)/J
bins = findInterval(t, seq(0, max(t), length.out=J))
K = 30
r = 2
Xi = rep(-3, (K+p)); Xi[(K+1):(K+p)] = 0

lambda = 100
gamma = 5

source("Functions.R")
Q = pen_mx(r, K)

Xi = NRA(lambda, gamma)$Xi_star      ## initial run, to get the starting Xi.

optim_log_hyper = optim(c(eta_tilde), log_hyper, control=list(fnscale=-1))
lambda = exp(optim_log_hyper$par[1])
gamma = exp(optim_log_hyper$par[2])

NRA_star = NRA(lambda, gamma)        ## final run, with the optimized \lambda and \gamma

eff = c("treatment", "female sex", "age [yrs]", "height [cm]", "weight [kg]")
for (xx in (K+1):(K+p)) {
  cat(eff[xx-K], "\t", round(NRA_star$Xi_star[xx], 5), "\t",
      round(sqrt(diag(-solve(NRA_star$Hessian_star))[xx]), 5), "\n")
}
print(round(gamma, 5))
print(round(lambda, 1))



############### Survival Curves with the CI Bands ################
theta=NRA_star$Xi_star[1:K, 1, drop=F]
beta=NRA_star$Xi_star[K+1]
Sigma_star=-solve(NRA_star$Hessian_star[1:(K+1), 1:(K+1)])
grid = seq(0, max(t), length=J)

############### Baseline hazard ###############
hazard=exp(B_sm %*% theta)

G = matrix(0, nrow=J, ncol=K+1)
denom = cumsum(hazard)
for (tp in 1:J) {
  for (wrt_theta in 1:K) {
    numer = sum(exp(B_sm[1:tp, ] %*% theta) * B_sm[1:tp, wrt_theta])
    G[tp, wrt_theta] = numer/denom[tp]
  }}

plot(1,1, xlim=c(0,400), ylim=c(0,1), type="n", xlab="Time [days]", ylab="Survival")
surv_curves(hazard, G, Sigma_star, Delta, "red")

############### Hazard with z = 1 ("treated") ###############
hazard=exp(B_sm %*% theta + beta)
G[ ,(K+1)] = 1
surv_curves(hazard, G, Sigma_star, Delta, "green")

km_fit = survfit(Surv(t, delta) ~ d$treat)
lines(km_fit, conf.int=T, col=c("red", "green"))

legend(0, 0.22, legend=c("IG treated", "untreated"), lwd=2, col=c("green", "red"), box.lty=0)


