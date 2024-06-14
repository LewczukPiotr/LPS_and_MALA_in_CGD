

############################### NRA iteration loop ###############################
NRA_one_iteration = function(Xi, Q_lambda, gamma, n_iter) {
  
  Theta=matrix(Xi[1:K], ncol=1); Beta=matrix(Xi[(K+1):(K+p)], ncol=1)
  
  #################### THE KERNEL ####################
  Gradient = kernel_g(Theta, Beta, gamma)
  Hessian = kernel_h(Theta, Beta, gamma)
  ll = logL(Xi, gamma)
  
  ##################### Updates #####################
  Gradient = Gradient - Q_lambda %*% Xi
  Hessian = Hessian - Q_lambda
  Xi = Xi - solve(Hessian) %*% Gradient
  
  return(list(Xi=Xi, Gradient=Gradient, Hessian=Hessian, ll=ll, n_iter=n_iter+1))
}
######################################################################################



###################################### Full NRA ######################################
NRA = function(lambda, gamma) {
  n_iter = 0
  Q_lambda = Q
  Q_lambda[1:K, 1:K] = lambda*Q_lambda[1:K, 1:K]
  not_converged=TRUE
  while (not_converged) {
    temp = NRA_one_iteration(Xi, Q_lambda, gamma, n_iter)
    Xi = temp$Xi
    n_iter = temp$n_iter
    not_converged = norm(temp$Gradient, "2") > 1e-6}
  
  return(list(Xi_star=temp$Xi, 
              Hessian_star=temp$Hessian,
              ll_star=temp$ll,
              n_iter_star=temp$n_iter))
}
################################################################################



############################## The penalty matrix ##############################
pen_mx = function(r, K) {
  Dr = diag(K)
  for(k in 1:r) {Dr = diff(Dr)}
  P = t(Dr) %*% Dr + diag(1e-6, K)
  Q = diag(1e-5, K+p)
  Q[1:K, 1:K] = P
  return(Q)
}
################################################################################



########################## Gradient for the MALA part ########################## 
gradient_for_MALA = function(Zeta_, Q_lambdaMALA_) {
  
  XiMALA_ = Zeta_[1:(K+p), drop=F]
  Theta = Zeta_[1:K, drop=F]
  Beta = Zeta_[(K+1):(K+p), drop=F]
  gammaMALA_ = exp(Zeta_[K+p+1])
  
  GradientMALA = kernel_g(Theta, Beta, gammaMALA_)
  
  S_aibi = sum(digamma(a_vv) - log(b_iv) - a_i_div_b_i)
  GradientMALA = GradientMALA - Q_lambdaMALA_ %*% XiMALA_
  deriv_gammaMALA = gammaMALA_*(I_clust*((log(gammaMALA_) + 1 - digamma(gammaMALA_)))
                                + S_aibi - b_gamma) + a_gamma
  GradientMALA = matrix(c(GradientMALA, deriv_gammaMALA), ncol=1)
  
  return(GradientMALA)
}
################################################################################



############################# The Kernel functions #############################
kernel_g = function(Theta, Beta, gamma) {
  
  assign("OMEGA_LIST", mapply(
    function(X,Y) {exp(B_sm %*% Theta %*% matrix(rep(1, nrow(X)), nrow=1)
                  + matrix(rep(1, J), ncol=1) %*% t((X) %*% Beta)) * Y},
                  Zx_LIST, DELTA_LIST, SIMPLIFY=F), envir=.GlobalEnv)
  assign("SS_b_LIST", lapply(OMEGA_LIST, ff1), envir=.GlobalEnv)
  assign("SS_zij_LIST", mapply(ff2, Zx_LIST, OMEGA_LIST, SIMPLIFY=F), envir=.GlobalEnv)
  assign("b_iv", ((unlist(lapply(OMEGA_LIST, sum))) + gamma), envir=.GlobalEnv)
  assign("a_vv", a_v + gamma, envir=.GlobalEnv)
  assign("a_i_div_b_i", a_vv/b_iv, envir=.GlobalEnv)
  assign("a_i_div_b_i2", a_vv/(b_iv)^2, envir=.GlobalEnv)
  
  G1 = matrix(rowSums(mapply(ff3, Bt_LIST, delta_LIST, a_i_div_b_i, SS_b_LIST)), ncol=1)
  G2 = matrix(rowSums(matrix( mapply(ff3, Zx_LIST, delta_LIST, a_i_div_b_i, SS_zij_LIST),
       nrow=p)), ncol=1)
  Gradient = matrix(c(G1, G2), ncol=1)
  
  return(Gradient)
}


kernel_h = function(Theta, Beta, gamma) {
  
  SS_b_b_LIST = lapply(OMEGA_LIST, function(X) {t(B_sm) %*% diag(rowSums(X), J) %*% B_sm})
  SS_zij_zij_LIST = mapply(function(X,Y) {t(X) %*% diag(colSums(Y), nrow(X)) %*% X},
           Zx_LIST, OMEGA_LIST, SIMPLIFY=F)
  SS_b_zij_LIST = mapply(function(X,Y) {t(B_sm) %*% X %*% Y},
           OMEGA_LIST, Zx_LIST, SIMPLIFY=F)
  
  H11 = Reduce("+", mapply(ff4, a_i_div_b_i2, SS_b_LIST, b_iv, SS_b_b_LIST, SIMPLIFY=F))
  H12 = Reduce("+", mapply(function(X,Y,Z,W,U) {X * (Y %*% t(Z) - W * U)},
               a_i_div_b_i2, SS_b_LIST, SS_zij_LIST, b_iv, SS_b_zij_LIST, 
               SIMPLIFY=F))
  H22 = Reduce("+", mapply(ff4, a_i_div_b_i2, SS_zij_LIST, b_iv, SS_zij_zij_LIST, SIMPLIFY=F))
  
  Hessian = rbind(cbind(H11, H12), cbind(t(H12), H22))
  
  return(Hessian)
}


logL = function(Xi_, g_) {
  Theta_ = Xi_[1:K, drop=F]
  Beta_ = Xi_[(K+1):(K+p), drop=F]
  
  temp_ = sum(mapply(function(X,Y,Z,W,U) 
  {(t(X) %*% (Y %*% Theta_ + Z %*% Beta_) - 
      W*log(U) + g_*log(g_) - lgamma(g_) + lgamma(W))},
  delta_LIST, Bt_LIST, Zx_LIST, a_vv, b_iv))
  
  return(temp_)
}
########################################################################################



#################################### Hyperparameter ####################################
nu = 2
a_delta = 1e-4
b_delta = 1e-4
a_gamma = 1e-4
b_gamma = 1e-4

v = log(lambda)
w = log(gamma)
eta_tilde = c(v,w)

log_hyper = function(eta_tilde) {
  v=eta_tilde[1]
  w=eta_tilde[2]
  
  NRA_con=NRA(exp(v), exp(w))
  H = NRA_con$Hessian_star
  ll = NRA_con$ll_star
  Xi = NRA_con$Xi_star
  
  output = 1/2*log(det(-solve(H))) + 
    (K+nu)*v/2 + a_gamma*w - 
    (nu/2+a_delta)*log(nu/2*exp(v)+b_delta) +
    ll - 1/2*exp(v)*t(Xi) %*% Q %*% Xi - b_gamma*exp(w)
  
  return(output)
}
#################################################################################



####################### Log-posterior for the MALA part #######################
log_posterior = function(Zeta_, Q_) {
  xi__ = Zeta_[1:(K+p)]
  g__ = exp(Zeta_[(K+p+1)])
  out = logL(xi__, g__) - 1/2 * t(xi__) %*% Q_ %*% xi__ - g__*b_gamma + a_gamma*log(g__)
  return (out)}
###############################################################################



################################ Matrices and Lists ############################
mx_slice = function(from, to, X) {
  temp=list(X[ ,from[1]:to[1], drop=F])
  for (i in 2:length(from)) {
    temp[[i]] = X[ ,from[i]:to[i], drop=F]}
  return(temp)}

from_v = cumsum(ni)-ni+1
to_v = cumsum(ni)

B_sm = matrix(0, nrow=J, ncol=K)
for(m in 1:J) {B_sm[m, ] =  cubicbs((m - 1/2)*Delta, lower=0, upper=max(t), K=K)$Bmatrix}

Bt = matrix(0, nrow=sum(ni), ncol=K)
for (j in 1: sum(ni)) {Bt[j, ] = cubicbs(t[j], lower=0, upper=max(t), K=K)$Bmatrix}

a_v=NULL; for (i in 1:length(ni)) a_v[i] = sum(delta[from_v[i]:to_v[i]])
Bt_LIST = lapply(mx_slice(from_v, to_v, t(Bt)), "t")
Zx_LIST = lapply(mx_slice(from_v, to_v, t(Zx)), "t")

DELTA = matrix(0, nrow=J, ncol=sum(ni))
for (n in (1:sum(ni))){DELTA[1:bins[n] ,n] = Delta}
DELTA_LIST = mx_slice(from_v, to_v, DELTA)

delta_LIST = lapply(mx_slice(from_v, to_v, matrix(delta, nrow=1)), "t")
########################################################################################



##################################### Helpers ########################################## 
ff1 = function(X) {matrix(rowSums(t(B_sm) %*% X), nrow=K)}

ff2 = function(X,Y) {matrix(rowSums(t(X) %*% t(Y)), nrow=p)}

ff3 = function(X,Y,Z,W) {matrix((t(X) %*% Y - Z * W), ncol=1)}

ff4 = function(X,Y,Z,W) {X * (Y %*% t(Y) - Z * W)}
########################################################################################

