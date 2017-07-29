
colors = sample(palette())[1:7]

get_table <- function(values){
  df = data.frame(rbind(t(values$mean),
                        t(values$var),
                        t(values$b_var),
                        t(values$gamma_con),
                        t(values$t_eff)))
  row.names(df) = c('Mean', 'Var','Tau_B','Gamma_G','T_eff')
  table = xtable(df, digits = 4)
  return(table)
}

log_l <-function(mu, tau, x, w){
  C = length(mu)
  l = 0
  for (c in 1:C){
    l = l + w[c] * dnorm(x, mu[c], sqrt(1/tau[c]))
  }
  
  return(sum(log(l)))
}

l <-function(mu, w, x, tau){
  return(exp(log_l(mu,tau, x, w)))
}

model_checking <-function(MC_, y, get_dic = TRUE){
  ncol_mc = ncol(MC_)
  
  if (which(colnames(MC_) == 'x_new')){
    par(mfrow=c(1,1))
    
    hist(y, freq = FALSE, xlab = 'x', ylab = 'Density', main = NA, ylim = c(0, 0.15))
    
    col_new_x  = which(colnames(MC_) == 'x_new')
    lines(density(MC_[,col_new_x]), col='red')
  }
  if(get_dic){
    dic = DIC(MC_[,c(1:ncol_mc-1)], y)  
    return(dic)
  }
}


new_x_sim <-function(MC_){
  
  new_x = rep(NA, nrow(MC_))
  for (i in 1:nrow(MC_)){
    comp = rmultinom(1,1,prob = c(MC_[1,1], MC_[1,2]))
    new_x[i] = rnorm(1, mean = c(MC_[1,3], MC_[1,4]), sd = sqrt(1/ MC_[1,5]))
  }
  return(cbind(MC_, new_x))
  
}

DIC <- function(Sample, y) {
  S = dim(Sample)[1]
  lik = rep(NA, S)
  
  for (i in 1:S) {
    w = Sample[i,1:2]; mu = Sample[i, 3:4]; tau = Sample[i, 5:ncol(Sample)]
    lik[i] = log_l(mu, tau, y, w)
  }
  theta_hat = colMeans(Sample)
  w = theta_hat[1:2]; mu = theta_hat[3:4]; tau = theta_hat[5:ncol(Sample)]
  
  D_bar = -2 * mean(lik)
  D_hat = -2 * log_l(mu, tau, y, w)
  pD = D_bar - D_hat
  pV = var(-2 * lik) / 2
  return (list(
    DIC = pD + D_bar,
    IC = 2 * pD + D_bar,
    pD = pD,
    pV = pV,
    Dbar = D_bar,
    Dhat = D_hat
  ))
}

show_mcmc <-function(MC_){
  
  J = dim(MC_)[2]
  
  if (any(colnames(MC_) == 'lambda[1]')){
    colnames(MC_)<- c('w1', 'w2', 'mu1', 'mu2', 'tau', 'x_new')
  }else if (J == 6 ){
    colnames(MC_)<- c('w1', 'w2', 'mu1', 'mu2', 'tau1', 'tau2')  
  }else{
    colnames(MC_)<- c('w1', 'w2', 'mu1', 'mu2', 'tau1', 'tau2', 'x_new')  
  }
  
  par(mfrow=c(2,3))
  for (i in 1:J){
    plot(cumsum(MC_[,i])/(1:nrow(MC_)), type = 'l', col = colors[i], xlab= 'Iteration', ylab='Running Mean') 
  }
  par(mfrow=c(2,3))
  for (i in 1:J){
    plot(MC_[,i], type='l', col=colors[i],xlab = 'Iteration', ylab = colnames(MC_)[i])
  }
  par(mfrow=c(2,3))
  for (i in 1:J){
    acf(MC_[,i], main = colnames(MC_)[i], col=colors[i])
  }
  par(mfrow=c(2,3))
  for (i in 1:J){
    densplot(as.mcmc(MC_[,i]), xlab = colnames(MC_)[i], ylab = 'Density', col = colors[i])
  }
  
  res = list(
    'Sample' = as.mcmc(MC_),
    'mean' = colMeans(MC_),
    'var' = apply(MC_, 2, var),
    'b_var' = batch_var(MC_, 100),
    'gamma_con' = gamma_est(MC_),
    'HPD' = HPDinterval(as.mcmc(MC_)),
    't_eff' = effectiveSize(as.mcmc(MC_))
  )
  return(res)
}

gamma_est <-function(MC_){
  s = ncol(MC_)
  var_est = rep(NA, s)
  for (i in 1:s){
    var_est[i] =  mcmc::initseq(MC_[,i])$var.pos
  }
    return(var_est)
}

batch_var <- function(MC_, batch_size = 200){
  n_col = ncol(MC_)
  group_split = ceiling(seq_along(MC_[,1])/batch_size)
  B = max(group_split)
  b_var = rep(NA, n_col)
  
  for (i in 1:n_col){
    
    samples = split(MC_[,i], group_split)
    b_means = lapply(samples, mean)
    i_hat = mean(MC_[,i])
    b_var[i] = B * 1 / (batch_size - 1)  * sum((unlist(b_means) - i_hat)^2)
  }
  return (b_var)  
}


fix_label_switching <-function(MC_, C, J, mc_clust){
  
  L = nrow(MC_)
  
  mcmc_samples = array(MC_, dim = c(L, C, J))
  
  pra_out = pra(
    mcmc.pars = mcmc_samples,
    pivot = cbind(
      mc_clust$parameters$mean,
      1 / mc_clust$parameters$variance$sigmasq,
      mc_clust$parameters$pro
    )
  )
  
  permuted_mcmc = matrix(permute.mcmc(mcmc =  mcmc_samples, permutations = pra_out$permutations)$output, nrow=L)
  return(permuted_mcmc)
 
}

gibbs_mix <- function(x, init, hyper, n_iter = 5000, burn_in = 2000, n_thin = 10,predict_x = FALSE){
  
  gr(alpha_0, mu_0, tau_0, a_0, b_0 )%=% hyper
  n = length(x)
  C = length(init$mu)
  
  Mu  = W = Tau = matrix(NA, nrow = n_iter + 1, ncol = C)
  
  if (predict_x == TRUE){
    X_new = rep(NA, n_iter +1)
  }
  
  gr(Mu[1, ], Z) %=% init
  if (is.na(Z)){
    w = rdirichlet(n = 1, alpha_0)
    Z = t(rmultinom(n, 1, w ))
  }
  
  for (k in 1:n_iter) {
    N = colSums(Z)
    
    W[k + 1, ] = rdirichlet(n = 1, alpha = alpha_0 + N)
    
    b_n = colSums(Z * t(matrix(
      x,
      nrow = C,
      ncol = n,
      byrow = TRUE
    ) - Mu[k, ]) ^ 2) / 2
    
    # Same tau for each group
    
    Tau[k + 1,] = rgamma(n = 2,
                        shape = a_0 + N / 2,
                        rate = b_0 + b_n)
    
    prec = tau_0 + N * Tau[k + 1, ]
    delta = N * Tau[k + 1,] / prec
    
    Mu[k + 1, ] = rnorm(
      n = C,
      mean = delta * colSums(Z * x) / N + (1 - delta) * mu_0,
      sd = sqrt(1 / prec)
    )

    if (anyNA(Mu[k+1, ])){
      print('problem')
      return(list(W, Mu, Tau))
    }
    
    P = matrix(0, nrow = n, ncol = C)
    
    for (i in 1:n) {
      p = W[k + 1,] * dnorm(x[i], mean = Mu[k + 1,], sd = sqrt(1 / Tau[k+1])) 
      P[i,] = p / sum(p)
      Z[i,] = t(rmultinom(1, 1, p = P[i,]))  
    }
    
    if (predict_x == TRUE){
      comp = sample(1:C, prob = W[k+1, ], size = 1, replace = TRUE)
      X_new[k] = rnorm(1, mean = Mu[k+1, comp], sd = sqrt(1/Tau[k+1]))
    }
    
  }
  if (predict_x == TRUE){
    Sample_mc = cbind(W, Mu, Tau, X_new)[-c(1:burn_in),]
    
  }else{
    Sample_mc = cbind(W, Mu, Tau)[-c(1:burn_in),]  
  }
  if (n_thin >0){
    mask = seq(1, nrow(Sample_mc), n_thin)
    Sample_mc = Sample_mc[mask-1, ]
  }
  return(Sample_mc)   
  
}

norm_mix <- function(x, mu, tau, w) {
  sum(w * dnorm(x , mu, (1 / tau) ^ 1 / 2))
}
d_norm_mix = Vectorize(FUN = norm_mix, vectorize.args = 'x')

# try 2D

mix_gauss_2D <-function(mu, tau, z){
  
  mus = mu * matrix(c(1,1,1,1), ncol=2)
  Sigma = sqrt(1/tau) * diag(2)
  N = colSums(z)
  
  all_S = c()
  for (l in 1:C){
    S = mvrnorm(N[l], mu = mus[l, ], Sigma = Sigma)
    S = cbind(S, rep(l, N[l]))
    all_S = rbind(all_S, S)
  }
  
  return(all_S)
  
}


# ==== Group Function ====

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
gr = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

