
colors = sample(palette())[1:6]

log_l <- function(mu,w,x, tau){
  return(sum(log(w%*%(sapply(x, FUN = dnorm, mean = mu, sd = sqrt(1/tau))))))
}

l <-function(mu, w, x, tau){
  return(exp(log_l(mu,w, x, tau)))
}

eval_l <- function(mu1,mu2,w, x, tau, log=TRUE) {
  if (log == FALSE) (likelihood = l) else (likelihood = log_l)
  return (apply(cbind(mu1, mu2),1, likelihood, w = w, x=x, tau=tau))
}

model_checking <-function(Sample, y){
  
  if(dim(Sample)[2]>5){
    par(mfrow=c(1,1))
    hist(y, freq = FALSE, xlab = 'x', ylab = 'Density', main = 'P.P. Check')
    lines(density(Sample[,6]), col='red')
    
  }
  dic = DIC(Sample[,c(1:5)], y)
  return(dic)
}

DIC <- function(Sample, y) {
  S = dim(Sample)[1]
  lik = rep(NA, S)
  
  for (i in 1:S) {
    gr(w1, w2, mu1, mu2, tau) %=% Sample[i, ]
    lik[i] = eval_l(mu1, mu2, c(w1, w2), y, tau, log = TRUE)
  }
  gr(w1, w2, mu1, mu2, tau) %=% colMeans(Sample)
  
  D_bar = -2 * mean(lik)
  D_hat = -2 * eval_l(mu1, mu2, c(w1, w2), y, tau, log = TRUE)
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
  if(J>6){
    colnames(MC_)<- c('w1', 'w2', 'mu1', 'mu2', 'tau_1', 'tau_2','X_new')  
  } else{
    colnames(MC_)<- c('w1', 'w2', 'mu1', 'mu2', 'tau_1','tau_2')
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
    'Mean' = colMeans(MC_),
    'sd' = apply(MC_, 2, sd),
    'var' = apply(MC_, 2, var),
    'HPD' = HPDinterval(as.mcmc(MC_)),
    'effective_size' = effectiveSize(as.mcmc(MC_))
  )
  return(res)
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
  # if (dim(permuted_mcmc)[2]>J+C){
  #   permuted_mcmc = permuted_mcmc[,-(J+C)]
  # }
  return(permuted_mcmc)
  # return(permuted_mcmc[,c(2,1,4,3,5)])
}

gibbs_mix <- function(x, init, hyper, n_iter = 5000, burn_in = 2000, n_thin = 10, predict_x = FALSE){
  
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
    
    # mu_prop = Mu[k, ] + runif(2 , min=-100  ,max=100)
    # 
    # omega = runif(1, min = 0, max = 1)
    # cond = c( target(mu_prop, tau=Tau[k+1]) / target (Mu[k,], tau = Tau[k+1,]))
    
    
    # if (anyNA(cond)){
    #   return(Mu[k,])
    # }
    # ACCEPT = (omega < min( cond, 1 ))
    # Mu[k+1, ] = Mu[k,]
    # 
    # if(ACCEPT){
    #   Mu[k+1, ] = mu_prop
    # }
    

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
    Sample_mc = Sample_mc[-mask, ]
  }
  return(Sample_mc)   
  
}

norm_mix <- function(x, mu, tau, w) {
  sum(w * dnorm(x , mu, (1 / tau) ^ 1 / 2))
}
d_norm_mix = Vectorize(FUN = norm_mix, vectorize.args = 'x')

# try 2D

mix_gauss_2D <-function(mu, tau, z){
  
  mus = mu * diag(2)
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

