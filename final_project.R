library(R2jags)
library(MASS)
library(PLMIX)
library(gtools)
library(mclust)
library(label.switching)
setwd("~/Desktop/Uni/computational/final_project")
source('helper.R')

# relax assumption variance on classes

set.seed(1234)

# Question:
#   - come conviene inizializzare, dando Z o dando i pesi? 
#   - Quali sono le prior su Z o su w?
#   - Mi sa dare dei parametri per la dist che metta in evidenza il label switching nella MC?
#   - Perche nell esempio di bug l output sono 50 medie?
#   - Ma il label switching emerge dalla parametrizzazione con l'indicatore Zi o e' proprio del mixture?
#   - Come funziona il metodo del pivot? Perche devo inizializzare la chain alla stima di MV e con le gisute proporzioni?

# === Gaussian Mixture ====


N = 100
w = c(0.2, 0.8)
C = 2

comp = sample(1:C,prob = w,size = N, replace = TRUE)

mu = rnorm(C, 0, 20)
tau = rgamma(1, 2, 6)

x = rnorm(n = N,
          mean = mu[comp],
          sd = sqrt(1 / tau))

z = binary_group_ind(comp, G = C)


curve(d_norm_mix(x, mu = mu, tau = tau, w = w), from = -10 + min(x), to = max(x) + 10, xlab = 'x', ylab='Density')

all_S = mix_gauss_2D(mu, tau, z)

plot(
  x = all_S[, 1],
  y = all_S[, 2],
  col = all_S[, 3],
  pch = 16,
  cex = .7,
  xlab = 'x2',
  ylab = 'x1',
)

mu_1 = seq(-50,50,length.out=length(x))  
mu_2 = seq(-50,50,length.out=length(x))  

posterior_eval = outer(mu_1, mu_2, eval_l, w=w, x=x, tau = tau)
# please double check on how to interpret this 
image(mu_1, mu_2, posterior_eval, xlab = expression(mu1), ylab = expression(mu2))
contour(mu_1, mu_2, posterior_eval, add = TRUE)
points(mu*diag(2))

perm = gtools::permutations(C,C)

# proof of label switching
log_l(mu[perm[1,]], w= w[perm[1,]], x= x, tau = tau)
log_l(mu[perm[2,]], w= w[perm[2,]], x= x, tau = tau)


# === Gibbs sampling ====

init = list(mu = mu, z = z)

hyper = list(
  alpha_0 = rep(1, C),
  mu_0 = 0,
  tau_0 = 10 ^ (-6),
  shape_0 = 10 ^ (-3),
  rate_0 = 10 ^ (-3)
)

mc_sample = gibbs_mix(x, init, hyper, predict_x = FALSE, n_iter = 10000)
values = show_mcmc(mc_sample)

# === Gibbs sampling with MLE  ====

# something doesn't wokr hear i have no clue
mc_clust = Mclust(x,G=C)

Z = t(rmultinom(length(x),1, mc_clust$parameters$pro[c(2,1)]))

mc_sample_mle = gibbs_mix(
  x,
  init = list(mu = mc_clust$parameters$mean, Z = Z),
  hyper,
  predict_x = FALSE,
  n_iter = 10000
)

MC_mle = show_mcmc(mc_sample_mle)
permuted_mcmc = fix_label_switching(MC_mle$Sample, 2, 3, mc_clust)
permuted_values = show_mcmc(permuted_mcmc)
permuted_values$Mean

# === EYES example vanilla  ====

# Microspectrophotometric measurements of retinal receptors are reported for
# eight species of Old World monkey
# http://vision.psychol.cam.ac.uk/jdmollon/papers/SaimiriPolymorphism1985.pdf
#
# Abstract-Microspectrophotometric measurements have been obtained for individual photoreceptors
# from four female squirrel monkeys (Suinriri sciureus) that had been shown behaviourally to be
# trichromatic. Relative to a normal human observer, two of the monkeys required more red light for a
# Rayleigh match; the other two required more green light than a normal human observer. In the red-green
# spectral region, the first type of monkey was found to have two cone pigments with peak sensitivities
# at approximately 536 and MQnm, whereas the second type was found to have pigments with peak
# sensitivities at approximately 549 and 564 nm. By maximum likelihood estimation it was shown that the
# microspectrophotometric data could be described by a model that assumed only three underlying
# distributions, two of which were present in each type of monkey. The tit of this model was as good as
# one in which a “double normal” distribution was fitted individually to the data for each animal.

monkeys_data = list(
  "N" = 48,
  "y" = c(529.0,530.0,532.0,533.1,533.4,533.6,533.7,534.1,534.8,535.3,
           535.4,535.9,536.1,536.3,536.4,536.6,537.0,537.4,537.5,538.3,
           538.5,538.6,539.4,539.6,540.4,540.8,542.0,542.8,543.0,543.5,
           543.8,543.9,545.3,546.2,548.8,548.7,548.9,549.0,549.4,549.9,
           550.6,551.2,551.4,551.5,551.6,552.8,552.9,553.2),
  "Itot" = c(1,1)
)

monkeys_inits = list(
  "lambda" = c(535,NA),
  "theta" =  5,
  "tau" =  0.1,
  "T" = c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,
           1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)
)

z_m = binary_group_ind(monkeys_inits$T, 2)
mu_m = c(500, 600)

# a week prior will turn to move all sampels in one cluster
hyper = list(
  alpha_0 = c(4,4),
  mu_0 = 0,
  tau_0 = 10 ^ (-6),
  shape_0 = 10 ^ (-3),
  rate_0 = 10 ^ (-3)
)

mc_monkey = gibbs_mix(
  monkeys_data$y,
  init = list(mu = mu_m, z = z_m),
  hyper = hyper,
  n_iter = 10000,
  n_thin = 20
)

values_m = show_mcmc(mc_monkey)

mc_clust_m = Mclust(data = monkeys_data$y)
# print(mc_clust_m)
# notice mc_clust_m scazza violentmente la proorzione

mc_monkey_mle = gibbs_mix(
  monkeys_data$y,
  init = list(mu = mc_clust_m$parameters$mean, z = t(
    rmultinom(monkeys_data$N, 1, mc_clust_m$parameters$pro)
  )),
  hyper = hyper,
  n_iter = 10000,
  n_thin = 20,
  predict_x = TRUE
)

values_m_mle = show_mcmc(mc_monkey_mle)
permuted_m = fix_label_switching(values_m_mle$Sample[,c(1:6)], 2, 3, mc_clust_m)
values_m_permuted = show_mcmc(permuted_m)

mle_dic = model_checking(mc_monkey_mle,monkeys_data$y)

# === EYES example JAGS  ====
params = c('P', 'lambda', 'tau')


# the bad thing about label swithcihng is that we understimate the variance
dir_file = 'examples/classic-bugs/vol2/eyes/eyes2.bug'

jags_model = jags.model(dir_file, monkeys_data, monkeys_inits, n.chains = 1)
mc_sample_jags = jags(
  data = monkeys_data,
  inits = list(init_1 = monkeys_inits),
  parameters.to.save = params,
  model.file = dir_file,
  n.chains = 1,
  n.iter = 10000,
  n.thin = 20,
  n.burnin = 2000
)
jags_m = show_mcmc(mc_sample_jags$BUGSoutput$sims.matrix[,-3])


# 
# mh1.rw.unif = function(nsim,initialstate,a,greekpi){
#   
#   # theta=rep(NA,nsim)
#   
#   theta = matrix(0, nrow = nsim, ncol = 2)
#   theta[1,] = initialstate
#   # theta[1]=initialstate
#   
#   for(t in 2:nsim){
#     
#     z=theta[t-1, ] 
#     
#     thetaprop = theta[t-1, ] + runif(2 ,min=-a,max=a)
#     
#     omega=runif(1,min=0,max=1)
#     
#     ACCEPT=(omega<min(c(greekpi(thetaprop)/greekpi(z),1)))
#     
#     theta[t, ]=z
#     
#     if(ACCEPT){
#       theta[t, ]=thetaprop
#     }
#     
#   }
#   return(theta)
# }
# 
# nsim = 10000
# a = 200
# initial_state = c(500, 600)
# initial_state_2 = c(550, 650)
# 
# target <- function(x, tau= 0.001){
#   
#   theta_0 = 550
#   # theta = 500
#   theta_seq = seq(-10000, 10000)
#   
#   num = (x - theta_0)^2 * dnorm(x, theta_0, sd = sqrt(1/tau))
#   denom = sum((theta_seq-theta_0)^2 * dnorm(theta_seq, theta_0, sd = sqrt(1/tau)))
#   return(num/denom)
#   # return(dmom(x, tau, baseDensity = 'normal'))
# }
# res = mh1.rw.unif(nsim, initial_state, a, target)
# res2 = mh1.rw.unif(nsim, initial_state, a, target)
# 
# res = as.mcmc(res)
# res2 = as.mcmc(res2)
# curve(target(x), to=700)
# acf(res)
# acf(res2)
# hist(res2)
# hist(res)
# curve(dmom(x, tau = 0.5), from = -3, to =3)



