library(R2jags)
library(MASS)
library(PLMIX)
library(mcmc)
library(mclust)
library(gtools)
library(label.switching)
library(xtable)

setwd("~/Desktop/Uni/computational/final_project/model")
source('helper.R')

set.seed(1234)

# === Gaussian Mixture ====

N = 50
w = c(0.6, 0.4)
C = 2

comp = sample(1:C,prob = w,size = N, replace = TRUE)

mu = c(530, 550)
tau = c(0.07,0.065)
# 
x = rnorm(n = N,
          mean = mu[comp],
          sd = sqrt(1 / tau[comp]))

z = binary_group_ind(comp, G = C)

par(mfrow=c(1,1))
curve(d_norm_mix(x, mu = mu, tau = tau, w = w), from = -10 + min(x), to = max(x) + 10, xlab = 'x', ylab='Density')
dev.copy(png, 'plot_1.png')
dev.off()

all_S = mix_gauss_2D(mu, tau, z)

plot(
  x = all_S[, 1],
  y = all_S[, 2],
  col = all_S[, 3],
  pch = 16,
  cex = .7,
  xlab = 'x2',
  ylab = 'x1'
)
dev.copy(png, 'plot_2.png')
dev.off()

# === Gibbs sampling ====

mu_init = c(500, 600)
z_init = t(rmultinom(N, 1, prob = c(0.5,0.5)))
init = list(mu = mu_init, z = z_init)

hyper = list(
  alpha_0 = c(2,2),
  mu_0 = 0,
  tau_0 = 10 ^ (-6),
  shape_0 = 10 ^ (-3),
  rate_0 = 10 ^ (-3)
)

mc_sample = gibbs_mix(x, init, hyper, predict_x = FALSE, n_iter = 10000, n_thin = 10)
values = show_mcmc(mc_sample)
dev.copy(png, 'plot_3.png')
dev.off()

get_table(values)

# === Gibbs sampling with MLE  ====

mc_clust = Mclust(x,G=C)

Z = t(rmultinom(length(x),1, mc_clust$parameters$pro[c(2,1)]))

mc_sample_mle = gibbs_mix(
  x,
  init = list(mu = mc_clust$parameters$mean, Z = Z),
  hyper,
  predict_x = TRUE,
  n_iter = 10000
)

MC_mle = show_mcmc(mc_sample_mle)
permuted_mcmc = fix_label_switching(MC_mle$Sample, 2, 3, mc_clust)
permuted_values = show_mcmc(permuted_mcmc)
dev.copy(png, 'plot_4.png')
dev.off()
tbl = get_table(permuted_values)
model_checking(MC_mle$Sample, x)
tbl

# === EYES example vanilla  ====

# http://vision.psychol.cam.ac.uk/jdmollon/papers/SaimiriPolymorphism1985.pdf


par(mfrow=c(1,1))
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

par(mfrow=c(1,1))
hist(monkeys_data$y, freq = F, xlab = 'x', ylab = 'Density', main =NA)
dev.copy(png, 'plot_5.png')
dev.off()
xtable(t(summary(monkeys_data$y)))

z_m = t(rmultinom(48, 1, c(0.5,0.5)))
mu_m = c(500, 600)

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
  n_iter = 20000,
  n_thin = 10
)

values_m = show_mcmc(mc_monkey)
dev.copy(png, 'plot_6.png')
dev.off()

tbl_m = get_table(values_m)

mc_clust_m = Mclust(data = monkeys_data$y)

mc_monkey_mle = gibbs_mix(
  monkeys_data$y,
  init = list(mu = mc_clust_m$parameters$mean, z = t(
    rmultinom(monkeys_data$N, 1, mc_clust_m$parameters$pro)
  )),
  hyper = hyper,
  n_iter = 20000,
  n_thin = 10,
  predict_x = TRUE
)
    
new_x = mc_monkey_mle[,7]
values_m_mle = show_mcmc(mc_monkey_mle)
permuted_m = fix_label_switching(values_m_mle$Sample[,c(1:6)], 2, 3, mc_clust_m)
values_m_permuted = show_mcmc(permuted_m)

tbl_m_permuted = get_table(values_m_permuted)
dev.copy(png, 'plot_7.png')
dev.off()


mle_dic = model_checking(values_m_mle$Sample,monkeys_data$y)

dev.copy(png, 'plot_8.png')
dev.off()

# === EYES example JAGS  ====
params = c('P', 'lambda', 'tau')
dir_file = 'eyes/eyes2.bug'

jags_model = jags.model(dir_file, monkeys_data, monkeys_inits, n.chains = 1)
mc_sample_jags = jags(
  data = monkeys_data,
  inits = list(init_1 = monkeys_inits),
  parameters.to.save = params,
  model.file = dir_file,
  n.chains = 1,
  n.iter = 20000,
  n.thin = 10,
  n.burnin = 2000
)

jags_m = new_x_sim(mc_sample_jags$BUGSoutput$sims.matrix[,-3])
jags_m = show_mcmc(jags_m)

dev.copy(png, 'plot_9.png')
dev.off()

tbl_jags = get_table(jags_m)

model_checking(jags_m$Sample, monkeys_data$y, get_dic = FALSE)

dev.copy(png, 'plot_10.png')
dev.off()

