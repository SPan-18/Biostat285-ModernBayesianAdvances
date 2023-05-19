rm(list = ls())

library(lattice)
library(latticeExtra)
library(ggplot2)

K.trunc = 20
Sethu_jump = function(alpha, K = K.trunc){
  v = w = array(dim = K)
  v[K] = 1
  v[1:(K - 1)] = rbeta(n = K - 1, 1, alpha)
  w[1] = v[1]
  for(h in 2:K){
    w[h] = v[h] * prod(1 - v[1:(h - 1)])
  }
  return(w)
}

generate_DPH = function(alpha, K, N.sample, 
                        base, baseparams, 
                        jump = Sethu_jump){
  if(base == 'norm' && length(baseparams) == 2){
    m = baseparams[1]
    s = baseparams[2]
    atoms = rnorm(n = K, mean = m, sd = s)
    w = jump(alpha, K)
    dph.sample = sample(x = atoms, replace = T, 
                        size = N.sample, prob = w)
  }else stop("Incorrect base measure/parameters.")
}

set.seed(1729)
N.S = 50
sample.size = 100
alphaseq = c(0.1, 0.5, 1, 10)
dph.samples = NULL
for(alpha in alphaseq){
  dph.samples.temp = array(dim = c(sample.size, N.S))
  for(i in 1:N.S){
    dph.samples.temp[, i] = generate_DPH(alpha = alpha, K = 20, 
                                         N.sample = sample.size,
                                         base = "norm", baseparams = c(0, 1))
  }
  xnam = paste("x", 1:N.S, sep = "")
  colnames(dph.samples.temp) = xnam
  dph.samples.temp = as.data.frame(dph.samples.temp)
  dph.samples = rbind(dph.samples, dph.samples.temp)
}
alpha.labels = paste("\U03B1 =",
                     as.character(alphaseq), sep = " ")
dph.samples$alpha = rep(alpha.labels, each = sample.size)

fmla1 = as.formula(paste("~ ", paste(xnam, collapse= "+"), " | alpha"))
ecdfplot(fmla1, data = dph.samples,
         xlab = "", ylab = "", xaxt = "n", asp = 1.5)

mv_functional = function(n.G, alpha, 
                         sample.size = 100,
                         K = 20, base = "norm",
                         baseparams = c(0, 1)){
  dph.samples.temp = array(dim = c(sample.size, n.G))
  for(i in 1:n.G){
    dph.samples.temp[, i] = generate_DPH(alpha = alpha, K = K, 
                                         N.sample = sample.size,
                                         base = base, 
                                         baseparams = baseparams)
  }
  m1.G = apply(dph.samples.temp, 2, mean)                   # first moment
  m2.G = apply(dph.samples.temp, 2, function(x) mean(x^2))  # second moment
  v.G = m2.G - (m1.G)^2                                     # variance functional
  return(cbind(m1.G, v.G))
}

n.G = 1000
res = mv_functional(n.G = n.G, alpha = 0.1)
hist(res[,1])

mv.df = data.frame(mean = as.numeric(),
                   var = as.numeric())
for(alpha in alphaseq){
  temp = mv_functional(n.G = n.G, alpha = alpha)
  mv.df = rbind(mv.df, temp)
}

names(mv.df) = c("mean", "var")
mv.df$alpha = rep(alpha.labels, each = n.G)

plot.mu = ggplot(mv.df) +
  geom_histogram(aes(x = mean, y = ..density..), 
                 binwidth = 0.1) +
  geom_density(aes(x = mean), lwd = 0.25, colour = 4,
               fill = 4, alpha = 0.25) +
  facet_wrap(~ alpha, ncol = 2) +
  xlab(expression(mu(G))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

plot.var = ggplot(mv.df) +
  geom_histogram(aes(x = var, y = ..density..), 
                 binwidth = 0.1) +
  facet_wrap(~ alpha, ncol = 2) +
  xlab(expression(mu(G))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

gridExtra::grid.arrange(plot.mu, plot.var, ncol = 2)



set.seed(1729)
N.S = 50
sample.size = 100
alphaseq = c(0.1, 0.5, 1, 10)
E.M.dph = array(dim = c(N.S, length(alphaseq)))
for(j in 1:length(alphaseq)){
  alpha = alphaseq[j]
  dph.samples.temp = array(dim = c(sample.size, N.S))
  for(i in 1:N.S){
    dph.samples.temp[, i] = generate_DPH(alpha = alpha, K = 20, 
                                         N.sample = sample.size,
                                         base = "norm", baseparams = c(0, 1))
  }
  m = apply(dph.samples.temp, 2, function(x) length(unique(x)))
  E.M.dph[, j] = m
}
E.M.theo = alphaseq * log((alphaseq + sample.size) / alphaseq)
alpha.grid = c(seq(0, 0.5, length.out = 20),
               seq(0.5, 1, length.out = 20),
               seq(1, 10, length.out = 20))
E.M.theo = alpha.grid * log((alpha.grid + sample.size) / alpha.grid)

plot(alpha.grid, E.M.theo, col = "red", xlim = c(0, 10.5),
     xlab = expression(paste(alpha)), ylab = expression(paste(E(M))), type = "l")
boxplot(E.M.dph, boxwex = c(0.25, 0.25, 0.25, 0.75), xaxt = "n", 
        whisklty = 1, whisklwd = 0.5, 
        staplelty = 1, staplelwd = 0.5,
        outpch = ".", outcol = "grey", outcex = 0.5, 
        medcol = "darkblue", col = "lightblue",
        at = alphaseq, add = T)



