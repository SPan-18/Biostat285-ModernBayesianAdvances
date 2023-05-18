rm(list = ls())

library(lattice)
library(latticeExtra)

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

N.S = 10
sample.size = 100
dph.samples = NULL
for(alpha in c(0.1, 0.5, 1, 10)){
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
                     c("0.1", "0.5", "1", "10"), sep = " ")
dph.samples$alpha = rep(alpha.labels, each = sample.size)

fmla1 = as.formula(paste("~ ", paste(xnam, collapse= "+"), " | alpha"))
ecdfplot(fmla1, data = dph.samples,
         xlab = "", ylab = "", xaxt = "n", asp = 1.5)
