rm(list = ls())
setwd("~/Documents/UCLA Study Stuff/Spring 23/Biostat 285/Biostat285HW/HW1")
# load library
library(readr)
library(MASS)
library(common)

dat <- readr::read_csv("dataHW1.csv", col_types = cols(...1 = col_skip()))
# View(dat)

Y <- dat$y
X <- as.matrix(dat[,-1])
# apply(X, 2, sd)

n = length(Y)
p = dim(X)[2]
tausq = 100
nu = 2
lambda = 1
a = 2
b = 2

# functions for gibbs

gamma_logmarginal = function(g){
  sub = which(g == 1)
  if(length(sub) > 1){
    A.g = solve((t(X[, sub]) %*% X[, sub]) + diag(1/tausq, length(sub)))
    XgtY = t(X[, sub]) %*% Y 
    S.g.sq = sum(Y^2) - t(XgtY) %*% A.g %*% XgtY
    res = 0.5 * log(det(A.g)) - (nu + 0.5 * n + 1) * log(lambda + 0.5 * S.g.sq)
  }else{
    res = - (nu + 0.5 * n + 1) * log(lambda + 0.5 * sum(Y^2))
  }
  return(res)
}

# g.trial = sample(0:1, 50, replace = T)
# gamma_logmarginal(g.trial)

alpha_j = function(g, j){
  if(g[j] == 0){
    res0 = gamma_logmarginal(g)
    g[j] = 1
    res1 = gamma_logmarginal(g)
  }else{
    res1 = gamma_logmarginal(g)
    g[j] = 0
    res0 = gamma_logmarginal(g)
  }
  return(exp(res0 - res1))
}

# alpha_j(g.trial, 4)

ci_post = function(mat, varname){
  ncol = dim(mat)[2]
  out = t(apply(mat, 2, function(x){
    c(mean(x, na.rm = T), sd(x, na.rm = T), 
      quantile(x, c(0.025, 0.5, 0.975), na.rm = T))}))
  rownames(out) = paste(varname, 1:ncol, sep = "")
  colnames(out) = c("Mean", "SD", "2.5%", "50%", "97.5%")
  return(out)
}

# find_best_gamma = function(mat){
#   
# }

n.iter = 10000
n.burnin = 1000
post.samples = vector(mode = "list", length = 4)
names(post.samples) = c("gamma", "sigmasq", "pi0", "beta")
post.samples$gamma = array(dim = c(n.iter + 1, p))
post.samples$beta = array(dim = c(n.iter + 1, p))
colnames(post.samples$gamma) = paste("gamma", 1:p, sep = "")
colnames(post.samples$beta) = paste("beta", 1:p, sep = "")
post.samples$sigmasq = array(dim = n.iter + 1)
post.samples$pi0 = array(dim = n.iter + 1)

post.samples$gamma[1,] = rep(1, p)
post.samples$pi0[1] = 0.5
A.g = solve((t(X) %*% X) + diag(1/tausq, p))
cat("Running SSVS Gibbs sampler with options:\n",
    "n.iter = ", n.iter, ", Burn-in = ", n.burnin, "\n",
    "Hyperpriors:\t\U1D6D1 ~ Beta (", a, ",", b, ")\n ", 
    "\t\t\U1D70E" %p% supsc("2"), " ~ InvGamma (", nu, ",", lambda, ")\n",
    "========================================\nProgress:\n")
t0 = Sys.time()
for(i in 1:n.iter){
  # sample gamma
  gamma.update = sample(1:p)
  post.samples$gamma[i + 1, ] = post.samples$gamma[i, ]
  for(j in gamma.update){
    alpha = alpha_j(g = post.samples$gamma[i + 1, ], j = j)
    pr = 1 / (1 + ((1 - post.samples$pi0[i]) / post.samples$pi0[i]) * alpha)
    post.samples$gamma[i + 1, j] = rbinom(1, 1, pr)
  }
  # sample sigmasq
  sub = which(post.samples$gamma[i+1, ] == 1)
  A.g = solve((t(X[, sub]) %*% X[, sub]) + diag(1/tausq, length(sub)))
  XgtY = t(X[, sub]) %*% Y 
  S.g.sq = sum(Y^2) - t(XgtY) %*% A.g %*% XgtY
  post.samples$sigmasq[i + 1] = 1 / rgamma(1, nu + 0.5 * n, lambda + 0.5 * S.g.sq)
  # sample pi0
  p1 = sum(post.samples$gamma[i+1, ])
  post.samples$pi0[i + 1] = rbeta(1, a + p1, b + p - p1)
  # sample beta
  g.zero = which(post.samples$gamma[i+1, ] == 0)
  post.samples$beta[i + 1, g.zero] = 0
  post.samples$beta[i + 1, sub] = mvrnorm(n = 1, mu = A.g %*% XgtY,
                                            Sigma = post.samples$sigmasq[i + 1] * A.g)
  t.elapsed = Sys.time() - t0
  if(i %% (0.1*n.iter) == 0) cat("Elapsed ", i, " iterations. (", 
                           round(t.elapsed,2), units(t.elapsed), ")\n")
}

post.gamma = post.samples$gamma[n.burnin:n.iter+1, ]
post.beta = post.samples$beta[n.burnin:n.iter+1, ]
colnames(post.beta) = paste("beta", 1:p, sep = "")
post.sigmasq = post.samples$sigmasq[n.burnin:n.iter+1]
post.pi0 = post.samples$pi0[n.burnin:n.iter+1]

plot(apply(post.samples$gamma, 1, sum), type = "l")
plot(post.samples$beta[,1], type = "l")
lines(post.samples$beta[,6], col = "red", lty = 2)

df.beta.summary = as.data.frame(ci_post(post.beta, "beta"))
acf(post.samples$sigmasq[-1])

require(reshape2)
df.beta = melt(post.beta)
boxplot(df.beta$value ~ df.beta$Var2,
        range = 1.96,
        xlab = "Posterior credible intervals", ylab = "",
        whisklty = 1, whisklwd = 0.5,
        staplelty = 1, staplelwd = 0.5,
        outpch = ".", outcol = "grey", outcex = 0.5,
        medcol = rep(c("darkgreen", "red"), times = c(5, 45)),
        col = rep(c("darkseagreen1", "coral"), times = c(5, 45)),
        yaxt = "n",
        horizontal = T, asp = 1)
axis(2, tick = T, at = 1:50, labels = paste("beta", 1:50, sep = ""),
     las = 2, cex.axis = 0.5)
abline(h = 5.5, lty = 2, lwd = 0.5, col = "orange")

library(plot.matrix)
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(post.samples$gamma[1:100,]), border = NA, las = 1,
     col = c("black", "lightyellow"),
     breaks = 2, cex.axis = 0.6, 
     xlab = "", ylab = "", key = NULL,
     main = expression(paste(gamma)))

acf(post.samples$beta[-1, 1])

ppi = data.frame(var = paste("X", 1:p, sep = ""),
                 prob = as.numeric(apply(post.gamma, 2, mean)))

ppi.ordered = ppi[order(ppi$prob), ]

dotchart(ppi.ordered$prob, labels = ppi.ordered$var, pch = 19,
         xlim = c(0, 1), pt.cex = 0.8,
         color = rep(c("lightblue4", "magenta4"),
                     times = c(45, 5)),
         cex = 0.5)

plot(post.samples$pi0, type = "l")

BFDR = function(lambda, ppi){
  decision = as.numeric(ppi > (1 / (1 + lambda)))
  return(sum((1 - ppi) * decision) / sum(decision))
}

BFDR(0.1826364, ppi = ppi$prob)

alpha = 0.05
lambdaseq = seq(0.001, 0.4, length.out = 100)
bfdr.vec = sapply(lambdaseq, function(x) BFDR(x, ppi = ppi$prob))

plot(lambdaseq, bfdr.vec, type = "l", las = 1,
     # ylim = c(0, 0.06),
     xlab = expression(paste(lambda)),
     ylab = "BFDR", cex.axis = 0.8, xaxt = "n")
abline(h = alpha, lty = 2, col = "brown1")
lambda.opt = lambdaseq[min(which(bfdr.vec > 0.05))]
arrows(x0 = lambda.opt, y0 =  bfdr.vec[min(which(bfdr.vec > 0.05))],
       x1 = lambda.opt, y1 = 0,
       length = 0.08, angle = 20,
       code = 2, lwd = 0.5,
       col = "forestgreen")
axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4), cex.axis = 0.8)
axis(1, at = lambda.opt, labels = NA, 
     col.ticks = "forestgreen", lwd = 4)

cutoff = 1/(1+lambda.opt)

library(horseshoe)

hs.object <- horseshoe(Y, X, nmc = 10000,
                       method.tau = "halfCauchy", 
                       method.sigma ="Jeffreys")
df <- data.frame(index = paste("beta", 1:p, sep = ""),
                 "Mean" = hs.object$BetaHat,
                 "SD" = apply(hs.object$BetaSamples, 1, sd),
                 "2.5%" = as.numeric(hs.object$LeftCI),
                 "50%" = as.numeric(hs.object$BetaMedian),
                 "97.5%" = as.numeric(hs.object$RightCI)
                 )
df$selected.CI <- HS.var.select(hs.object, Y, method = "threshold",
                                threshold = 0.3)

hs.beta = t(hs.object$BetaSamples)
hs.beta = melt(hs.beta)
boxplot(hs.beta$value ~ hs.beta$Var2,
        range = 1.96,
        xlab = "Horseshoe HPD intervals", ylab = "",
        whisklty = 1, whisklwd = 0.5,
        staplelty = 1, staplelwd = 0.5,
        outpch = ".", outcol = "grey", outcex = 0.5,
        medcol = rep(c("darkgreen", "red"), times = c(5, 45)),
        col = rep(c("darkseagreen1", "coral"), times = c(5, 45)),
        yaxt = "n",
        horizontal = T, asp = 1)
axis(2, tick = T, at = 1:50, labels = paste("beta", 1:50, sep = ""),
     las = 2, cex.axis = 0.5)
abline(h = 5.5, lty = 2, lwd = 0.5, col = "orange")

df$selected.CI <- HS.var.select(hs.object, Y, method = "threshold",
                                threshold = cutoff)


######################################################################

library(BoomSpikeSlab)

generate.X = function(n, p){
  return(matrix(rnorm(n * p), ncol = p))
}

n <- 100
p <- 10000
p.star <- 10
sigma <- 0.5
niter <- 1000

X <- generate.X(n, p)
beta.star <- rep(0, p)
beta.star[1:10] <- c(-5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5)
Y <- as.numeric(X %*% beta.star + rnorm(n = n, mean = 0, sd = sigma))

prior <- SpikeSlabPrior(X, Y,
                        optional.coefficient.estimate = rep(0, p))
model <- lm.spike(Y ~ 0 + X, niter = niter)

# plot.ts(model$beta)
hist(model$sigma)
mean(model$sigma)
median(model$sigma)
plot(model)
res = summary(model)
plot(model, "scaled.coefficients", whisklty = 1, whisklwd = 0.5,
     staplelty = 1, staplelwd = 0.5,
     outpch = ".", outcol = "grey", outcex = 0.5)
relev.true = paste("X", 1:10, sep = "")
relev.est = names(which(res$coefficients[, "inc.prob"] > 0.5))
false.disc = setdiff(relev.est, relev.true)
false.neg = setdiff(relev.true, relev.est)
res$coefficients[false.neg, "inc.prob"]
beta.star[as.numeric(sub('X', "", false.neg))]


relev.true = paste("X", 1:10, sep = "")
out2 <- data.frame(p = as.numeric(),
                   s = as.numeric(),
                   mean.s = as.numeric(),
                   median.s = as.numeric(),
                   FN = as.character(),
                   FP = as.character())
pseq = rep(c(100, 1000, 10000), each = 2)
sseq = rep(c(0.5, 3), times = 3)
ps.opt = cbind(pseq, sseq)
niter = 1000
burnin = 100
for(i in 1:6){
  p <- ps.opt[i, 1]
  sigma <- ps.opt[i, 2]
  out2[i, 1] <- p
  out2[i, 2] <- sigma 
  X <- generate.X(n, p)
  beta.star <- rep(0, p)
  beta.star[1:10] <- c(-5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5)
  Y <- as.numeric(X %*% beta.star + rnorm(n = n, mean = 0, sd = sigma))
  prior <- SpikeSlabPrior(X, Y)
  model <- lm.spike(Y ~ 0 + X, niter = niter)
  res <- summary(model)
  out2[i, 3] <- mean(model$sigma[-(1:burnin)])
  out2[i, 4] <- median(model$sigma[-(1:burnin)])
  relev.est = names(which(res$coefficients[, "inc.prob"] > 0.5))
  false.disc = setdiff(relev.est, relev.true)
  false.neg = setdiff(relev.true, relev.est)
  out2[i, 5] <- paste(false.disc, collapse = ", ")
  out2[i, 6] <- paste(false.neg, collapse = ", ")
}

print(out2)

library(mombf)
n <- 100
p <- 100
sigma <- 0.5
X <- generate.X(n, p)
beta.star <- rep(0, p)
beta.star[1:10] <- c(-5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5)
Y <- as.numeric(X %*% beta.star + rnorm(n = n, mean = 0, sd = sigma))

fit <- modelSelection(Y, X)
head(postProb(fit))

ps.options = cbind(c(rep(100, 4), 500),
                   c(100, 100, 1000, 1000, 1000),
                   c(0.5, 3, 0.5, 3, 3))
out.nlp = vector(mode = "list", length = 5)
for(i in 1:5){
  n <- ps.options[i, 1]
  p <- ps.options[i, 2]
  sigma <- ps.options[i, 3]
  X <- generate.X(n, p)
  beta.star <- rep(0, p)
  beta.star[1:10] <- c(-5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5)
  Y <- as.numeric(X %*% beta.star + rnorm(n = n, mean = 0, sd = sigma))
  fit <- invisible(modelSelection(Y, X, verbose = F))
  out.nlp[[i]] = postProb(fit)[1:6, ]
}

for(i in 1:5){
  cat("n = ", ps.options[i, 1], 
      ", p = ", ps.options[i, 2], 
      ", s = ", ps.options[i, 3], "\n")
  print(out.nlp[[i]])
  cat("\n")
}


