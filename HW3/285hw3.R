rm(list = ls())

##### Problem 2 ######
pvec <- c(0.1, 0.5, 0.4)
muvec <- c(-1, 0, 1)
sdvec <- c(0.2, 1, 0.4)

normmixture <- function(N, p, m, s){
  if(!length(p) == length(m) | !length(m) == length(s)){
    stop("Incorrect dimensions of input parameters.")
  }else{
    out <- NULL
    count <- rmultinom(1, N, prob = p)
    for(i in 1:length(count)){
      temp <- rnorm(n = count[i], mean = m[i], sd = s[i])
      out <- c(out, temp)
    }
  }
  return(out)
}

truemixnorm <- function(x, p, m, s){
  if(!length(p) == length(m) | !length(m) == length(s)){
    stop("Incorrect dimensions of input parameters.")
  }else{
    temp <- 0
    for(i in 1:length(p)){
      temp <- temp + p[i] * dnorm(x, mean = m[i], sd = s[i])
    } 
  }
  return(temp)
}

curveplotter <- function(x){
  return(truemixnorm(x, p = pvec, m = muvec, s = sdvec))
}

dat = normmixture(N = 50, p = pvec, m = muvec, s = sdvec)

fhat_dat = density(dat)
plotmax = max(fhat_dat$y, curveplotter(seq(-2, 2, length.out = 100)))
plot(fhat_dat, lty = 3, ylim = c(0, 1.2 * plotmax))
curve(curveplotter, n = 1000, add = T, col = "red")

Sethu_post = function(alpha, A, B, H){
  v = w = array(dim = H)
  v[H] = 1
  for(i in 1:(H - 1)){
    v[i] = rbeta(n = 1, A[i] + 1, B[i] + alpha)
  }
  w[1] = v[1]
  for(h in 2:H){
    w[h] = v[h] * prod(1 - v[1:(h - 1)])
  }
  return(w)
}

bnp_blockgibbs <- function(y, H = 20, chainparams, hyparams){
  # Set hyper-parametrs
  alpha <- hyparams$alpha
  mu0 <- hyparams$mu0
  kappa <- hyparams$kappa
  a_tau <- hyparams$a
  b_tau <- hyparams$b
  
  # Chain size
  N.samples <- round((chainparams$n.samples * 
                  chainparams$n.thin) / (1 - chainparams$burnin))
  
  # Pre-allocation
  n <- length(y)
  r <- array(dim = c(N.samples, n))
  w <- array(dim = c(N.samples + 1, H))
  m <- array(dim = c(N.samples + 1, H))
  tausq <- array(dim = c(N.samples + 1, H))
  library(common)
  cat("Running Blocked Gibbs sampler with options:\n",
      "n.iter = ", N.samples, 
      ", Burn-in = ", 0.1 * N.samples, "\n",
      "Prior: DPM(\U03B1, G"  %p% subsc("0"), ")\n",
      "Hyperpriors: G" %p% subsc("0"), 
      "= N(\U03BC"%p% subsc("0"), ",\U03BA)\n", 
      "\t\t\U03BC"%p% subsc("0"), " = ", mu0, ", \U03BA = ", kappa, ";\n",
      "\t\t \U03B1 = ", alpha, "; a = ", a_tau, ", b = ", b_tau, ".\n",
      "========================================\n")
  
  # Initialization
  w[1, ] = rep(1/H, H)
  m[1, ] = runif(H, min = min(y), max = max(y))
  tausq[1, ] = rep(1, H)
  
  # Chain
  t0 = Sys.time()
  for(t in 1:N.samples){
    for(i in 1:n){
      samp_prob_i = w[t, ] * dnorm(y[i], mean = m[t, ], sd = sqrt(tausq[t, ]))
      r[t, i] = sample(1:H, size = 1, prob = samp_prob_i)
    }
    A = B = array(dim = H)
    for(h in 1:H){
      A[h] = sum(r[t, ] == h)
      B[h] = sum(r[t, ] > h)
    }
    w[t+1, ] = Sethu_post(alpha = alpha, A = A, B = B, H = H)
    for(h in 1:H){
      if(sum(r[t, ] == h) == 0){
        m[t+1, h] = rnorm(n = 1, mean = mu0, sd = sqrt(kappa))
        tausq[t+1, h] = 1/rgamma(n = 1, shape = a_tau, rate = b_tau)
      }else{
        S_h = which(r[t, ] == h)
        atom_h_prec = (length(S_h) / tausq[t, h]) + (1 / kappa)
        atom_h_mean = ((length(S_h) * mean(y[S_h]) / tausq[t, h]) + 
                         (mu0 / kappa)) / atom_h_prec
        m[t+1, h] = rnorm(n = 1, mean = atom_h_mean, sd = sqrt(1/atom_h_prec))
        sse_h = sum((y[S_h] - m[t, h])^2)
        tausq[t+1, h] = 1/rgamma(n = 1, shape = a_tau + 0.5*length(S_h), 
                                 rate = b_tau + 0.5*sse_h)
      }
    }
  }
  t.elapsed = Sys.time() - t0
  cat("Elapsed", round(t.elapsed,2), units(t.elapsed), "\n")
  
  # Chain cleaning
  w = w[-1, ]
  m = m[-1, ]
  tausq = tausq[-1, ]
  indices = chainparams$burnin * N.samples + 
    (1:chainparams$n.samples) * chainparams$n.thin
  return(list(weights = w[indices, ], 
              atom_m = m[indices, ], 
              atom_v = tausq[indices, ]))
}

hyperparams <- list(alpha = 1, mu0 = 0, kappa = 1, a = 1, b = 1)
chain_params <- list(n.samples = 1000, burnin = 0.1, n.thin = 10)
out <- bnp_blockgibbs(y = dat, 
                      H = 20, 
                      chainparams = chain_params, 
                      hyparams = hyperparams)

# Plot: Blocked Gibbs
plot(fhat_dat, lty = 3, lwd = 1.5, ylim = c(0, 1.25 * plotmax),
     main = "", xlab = expression(paste(y)))
curve(curveplotter, lwd = 1, n = 1000, add = T, col = "red")
for(i in 950:1000){
  post_curveplotter <- function(x){
    return(truemixnorm(x, p = out$weights[i, ], 
                       m = out$atom_m[i, ], 
                       s = sqrt(out$atom_v[i, ])))
  }
  curve(post_curveplotter, n = 1000, add = T, 
        lwd = 0.5, lty = 1,
        col = adjustcolor("green", alpha.f = 0.3))
}
legend("topright", legend = c("Frequentist", "BNP"),
       bty = "n", lty = c(2, 1),
       col = c("black", "green"))


##### Slice sampler #####

v_to_w <- function(v){
  w = array(dim = length(v))
  w[1] = v[1]
  for(h in 2:length(v)){
    w[h] = v[h] * prod(1 - v[1:(h - 1)])
  }
  return(w)
}

y <- dat
n <- length(y)
N.samples <- 10
alpha <- 1
mu0 <- 1
kappa <- 1
a_tau <- 1
b_tau <- 1
nC <- 20
r <- array(dim = c(N.samples + 1, n))
u <- array(dim = n)
w <- vector(mode = "list", length = N.samples + 1)
m <- vector(mode = "list", length = N.samples + 1)
tausq <- vector(mode = "list", length = N.samples + 1)

r[1, ] <- sample(1:nC, size = n, replace = TRUE)
v <- rbeta(nC, 1, alpha)
w[[1]] <- v_to_w(v)
m[[1]] <- runif(nC, min = min(y), max = max(y))
tausq[[1]] <- rep(1, nC)

for(t in 1:N.samples){
  A = B = array(dim = nC)
  for(h in 1:nC){
    A[h] = sum(r[t, ] == h)
    B[h] = sum(r[t, ] > h)
  }
  for(h in 1:nC){
    if(sum(r[t, ] == h) == 0){
      m[[t+1]][h] = rnorm(n = 1, mean = mu0, sd = sqrt(kappa))
      tausq[[t+1]][h] = 1/rgamma(n = 1, shape = a_tau, rate = b_tau)
    }else{
      S_h = which(r[t, ] == h)
      atom_h_prec = (length(S_h) / tausq[[t]][h]) + (1 / kappa)
      atom_h_mean = ((length(S_h) * mean(y[S_h]) / tausq[[t]][h]) + 
                       (mu0 / kappa)) / atom_h_prec
      m[[t+1]][h] = rnorm(n = 1, mean = atom_h_mean, sd = sqrt(1/atom_h_prec))
      sse_h = sum((y[S_h] - m[[t]][h])^2)
      tausq[[t+1]][h] = 1/rgamma(n = 1, shape = a_tau + 0.5*length(S_h), 
                               rate = b_tau + 0.5*sse_h)
    }
  }
  
  r[t+1, ] <- sample(1:nC, size = n, replace = TRUE)
}







