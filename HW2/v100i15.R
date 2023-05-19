#----------------------------------
#load the library and the data
library("BNPmix")
data(CPP, package = "BNPmix")

#----------------------------------
### Univariate (Sec. 5.1)
#----------------------------------

y <-  CPP[CPP$hosp == 11,]

# Prior elicitation assuming 3 clusters
DPprior <- PYcalibrate(Ek = 3, n = nrow(y), discount = 0)

prior <- list(strength = DPprior$strength, discount = 0)

mcmc <- list(niter = 5000, nburn = 4000)
output <- list(grid = seq(30, 46, length.out = 100), out_type = "FULL")

# Run the estimation procedure and plot the results
set.seed(42)
fit1 <- PYdensity(y = y$gest, mcmc = mcmc, prior = prior, output = output)
print(fit1)

plot(fit1, show_hist = TRUE,  xlab = "gestational age")
plot(fit1, band = FALSE, show_clust = TRUE,  xlab = "gestational age")

#-----------------------------------
# Multivariate example 1 (Sec. 5.2a)
#-----------------------------------

# Define y and tuning parameters
y <-  cbind(CPP$gest, log(CPP$dde), CPP$weight)
prior <- list(strength = 0.05, discount = 0.1)

grid <- expand.grid(seq(30, 46, length.out = 25),
                    seq(0.7, 5.2, length.out = 25), 
                    seq(0, 130, length.out = 25))
output = list(grid = grid, out_type = "FULL")
mcmc <- list(niter = 2000, nburn = 1000, m_imp = 100)

# Run the estimation procedure
set.seed(42)
fit2 <- PYdensity(y = y, mcmc = mcmc, prior = prior, output = output)

# Plot the results (intensive computation, be patient)
p12 <- plot(fit2, dim = c(1,2), show_clust =  TRUE,
            xlab="gestational age", ylab="log(DDE)")
p21 <- plot(fit2, dim = c(2,1), show_clust = TRUE,
            ylab="gestational age", xlab="log(DDE)")
p13 <- plot(fit2, dim = c(1,3), show_clust = TRUE,
            xlab="gestational age", ylab="weight")
p23 <- plot(fit2, dim = c(2,3), show_clust = TRUE,
            xlab="log(DDE)", ylab="weight")
p32 <- plot(fit2, dim = c(3,2), show_clust = TRUE,
            ylab="log(DDE)", xlab="weight")
p31 <- plot(fit2, dim = c(3,1), show_clust = TRUE,
            ylab="gestational age", xlab="weight")
p1 <- plot(fit2, dim = c(1,1), show_clust = TRUE,
           xlab="gestational age", ylab="density")
p2 <- plot(fit2, dim = c(2,2), show_clust = TRUE,
           xlab="log(DDE)", ylab="density")
p3 <- plot(fit2, dim = c(3,3), show_clust = TRUE,
           xlab="weight", ylab="density")

gridExtra::grid.arrange(
  p1, p12, p13,
  p21, p2, p23,
  p31, p32, p3,
  layout_matrix = matrix(1:9,3,3))

#----------------------------------
# Multivariate example 2 (Sec. 5.2b)
#----------------------------------

output <- list(out_type = "CLUST")
mcmc <- list(niter = 5000, nburn = 4000, model = "DLS", method = "MAR")

# Simulate the synthetic data set
p <- 25
set.seed(42)
ysim <- rbind(
  mnormt::rmnorm(50, mean = rep(1, p), varcov = diag(1, p)),
  mnormt::rmnorm(40, mean = sample(1:5, p, rep=TRUE), varcov = diag(1, p)),
  matrix(rt(10*p, df = 2), 10, p) * 2)
ysim <- scale(ysim)

# Fit the PY mixture
prior <- list(strength = 1, discount = 0.1, hyper = TRUE)
fit.sim <- PYdensity(y = ysim, mcmc = mcmc,
                     prior = prior, output = output)

# Explore the clustering structure
fit.sim.part <- partition(fit.sim)
sort(ftable(fit.sim.part$partitions[3, ]), decreasing = T)

dissmat <- dist(1-fit.sim.part$psm)
clus <- hclust(dissmat)
heatmap(as.matrix(dissmat), Rowv = as.dendrogram(hclust(dissmat)),
        Colv = NA, labRow = FALSE, labCol = FALSE)

pdf("dendro.pdf",7,5)
heatmap(as.matrix(dissmat), Rowv = as.dendrogram(hclust(dissmat)),
        Colv = NA, labRow = FALSE, labCol = FALSE)
dev.off()
#---------------------------------------
# Density regression example (Sec. 5.3)
#---------------------------------------

# Prior initialization and parameter setting
PYpar <- PYcalibrate(Ek = 3, n = nrow(CPP), discount = 0.25)
prior <- list(strength = PYpar$strength, discount = PYpar$discount)

grid_y <- seq(26, 47, length.out = 100)
grid_x <- round(quantile(CPP$dde, c(0,0.25,.5,.75,.95,.99)))
mcmc <- list(niter = 5000, nburn = 4000)
output = list(grid_x = grid_x, grid_y = grid_y, out_type = "FULL", out_param = TRUE)

# Fit the model
set.seed(42)
fit.reg <- PYregression(y = CPP$gest, x = CPP$dde, prior = prior,
                        mcmc = mcmc, output = output)

# Plot the results
regplot <- data.frame(
  dens = as.vector(apply(fit.reg$density, c(1,2), mean)),
  qlow = as.vector(apply(fit.reg$density, c(1,2),
                         quantile, probs = 0.025)),
  qupp = as.vector(apply(fit.reg$density, c(1,2),
                         quantile, probs = 0.975)),
  grid = rep(grid_y, 6),
  label = factor(rep(paste("DDE = ", grid_x), each = length(grid_y)),
                 level = rep(paste("DDE = ", grid_x))))
library(ggplot2)
ggplot(regplot) +
  theme_bw() +
  geom_line(data = regplot, map = aes(x = grid, y = dens)) +
  geom_ribbon(data = regplot, map = aes(x = grid, ymin = qlow,
                                      ymax = qupp), fill = "blue", alpha = 0.3) +
  facet_wrap(~label, ncol = 3, nrow = 2) +
  labs(x = "gestational age", y = "density")

#------------------------------------------------------
# Density estimation for correlated samples (Sec. 5.4)
#------------------------------------------------------

# Parameter setting
mcmc <- list(niter = 5000, nburn = 4000)
output <- list(grid = seq(30, 46, length.out = 100), out_type = "FULL")

# Estimate the model
fit.ddp <- DDPdensity(y = CPP$gest, group = CPP$hosp,
                        mcmc = mcmc, output = output)
fit.ddp$group <- factor(CPP$hosp, labels = paste("Hospital ", levels(CPP$hosp)))
plot(fit.ddp, wrap_dim = c(3,4), ylab = "density", xlab = "gestiational age")
