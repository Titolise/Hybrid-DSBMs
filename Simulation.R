rm(list = ls())

require(mclust)
require(gtools)
require(Rcpp)
require(RcppArmadillo)

source("./src/est_hyb_sbm_dyn_dec.R")
source("./src/draw_sn_dyn.R")
source("./src/best_perm.R")

# Compile and load the C++ code
sourceCpp("./Cpp/Core_hyb.cpp")

# Simulation settings
mod     = 1    
persist = "high"   
n       = 20       
B       = 2       
TT      = 6               
k       = 3
nrep    = 15
maxit   = 5000

# Initial Probabilities
piv0 = rep(1/k, k)

if (persist == "high") {
  rho = 0.1
} else {
  rho = 0.3
}

Pi0 = matrix(0, k, k)
for (u in 1:k) {
  for (v in 1:k) {
    Pi0[u, v] = rho^abs(v - u)
  }
}
Pi0 = diag(1 / rowSums(Pi0)) %*% Pi0

# Model Selection
if (mod == 1) {
  set.seed(6) 
  # High intra-groups / Low inter-groups
  Psi0 = matrix(0.03, k, k)
  grid = 0.3 * runif(k, 0.5, 1.5)
  diag(Psi0) = grid
} else if (mod == 2) {
  set.seed(6)
  # Low intra-groups / High inter-groups
  Psi0 = matrix(0.3, k, k)
  grid = 0.03 * runif(k, 0.5, 1.5)
  diag(Psi0) = grid
} 

filename = sprintf(
  ".Res_mod%g_n%g_k%g_T%g_rep%g_perist_%s.RData",
  mod,
  n,
  k,
  TT,
  nrep,
  persist
)


#---- Simulation ----
ari_pred       = rep(0, B)
ari_predv      = rep(0, B)
sum_diag_pred  = rep(0, B)
sum_diag_predv = rep(0, B)

out = list()

for (b in 1:B) {
  try({
    print("%%%")
    cat(sprintf("\nIteration: %d / %d \n", b, B))
    set.seed(b + 136532)
    
    res   = draw_sn_dyn(n, k, TT, piv0, Pi0, Psi0)
    Y     = res$Y
    Utrue = res$U + 1
    
    # Initialization via K-means across flattened slices
    Tau = array(0, c(k, n, TT))
    YY  = Y[,, 1]
    for (t in 2:TT) {
      YY = rbind(YY, Y[,, t])
    }
    
    U = matrix(kmeans(YY, k, nstart = 100)$cl, n, TT)
    for (t in 1:TT) {
      for (i in 1:n) {
        Tau[U[i, t], i, t] = 1
      }
    }
    
    # Diagonal set to NA for network adjacencies
    YY = Y
    for (t in 1:TT) {
      diag(YY[,, t]) = NA
    }
    
    out[[b]] = list()
    out[[b]]$sim   = res
    out[[b]]$Utrue = Utrue
    
    # Hybrid SBM - Initial Fit
    print("Fitting Hybrid DSBM...")
    pred = est_pred_sbm_dyn_decoding(
      YY,
      k,
      start = 2,
      tol   = 10^-10,
      Tau   = Tau,
      maxit = maxit
    )
    
    ari_pred[b]  = adjustedRandIndex(pred$cl, Utrue)
    ari_predv[b] = adjustedRandIndex(pred$clv, Utrue)
    
    tmp = best_perm(Utrue, pred$cl)
    sum_diag_pred[b] = sum(diag(tmp$Tab1))
    
    tmp = best_perm(Utrue, pred$clv)
    sum_diag_predv[b] = sum(diag(tmp$Tab1))
    
    out[[b]]$pred         = pred
    out[[b]]$pred_lktrace = pred$lk
    
    # Random initializations via perturbation
    if (nrep > 0) {
      # Base clustering via spectral decomposition of mean adjacency
      Y_mean = apply(YY, c(1, 2), mean, na.rm = TRUE)
      Y_mean[is.na(Y_mean)] = 0
      vecs    = svd(Y_mean)$u[, 1:k]
      cl_base = kmeans(vecs, k, nstart = 50)$cluster
      
      for (h in 1:nrep) {
        perturb_rate = runif(1, 0.05, 0.35)
        Taur = array(0.01, c(k, n, TT))
        
        for (t in 1:TT) {
          for (i in 1:n) {
            assigned_class = cl_base[i]
            if (runif(1) < perturb_rate) {
              alt_class = setdiff(1:k, cl_base[i])
              if (length(alt_class) > 0) {
                assigned_class = sample(alt_class, 1)
              }
            }
            Taur[assigned_class, i, t] = 0.98
          }
        }
        
        for (t in 1:TT) {
          Taur[,, t] = sweep(Taur[,, t], 2, colSums(Taur[,, t]), "/")
        }
        
        cat(sprintf(
          "Random init %d/%d : Hybrid SBM (Perturb: %.0f%%)\n",
          h, nrep, perturb_rate * 100
        ))
        
        predh = est_pred_sbm_dyn_decoding(
          YY,
          k,
          start = 2,
          tol   = 10^-10,
          Tau   = Taur,
          maxit = maxit
        )
        
        out[[b]]$pred_lktrace = c(out[[b]]$pred_lktrace, predh$lk)
        out[[b]]$pred_arivtrace = c(
          out[[b]]$pred_arivtrace,
          adjustedRandIndex(predh$clv, Utrue)
        )
        out[[b]]$pred_aritrace = c(
          out[[b]]$pred_aritrace,
          adjustedRandIndex(predh$cl, Utrue)
        )
        
        # Keep best fit based on log-likelihood
        if (predh$lk > out[[b]]$pred$lk) {
          out[[b]]$pred = predh
          ari_pred[b]   = adjustedRandIndex(predh$cl, Utrue)
          ari_predv[b]  = adjustedRandIndex(predh$clv, Utrue)
          
          tmp = best_perm(Utrue, predh$cl)
          sum_diag_pred[b] = sum(diag(tmp$Tab1))
          
          tmp = best_perm(Utrue, predh$clv)
          sum_diag_predv[b] = sum(diag(tmp$Tab1))
        }
      }
    }
  })
}

save.image(filename)