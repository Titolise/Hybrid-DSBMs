#' Estimation of the Hybrid DSBM model via the CEM algorithm
#'
#' @param Y 3D binary array (n x n x TT) representing the temporal adjacency matrix.
#' @param k Number of latent blocks/classes.
#' @param start Initialization strategy (0 = K-means, 1 = random, 2 = custom Tau).
#' @param tol Convergence tolerance for the log-likelihood.
#' @param maxit Maximum number of CEM iterations.
#' @param Tau Initial membership array (required if start = 2).
#' @return List containing parameter estimates, likelihood, and classifications (cl and clv).
#' @importFrom stats kmeans
#' @export

est_pred_sbm_dyn_decoding <- function(Y, k, start=0, tol=10^-10, maxit=100, Tau=NULL){
  n = dim(Y)[1]
  TT = dim(Y)[3]
  
  
  Y[is.na(Y)] = 0
  
  # Starting
  if(start==0){
    Tau = array(0, c(k,n,TT))
    YY = Y[,,1]
    for(t in 2:TT) YY = rbind(YY,Y[,,t])
    U = matrix(kmeans(YY, k, nstart=100)$cl,n,TT)
    for(t in 1:TT) for(i in 1:n) Tau[U[i,t],i,t]  = 1 
  } else if(start==1){
    U = matrix(sample(1:k,n*TT,replace=TRUE),n,TT)
    Tau = array(0, c(k,n,TT))
    for(t in 1:TT) for(i in 1:n) Tau[U[i,t],i,t]  = 1
  } else if(start==2){
    if(is.null(Tau)) stop("initial value of Tau must be given in input")
  }		
  
  Ub = matrix(0,n,TT)
  for(t in 1:TT) Ub[,t] = apply(Tau[,,t],2,which.max)
  
  # Compute the parameters
  Psi = matrix(NA, k,k)
  for(u in 1:k) for(v in 1:k){
    num = den = 0
    for(i in 1:(n-1)){ 
      ind = (i+1):n 
      for(t in 1:TT){
        num = num + Tau[u,i,t] * sum(Tau[v,ind,t] * Y[i,ind,t]) + Tau[v,i,t] * sum(Tau[u,ind,t] * Y[ind,i,t])
        den = den + Tau[u,i,t] * sum(Tau[v,ind,t]) + Tau[v,i,t] * sum(Tau[u,ind,t])
      }
    }
    Psi[u,v] = num/den
  }
  
  piv = rep(0,k)
  for(u in 1:k) piv[u] = sum(Tau[u,,1])/n
  
  Pi = matrix(0,k,k)
  for(u in 1:k){
    for(v in 1:k) for(t in 2:TT) Pi[u,v] = Pi[u,v]+sum(Tau[u,,t-1]*Tau[v,,t])
    Pi[u,] = Pi[u,]/sum(Pi[u,])
  }
  
  out = complk_dyn_cpp(Y, piv, Pi, Psi, Ub, k, n, TT)
  lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  
  cat("------------|-------------|-------------|-------------|\n")
  cat("  iteration |   classes   |      lk     |    lk-lko   |\n")
  cat("------------|-------------|-------------|-------------|\n")
  cat(sprintf("%11g", c(0, k, lk, NA)), "\n", sep = " | ")
  
  lk0 = lk; it = 0
  
  # Staring CEM alg
  while((abs(lk-lk0)/abs(lk0) > tol | it==0) & it < maxit){
    lk0 = lk; it = it+1
    
    # E-step 
    V = array(0,c(n,k,TT)); U = array(0,c(k,k,TT))
    Yvp = matrix(1/pv,n,k)
    M = matrix(1,n,k)
    V[,,TT] = Yvp*L[,,TT]
    U[,,TT] = (t(L[,,TT-1])%*%(Yvp*Phi[,,TT]))*Pi
    if(TT>2){
      for(t in seq(TT-1,2,-1)){
        M = (Phi[,,t+1]*M)%*%t(Pi);
        V[,,t] = Yvp*L[,,t]*M
        U[,,t] = (t(L[,,t-1])%*%(Yvp*Phi[,,t]*M))*Pi
      }
    }		
    M = (Phi[,,2]*M)%*%t(Pi)
    V[,,1] = Yvp*L[,,1]*M
    
    # M-step
    piv = colSums(V[,,1])/n
    Ut = apply(U[,,2:TT],c(1,2),sum)
    Pi = diag(1/rowSums(Ut))%*%Ut
    
    # Asymmetric Update
    Num = matrix(0,k,k)
    Den = matrix(0,k,k)
    
    for(t in 1:TT){
      Vt = V[,,t]           
      Ut_hard = matrix(0, n, k)
      for(i in 1:n) Ut_hard[i, Ub[i,t]] = 1 
      
      Num = Num + t(Vt) %*% Y[,,t] %*% Ut_hard
      Den = Den + outer(colSums(Vt), colSums(Ut_hard)) - t(Vt) %*% Ut_hard
    }
    
    Num = Num+t(Num)
    Den = pmax(Den+t(Den),10^-300)
    Psi = Num/Den
    
    # Improve classification  C++
    out = greedy_search_cpp_original(Y, Ub, Psi, Pi, piv)
    lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv; Ub = out$Ub
    
    if(it/1==floor(it/1)) cat(sprintf("%11g", c(it, k, lk, lk -lk0)), "\n", sep = " | ")
  }
  if(it/10>floor(it/10)) cat(sprintf("%11g", c(it, k, lk, lk -lk0)), "\n", sep = " | ")
  cat("------------|-------------|-------------|-------------|\n")
  
  # Improve classification with Viterbi
  R = L; Ubv = matrix(0,n,TT)
  for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*Pi[,u])
  if(n==1) Ubv[,TT] = which.max(R[,,TT])
  else Ubv[,TT] = apply(R[,,TT],1,which.max)
  for(i in 1:n) for(t in seq(TT-1,1,-1)) Ubv[i,t] = which.max(R[i,,t]*Pi[,Ubv[i,t+1]])
  if(n==1) Ubv = as.vector(Ubv)
  
  out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,cl=Ub,clv=Ubv,it=it)
  return(out)
}	