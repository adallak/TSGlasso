# TS GLASSO as described in Dallakyan et. al 2021 (https://arxiv.org/pdf/2107.01659.pdf) 
# Tugnait 2018 (https://ieeexplore.ieee.org/document/8645324)
# Jung et. al 2015 (https://research.cs.aalto.fi/MLBigDat/papers/AJung_LSP15_Graphical_LASSO_Model_Selection.pdf)
## Graphical modeling of high-dimensional time series
##@dta               - Time Series Data
##@rho               - ADMM convergence tuning parameter
##@halfWindowLength  - length for halfwindow 
##@lambda            - ADMM penalty parameter
##@alpha             - additional relaxation for ADMM
##@MAX_ITER          - Max iteration
##@rho.flex          - If TRUE, rho is chosen adaptivaly 
##@smooth            - If TRUE, smoothes the sample estimate of periodagram
##@diag              - If TRUE, penalizes diagonal elements
##@thr               - Thresholding the inverse spectral density
##@demean            - If TRUE, demeanizes data
##@RELTOL            - Tolerance for convergence
##@ABSTOL            - Tolerance for convergence
glasso_TimeSeries<-function(dta, lambda, halfWindowLength , 
                            rho = 100, alpha = 1, MAX_ITER = 100, 
                            thresh=TRUE, thr = 1e-5, ABSTOL   = 1e-5, 
                            demean = FALSE, RELTOL   = 1e-4, rho.flex = FALSE, 
                            diag = FALSE, smooth = FALSE)
{
      x       <- as.matrix(dta)
      p = dim(x)[2]
      x <-as.ts(x)
      if (demean)
      {
        x <- scale(x, center = TRUE, scale = FALSE)
      }
      xfft  <- mvfft(x)/ sqrt(dim(x)[1]) ## Convert to Fourier 
      xfft  <- xfft[2:dim(x)[1],]
      K = 2 * halfWindowLength + 1       
      M = floor((nrow(x)/ 2 - halfWindowLength - 1) / K)
      S = array(0, dim = c(p, p, M))
      if( M <= 1)
      {
            stop("halfWindowLength is too large")
      }
      
      for(f in 1 : M)
      {
            d_freq = (f - 1) * K + halfWindowLength + 1  + (
               -halfWindowLength : halfWindowLength )
            for(wind in -halfWindowLength : halfWindowLength)
            {
              if(smooth)
              {
                kernel <- kernel("modified.daniell",halfWindowLength)
                S[, , f] = S[, , f] + kernel[halfWindowLength] * (xfft[((f - 1) * K 
                                          + halfWindowLength + 1 + wind), ] %*% 
                                         Conj(t(xfft[((f - 1) * K + halfWindowLength 
                                                      + 1 + wind), ]))) / K
              }else{
                S[, , f] = S[, , f] + (xfft[((f - 1) * K + halfWindowLength
                                             + 1 + wind), ] %*% Conj(
                                                t(xfft[((f - 1) * K + halfWindowLength 
                                                + 1 + wind), ]))) / K
              }
                 
            }

      }
      ### Parameters for faster convergence of rho
      tau.iner = tau.decr = 2
      mu <- 10
      
      MAX_ITER = MAX_ITER
      
      p = dim(S)[1]
      F = dim(S)[3] ; 
      ## ADMM solver
      history.objval = history.r_norm = history.s_norm = 
      history.eps_pri = history.eps_dual = numeric(MAX_ITER)
      X = array(0,dim=c(p,p,F));
      tilde_S = array(0, dim = c(p, p, F))
      Z = array(0,dim=c(p,p,F))
      U = array(0,dim=c(p,p,F))
      A =array(0,dim=c(p,p,F))
      
      dmy = array(0,dim=c(F,1))
      
      
      for (k in 1:MAX_ITER)
      {
            for (iter_f in 1:F)
            {
                  # x-update
                  ev = eigen(rho/ K * (Z[,,iter_f] -
                                          U[,,iter_f]) - S[,,iter_f]) 
                  Q <- ev$vectors
                  L <- ev$values
                  ind = order(L)
                  es = L[ind]
                  Q <- Q[,ind]
                  xi = K * (es + sqrt(es^2 + 4*rho / K)) / (2*rho)
                  X[,,iter_f] = Q %*% diag(xi) %*% Conj(t(Q))
                  tilde_S[,,iter_f] = (Q) %*% diag(1/xi) %*% Conj(t(Q))
            }
            
            weight = array(1,dim=c(p,p))
            # z-update 
            Zold = Z
            X_hat = alpha * X + (1 - alpha) * Zold
            A = X_hat + U 
            for (iter_idx1 in 1:p)
            {
                  for (iter_idx2 in 1:p)
                  {
                        if(diag)
                        {
                              dmy = A[iter_idx1, iter_idx2,]; 
                              weight[iter_idx1, iter_idx2] = max(1 - 
                                       lambda / (rho * norm(dmy, type="2")),0); 
                        }else{
                              if(iter_idx1 != iter_idx2)
                              { 
                                    dmy = A[iter_idx1, iter_idx2,]
                                    weight[iter_idx1, iter_idx2] = max(
                                       1 - lambda / (rho * sqrt(sum(Mod(dmy)))),0); 
                              }else{
                                    A[iter_idx1, iter_idx2, ] = X[iter_idx1, iter_idx2, ]
                              }
                        }
                  }
            }
            for (iter_f in 1:F) {
                  Z[,,iter_f] = A[,,iter_f] * weight  
                  U[,,iter_f] = U[,,iter_f] + (X_hat[,,iter_f] - Z[,,iter_f])
                  
                  history.r_norm[k] = history.r_norm[k] + 
                     sqrt(sum(Mod(X[, ,iter_f] - Z[, , iter_f])^2))
                  
                  history.s_norm[k]  = rho * sqrt(
                     sum(abs((Z[, , iter_f] - Zold[,,iter_f]))^2))
                  
                  history.eps_pri[k] = history.eps_pri[k] + max(
                     sqrt(sum(abs(X[,,iter_f])^2)), sqrt(sum(abs(Z[,,iter_f])^2)),0)
                  
                  history.eps_dual[k] = history.eps_dual[k] +  sqrt(
                     sum(abs(rho * U[,,iter_f])^2))
            }
            
            history.objval[k]  = objective.n(S, X, Z, lambda,p,K);
            history.eps_pri[k] =F * sqrt(p) * ABSTOL + RELTOL * history.eps_pri[k]
            history.eps_dual[k] = F * sqrt(p) * ABSTOL + RELTOL * history.eps_dual[k];
            
            if ((history.r_norm[k] < history.eps_pri[k]) & (
               history.s_norm[k] < history.eps_dual[k]))
                  break
            
            if (rho.flex)
            {
                  if (history.r_norm[k] > mu * history.s_norm[k])
                  {
                        rho = tau.iner * rho
                        U   = U / tau.iner
                  }else if (history.r_norm[k] * mu < history.s_norm[k])
                  {
                        rho = rho / tau.decr
                        U   = U * tau.decr
                  } else {
                        rho = rho
                  }
            } 
            
      }
      if(thresh)
      {
            for(f in 1: F)
            {
              X[,,f] <- threshold(X[,,f], th = thr)
            }
      }
      return <- list("history.objval" = history.objval[1:k],
                     "history.r_norm" = history.r_norm[1:k] ,
                     "history.s_norm" = history.s_norm[1:k],
                     "K" = K, "history.eps_dual" = history.eps_dual[1:k], 
                     "history.eps_pri" = history.eps_pri[1:k], 
                     "X"= X, "Z" = Z, 
                     "est_S" = tilde_S, "hat_S" = S)
      return(return)
}


objective.n = function(S, X, Z, lambda, p, K)
{
      F = dim(X)[3]
      sum_var = 0 ; 
      mtx1 = array(0,dim=c(p,p)) 
      mtx2 = array(0,dim=c(p,p))  
      dmy = array(0,dim=c(F,1)) 
      sum_penalty =0 ; 
      
      for (iter_i in 1:p) 
      {
            for (iter_j in 1:p)
            {
                  if (iter_i != iter_j) 
                  {   
                        dmy = X[iter_i,iter_j,]  
                        sum_penalty = sum_penalty + norm(dmy, type="2") 
                  }
            }
      }  
      
      for (iter_f in 1:F)
      {
            mtx1 = X[,,iter_f]  
            mtx2 = S[,,iter_f]  
            sum_var = sum_var + K * (sum(diag(mtx1 %*% mtx2)) - log(determinant(mtx1)))  
      }
      
      obj = sum_var + lambda * sum_penalty 
      return(Re(obj))
}

shrinkage<-function(a, kappa)
{
      y = pmax((a-kappa),0) - pmax( (-a-kappa),0)
      return(y)
}

determinant<-function(x)
{
      values <- eigen(x,only.values = TRUE)$values
      det <- prod(as.matrix(values))
      return(det)
}


threshold<-function(x, th = 1e-6)
 {
   x[abs(Re(x)) <= th] = 0
   return(x)
 }


##############################################
########## COMPARE time series graph #######

comparetsg <- function (esttsg, truetsg) 
{
   #### This function compares esitmated
   ### adjacency matrix of A with the true Adj matrix 
   ### of coefficient matrix A
   
   p = dim(esttsg)[2]
   f = dim(esttsg)[3]
   
   #### Construct adjacency matrices
   estAdj = trueAdj = matrix(0, p, p)
   for (iter_i in 1:p) 
   {
      for (iter_j in 1:p)
      {
         if (iter_i != iter_j) 
         {   
            dmy_est = esttsg[iter_i,iter_j,] 
            dmy_true = truetsg[iter_i,iter_j,] 
            if(all(Mod(dmy_est) != 0))
            {
               estAdj[iter_i, iter_j] = 1
            }
            if(all(Mod(dmy_true) != 0))
            {
               trueAdj[iter_i, iter_j] = 1
            }
            
         }
      }
   }  
   
   
   ml <- estAdj
   mt <- trueAdj
   p <- dim(ml)[2]
   mt[mt != 0] <- rep(1, sum(mt != 0))
   ml[ml != 0] <- rep(1, sum(ml != 0))
   diffm <- ml - mt
   nmbTrueGaps <- (sum(mt == 0) - p)/2
   fpr <- if (nmbTrueGaps == 0) 
      1
   else (sum(diffm > 0)/2)/nmbTrueGaps
   diffm2 <- mt - ml
   nmbTrueEdges <- (sum(mt == 1)/2)
   tpr <- if (nmbTrueEdges == 0) 
      0
   else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
   trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
   tdr <- if (sum(ml == 1) == 0) {
      if (trueEstEdges == 0) 
         1
      else 0
   }
   else trueEstEdges/(sum(ml == 1)/2)
   c(tpr = tpr, fpr = fpr, tdr = tdr)
}


