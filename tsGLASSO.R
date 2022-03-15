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
##@thresh              - Thresholding the inverse spectral density
##@demean            - If TRUE, demeanizes data
##@RELTOL            - Tolerance for convergence
##@ABSTOL            - Tolerance for convergence
glasso_TimeSeries<-function(dta, lambda, halfWindowLength, L_init = NULL, 
                            rho = 100, alpha = 1, MAX_ITER = 100, 
                            thresh = 1e-5, ABSTOL   = 1e-5, 
                            demean = FALSE, RELTOL   = 1e-4, rho.flex = FALSE, 
                            diag = FALSE, smooth = FALSE)
{
   x       <- as.matrix(dta)
   p = dim(x)[2]
   n = dim(x)[1]
   x <-as.ts(x)
   K = 2 * halfWindowLength + 1       
   M = floor((n/ 2 - halfWindowLength - 1) / K)
   S = mysdf(dta = x, halfWindowLength = halfWindowLength, 
             smooth = smooth, demean = demean)
   ### Parameters for faster convergence of rho
   tau.iner = tau.decr = 2
   mu <- 10
   
   MAX_ITER = MAX_ITER
   
   p = dim(S)[1]
   F = dim(S)[3] ; 
   ## ADMM solver
   if (is.null(L_init)) {
      X = array(0,dim=c(p,p,F));
   } else{
      X = L_init
   }
   history.objval = history.r_norm = history.s_norm = 
   history.eps_pri = history.eps_dual = numeric(MAX_ITER)

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
         ev = eigen(rho * (Z[,,iter_f] -
                              U[,,iter_f]) - S[,,iter_f]) 
         Q <- ev$vectors
         L <- ev$values
         ind = order(L)
         es = L[ind]
         Q <- Q[,ind]
         xi = (es + sqrt(es^2 + 4*rho)) / (2*rho)
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
                     1 - lambda / (rho * norm(dmy, type="2")),0); 
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
            sum(abs(U[,,iter_f])^2))
      }
      
      history.objval[k]  = objective.n(S, X, Z, lambda,p,K);
      history.eps_pri[k] = F * p * ABSTOL + RELTOL * history.eps_pri[k]
      history.eps_dual[k] = F * p * ABSTOL + RELTOL * history.eps_dual[k]
      
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
   if(!is.null(thresh))
   {
      for(f in 1: F)
      {
         X[,,f] <- threshold(X[,,f], th = thresh)
      }
   }
   return <- list("history.objval" = history.objval[1:k],
                  "history.r_norm" = history.r_norm[1:k] ,
                  "history.s_norm" = history.s_norm[1:k],
                  "K" = K, "history.eps_dual" = history.eps_dual[1:k], 
                  "history.eps_pri" = history.eps_pri[1:k], 
                  "X"= X, "Z" = Z, 
                  "est_S" = tilde_S, 
                  "hat_S" = S,
                  "M" = M)
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
      sum_var = sum_var + (sum(diag(mtx1 %*% mtx2)) -  
                              log(determinant(mtx1)))  
   }
   
   obj = sum_var + lambda * sum_penalty 
   return(Re(obj))
}

shrinkage<-function(a, kappa)
{
   y = pmax((a-kappa),0) - pmax( (-a-kappa),0)
   return(y)
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


selectlambda <- function(X, 
                         lambdalist = NULL, nlam = 40, trim_min = 0.01,
                         nfreq = NULL, thresh = 1e-5,  trim_max = 1,
                         halfWindowLength , rho = 100, 
                         Max_Iter = 20, diag = FALSE, 
                         rho.flex = FALSE, alpha = 1.5,
                         freq = NULL, gam = 0.25, 
                         smooth = FALSE, demean = FALSE,
                         ABSTOL = 1e-7, RELTOL = 1e-6, 
                         criteria = c("BIC"), verbose = FALSE)
{
   crit = criteria
   xfreq   <- frequency(X)
   N       <- nrow(X) 
   nser    <- ncol(X)
   Nspec <- floor(N/2)   ## +1 for frequency start from 0 and and at 0.5
   freq  <- seq.int(from = (xfreq)/N, by = (xfreq)/N, length.out = Nspec)
   if (is.null(lambdalist)) {
      S_hat = mysdf(X, halfWindowLength)
      lambda_max = lammax(S_hat, trim_max)
      lambdalist = pathGen(nlam = nlam, lam_max = lambda_max,
                           flmin = trim_min)
   }
   if (criteria == "all") {
      ic =  matrix(0, nrow = length(lambdalist), ncol = 4)
      colnames(ic) = c("AIC", "BIC", "eBIC", "CV")
   } else {
      ic =  numeric(length(lambdalist))
   }
   L_init = NULL
   for (i in 1 : length(lambdalist))
   {
      gmodel <- glasso_TimeSeries(dta = dta, lambda = lambdalist[i], 
                                  rho = rho, alpha = alpha,
                                  MAX_ITER = Max_Iter, L_init = L_init,
                                  ABSTOL   = ABSTOL, 
                                  RELTOL   = RELTOL, 
                                  thresh = thresh,
                                  rho.flex = rho.flex, 
                                  diag = diag, demean = demean,
                                  halfWindowLength = halfWindowLength,
                                  smooth = smooth)
      #K = gmodel$interval
      est_isdf = gmodel$X  
      #S_tilde = gmodel$X
      L_init = est_isdf  
      #  fr = gmodel$fr
      K = gmodel$K
      fxx = gmodel$hat_S
      M   = gmodel$M
      #  print(isSymmetric(fxx[,,1]))
      for (f in 1: dim(est_isdf)[3])
      {
         nz = sum(est_isdf[,,f] != 0) - nser
         
         if (crit == "eBIC") {
            ic[i] = ic[i] +  
               2 * K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) -  
                          log(determinant(est_isdf[,,f]))) + 
               nz * log(2 * K * M) + 4 * nz * gam * log(nser)
         }
         else if (crit == "AIC") {
            ic[i] = ic[i] + 
               2* K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) - 
                         log(determinant(est_isdf[,,f]))) + 2*nz #dim(fxx)[3] *
         } else if (crit == "BIC") {
            ic[i] = ic[i] +  
               2 * K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) - 
                          log(determinant(est_isdf[,,f]))) + nz * log(2 * K * M)
         } else if (crit == "CV") {
            ic[i] = (ic[i] + cv_selection_f(S_hat = fxx,
                                            S_tilde = est_isdf, f = f)) / M
         } else if (crit == "all") { 
            ic[i,1] = ic[i,1] + 
               2* K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) - 
                         log(determinant(est_isdf[,,f]))) + 2*nz
            
            ic[i,2] =  ic[i,2] + 
               2* K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) - 
                         log(determinant(est_isdf[,,f]))) + nz * log(2 * K * M)
            
            ic[i,3] = ic[i,3] +  
               2 * K *(Re(sum(diag((fxx[,,f] %*% est_isdf[,,f])))) -  
                          log(determinant(est_isdf[,,f]))) + 
               nz * log(2 * K * M) + 4 * nz * gam * log(nser)
            
            ic[i, 4] = (ic[i,4] + cv_selection_f(S_hat = fxx,
                                            S_tilde = est_isdf, f = f)) / M
            
            
         }else {
            stop("Criteria should be one of AIC, BIC, eBIC, CV or all")
         }
      }
   }
   if(crit != "all") {
      ind = which.min(rev(Re(ic)))
      rev_ic = rev(Re(ic))
      lambda = rev(lambdalist)[ind]
      if(verbose) {
         plot(rev_ic~rev(lambdalist), xlab = "Threshold", 
              ylab = criteria, type = "l", col = "blue")
         abline(v = rev(lambdalist)[ind], col = "red", lwd = 3, lty = 2)
      }
   } else{
      rev_ic = apply(Re(ic), 2, rev)
      ind = apply(rev_ic, 2, which.min)
      lambda = cbind(rev(lambdalist)[ind[1]], rev(lambdalist)[ind[2]],
                 rev(lambdalist)[ind[3]],rev(lambdalist)[ind[4]])
      colnames(lambda) = c("AIC", "BIC", "eBIC", "CV")
      if (verbose) {
         warning("verbose is ignored when criteria is all")
      }
   }

cv_selection_f <- function(S_hat, S_tilde, f) {
   p = dim(S)[1]
   M = dim(S)[3]
   S = S_hat[,, -f]
   X = S_tilde[,, -f]
   F = M - 1
   for (iter_f in 1:F)
   {
      mtx1 = X[,,iter_f]  
      mtx2 = S[,,iter_f]  
      sum_var = sum_var + (sum(diag(mtx1 %*% mtx2)) -  
                              log(determinant(mtx1)))  
   }
   
   return(sum_var)
   
}
   
   # icresult <- glasso_TimeSeries(dta = dta, lambda = lambdalist[ind], 
   #                                  rho = rho, alpha = alpha, MAX_ITER = Max_Iter, 
   #                                  ABSTOL   = ABSTOL, RELTOL   = RELTOL,
   #                                  thr = thr, thresh = thresh,
   #                                  rho.flex = rho.flex, diag = diag, 
   #                                  halfWindowLength = halfWindowLength,
   #                                  smooth = smooth)
   # 
   return(list("lambda" = lambda,
               "lamlist" = rev(lambdalist),
               "IC" = rev_ic))
}

lammax <- function(S, trim_max = 1){
   #### This function calculates the max value in the tuning parameter list
   # such that the estimator L_{\lambda} is a diagonal matrix
   # NOTE: this is not necessarily true, but generally
   # a upper bound of the value we are looking for.
   
   # Args:
   #     S: the p-by-p-by M sample spectral density matrix
   
   p <- dim(S)[1]
   M <- dim(S)[3]
   sighat <- matrix(0, p-1, p - 1)
   for (i in seq(2, p)){
      for (j in seq(1,p - 1)){
         if (i != j){
            max_ind = which.max(abs(S[i, j, ]))
            sighat[i - 1, j] <- abs(S[i, j,max_ind ])/ sqrt(abs(S[j, j,max_ind]))
         }
      }
   }
   max(sighat) * trim_max
}

pathGen <- function(nlam, lam_max, flmin = 0.01){
   # Generate a path of lambda, with
   # nlam/2 decreasing exponentially
   # nlam/2 decreasing linearly
   # lam_max <- lammax(S)
   lamlist_lin <- lam_max * exp(seq(0, log(flmin), length = nlam/2))
   lamlist_exp <- seq(lam_max, lam_max*flmin, length.out = nlam/2)
   return(sort(unique(c(lamlist_lin, lamlist_exp)), decreasing = TRUE))
}


mysdf <- function(dta, halfWindowLength, smooth = FALSE, demean= TRUE) {
   x       <- as.matrix(dta)
   p = dim(x)[2]
   x <-as.ts(x)
   n = nrow(x)
   if (demean)
   {
      x <- scale(x, center = TRUE, scale = FALSE)
   }
   xfft  <- mvfft(x)/ sqrt(n) ## Convert to Fourier 
   xfft  <- xfft[2:(floor(n/2) - 1),]
   K = 2 * halfWindowLength + 1       
   M = floor((n/ 2 - halfWindowLength - 1) / K)
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
   return(S)
   
}

determinant<-function(x)
{
   values <- eigen(x,only.values = TRUE)$values
   det <- prod(as.matrix(values))
   return(det)
}

cv_selection_f <- function(S_hat, S_tilde, f) {
   p = dim(S_hat)[1]
   M = dim(S_hat)[3]
   S = S_hat[,, -f]
   X = S_tilde[,, -f]
   sum_var = 0

   for (iter_f in 1:(M - 1))
   {
      mtx1 = X[,,iter_f]  
      mtx2 = S[,,iter_f]  
      sum_var = sum_var + (sum(diag(mtx1 %*% mtx2)) -  
                              log(determinant(mtx1)))  
   }
   
   return(sum_var)
   
}
