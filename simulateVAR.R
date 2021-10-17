## This function is borrowed from Davis et.al 2015
## https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1092978


simulateVAR = function(coefMatr=NULL,intercept,size,burn=100,Sigma,error.dist=c("normal","t"),df=NULL,checkCausal=TRUE){
      error.dist  = match.arg(error.dist)
      if(!is.null(coefMatr)) {coefMatr = as.matrix(coefMatr)}
      inteVect    = intercept
      #	inteVect    = matrix(inteVect,ncol=1) Command was in actual code
      K           = nrow(Sigma)
      #  print(K)
      if(error.dist=="t"){
            if(is.null(df)|df<=2) {stop("error: the d.f. of t-distribution is NULL or <=2")}
      }
      if (min(eigen(Sigma)$val)<0){stop("error: the covariance matrix is not positive-definite")}
      if (!is.null(coefMatr)){ 
            if(checkCausal & !companionVAR(coefMatr)$isCausal){stop("error: the VAR model is non-causal")}
      }
      S.scale        = diag(NA,K) 
      if(error.dist=="t") {S.scale = Sigma*(df-2)/df}
      if (!is.null(coefMatr)) {p = ncol(coefMatr)/K} else {p=0}
      simData        = matrix(NA,nrow=K,ncol=size+burn)
      cholSigma      = t(chol(Sigma))
      if(p!=0){
            simData[,1:p]  = matrix(rnorm(K*p),nrow=K,ncol=p)
            for (i in (p+1):(size+burn)){
                  if (error.dist=="normal") {error = cholSigma %*% matrix(rnorm(K),ncol=1)}
                  if (error.dist=="t")      {error = matrix(rmt(n=1,rep(0,K),S.scale,df),ncol=1)}
                  simData[,i] = inteVect + coefMatr %*% matrix(c(simData[,(i-1):(i-p)]),ncol=1) + error
            }}
      if(p==0){
            for (i in (p+1):(size+burn)){
                  if (error.dist=="normal") {error = cholSigma %*% matrix(rnorm(K),ncol=1)}
                  if (error.dist=="t")      {error = matrix(rmt(n=1,rep(0,K),S.scale,df),ncol=1)}
                  simData[,i] = inteVect + error
            }}
      simData  = t(simData[,(1+burn):(size+burn)])
      result   = list(simData=simData,coefMatr=coefMatr,inteVect=inteVect,K=K,p=p,size=size,checkCausal=checkCausal,Sigma=Sigma,error.dist=error.dist,df=df)
      return(result)
}

#------------------------------------------------------------------------------------------------------------------------
# companionVAR(A)
# compute the companion matrix of the VAR model's coefficient matrix A and decide whether or not the VAR model is causal.
#------------------------------------------------------------------------------------------------------------------------
companionVAR = function(A){ 
      A         = as.matrix(A)
      K         = nrow(A)
      p         = ncol(A)/K
      companion = matrix(0, nrow=K*p, ncol=K*p)
      companion[1:K,1:(K*p)] = A
      if (p > 1) {companion[(K+1):(K*p),1:(K*p-K)] = diag(1,K*p-K)}
      eigens    = eigen(companion)$val
      mods      = Mod(eigens)
      result    = list(coefMatr=A, isCausal = max(mods) < 1, companion=companion, eigens=eigens, mods=mods)
      return(result)
}

