source("simulateVAR.R")
source("tsGLASSO.R")
K = 6
size = 128
burn = 100
halfWindowLength = 6
A_exp = matrix(c(0, 0.5, 0.5,
                 0, 0, 0.3,
                 0, 0.25, 0.5), 3, 3, byrow = TRUE)



sigma_u = matrix(c(18, 0, 6,
                   0, 1, 0,
                   6, 0 ,3), 3 ,3)

sigma = matrix(0, 6, 6)
sigma[1:3,1:3] = sigma_u
sigma[4:6, 4:6] = sigma_u / 1.5

A_6 = matrix(0, 6, 6)
A_6[1:3, 1:3] = A_exp
A_6[4:6, 4:6] = A_exp / 1.5
A_6[6,2]= 0.5
A_6[1,4] = 0.2
th_u = solve(sigma_u)
A = A_6
## simulate data from VAR process
set.seed(1234)
dta = simulateVAR(coefMatr = A, intercept = rep(0,K), size=size, 
                  burn = burn, Sigma = sigma, error.dist="normal")$simData
## run time series graphical lasso
TS = glasso_TimeSeries(dta, lambda = 2.5, halfWindowLength , 
                       rho = 100, alpha = 1.5, MAX_ITER = 100, 
                       thresh=TRUE, thr = 1e-5, rho.flex = TRUE, 
                       smooth = TRUE)
plot(TS$history.objval, type = 'l', main = "Objective value")
