## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=6)

## ----load_libraries, warning=FALSE, message=FALSE-----------------------------
require(sparsevar)
require(fsMTS)
require(plot.matrix)
require(svMisc)
require(MTS)

## ----generate_structures------------------------------------------------------
k <- 5
max.lag <- 3
N <- 512
show.progress = F
genM <- function(k, sparsity_lag = 0.2){
  nonzero <- round(sparsity_lag*k*k)
  m = diag(k)*runif(k*k, min = 0, max=0.1)
  size <- k
  samples <- nonzero
  vals <- sample.int(size ^ 2, samples)
  pair<-cbind(vals %/% size + 1, vals %% size+1)
  pair[pair>k]<-k
  m[pair] <- runif(nonzero, min = -0.2, max=0.2)
  return(m)
}
set.seed(9)
pars <- list()
for (i in 1:max.lag){
  pars[[paste0("l",i)]] <- genM(k) 
}
res<-data.frame()
for (i in pars){
  res <-rbind(res,i)
}
mreal<-as.matrix(res)
plot(abs(mreal), digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Real structure of dependencies",
     xlab="Variables")

## ----generate_MTS-------------------------------------------------------------
mts.sim <- simulateVAR(N=k, p=max.lag, nobs=N,fixedMat=pars)
data <- mts.sim$series
colnames(data)<-paste0("V",1:ncol(data))

## ----generate_distances-------------------------------------------------------
shortest <- matrix(rexp(k*k, rate = 0.2), nrow=k)
shortest <- shortest-diag(k)*shortest
colnames(shortest)<-colnames(data)
rownames(shortest)<-colnames(data)
plot(shortest, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Shortest distances between nodes")

data <- scale(data)
data[5,]<-NA

## ----fsOwnlags----------------------------------------------------------------
mIndep<-fsMTS(data, max.lag=max.lag, method="ownlags")
plot(mIndep, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Only own lags")

## ----fsCCF--------------------------------------------------------------------
mCCF<-fsMTS(data, max.lag=max.lag, method="CCF",show.progress=show.progress)
plot(mCCF, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Cross-correlations")

## ----fsDistance---------------------------------------------------------------
mDistance<-fsMTS(data, max.lag=max.lag, method="distance", shortest = shortest, step = 1)
plot(mDistance, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Distance-based feature selection")

## ----fsGLASSO-----------------------------------------------------------------
mGLASSO.localized<-fsMTS(data, max.lag=max.lag,method="GLASSO", rho = 0.1,show.progress=show.progress, localized= TRUE)
plot(mGLASSO.localized, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Graphical LASSO-based (localized) feature selection")
mGLASSO.global<-fsMTS(data, max.lag=max.lag,method="GLASSO", rho = 0.1,show.progress=show.progress, localized = FALSE)
plot(mGLASSO.global, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Graphical LASSO-based (global) feature selection")

## ----fsLARS-------------------------------------------------------------------
mLARS<-fsMTS(data, max.lag=max.lag,method="LARS",show.progress=show.progress)
plot(mLARS, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Least angle  regression-based feature selection")

## ----fsRF---------------------------------------------------------------------
mRF.localized<-fsMTS(data, max.lag=max.lag,method="RF",show.progress=show.progress, localized = TRUE)
plot(mRF.localized, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Random forest-based (localized) feature selection")
mRF.global<-fsMTS(data, max.lag=max.lag,method="RF",show.progress=show.progress, localized = FALSE)
plot(mRF.global, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Random forest-based (global) feature selection")

## ----fsMI---------------------------------------------------------------------
mMI.localized<-fsMTS(data, max.lag=max.lag,method="MI",show.progress=show.progress, localized= TRUE)
plot(mMI.localized, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Mutual information-based (localized) feature selection")
mMI.global<-fsMTS(data, max.lag=max.lag,method="MI",show.progress=show.progress, localized= FALSE)
plot(mMI.global, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Mutual information-based (global) feature selection")

## ----fsPSC--------------------------------------------------------------------
mPSC<-fsMTS(data, max.lag=max.lag,method="PSC",show.progress=show.progress)
plot(mPSC, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="PSC-based (global) feature selection")

## ----fsEnsemble---------------------------------------------------------------
mlist <- list(Independent = mIndep,
              Distance = mDistance,
              CCF = mCCF,
              GLASSO.localized = mGLASSO.localized,
              GLASSO.global=mGLASSO.global,
              LARS = mLARS,
              RF.localized = mRF.localized,
              RF.global = mRF.global,
              MI.localized = mMI.localized,
              MI.global = mMI.global,
              PSC=mPSC)

th<-0.1
mE1 <- fsEnsemble(mlist, threshold = th, method="ranking")
plot(mE1, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Ensemble feature selection  using Ranking")
mlist[["EnsembleRank"]] <- mE1


mE2 <- fsEnsemble(mlist, threshold = th, method="majority")
plot(mE2, digits=2, col=rev(heat.colors(10)), key=NULL, 
     main="Ensemble feature selection  using Majority Voting")
mlist[["EnsembleMajV"]] <- mE2

## ----comparison, fig.width=9, fig.height=9------------------------------------
mlist[["Real"]] <- ifelse(mreal!=0,1,0)
msimilarity <- fsSimilarityMatrix(mlist, threshold=th, method="Kuncheva")
plot(msimilarity, digits=2, col=rev(heat.colors(ncol(msimilarity))), key=NULL,
     main="Pairwise comparison of feature sets", cex.axis=0.7)

## ----predict------------------------------------------------------------------
dat <- stats::na.omit(data)
model <- VAR(dat, p=max.lag, include.mean = F, fixed = mE2)
print(model$coef)
VARpred(model,1)

