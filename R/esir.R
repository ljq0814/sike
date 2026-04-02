#' Spatial Kendall's Tau Matrix
#'
#' Computes the spatial Kendall's tau matrix of a multivariate dataset,
#' with an optional rescaling by the number of variables.
#'
#' @param X A numeric matrix of observations with dimension \code{n x p}.
#' @param type A character string, either \code{"redefined"} (default) or
#'   \code{"standard"}. If \code{"redefined"}, the result is multiplied by
#'   \code{p}.
#'
#' @return A \code{p x p} symmetric matrix of spatial Kendall's tau estimates.
#'
#' @keywords internal
spkendalltau <- function(X,type = "redefined"){
    p <- ncol(X)
    res <- kendalltau(X)
    if(type == "redefined")res = res*p
    return(res)
}
#' Matrix Power
#'
#' Computes an arbitrary real power of a symmetric positive definite matrix
#' via eigendecomposition.
#'
#' @param mat A \code{p x p} symmetric positive definite matrix.
#' @param power A numeric scalar specifying the power to raise \code{mat} to.
#'
#' @return A \code{p x p} matrix equal to \eqn{A^{\text{power}}}.
#'
#' @keywords internal
matpower <- function(mat,power){
    eigmat <- eigen(mat)
    return(eigmat$vectors%*%diag((eigmat$values)^power)%*%t(eigmat$vectors))
}
#' Sort Data by Response
#'
#' Sorts the response vector and the corresponding rows of the predictor
#' matrix in ascending order of the response.
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{y}{The sorted response vector.}
#'     \item{X}{The predictor matrix with rows reordered accordingly.}
#'   }
#'
#' @keywords internal
sortdata <- function(y,X){
    n <- nrow(X)
    ord <- order(y)
    return(list(
        y = sort(y),
        X = X[ord,]
    ))
}
#' Update Location Estimate
#'
#' Updates the location vector in the HR divergence minimization algorithm
#' using the spatial sign function.
#'
#' @param resi A numeric matrix of residuals with dimension \code{n x p}.
#' @param sqB A \code{p x p} matrix representing the square root of the
#'   current scatter matrix estimate.
#'
#' @return A numeric vector of length \code{p} giving the updated location
#'   estimate.
#'
#' @keywords internal
updtloc <- function(resi, sqB){
    n = nrow(resi);p = ncol(resi)
    res = matrix(0,p,p)
    nume = rep(0,p)
    deno = 0.0;
    for (i in 1:n){
        nume = nume + spsign(resi[i,]);
        deno = deno + 1/Matrix::norm(resi[i,],"2")
    }
    return(sqB%*%nume/deno)
}
#' Update Scatter Matrix Estimate
#'
#' Updates the scatter matrix in the HR divergence minimization algorithm
#' using the spatial covariance of the residuals.
#'
#' @param resi A numeric matrix of residuals with dimension \code{n x p}.
#' @param sqB A \code{p x p} matrix representing the square root of the
#'   current scatter matrix estimate.
#'
#' @return A \code{p x p} symmetric positive definite matrix giving the
#'   updated scatter matrix estimate.
#'
#' @keywords internal
updtshape <- function(resi, sqB){
    n = nrow(resi);p = ncol(resi)
    return(sqB%*%spcov(resi)%*%sqB*p)
}

#' Scatter Matrix Estimation by Hettmansperger-Randles estimate
#'
#' Estimates the location vector and scatter matrix of the covariate matrix
#' by Hettmansperger-Randles iterative estimate. The algorithm initializes with
#' a spatial Kendall's tau covariance estimate and iterates until convergence.
#'
#' @param X A numeric matrix of observations with dimension \code{n x p}.
#' @param epsilon A positive numeric scalar specifying the convergence
#'   tolerance. Defaults to \code{5e-3}.
#' @param maxiter A positive integer specifying the maximum number of
#'   iterations. Defaults to \code{1000}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{location}{A numeric vector of length \code{p} giving the
#'       estimated location.}
#'     \item{scatter}{A \code{p x p} symmetric positive definite matrix
#'       giving the estimated scatter matrix.}
#'   }
#'
#' @export
hrstd <- function(X, epsilon = 5e-3, maxiter = 1000){
    n = nrow(X);p = ncol(X)
    initcov = spkendalltau(X)
    initloc = rep(0,p)
    criteria = 1E7
    i = 0
    resi = matrix(0,p,p)

    while (criteria > epsilon && i < maxiter){
        sqB = matpower(initcov,0.5)
        invsqB = matpower(initcov,-0.5)

        resi = (X - rep(1,n)%*%t(initloc))%*%invsqB
        criteria = Matrix::norm(spcov(resi)*p - diag(p), "f")
        initloc = updtloc(resi,sqB)
        initcov = updtshape(resi,sqB)
        initcov = initcov*p/sum(diag(initcov))
        i = i+1
    }

    return(list(
        location = initloc,
        scatter = initcov
        )
    )
}

#' Scatter Matrix Estimation for Elliptical Distributions by the variant of Hettmansperger-Randles estimate
#'
#' Estimates the scatter matrix of the covariate matrix by iteratively updating eigenvalues 
#' of the spatial Kendall's tau covariance matrix until convergence. See Li et al. (2026). 
#'
#' @param X A numeric matrix of observations with dimension \code{n x p}.
#' @param epsilon A positive numeric scalar specifying the convergence
#'   tolerance. Defaults to \code{5e-3}.
#' @param maxiter A positive integer specifying the maximum number of
#'   iterations. Defaults to \code{1000}.
#'
#' @return A \code{p x p} symmetric positive definite scatter matrix.
#'
#' @export
elptstd <- function(X, epsilon = 5e-3, maxiter = 1000){
  n = nrow(X);p = ncol(X)
  initcov = spkendalltau(X)
  eigmat <- eigen(initcov)
  initeigval <- eigmat$values

  ind <- which(initeigval <= 0.05)
  initeigval[ind] <- 0.5
  Z <- X%*%eigmat$vectors

  criteria = 1E7
  i = 0

  resi = diag(spkendalltau(Z%*%diag(1/sqrt(initeigval))))
  criteria = norm(resi - rep(1,p), "2")
  while (criteria > epsilon && i < maxiter){
    initeigval = sqrt(initeigval)*resi*sqrt(initeigval)
    initeigval = initeigval*p/sum(initeigval)
    resi = diag(spkendalltau(Z%*%diag(1/sqrt(initeigval))))
    criteria = norm(resi - rep(1,p), "2")
    i = i+1
  }

  return(
    eigmat$vectors%*%diag(initeigval)%*%t(eigmat$vectors)
  )
}

#' Sliced Inverse Kendall's Tau Estimation (SIKE)
#'
#' Performs sufficient dimension reduction by estimating the central subspace
#' using sliced inverse Kendall's tau estimation. Supports a standard mode and
#' a refinement mode that selects directions via distance correlation.
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param type A character string specifying the estimation method. Either
#'   \code{"standard"} (default) or \code{"refinement"}. Partial matching is
#'   supported.
#' @param eigenmed Either \code{NULL} (default) or a numeric scalar. If
#'   \code{NULL}, the median eigenvalue is computed separately within each
#'   slice. If a numeric value is provided, it is used as a fixed median
#'   eigenvalue across all slices.
#' @param Sig A \code{p x p} symmetric positive definite matrix used for
#'   standardization. Defaults to the identity matrix.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
sike <- function(y, X, K, type = c("standard", "refinement"), eigenmed = NULL, Sig = diag(ncol(X)), slices = 10){
  n <- nrow(X);p <- ncol(X)
  type = match.arg(type)
  stdZ <- X%*%matpower(Sig, -0.5)
  sortdt <- sortdata(y, stdZ)
  l = floor(n/slices)
  savemat <- matrix(0,p,p)
  if (tolower(type) == "refinement"){
    for (h in 1:slices){
      ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
      varZ <- spkendalltau(sortdt$X[ind,])
      savemat <- savemat + 1/slices*varZ
    }
    eigensave <- eigen(savemat)
    choose_arr <- c(1:(2*K))
    betak <- matpower(Sig, -0.5)%*%eigensave$vectors[,c(1:K,(p-K+1):p)]
    refine_res <- rep(0,2*K)
    for (j in 1:(2*K)){
      refine_res[j]<-abs(dcov(y,X%*%betak[,j]))
    }
    arr_chosen <- choose_arr[order(refine_res,decreasing = T)[1:K]]
    return(betak[,arr_chosen])
  } else {
    for (h in 1:slices){
      ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
      varZ <- spkendalltau(sortdt$X[ind,])
      eigenvarZ <- eigen(varZ)
      med <- 1
      if (is.null(eigenmed)){
        med <- median(eigenvarZ$values)
      } else {
        med <- eigenmed
      }
      savemat <- savemat + 1/slices*(med*diag(p) - varZ)%*%(med*diag(p) - varZ)
    }
    eigensave <- eigen(savemat)
    betak <- matpower(Sig, -0.5)%*%eigensave$vectors[,1:K]
  }
  return(betak)
}
#' Elliptical Sliced Inverse Regression (ESIR)
#'
#' Performs sufficient dimension reduction via sliced inverse regression under
#' an elliptical distribution assumption, using the spatial Kendall's tau
#' matrix for standardization instead of the sample covariance matrix. 
#' See Chen et al. (2022).
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
esir <- function(y, X, K, slices = 10){
    n <- nrow(X);p <- ncol(X)
    stdX <- X - rep(1,n)%*%t(apply(X,2,mean))
    M <- spkendalltau(X,type = "redefined")
    halfM <- matpower(M,-0.5)
    stdZ <- stdX%*%halfM
    sortdt <- sortdata(y, stdZ)
    l = floor(n/slices)
    mh = matrix(0,slices,p)
    for (h in 1:slices){
        ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
        mh[h,] <- apply(sortdt$X[ind,], 2, mean)
    }
    Mm <- spkendalltau(mh,type = "redefined")
    eigenMm <- eigen(Mm)
    betak <- halfM%*%eigenMm$vectors[,1:K]
    return(betak)
}

#' Sliced Inverse Regression (SIR)
#'
#' Performs sufficient dimension reduction via the classical sliced inverse
#' regression method using the sample covariance matrix for standardization.
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
sir <- function(y,X,K,slices = 10){
    n <- nrow(X);p <- ncol(X)
    stdX <- X - rep(1,n)%*%t(apply(X,2,mean))
    M <- cov(X)
    halfM <- matpower(M,-0.5)
    stdZ <- stdX%*%halfM
    sortdt <- sortdata(y, stdZ)
    l = floor(n/slices)
    Mm <- matrix(0,p,p)
    for (h in 1:slices){
        ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
        mh <- apply(sortdt$X[ind,], 2, mean)
        Mm <- Mm + mh%*%t(mh)/l
    }
    eigenMm <- eigen(Mm)
    betak <- halfM%*%eigenMm$vectors[,1:K]
    return(betak)
}
#' Sliced Average Variance Estimation (SAVE)
#'
#' Performs sufficient dimension reduction via the sliced average variance
#' estimation method, which captures symmetric and nonlinear dependencies
#' between the response and predictors.
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
slicedave <- function(y,X,K,slices = 10){
    n <- nrow(X);p <- ncol(X)
    stdX <- X - rep(1,n)%*%t(apply(X,2,mean))
    M <- cov(X)
    halfM <- matpower(M,-0.5)
    stdZ <- stdX%*%halfM
    sortdt <- sortdata(y, stdZ)
    l = floor(n/slices)
    savemat <- matrix(0,p,p)
    for (h in 1:slices){
        ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
        varZ <- var(sortdt$X[ind,])
        savemat <- savemat + 1/slices*1/l*(diag(p) - varZ)%*%(diag(p) - varZ)
    }
    eigensave <- eigen(savemat)
    betak <- halfM%*%eigensave$vectors[,1:K]
    return(betak)
}
#' Kernel Weight Function
#'
#' Computes a scalar kernel weight for a given vector, defined as
#' \eqn{w(x) = \|x\|_2 / (1 + \|x\|_2^2)}.
#'
#' @param x A numeric vector.
#'
#' @return A non-negative numeric scalar giving the kernel weight.
#'
#' @keywords internal
k_fun <- function(x)
{
  z <- norm(x,"2")
  return(z/(1 + z^2))
}
#' Weighted Inverse Regression Estimation (WIRE)
#'
#' Performs sufficient dimension reduction via weighted inverse regression,
#' See Dong et al.(2015).
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
wire <- function(y,X,K,slices = 10){
  n <- nrow(X);p <- ncol(X)
  res <- robustbase::covMcd(X)
  mu_hat <- res$center
  Gamma_hat <- res$cov
  S <- matpower(Gamma_hat,-0.5)
  Xtilde <- (X - rep(1,n)%*%t(mu_hat))%*%S
  sortdt <- sortdata(y, Xtilde)
  Xtilde <- sortdt$X

  #2.
  mhat <- matrix(0,slices,p)
  multi <- floor(n/slices)
  for(h in 1:slices)
  {
    a <- (h-1) * multi + 1
    b <- h * multi
    Xtilde_h <- Xtilde[a:b,]
    weight <- apply(Xtilde_h,1,k_fun)
    mhat[h,] <- t(weight) %*% Xtilde_h / multi
  }
  #3.
  Vhat<-cov(mhat)
  #4.
  res<-eigen(Vhat)
  Beta <- S %*% res$vectors[,1:K]
  return(Beta)
}
#' Spatial Inverse Median Estimation (SIME)
#'
#' Performs sufficient dimension reduction via sliced inverse regression with
#' slice means replaced by L1 medians. See Dong et al.(2015).
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
sime <- function(y,X,K,slices = 10){
  n <- nrow(X);p <- ncol(X)
  res <- robustbase::covMcd(X)
  mu_hat <- res$center
  Gamma_hat <- res$cov
  S <- matpower(Gamma_hat,-0.5)
  Xtilde <- (X - rep(1,n)%*%t(mu_hat))%*%S
  sortdt <- sortdata(y, Xtilde)
  Xtilde <- sortdt$X

  #2.
  mhat <- matrix(0,slices,p)
  multi <- floor(n/slices)
  for(h in 1:slices)
  {
    a <- (h-1) * multi + 1
    b <- h * multi
    Xtilde_h <- Xtilde[a:b,]
    mhat[h,] <- pcaPP::l1median(Xtilde_h)
  }
  #3.
  Vhat<-cov(mhat)
  #4.
  res<-eigen(Vhat)
  Beta <- S%*%res$vectors[,1:K]
  return(Beta)
}
#' Contour Sliced Inverse Regression (CSIR)
#'
#' Performs sufficient dimension reduction via contour sliced inverse
#' regression, which projects predictors onto the unit sphere prior to
#' slicing. Location and scatter are estimated via HR divergence
#' minimization using \code{hrstd()}. See Luo et al.(2009).
#'
#' @param y A numeric vector of response values of length \code{n}.
#' @param X A numeric matrix of predictors with dimension \code{n x p}.
#' @param K A positive integer specifying the dimension of the central subspace.
#' @param slices A positive integer specifying the number of slices. Defaults
#'   to \code{10}.
#'
#' @return A \code{p x K} matrix whose columns are the estimated basis vectors
#'   of the central subspace.
#'
#' @export
csir <- function(y,X,K,slices = 10){
  n <- nrow(X);p <- ncol(X)
  hrls <- hrstd(X)
  loc_es <- hrls$location
  sig_es <- hrls$scatter*p/sum(diag(hrls$scatter))
  stdX <- (X - rep(1,n)%*%t(loc_es))%*%matpower(sig_es,-0.5)
  vector_len <- apply(stdX,1,norm,type = "2")
  contourX <- stdX/(vector_len%*%t(rep(1,p)))
  sortdt <- sortdata(y, contourX)
  l = floor(n/slices)
  Mm <- matrix(0,p,p)
  overmean <- apply(sortdt$X, 2, mean)
  for (h in 1:slices){
    ind <- if (h < slices) ((h-1)*l+1):(h*l) else ((h-1)*l+1):n
    mh <- apply(sortdt$X[ind,], 2, mean)
    Mm <- Mm + (mh-overmean)%*%t(mh-overmean)/l
  }
  eigenMm <- eigen(Mm)
  betak <- eigenMm$vectors[,1:K]
  return(betak)
}
