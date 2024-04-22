
innerprod_q2 <- function(
  q1,
  q2
){
  T1 = ncol(q1)
  val = sum(q1*q2)/T1
  return(val)
}


find_rotation_seed_unique <- function(
  q1, q2,
  mode = "O",
  alignment = TRUE,
  rotation = FALSE,
  scale = FALSE,
  norm_ratio = 1.0,
  lambda = 0.0
) {
  L <- nrow(q1)
  M <- ncol(q1)

  alignment = TRUE
  rotation = FALSE
  scale = FALSE

  # Variables for DPQ2 algorithm
  grd <- seq(0, 1, length.out = M)
  nbhd_dim <- 7L
  Gvec <- rep(0, M)
  Tvec <- rep(0, M)
  size <- 0

  scl <- 4
  minE <- Inf

  end_idx <- 0

  for (ctr in 0:end_idx) {
    q2n <- q2

    Rbest <- diag(nrow(q2n))

    if (norm(q1 - q2n, 'F') > 0.0001) {

      q1i <- q1 / sqrt(innerprod_q2(q1, q1))
      q2ni <- q2n / sqrt(innerprod_q2(q2n, q2n))

      dim(q1i) <- M * L
      dim(q2ni) <- M * L

      ret <- .Call(
        "DPQ2", PACKAGE = "fdasrvf",
        q1i, grd, q2ni, grd, L, M, M, grd, grd, M,
        M, Gvec, Tvec, size, lambda, nbhd_dim
      )
      Gvec <- ret$G[1:ret$size]
      Tvec <- ret$T[1:ret$size]
      gamI <- stats::approx(Tvec, Gvec, xout = grd)$y

      gam <- (gamI - gamI[1]) / (gamI[length(gamI)] - gamI[1])
      q2new <- group_action_by_gamma(q2n, gam, scale = scale)

    } else {
      q2new <- q2n
      gam <- seq(0, 1, length.out = M)
    }

    Ec <- amplitude_distance(q1, q2new, scale = scale, norm_ratio = norm_ratio)

    if (Ec < minE) {
      Rbest <- Rbest
      q2best <- q2new
      gambest <- gam
      minE <- Ec
      tau <- scl * ctr
    }
  }

  list(
    q2best = q2best,
    Rbest = Rbest,
    gambest = gambest,
    tau = tau,
    d = minE
  )
}

amplitude_distance <- function(
  q1,
  q2,
  scale = FALSE,
  norm_ratio = 1
) {
  if (scale) {
    q1dotq2 <- innerprod_q2(q1, q2)
    if (q1dotq2 >  1) q1dotq2 <-  1
    if (q1dotq2 < -1) q1dotq2 <- -1
    d <- sqrt(acos(q1dotq2) ^ 2 + log(norm_ratio) ^ 2)
  } else {
    v <- q1 - q2
    d <- sqrt(innerprod_q2(v, v))
  }
  d
}


group_action_by_gamma <- function(q, gamma, scale = TRUE) {
  L <- nrow(q)
  M <- ncol(q)
  grd <- seq(0, 1, length.out = M)
  gammadot <- fdasrvf::gradient(gamma, 1.0 / M)
  qn <- matrix(nrow = L, ncol = M)

  for (l in 1:L)
    qn[l, ] <- stats::spline(grd, q[l, ], xout = gamma)$y * sqrt(gammadot)

  if (scale)
    qn <- qn / sqrt(innerprod_q2(qn, qn))

  qn
}


group_action_by_gamma_coord <- function(
    f,  # matrix of sizes LxM: L-dimensional curve, sampled in M points
    gamma  # vector of length M: warping
) {
  L <- nrow(f)
  M <- ncol(f)
  fn <- matrix(nrow = L, ncol = M)
  grd <- seq(0, 1, length.out = M)

  for (l in 1:L)
    fn[l, ] <- stats::spline(grd, f[l, ], xout = gamma)$y

  fn
}


pvecnorm = function ( ###########
    v,
    p = 2
) {
  sum(abs(v)^p)^(1/p)
}


curve_to_srvf <- function(
  beta,
  scale = TRUE
) {
  centroid <- calculatecentroid(beta)
  beta <- beta - centroid
  out <- fdasrvf::curve_to_q(beta, scale = scale)
  list(q = out$q, qnorm = out$lenq, centroid = centroid)
}


calculatecentroid <- function(
    beta,
    returnlength = F
){
  n = nrow(beta)
  T1 = ncol(beta)

  betadot = apply(beta,1,fdasrvf::gradient,1.0/(T1-1))
  betadot = t(betadot)

  normbetadot = apply(betadot,2,pvecnorm,2)
  integrand = matrix(0, n, T1)
  for (i in 1:T1){
    integrand[,i] = beta[,i] * normbetadot[i]
  }

  scale = integrate(seq(0,1,length.out=T1), normbetadot)
  centroid = apply(integrand,1,integrate,x = seq(0,1,length.out=T1))/scale
  if(returnlength)  return(list("length" = scale,"centroid" = centroid))
  return(centroid)
}


## R equivalent of repmat (matlab)
repmat <- function( #############
    X,  # matrix
    m,  # nrow multiply
    n  # ncol multiply
) {
  mx = dim(X)[1]
  if (is.null(mx)){
    mx = 1
    nx = length(X)
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else {
    nx = dim(X)[2]
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }

  return(mat)
}


# Utilizzate per sqrt_weighted_mean_inverse :

exp_map <- function(
    psi,
    v,
    wnorm = l2_norm
){
  v_norm <- wnorm(v)
  expgam <- cos(v_norm) * psi + sin(v_norm) * v / v_norm
  return(expgam)
}

l2_norm <- function(
    psi,
    time = seq(0, 1, length.out = length(psi))
){
  l2norm <- sqrt(integrate(time, psi*psi))
  return(l2norm)
}


inner_product <- function(
    psi1,
    psi2,
    time = seq(0, 1, length.out = length(psi1))
){
  ip <- integrate(time,psi1*psi2)
  return(ip)
}

inv_exp_map <- function(
    Psi,
    psi
){
  ip <- inner_product(Psi, psi)
  if(ip < -1){
    ip = -1
  }else if(ip > 1){
    ip = 1
  }
  theta <- acos(ip)

  if (theta < 1e-10){
    exp_inv = rep(0, length(psi))
  } else {
    exp_inv = theta / sin(theta) * (psi-cos(theta)*Psi)
  }
  return(exp_inv)
}
