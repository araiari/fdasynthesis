#' Karcher (weighted) Mean of multidimensional functions
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic
#' square-root velocity (srvf) framework.
#'
#' @param beta Array of sizes \eqn{M \times L \times N} describing \eqn{N}
#' curves of dimension \eqn{L} evaluated on \eqn{M} points
#' @param maxit Maximum number of iterations
#' @param wts Vector of weights for computing the weighted Karcher
#' mean (default = `NULL`)
#'
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{beta}{centered data}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' \item{cent}{centers of original curves}
#' \item{len}{length of curves}
#' \item{len_q}{length of srvfs}
#' \item{mean_scale}{mean length}
#' \item{mean_scale_q}{mean length srvf}
#' \item{E}{energy}
#' \item{qun}{cost function}
weighted_karcher_mean <- function (beta, lambda = 0.0, maxit = 20, wts = NULL)
{

  if (is.null(wts))
    wts = rep(1, N)

  tmp = dim(beta)
  L = tmp[1]
  M = tmp[2]
  N = tmp[3]
  q = array(0, dim = tmp)
  len = rep(0, N)
  len_q = rep(0, N)
  cent = matrix(0, nrow = L, ncol = N)
  for (ii in 1:N) {
    beta1 = beta[ , , ii]
    centroid1 = calculate_centroid(beta = beta1)
    cent[ , ii] = -1 * centroid1
    dim(centroid1) = c(length(centroid1), 1)
    beta1 = beta1 - repmat(centroid1, 1, M)
    beta[ , , ii] = beta1
    out = fdasrvf::curve_to_q(beta = beta1)
    q[, , ii] = out$q
    len[ii] = out$len
    len_q[ii] = out$lenq
  }

  mu = q[, , 1] # initialize centroid as the first observation
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  v = array(0, c(L, M, N))
  normvbar = rep(0, maxit + 1)


  cat("\nInitializing...\n")
  gam = matrix(0, M, N)
  for (k in 1:N) {
    out = find_rotation_seed_unqiue(mu, q[, , k], mode = "O", rotated = FALSE,
                                    lambda)
    gam[, k] = out$gambest
  }

  gam = t(gam)
  gamI = fdasrvf::SqrtMeanInverse(t(gam))
  bmu = group_action_by_gamma_coord(bmu, gamI) #in utils_curve -> beta * warp
  mu = fdasrvf::curve_to_q(bmu)$q
  mu[is.nan(mu)] = 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu / sqrt(innerprod_q2(mu, mu)) #in utils_curve

    for (i in 1:N) {
      q1 = q[ , , i]

      out = find_rotation_seed_unqiue(mu, q1, mode = "O", rotated = FALSE,
                                      lambda)
      qn_t = out$q2best / sqrt(innerprod_q2(out$q2best, out$q2best)) #in utils_curve

      q1dotq2 = innerprod_q2(mu, qn_t) #in utils_curve

      if (q1dotq2 > 1){
        q1dotq2 = 1
      }
      if (q1dotq2 < -1){
        q1dotq2 = -1
      }

      dist = acos(q1dotq2)

      u = qn_t - q1dotq2 * q1 # shooting vector di q1 rispetto a qn_t
      normu = sqrt(innerprod_q2(u, u))
      if (normu > 1e-4){
        w = u * acos(q1dotq2)/normu # perché non c'è il seno a denominatore?
      } else {
        w = matrix(0, L, M)
      }

      v[, , i] = w

      sumd[itr + 1] = sumd[itr + 1] + dist^2
    }

    # if(ms == "median"){#run for median only
    #   sumv = rowSums(v_d, dims = 2)
    #   sum_dinv = sum(1/d_i)
    #   vbar = sumv/sum_dinv
    # }
    # else{ #run for mean only
    # sumv = rowSums(v, dims = 2) ##* reso weighted
    # vbar = sumv/N ##*
    v_w = array(sapply(1:2, function(n) v[, , n] * wts[n]), dim = dim(v))
    vbar = rowSums(v_w, dims = 2) / sum(wts)
    # }

    normvbar[itr] = sqrt(innerprod_q2(vbar, vbar)) #in utils_curve
    normv = normvbar[itr]
    if ((sumd[itr]-sumd[itr+1]) < 0){
      break
    } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
      mu = cos(delta * normvbar[itr]) * mu +
        sin(delta *  normvbar[itr]) * vbar/normvbar[itr] # dallo shooting vector alla media
      x = fdasrvf::q_to_curve(mu)
      a = -1 * calculate_centroid(x)
      dim(a) = c(length(a), 1)
      betamean = x + repmat(a, 1, M) #non ci andrebbe un -??
    }
    else {
      break
    }
    itr = itr + 1
  }

  if (scale){
    mean_scale = prod(len)^(1/length(len))
    mean_scale_q = prod(len_q)^(1/length(len))
    betamean = mean_scale*betamean
  }

  ifelse(weighted, type <- "Weighted Karcher mean", type <- "Karcher mean")
  return(list(beta = beta, mu = mu, type = type, betamean = betamean,
              v = v, q = q, E = normvbar[1:itr], cent = cent, len = len,
              len_q = len_q, qun = sumd[1:itr]))
}


#' Karcher Mean of Curves
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic
#' square-root velocity (srvf) framework.
#'
#' @param beta Array of sizes \eqn{n \times T \times N} describing \eqn{N}
#' curves of dimension \eqn{n} evaluated on \eqn{T} points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotated Optimize over rotation (default = `TRUE`)
#' @param scale Include scale (default = `FALSE`)
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit maximum number of iterations
#' @param wts Vector of weights for computing the weighted Karcher
#' mean (default = `NULL`)
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{beta}{centered data}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' \item{cent}{centers of original curves}
#' \item{len}{length of curves}
#' \item{len_q}{length of srvfs}
#' \item{mean_scale}{mean length}
#' \item{mean_scale_q}{mean length srvf}
#' \item{E}{energy}
#' \item{qun}{cost function}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
curve_weighted_karcher_mean <- function (beta, mode = "O", rotated = TRUE,
                                         scale = FALSE, lambda = 0.0,
                                         maxit = 20, wts = NULL)
  {
  mean_scale = NA
  mean_scale_q = NA
  tmp = dim(beta)
  n = tmp[1] # dimension of the curves/functions
  T1 = tmp[2] # n. evaluation points of the curves/functions
  N = tmp[3] # n. curves/functions
  q = array(0, c(n, T1, N))
  len = rep(0, N)
  len_q = rep(0, N)
  cent = matrix(0, n, N)
  for (ii in 1:N) {
    beta1 = beta[ , , ii]
    centroid1 = calculatecentroid(beta1)
    cent[ , ii] = -1 * centroid1
    dim(centroid1) = c(length(centroid1), 1)
    beta1 = beta1 - repmat(centroid1, 1, T1)
    beta[ , , ii] = beta1
    out = fdasrvf::curve_to_q(beta1)
    q[, , ii] = out$q
    len[ii] = out$len
    len_q[ii] = out$lenq
  }

  mu = q[, , 1]
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  v = array(0, c(n, T1, N))
  normvbar = rep(0, maxit + 1)


  cat("\nInitializing...\n")
  gam = matrix(0,T1,N)
  for (k in 1:N) {
    out = find_rotation_seed_unqiue(mu,q[, , k],mode,rotated,lambda)
    gam[, k] = out$gambest
  }

  gam = t(gam)
  gamI = fdasrvf::SqrtMeanInverse(t(gam))
  bmu = group_action_by_gamma_coord(bmu, gamI)
  mu = fdasrvf::curve_to_q(bmu)$q
  mu[is.nan(mu)] <- 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu/sqrt(innerprod_q2(mu, mu))

    if (mode == "C"){
      basis = find_basis_normal(mu)
    }

    for (i in 1:N) {
      q1 = q[, , i]

      out = find_rotation_seed_unqiue(mu, q1, mode, rotated, lambda)
      qn_t = out$q2best / sqrt(innerprod_q2(out$q2best, out$q2best))

      q1dotq2 = innerprod_q2(mu,qn_t)

      if (q1dotq2 > 1){
        q1dotq2 = 1
      }
      if (q1dotq2 < -1){
        q1dotq2 = -1
      }

      dist = acos(q1dotq2)

      u = qn_t - q1dotq2 * q1 # shooting vector di q1 rispetto a qn_t
      normu = sqrt(innerprod_q2(u, u))
      if (normu > 1e-4){
        w = u * acos(q1dotq2) / normu # perché non c'è il seno a denominatore?
      } else {
        w = matrix(0, nrow(beta1), T1)
      }

      if (mode=="O"){
        v[, , i] = w
      } else {
        v[, , i] = project_tangent(w, q1, basis)
      }

      sumd[itr + 1] = sumd[itr + 1] + dist^2
    }


    if (is.null(wts)){
      sumv = rowSums(v, dims = 2)
      vbar = sumv/N
    }
    else {
      v_w = array(sapply(1:2, \(n) v[, , n] * wts[n]), dim = dim(v))
      vbar = rowSums(v_w, dims = 2) / sum(wts) #weighted mean
    }


    normvbar[itr] = sqrt(innerprod_q2(vbar, vbar))
    normv = normvbar[itr]
    if ((sumd[itr]-sumd[itr+1]) < 0){
      break
    } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
      mu = cos(delta * normvbar[itr]) * mu +
        sin(delta * normvbar[itr]) * vbar/normvbar[itr] # dallo shooting vector alla media
      if (mode == "C") {
        mu = project_curve(mu)
      }
      x = fdasrvf::q_to_curve(mu)
      a = -1 * calculatecentroid(x)
      dim(a) = c(length(a), 1)
      betamean = x + repmat(a, 1, T1)  #non ci andrebbe un -??
    }
    else {
      break
    }
    itr = itr + 1
  }

  if (scale){
    mean_scale = prod(len)^(1/length(len))
    mean_scale_q = prod(len_q)^(1/length(len))
    betamean = mean_scale*betamean
  }

  ifelse(weighted, type <- "Weighted Karcher mean", type <- "Karcher mean")
  return(list(beta = beta, mu = mu, type = type, betamean = betamean,
              v = v, q = q, E = normvbar[1:itr], cent = cent, len = len,
              len_q = len_q, qun = sumd[1:itr]))
}
