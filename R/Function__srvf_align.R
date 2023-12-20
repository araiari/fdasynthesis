#' Alignment
#'
#' This function aligns functions using the elastic square-root
#' slope function (SRSF) framework.
#'
#' @param f Either a numeric matrix or a numeric 3D array specifying the
#'   functions that need to be jointly clustered and aligned.
#'   - If a matrix, it must be of shape \eqn{M \times N}. In this case, it is
#'   interpreted as a sample of \eqn{N} curves observed on a grid of size
#'   \eqn{M}.
#'   - If a 3D array, it must be of shape \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   the curves are evaluated.
#' @param wgt A numeric vector of length \eqn{N} specifying the weights for
#'   each curve to be considered during the averaging process. If NULL, then
#'   all curves are set to have the same weight = 1.
#' @param seeds An integer vector of length `K` specifying the indices of the
#'   curves in `f` which will be chosen as initial centroids. Defaults to `NULL`
#'   in which case such indices are randomly chosen.
#' @param centroid_type A string specifying the type of centroid to compute.
#'   Choices are `"mean"` or `"medoid"`. Defaults to `"mean"`.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param smooth_data A boolean specifying whether to smooth data using a box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of box filters applied.
#'   Defaults to `25L`.
#' @param parallel A boolean specifying whether parallel mode (using
#'   [foreach::foreach()] and the **doParallel** package) should be activated.
#'   Defaults to `FALSE`.
#' @param alignment A boolean specifying whether to perform alignment. Defaults
#'   to `TRUE`.
#' @param omethod A string specifying which method should be used to solve the
#'   optimization problem that provides estimated warping functions. Choices are
#'   `"DP"`, `"DP2"` or `"RBFGS"`. Defaults to `"DP"`.
#' @param use_verbose A boolean specifying whether to display information about
#'   the calculations in the console. Defaults to `FALSE`.
#'
#' @return An object of class `fdakma` which is a list containing:
#'
#' - `f0`: the original functions;
#' - `q0`: the original SRSFs;
#' - `fn`: the aligned functions as matrices or a 3D arrays of the same shape
#' than `f0` by clusters in a list;
#' - `qn`: the aligned SRSFs as matrices or a 3D arrays of the same shape
#' than `f0` separated in clusters in a list;
#' - `template`: the centroids in the original functional space;
#' - `template.q`: the centroids in SRSF space;
#' - `template_start` : the initial point `template(0)`
#' - `gam`: the warping functions as matrices or a 3D arrays of the same shape
#' than `f0` by clusters in a list;
#' - `qun`: cost function value.
#'
#' @importFrom foreach %dopar%

srvf_align = function (
    f, time,
    wgt = NULL,
    seeds = NULL,
    centroid_type = c("mean","medoid"),
    lambda = 0,
    smooth_data = FALSE,
    sparam = 25L,
    parallel = FALSE,
    alignment = TRUE,
    omethod = c("DP", "DP2", "RBFGS"),
    use_verbose = FALSE
){

  #Initialization of all things
  omethod <- rlang::arg_match(omethod)
  centroid_type <- rlang::arg_match(centroid_type)
  if (parallel) {
    cores <- max(parallel::detectCores() - 1, 1)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  else {
    foreach::registerDoSEQ()
  }
  dims <- dim(f)
  if (length(dims) == 2) {
    dim(f) <- c(1, dims)
    dims <- dim(f)
  }
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]
  if (is.null(wgt))
    wgt = rep(1, N)


  #Initialization of the template of f
  if (is.null(seeds))
    template.ind <- sample(1:N,1)
  else template.ind <- seeds
  template <- matrix(f[, , template.ind], nrow=L, ncol=M) #matrix of size L x M

  #Initialization: functions f and srsf q
  if (smooth_data) {
    for (l in 1:L) f[l, , ] <- fdasrvf::smooth.data(f[l, , ], sparam = sparam)
  }
  q <- fdasrvf::f_to_srvf(f, time, multidimensional = (L > 1))

  #Initialization of the template of q
  template.q <- matrix(q[, , template.ind], nrow=L, ncol=M)


    if (use_verbose)
      cli::cli_alert_info("Running iteration {itr}...")
    if (use_verbose)
      cli::cli_alert_info("----> Alignment step")
    gam <- matrix(0, nrow = M, ncol = N)
    qn <- array(0, dim = c(L, M, N))
    fn <- array(0, dim = c(L, M, N))
    fw <- matrix(nrow = L, ncol = M)
    template_start <- matrix(0, nrow = L, ncol = 1)


    # Find the optimal warping for each function wrt the template
    # Warp the functions to fn, get their srsf qn
    n = NULL
    outfor <- foreach::foreach(n = 1:N, .combine = cbind,
                               .packages = "fdasrvf") %dopar% {
                                 if (alignment) {
                                   gam_tmp <- fdasrvf::optimum.reparam(
                                     Q1 = template.q, T1 = time, Q2 = matrix(q[, , n], nrow=L, ncol=M),
                                     T2 = time,
                                      lambda = lambda, method = omethod, w = 0,
                                      f1o = template[, 1], f2o = f[, 1, n])
                                 }
                                 else gam_tmp <- seq(0, 1, length.out = M)
                                 for (l in 1:L) {
                                   fw[l, ] <- stats::approx(x = time, y = f[l,
                                                                            , n], xout = (time[M] - time[1]) * gam_tmp +
                                                              time[1])$y
                                 }
                                 qw <- fdasrvf::f_to_srvf(fw, time, multidimensional = (L >
                                                                                 1))
                                 dist <- 0
                                 for (l in 1:L) {
                                   dist <- dist + integrate(time, (qw[l, ] - template.q[l,])^2)
                                 }
                                 dist <- sqrt(dist)

                                 list(gam_tmp, fw, qw, dist)
                               }
      gam <- do.call(cbind, outfor[1, ]) #output1: optimal gamma
      f_temp <- unlist(outfor[2, ]) #output2 : aligned f
      dim(f_temp) <- c(L, M, N)
      q_temp <- unlist(outfor[3, ]) #output3: (q,gamma)
      dim(q_temp) <- c(L, M, N)
      qn <- q_temp
      fn <- f_temp
      Dy <- unlist(outfor[4, ]) #output4 : dist


    # Among the equivalent class of warpings, use the optimal one
    if (use_verbose)
      cli::cli_alert_info("----> Normalisation step")
    ftmp <- fn[, , , drop = FALSE]
    gamtmp <- gam[, , drop = FALSE]
    gamI <- fdasrvf::SqrtMeanInverse(gamtmp)
    fw <- matrix(nrow = L, ncol = M)
    outfor <- foreach::foreach(n = 1:N, .combine = cbind,
                               .packages = "fdasrvf") %dopar% {
                                 for (l in 1:L) {
                                   fw[l, ] <- stats::approx(x = time, y = ftmp[l,
                                                                               , n], xout = (time[M] - time[1]) * gamI +
                                                              time[1])$y
                                 }
                                 qw <- fdasrvf::f_to_srvf(fw, time, multidimensional = (L >
                                                                                 1))
                                 gamt1 <- stats::approx(x = time, y = gamtmp[,
                                                                             n], xout = (time[M] - time[1]) * gamI + time[1])$y
                                 list(gamt1, fw, qw)
                               }
    gam <- do.call(cbind, outfor[1, ])
    f_temp <- unlist(outfor[2, ])
    dim(f_temp) <- c(L, M, N)
    q_temp <- unlist(outfor[3, ])
    dim(q_temp) <- c(L, M, N)
    qn <- q_temp
    fn <- f_temp

    #Identify the template
    if (use_verbose)
      cli::cli_alert_info("----> Template identification step")
    qun.t <- 0
    old.template.q <- template.q


      if (centroid_type == "mean") {

        wgt_temp = matrix(rep(wgt, M), nrow=M, byrow=T)
        wgt_temp2 = matrix(rep(wgt, L), nrow=L, byrow=T)
        for (l in 1:L) {
          template.q[l, ] <- rowSums(wgt_temp * qn[l, , ]) / sum(wgt)
        }

        template_start = rowSums(wgt_temp2 * f[ , 1, ]) / sum(wgt)
        template = srvf_to_f(
          template.q, time, template_start, multidimensional = T
          )
      }
      else if (centroid_type == "medoid") {
        idx <- which.min(Dy)
        for (l in 1:L) {
          template.q[l, ] <- qn[l, , idx]
          template[l, ] <- fn[l, , idx]
        }
      }

  if (use_verbose)
    cli::cli_alert_info("Consolidating output...")
  ftmp <- qtmp <- gamtmp <- list()
  out <- list(f0 = f, q0 = q, time = time, fn = fn,
              qn = qn, gam = gam,  template = template,
              template.q = template.q, template_start = template_start,
              lambda = lambda, omethod = omethod)
  if (parallel)
    parallel::stopCluster(cl)
  out
}
