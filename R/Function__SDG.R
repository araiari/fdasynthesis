#' Generate a synthetic dataset
#'
#' @description
#' Given a set of original functions, generate a synthetic dataset.
#'
#' @details
#' It can handle multidimensional functions of \eqn{L} dimensions.
#'
#' `is_constrained` is used to set constraints for the values of the synthetic
#' functions. For each dimension of the function, it states whether values are
#' (strictly) positive or (strictly) monotone. Accepted values are: `pos` (positive functions),
#' `strict-pos` (non-negative functions), `mon` (for monotone increasing functions),
#' `strict-mon` (for strictly-monotone increasing functions)
#'
#' Remark: if `D_list` is not provided, the elastic distances are computed.
#' This may slow down the process of data generation.
#'
#' Remark: if `w_a` is not provided, it is computed as the parameter
#' allowing for the optimal compromise between amplitude and phase distances.
#' It is computed using the package `SPARTAAS`.
#' The computation may slow down the process of data generation.
#'
#'
#' @param fun_array array of sizes \eqn{L \times M \times N} : contains the values
#' of the \eqn{N} original functions of dimension \eqn{L} evaluated on \eqn{M} points.
#' @param i_try integer vector : subset of the curves to be synthetized
#' @param time vector of size \eqn{M} : specifies the grid on which the curves are
#' evaluated. If `NULL`, it is set as an uniform grid in \eqn{[0,1]}.
#' @param D_list list of two elements, each of which is either a dist objector a square
#' matrix, which must have the same dimension \eqn{N} : represent the elastic
#' distances of (respectively) phase and amplitude between the functions.
#' @param w_a scalar with values in \eqn{[0,1]} : weights of the
#' amplitude distance in the linear combination of the elastic distances.
#' @param K scalar integer : number of neighbors to be considered.
#' @param alpha_0 scalar, with \eqn{\alpha_0 >0} : concentation parameter of the
#' Dirichlet distribution from which the weights of neighbors are sampled.
#' @param beta_0 scalar, with \eqn{\beta_0 > 0} : rate parameter
#' weighting the inverse distances in the weight computation.
#' @param add_noise boolean : whether noise should be added to the weighted mean.
#'  Defaults to `FALSE`.
#' @param is_constrained vector of size \eqn{L} : each element indicates any constraint of
#' the dimension \eqn{l = 1,\dots,L}. It is considered only if `add_noise = TRUE` (see Details)
#' @param cv scalar with \eqn{cv > 0} : coefficient of variation of the covariance
#' amplitude with respect to the mean squared norm. It is the input of the function add_noise.
#' @param clust_labels vector of size \eqn{N} : specifies the labels
#' of the functions. If provided, the \eqn{K} neighbors are searched among the
#' subset of functions with the same label. Defaults to `NULL`.
#' @param use_verbose boolean : specifying whether to display information about
#'   the calculations in the console. Defaults to `FALSE`.
#'
#' @return a list of the following elements:
#' - `fun_s_array` array of sizes \eqn{L \times M \times N} : contains the values
#' of the \eqn{N} synthetic functions of dimension \eqn{L} evaluated on \eqn{M}
#' points.
#'
#' @export
SDG = function (
    fun_array, i_try = NULL, time = NULL, D_list = NULL, w_a = NULL, K = 5, alpha_0 = 1,
    beta_0 = 1, add_noise = F, is_constrained = NULL, cv = 0,
    clust_labels = NULL, use_verbose = F
) {


  dims = dim(fun_array)
  if (length(dims) == 2) { #if matrix -> unidimensional functions
    fun_array = array(fun_array, dim=c(1,dims))
    dims = dim(fun_array)
  }
  L = dims[1]
  M = dims[2]
  N = dims[3]

  if (is.null(time))
    time = seq(0,1,length.out=M)
  if (length(time) != M)
    cli::cli_abort("Error: time must be of length M")

  if (is.null(clust_labels))
    cluster_search = FALSE
  else {
    if (length(clust_labels) != N)
      cli::cli_abort("Error: cluster labels must be of length N")
    cluster_search = TRUE
  }
  if (!is.null(is_constrained) &
      sum(is_constrained %in% c("real", "pos","strict-neg","mon","strict-mon"))!= L)
    cli::cli_abort("Error: is_constrained accept only values among pos, strict-pos, real")
  do_positive = (sum(is_constrained %in% c("pos","strict-pos")) > 0)
  is_constrained_2 = factor(is_constrained,
                            labels = c("real", "real","real","pos","strict-pos"),
                            levels=c("real", "pos","strict-neg","mon","strict-mon"))

  if (is.null(i_try))
    i_try = 1:N


  if (do_positive) {
    for (l in which(is_constrained %in% c("pos","strict-pos"))){
      if (is_constrained[l] == "pos")
        MU[l,] = sqrt(MU[l,])
      if (is_constrained[l] == "strict-pos")
        MU[l,] = log(MU[l,])
    }
  }

  # Computation of the elastic distances
  if (is.null(D_list)) {
    if (use_verbose)
      cli::cli_alert_info("Starting the computation of elastic distances")
    D_list = compute_elastic_distance_N(fun_array, time)
    if (use_verbose)
      cli::cli_alert_info("End of the computation of elastic distances")
  }

  # Computation of the total distance as convex combination of elastic distances
  if (is.null(w_a) && requireNamespace("SPARTAAS", quietly = TRUE)){
    if (use_verbose)
      cli::cli_alert_info("Starting the computation of the optimal w_a")
    w_a = estimate_optimal_w_a (D1 = D_list$Dy_tot, D2 = D_list$Dx_tot)
    if (use_verbose)
      cli::cli_alert_info("End of the computation of the optimal w_a")
  }
  if (is.null(w_a) && !requireNamespace("SPARTAAS", quietly = TRUE))
    cli::cli_abort("Error: parameter {.arg w_a} not provided and {.pkg SPARTAAS} not installed.")
  Dtot = w_a * as.matrix(D_list$Dy_tot) + (1-w_a) * as.matrix(D_list$Dx_tot)

  # Synthetic data generation
  fun_s_array = array(0, dim=dims)

  if (use_verbose)
    cli::cli_alert_info("Starting SDG - iteration...")
  pb = utils::txtProgressBar(0,N,style=3)
  for (i in 1:N) {
    if (!i %in% i_try)
      next
    utils::setTxtProgressBar(pb, i)
    # Weight computation
    if (cluster_search)
      i_other = which(clust_labels == clust_labels[i])
    else
      i_other = 1:N
    out_weight = compute_weights(
      D_i = Dtot[i, i_other],
      K = K,
      alpha_0 = alpha_0,
      beta_0 = beta_0
      )

    #Alignment
    res = srvf_align(
      f = fun_array[ , , i_other[out_weight$index_k]],
      time = time,
      wgt = out_weight$p_k,
      centroid_type = 'mean'
    )

    #SDG as the mean
    if (!add_noise)
      fun_s_array[ , ,i] = res$template
    #SDG as the mean + noise
    else {
      q_s_temp = add_noise(
        qn = res$qn,
        template_q = res$template.q,
        time = time,
        cv = cv,
        is_constrained = is_constrained_2
        )
      f0_temp = res$template_start

      fun_s_array[ , ,i] = srvf_to_f(q_s_temp, time, f0_temp, multidimensional = T)
    }

  }
  close(pb)

  if (do_positive){
    for (l in which(is_constrained != "real")) {
      if (is_constrained[l] == "pos")
        fun_s_array[l,,] = (fun_s_array[l, , ])^2
      if (is_constrained[l] == "strict-pos")
        fun_s_array[l,,] = exp(fun_s_array[l, , ])
    }
  }

  return(fun_s_array)

}
