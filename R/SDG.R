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
#' `strict-mon` (for strictly-monotone increasing functions).
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
#' @param fun_array array of sizes \eqn{L \times M \times N} : contains the
#'   values of the \eqn{N} original functions of dimension \eqn{L} evaluated on
#'   \eqn{M} points.
#' @param i_synth integer vector of length \eqn{N_{synth} \leq N}: indices of
#'   the subset of the curves to be synthetized. If `NULL` all the curves are
#'   synthetized.
#' @param time vector of size \eqn{M} : specifies the grid on which the curves are
#'   evaluated. If `NULL`, it is set as an uniform grid in \eqn{[0,1]}.
#' @param D_list list of two elements `Dx_tot` and `Dy_tot` for (respectively)
#'   phase and amplitude distance matrices, provided as `matrix` objects.
#' - If `w_a` is not specified, list of two elements, each of which is either a
#' distance matrix of dimension \eqn{N \times N} : elastic distances
#' between all the functions.
#' - If `w_a` is specified, list of two elements, each of which is either a
#' distance matrix of dimension \eqn{N_{synth} \times N} : elastic
#' distances between functions to be synthetized and all the functions.
#' - If `NULL` it is computed (see Details).
#' @param w_a scalar with values in \eqn{[0,1]} : weights of the
#'   amplitude distance in the linear combination of the elastic distances.
#'   If `NULL` it is estimated (see Details).
#' @param K scalar integer : number of neighbors to be considered.
#' @param alpha_0 scalar, with \eqn{\alpha_0 >0} : concentation parameter of the
#'   Dirichlet distribution from which the weights of neighbors are sampled.
#' @param beta_0 scalar, with \eqn{\beta_0 > 0} : rate parameter
#'   weighting the inverse distances in the weight computation.
#' @param add_noise boolean : whether noise should be added to the weighted mean.
#'  Defaults to `FALSE`.
#' @param is_constrained vector of size \eqn{L} : each element indicates any constraint of
#'   the dimension \eqn{l = 1,\dots,L}. It is considered only if
#'   `add_noise = TRUE` (see Details).
#' @param cv scalar with \eqn{cv > 0} : coefficient of variation of the
#'   covariance amplitude with respect to the mean squared norm. It is the
#'   input of the function add_noise.
#' @param clust_labels vector of size \eqn{N} : specifies the labels
#' of the functions. If provided, the \eqn{K} neighbors are searched among the
#' subset of functions with the same label. Defaults to `NULL`.
#' @param use_verbose boolean : specifying whether to display information about
#'   the calculations in the console. Defaults to `FALSE`.
#'
#' @return A list of the following elements:
#' - `fun_s_array` array of sizes \eqn{L \times M \times N_{synth}} : contains the values
#' of the \eqn{N_{synth}} synthetic functions of dimension \eqn{L} evaluated on \eqn{M}
#' points.
#'
#' @export
SDG = function (
    fun_array,
    i_synth = NULL,
    time = NULL,
    D_list = NULL,
    w_a = NULL,
    K = 5,
    alpha_0 = 1,
    beta_0 = 1,
    add_noise = FALSE,
    is_constrained = NULL,
    cv = 0,
    clust_labels = NULL,
    use_verbose = FALSE
) {


  dims = dim(fun_array)

  L = dims[1]
  M = dims[2]
  N = dims[3]

  sub_synth = (!is.null(i_synth))
  if (!sub_synth)
    i_synth = 1:N
  N_synth = length(i_synth)

  #Check incompatibilities
  if (use_verbose)
    cli::cli_alert_info("Check possible incompatibilities...")
  if (is.null(time))
    time = seq(0,1,length.out=M)
  if (length(time) != M)
    cli::cli_abort("Error: `time` must be of length M")

  if (add_noise && !is.null(is_constrained)) {
    if (length(is_constrained) != L)
      cli::cli_abort("Error: is_constrained is of length {length(is_constrained)} but must be of length L")
  if (sum(is_constrained %in% c("real", "pos","strict-neg","mon","strict-mon"))!= L)
    cli::cli_abort("Error: is_constrained accepts only values between 'pos', 'strict-pos', 'mon', 'strict-mon', and 'real'")
  }
  do_positive = (sum(is_constrained %in% c("pos","strict-pos")) > 0)
  is_constrained_2 = is_constrained
  is_constrained_2[which(is_constrained_2 %in% c("pos", "strict-pos"))] = "real"
  is_constrained_2[which(is_constrained_2 %in% c("mon"))] = "pos"
  is_constrained_2[which(is_constrained_2 %in% c("strict-mon"))] = "strict-pos"

  if (is.null(w_a) && !requireNamespace("SPARTAAS", quietly = TRUE))
    cli::cli_abort("Error: parameter {.arg w_a} not provided and {.pkg SPARTAAS} not installed.")

  estimate_w_a = (is.null(w_a) && requireNamespace("SPARTAAS", quietly = TRUE))
  if (estimate_w_a &&
      !is.null(D_list) &&
      (nrow(D_list$Dx_tot)!= N || ncol(D_list$Dx_tot)!= N  ||
       nrow(D_list$Dy_tot)!= N || ncol(D_list$Dy_tot)!= N ))
    cli::cli_abort("Error: distance matrices must be of dimensions NxN = {N}x{N}")
  if (!estimate_w_a &&
      !is.null(D_list) &&
      (nrow(D_list$Dx_tot)!= N_synth || ncol(D_list$Dx_tot)!= N  ||
       nrow(D_list$Dy_tot)!= N_synth || ncol(D_list$Dy_tot)!= N ))
    cli::cli_abort("Error: distance matrices must be of dimensions N_synthxN = {N_synth}x{N}")
  if (is.null(clust_labels))
    cluster_search = FALSE
  else {
    if (length(clust_labels) != N)
      cli::cli_abort("Error: cluster labels must be of length N={N}")
    cluster_search = TRUE
  }
  if (use_verbose)
    cli::cli_alert_info("End: no incompatibilities found!")


  if (do_positive) {
    for (l in which(is_constrained %in% c("pos","strict-pos"))){
      if (is_constrained[l] == "pos")
        fun_array[l, , ] = sqrt(fun_array[l, , ])
      if (is_constrained[l] == "strict-pos")
        fun_array[l, , ] = log(fun_array[l, , ])
    }
  }

  # Computation of the elastic distances
  if (is.null(D_list)) {
    if (use_verbose)
      cli::cli_alert_info("Starting the computation of elastic distances...")
    if (!sub_synth || estimate_w_a)
      D_list = compute_elastic_distance_one_set(fun_array, time)
    else
      D_list = compute_elastic_distance_two_sets(fun_array[,,i_synth], fun_array, time)
    if (use_verbose)
      cli::cli_alert_info("End of the computation of elastic distances!")
  }


  # Computation of the total distance as convex combination of elastic distances
  if (estimate_w_a){
    if (use_verbose)
      cli::cli_alert_info("Starting the computation of the optimal {.arg w_a}...")
    w_a = estimate_optimal_w_a (D1 = D_list$Dy_tot, D2 = D_list$Dx_tot)
    if (use_verbose)
      cli::cli_alert_info("End of the computation of the optimal {.arg w_a}!")
    if (sub_synth) {
      D_list$Dy_tot = D_list$Dy_tot[i_synth,]
      D_list$Dx_tot = D_list$Dx_tot[i_synth,]
    }
  }
  Dtot = w_a * D_list$Dy_tot + (1-w_a) * D_list$Dx_tot



  # Synthetic data generation
  fun_s_array = array(0, dim=c(M, L, N_synth))

  if (use_verbose)
    cli::cli_alert_info("Starting SDG - iteration...")
  pb = utils::txtProgressBar(0, N_synth, style=3)
  for (i in 1:N_synth) {
    utils::setTxtProgressBar(pb, i)
    # Weight computation
    if (cluster_search)
      i_other = which(clust_labels == clust_labels[i_synth[i]])
    else
      i_other = 1:N
    out_weight = compute_weights(
      D_i = Dtot[i, i_other],
      K = K,
      alpha_0 = alpha_0,
      beta_0 = beta_0
      )

    #starting point
    f0_temp = apply(fun_array[, 1, i_other[out_weight$index_k]], 1, stats::weighted.mean, out_weight$p_k)

    #karcher #################################################### Da rifare
    res = multivariate_weighted_karcher_mean(
      beta = fun_array[ , , i_other[out_weight$index_k]],
      wts = out_weight$p_k
    )

    cli::cli_alert_info("Mean computed")

    #SDG as the mean
    if (!add_noise)
      fun_s_array[ , ,i] = res$betamean + f0_temp
    #SDG as the mean + noise
    else {
      q_s_temp = add_noise(
        qn = res$qn,
        template_q = res$mu,
        cv = cv,
        is_constrained = is_constrained_2,
        use_verbose = use_verbose
        )

      fun_s_array[ , , i] = fdasrvf::q_to_curve(q_s_temp, scale = 1) + f0_temp
    }

  }
  close(pb)

  if (do_positive){
    for (l in which(is_constrained != "real")) {
      if (is_constrained[l] == "pos")
        fun_s_array[l, , ] = (fun_s_array[l, , ])^2
      if (is_constrained[l] == "strict-pos")
        fun_s_array[l, , ] = exp(fun_s_array[l, , ])
    }
  }

  return(fun_s_array)

}
