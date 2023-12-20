#' Add Noise to the Template
#'
#' @description
#' Generate noise with zero mean and covariance function, and use it to perturb
#' the given template function
#'
#' @details
#' This function adds noise to the given template function `template_q`.
#'
#' The noise function is sampled from a gaussian distribution with null mean function
#' and exponential covariance function. The hyperparameter of the covariance function
#' are estimated by minimizing the squared Frobenius distance to the sample covariance
#' function of the data in `qn`.
#'
#' By setting the coefficient of variation `cv`, users can fix the amount of noise to be applied.
#' Specifically, the amplitude of the exponential covariance function would be
#' proportional to the squared norm of the `template_q` with proportionality constant
#' `cv`. If `cv` is not set, the optimal amplitude is used.
#'
#' `add_noise()` handles the case in which dimensions of the function are constrained
#' to be non-negative (or positive) by transforming them though a `sqrt` transformation
#'  (or `log` transformation). Such dimensions must be specified in the parameter
#'  `is_positive`
#'
#' `is_constrained` is used to set constraints for the values of the synthetic
#' functions. For each dimension of the function, it states whether values are
#' (strictly) positive. Accepted values are: `pos` (positive functions),
#' `strict-pos` (non-negative functions), and `real` (non-constrained functions).
#' Constraied dimensions are transformed through sqrt-transformation (if positive
#' functions) or log-transformation (strict-positive functions)
#'
#' @param template_q either array of dimension \eqn{L \times M \times 1} or matrix
#' of dimension \eqn{L \times M} : contains the template function to be perrurbed
#' by adding noise. The function is of dimension \eqn{L} and is evaluated in \eqn{M} points.
#' @param qn array of dimension \eqn{L \times M \times N} : contains the set of \eqn{N}
#' functions for which `template_q` is the template. Functions are of the same
#' dimension \eqn{L} and evaluated at the same points \eqn{M} as `template_q`
#' @param time vector of dimension \eqn{M} : specifies the grid on which the functions
#' are evaluated. If `NULL`, it is set as an uniform grid in \eqn{[0,1]}.
#' @param cv scalar, only non-negative values are accepted : coefficient of variation (see Details).
#' @param is_constrained vector of size \eqn{L} : each element indicates
#' whether the dimension \eqn{l = 1,\dots,L} is constrained and how (see Details)
#'
#' @return
#' - `template_q_new` matrix of dimensions \eqn{L \times M} : original function
#' `template_q` perturbed by noise.
#'
#' @export
add_noise = function(
    qn, template_q, time = NULL, cv = 0, is_constrained = NULL
    ) {

  dims = dim(qn)
  if (length(dims) == 2) {
    dim(qn) <- c(1, dims)
    dims <- dim(qn)
    cli::cli_alert_warning("Warning: qn must be provided as an array of dimensions L * M * N. Now: L=1")
  }
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]
  if (is.null(time))
    time = seq(0, 1, length.out = M)
  if (length(time) != M)
    cli::cli_abort("Error: time must be of length M")
  dims_template = dim(template_q)
  if (sum(dims_template[1:2] != dims[1:2])>0)
    cli::cli_abort("Error: dimensions of qn and template_qn must be coherent")
  if(length(dims_template)==2)
    template_q = array(template_q, dim=c(L,M,1))
  if (length(cv) > 1)
    cli::cli_alert_warning("Warning: cv has length > 1. Only the first element is used")
  cv = rep(cv[1], L)
  if (cv[1] > 1)
    cli::cli_abort("Error: cv must be non-negative")
  use_opt_ampl = ifelse(cv[1] == 0, TRUE, FALSE)
  if (is.null(is_constrained))
    is_constrained = rep("real", L)
  if (length(is_constrained != L))
    cli::cli_abort("Error: is_constrained must be of length L")
  if (sum(is_constrained == "real") + sum(is_constrained == "pos") + sum(is_constrained == "strict-pos") != L)
    cli::cli_abort("Error: is_constrained accept only values among pos, strict-pos, real")
  do_positive = (sum(is_constrained != "real") > 0)


  # Mean function
  MU = matrix(template_q[,,1], nrow=L, ncol=M)

  if (do_positive) {
    for (l in which(is_constrained != "real")){
      if (is_constrained[l] == "pos")
        MU[l,] = sqrt(MU[l,])
      if (is_constrained[l] == "strict-pos")
        MU[l,] = log(MU[l,])
    }
  }

  # Covariance operator
  norm_MU_2 = apply(MU, 1, function(y){integrate(time,y^2)})
  cv = cv * norm_MU_2

  temp_list = purrr::map(1:L, \(l){t(qn[l,,])})
  qn_mfData = roahd::mfData(
    grid=time,
    Data_list=temp_list
  )


  COV_true = roahd::cov_fun(qn_mfData)
  COV = purrr::map(1:L, \(l){
    COV_true_temp = COV_true[[paste0(l,"_",l)]]$values

    my_fun_min = function(par) {
      COV_ab = roahd::exp_cov_function(time, alpha = par[1], beta = par[2])
      sum(diag((COV_ab-COV_true_temp)%*%t(COV_ab-COV_true_temp))) # Squared Frob. distance
    }

    par_opt = stats::constrOptim(theta = c(norm_MU_2[l], 1), f = my_fun_min,
                       grad = NULL, ci = c(0, 0), ui = rbind(c(1,0), c(0,1)))

    if (use_opt_ampl)
      roahd::exp_cov_function(time, alpha = par_opt$par[1], beta = par_opt$par[2])
    else
      roahd::exp_cov_function(time, alpha = cv[l], beta = par_opt$par[2])
    })


  # Correlation between dimensions (only if L>1)
  if (L>1) {
    cor_temp = roahd::cor_spearman(qn_mfData)
    RHO = unlist(purrr::map(1:(L-1), \(l){cor_temp[l,(l+1):L]}))
  }



  # Adding gaussian noise to the mean
  if (L>1)
    q_s_roahd = roahd::generate_gauss_mfdata(
      N = 1,
      L = L,
      centerline = MU,
      correlations = RHO,
      listCov = COV
    )
  else
    q_s_roahd = roahd::generate_gauss_fdata(
      N = 1,
      centerline = c(MU),
      Cov = COV[[1]]
    )

  template_q_new = matrix(nrow=L, ncol=M)
  for (l in 1:L)
    template_q_new[l,] = q_s_roahd[[l]]

  if (do_positive){
    for (l in which(is_constrained != "real")) {
      if (is_constrained[l] == "pos")
        template_q_new[l,] = (template_q_new[l,])^2
      if (is_constrained[l] == "strict-pos")
        template_q_new[l,] = exp(template_q_new[l,])
    }
  }

  template_q_new
}
