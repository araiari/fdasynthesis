#' Elastic distances computation for 2 functions
#'
#' @description
#' Calculates the amplitude distance and the phase distance of two functions.
#'
#' @param f1 An \eqn{M}-dimensional vector (if unidimensional) or an \eqn{L
#'   \times M} matrix (if multidimensional) with the \eqn{M} measurements in the
#'   \eqn{L} dimensions of the function `f1`.
#' @param f2 An \eqn{M}-dimensional vector (if unidimensional) or an \eqn{L
#'   \times M} matrix (if multidimensional) with the \eqn{M} measurements in the
#'   \eqn{L} dimensions of the function `f2`.
#' @param time An \eqn{M}-dimensional vector, sample points of functions.
#' @param lambda A numeric value in \eqn{[0,1]}, controlling the amount of
#'   warping. Defaults to `0`.
#' @param pen A string specifying the type of alignment penalty. Choices are
#'   `"roughness"` (second derivative), `"geodesic"` (geodesic distance from
#'   `id`), `"norm"` (norm from `id`). Defaults to `"roughness"`.
#'
#' @return list of two elements:
#' - `Dx` scalar : phase distance
#' - `Dy` scalar : amplitude distance
compute_elastic_distance = function (f1, f2, time, lambda = 0, pen = "roughness") {
  q1 <- fdasrvf::f_to_srvf(f1, time, multidimensional = TRUE)
  q2 <- fdasrvf::f_to_srvf(f2, time, multidimensional = TRUE)

  gam <- fdasrvf::optimum.reparam(
    Q1 = q1,
    T1 = time,
    Q2 = q2,
    T2 = time,
    lambda = 0,
    pen = "roughness",
    f1o = c(f1[, 1]),
    f2o = c(f2[, 1])
  )

  fw <- qw <- NULL
  for (i in 1:nrow(f2)) {
    fw <- rbind(fw, fdasrvf::warp_f_gamma(f2[i, ], time, gam))
    qw <- rbind(qw, fdasrvf::warp_q_gamma(q2[i, ], time, gam))
  }

  Dy <- 0
  for (i in 1:nrow(f2))
    Dy <- Dy + pracma::trapz(time, (q1[i, ] - qw[i, ]) ^ 2)
  Dy <- sqrt(Dy) ## amplitude distance

  time1 <- seq(0, 1, length.out = length(time))
  binsize <- mean(diff(time1))
  psi <- sqrt(fdasrvf::gradient(gam, binsize))
  q1dotq2 <- pracma::trapz(time1, psi)
  if (q1dotq2 > 1)
    q1dotq2 <- 1
  else if (q1dotq2 < -1)
    q1dotq2 <- -1
  Dx <- acos(q1dotq2) ## phase distance

  list(Dy = Dy, Dx = Dx)
}

#' Elastic distances computation between a set of functions
#'
#' Calculates pairwise the amplitude distance and the phase distance of a set of
#' functions. Iteratively calls the function `compute_elastic_distance()`.
#'
#' @param f_array An array of dimensions \eqn{L \times M \times N} : contains
#'   the values of \eqn{N} functions of dimension \eqn{L} observed over \eqn{M}
#'   points in the domain.
#' @param time An \eqn{M}-dimensional vector, sample points of functions.
#'
#' @return A list of two elements:
#' - `Dx_tot` matrix \eqn{N \times N} : phase distance
#' - `Dy_tot` matrix \eqn{N \times N} : amplitude distance
#' For each matrix, the element in position `[i,j]` represents the distance
#' between the `i`-th and the `j`-th function in `f_array`.
compute_elastic_distance_one_set <- function(f_array, time) {
  call <- rlang::call_match()

  N <- dim(f_array)[3]

  if (is.null(labels))
    labels <- 1:N

  index_table <- linear_index(N)

  .pairwise_distances <- function(index_table) {
    pb <- progressr::progressor(steps = nrow(index_table))
    furrr::future_map2(index_table$i, index_table$j, \(i, j) {
      pb()
      compute_elastic_distance(
        f1 = f_array[, , i],
        f2 = f_array[, , j],
        time = time,
        lambda = 0,
        pen = "roughness"
      )
    }, .options = furrr::furrr_options(seed = TRUE))
  }

  Dlist <- .pairwise_distances(index_table)

  Dx_tot <- purrr::map_dbl(Dlist, "Dx")
  attributes(Dx_tot) <- NULL
  attr(Dx_tot, "Labels") <- labels
  attr(Dx_tot, "Size") <- N
  attr(Dx_tot, "Diag") <- FALSE
  attr(Dx_tot, "Upper") <- FALSE
  attr(Dx_tot, "call") <- call
  class(Dx_tot) <- "dist"

  Dy_tot <- purrr::map_dbl(Dlist, "Dy")
  attributes(Dy_tot) <- NULL
  attr(Dy_tot, "Labels") <- labels
  attr(Dy_tot, "Size") <- N
  attr(Dy_tot, "Diag") <- FALSE
  attr(Dy_tot, "Upper") <- FALSE
  attr(Dy_tot, "call") <- call
  class(Dy_tot) <- "dist"

  list(Dx_tot = Dx_tot, Dy_tot = Dy_tot)
}

#' Elastic distances computation between two set of functions
#'
#' Calculates pairwise the amplitude distance and the phase distance of the
#' items in a set of functions with the items in another set of functions.
#' Iteratively calls the function `compute_elastic_distance()`.
#'
#' @param f_array1 An array of dimensions \eqn{L \times M \times N1} : contains
#'   the values of \eqn{N1} functions of dimension \eqn{L} observed over \eqn{M}
#'   points.
#' @param f_array2 An array of dimensions \eqn{L \times M \times N2} : contains
#'   the values of \eqn{N2} functions of dimension \eqn{L} observed over \eqn{M}
#'   points in the domain.
#' @param time An \eqn{M}-dimensional vector, sample points of functions.
#'
#' @return A list of two elements:
#' - `Dx_tot` matrix \eqn{N1 \times N2} : phase distance
#' - `Dy_tot` matrix \eqn{N1 \times N2} : amplitude distance
#' For each matrix, the element in position `[i,j]` represents the distance
#' between the `i`-th function in `f_array1` and the `j`-th function in
#' `f_array_2`.
compute_elastic_distance_two_sets = function(f_array1, f_array2, time) {
  L <-  dim(f_array1)[1]
  M <- dim(f_array1)[2]
  N1 <- dim(f_array1)[3]
  N2 <- dim(f_array2)[3]

  if (dim(f_array2)[1] != L)
    cli::cli_abort("The functions in the two sets have different codomain dimensions.")
  if (dim(f_array2)[2] != M)
    cli::cli_abort("The functions in the two sets have different grid sizes.")

  Dx_tot <- Dy_tot <- matrix(nrow = N1, ncol = N2)
  .pairwise_distances <- function(f1, f2, time) {
    N1 <- dim(f1)[3]
    N2 <- dim(f2)[3]
    pb <- progressr::progressor(steps = N1)
    furrr::future_walk(1:N1, \(n1) {
      pb()
      out <- purrr::map(1:N2, \(n2) {
        compute_elastic_distance(
          f1 = f1[, , n1],
          f2 = f2[, , n2],
          time = time,
          lambda = 0,
          pen = "roughness"
        )
      })
      Dx_tot[n1, ] <- purrr::map_dbl(out, "Dx")
      Dy_tot[n1, ] <- purrr::map_dbl(out, "Dy")
    }, .options = furrr::furrr_options(seed = TRUE))
  }

  list(Dx_tot = Dx_tot, Dy_tot = Dy_tot)
}
