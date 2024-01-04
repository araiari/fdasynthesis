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

  L <- dim(f_array)[1]
  N <- dim(f_array)[3]

  index_table <- fdasynthesis:::linear_index(N)

  # Multidimensional case:
  if (L > 1) {
    .pairwise_distances <- function(index_table) {
      pb <- progressr::progressor(steps = nrow(index_table))
      furrr::future_map2(index_table$i, index_table$j, \(i, j) {
        pb()
        fdasrvf::calc_shape_dist(
          beta1 = f_array[, , i],
          beta2 = f_array[, , j],
          mode = "O",
          rotation = FALSE,
          scale = TRUE
        )
      }, .options = furrr::furrr_options(seed = TRUE))
    }

    Dlist <- .pairwise_distances(index_table)

    Dx_tot <- purrr::map_dbl(Dlist, "dx")
    Dy_tot <- purrr::map_dbl(Dlist, "d")
  }
  #Unidimensional case:
  else {
    .pairwise_distances <- function(index_table) {
      pb <- progressr::progressor(steps = nrow(index_table))
      furrr::future_map2(index_table$i, index_table$j, \(i, j) {
        pb()
        fdasrvf::elastic.distance(
          f1 = f_array[1, , i],
          f2 = f_array[1, , j],
          time = time
        )
      }, .options = furrr::furrr_options(seed = TRUE))
    }

    Dlist <- .pairwise_distances(index_table)

    Dx_tot <- purrr::map_dbl(Dlist, "Dx")
    Dy_tot <- purrr::map_dbl(Dlist, "Dy")
  }

  attributes(Dx_tot) <- NULL
  attr(Dx_tot, "Labels") <- 1:N
  attr(Dx_tot, "Size") <- N
  attr(Dx_tot, "Diag") <- FALSE
  attr(Dx_tot, "Upper") <- FALSE
  attr(Dx_tot, "call") <- call
  class(Dx_tot) <- "dist"

  attributes(Dy_tot) <- NULL
  attr(Dy_tot, "Labels") <- 1:N
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
#'   points. For unidimensional functions, provide an array of dimensions
#'   \eqn{1 \times M \times N1}.
#' @param f_array2 An array of dimensions \eqn{L \times M \times N2} : contains
#'   the values of \eqn{N2} functions of dimension \eqn{L} observed over \eqn{M}
#'   points in the domain. For unidimensional functions, provide an array of
#'   dimensions \eqn{1 \times M \times N2}.
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

  # Multidimensional case:
  if (L > 1) {
    .pairwise_distances <- function(f1, f2, time) {
      N1 <- dim(f1)[3]
      N2 <- dim(f2)[3]

      pb <- progressr::progressor(steps = N1)
      furrr::future_walk(1:N1, \(n1) {
        pb()
        out <- purrr::map(1:N2, \(n2) {
          fdasrvf::calc_shape_dist(
            beta1 = f1[, , n1],
            beta2 = f2[, , n2],
            mode = "O",
            rotation = FALSE,
            scale = TRUE
          )
        })

        Dx_tot[n1, ] <- purrr::map_dbl(out, "dx")
        Dy_tot[n1, ] <- purrr::map_dbl(out, "d")
      }, .options = furrr::furrr_options(seed = TRUE))
    }

  }

  # Unidimensional case:
  else {
    .pairwise_distances <- function(f_array1, f_array2, time) {
      N1 <- dim(f_array1)[3]
      N2 <- dim(f_array2)[3]
      pb <- progressr::progressor(steps = N1)
      furrr::future_walk(1:N1, \(n1) {
        pb()
        out <- purrr::map(1:N2, \(n2) {
          fdasrvf::elastic.distance(
            f1 = f_array1[1, , n1],
            f2 = f_array2[1, , n2],
            time = time
          )
        })
        Dx_tot[n1, ] <- purrr::map_dbl(out, "Dx")
        Dy_tot[n1, ] <- purrr::map_dbl(out, "Dy")
      }, .options = furrr::furrr_options(seed = TRUE))
    }
  }

  Dlist = .pairwise_distances(f_array1, f_array2, time)

  list(Dx_tot = Dx_tot, Dy_tot = Dy_tot)

} ############## NOT WORKING
