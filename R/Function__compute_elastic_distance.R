
#' Elastic distances computation for 2 functions
#'
#' @description
#' Calculates the amplitude distance and the phase distance of two functions
#'
#' @param f1 \eqn{M}-dimensional vector (if unidimensional) or  \eqn{L \times M} matrix (if multidimensional)
#'        with the \eqn{M} measurements in the \eqn{L} dimensions of the function `f1`
#' @param f2 \eqn{M}-dimensional vector (if unidimensional) or a \eqn{L \times M} matrix (if multidimensional)
#'        with the \eqn{M} measurements in the \eqn{L} dimensions of the function `f2`
#' @param time \eqn{M}-dimensional vector, sample points of functions
#' @param lambda numeric in \eqn{[0,1]}, controlling the amount of warping
#' @param pen alignment penalty, between "roughness" (second derivative), "geodesic"
#'        (geodesic distance from id), "norm" (norm from id)
#'
#' @return list of two elements:
#' - `Dx` scalar : phase distance
#' - `Dy` scalar : amplitude distance
compute_elastic_distance = function (
    f1, f2, time, lambda = 0, pen = "roughness"
    )  {

  q1 <- fdasrvf::f_to_srvf(f1, time, multidimensional = T) #aggiunto multidimensional=T
  q2 <- fdasrvf::f_to_srvf(f2, time, multidimensional = T) #aggiunto multidimensional=T

  gam <- fdasrvf::optimum.reparam(q1, time, q2, time, lambda=0, pen='roughness',
                         # f1o = rep(0, nrow(q1)), f2o = rep(0, nrow(q2))) #modificato -> errore inizializzazione delle fo
                         f1o = c(f1[,1]), f2o = c(f2[,1])) #modificato -> errore inizializzazione delle fo in 0

  fw <- qw <- NULL #modificato -> questo non è multidimensional e con un ciclo for lo rendo tale
  for (i in 1:nrow(f2)) {
    fw = rbind(fw, fdasrvf::warp_f_gamma(f2[i,], time, gam))
    qw <- rbind(qw, fdasrvf::warp_q_gamma(q2[i,], time, gam))
  }

  Dy <- 0
  for (i in 1:nrow(f2)) #modificato -> questo non è multidimensional e con un ciclo for lo rendo tale
    Dy = Dy + pracma::trapz(time, (q1[i,] - qw[i,])^2)
  Dy = sqrt(Dy) ## amplitude distance

  time1 <- seq(0, 1, length.out = length(time))
  binsize <- mean(diff(time1))
  psi <- sqrt(fdasrvf::gradient(gam, binsize))
  q1dotq2 = pracma::trapz(time1, psi)
  if (q1dotq2 > 1) {
    q1dotq2 = 1
  }
  else if (q1dotq2 < -1) {
    q1dotq2 = -1
  }
  Dx <- acos(q1dotq2) ## phase distance

  return(list(Dy = Dy, Dx = Dx))
}


#' Elastic distances computation between a set of functions
#'
#' Calculates pairwise the amplitude distance and the phase distance of a set of
#' functions. Iteratively calls the function elastic.distance.ari.
#'
#' @param f_array an array of dimensions \eqn{L \times M \times N} : contains the
#' values of \eqn{N} functions of dimension \eqn{L} observed over \eqn{M} points
#' in the domain
#' @param time \eqn{M}-dimensional vector, sample points of functions
#'
#' @return list of two elements:
#' - `Dx_tot` matrix \eqn{N \times N} : phase distance
#' - `Dy_tot` matrix \eqn{N \times N} : amplitude distance
compute_elastic_distance_N = function(
    f_array, time
    ) {

  L = dim(f_array)[1]
  M = dim(f_array)[2]
  N = dim(f_array)[3]

  Dx_tot = Dy_tot = matrix(data=0, nrow=N, ncol=N)
  pb = utils::txtProgressBar(min=0, max=N*(N-1)/2, style=3); pb_iter=0
  for (i_f1 in 1:(N-1)){
    for (i_f2 in (i_f1+1):N){
      pb_iter=pb_iter+1
      utils::setTxtProgressBar(pb, pb_iter)
      out = compute_elastic_distance(f_array[,,i_f1], f_array[,,i_f2], time, lambda=0, pen='roughness')
      Dx_tot[i_f1, i_f2] = out$Dx
      Dy_tot[i_f1, i_f2] = out$Dy
    }
  }
  close(pb)

  Dx_tot = Dx_tot + t(Dx_tot)
  Dy_tot = Dy_tot + t(Dy_tot)

  return(list(Dx_tot = Dx_tot, Dy_tot = Dy_tot))
}

