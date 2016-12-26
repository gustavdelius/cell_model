#' Evolve cell population density using the population balance equation
#'
#' see eq.(2.10)
#' @param r PlanktonSim object
#' @return list of two elements:
#'     Nt x N x Ns  array of population densities
#'     vector of length Nt containing the nutrient densities
evolve_cell_pop <- function(r) {
    N <- r@N
    Ns <- r@Ns

    # Predation
    mkernel <- fft(r@S(-r@xa)*exp(r@xa*(-r@xi)))
    gkernel <- fft(r@epsilon*r@S(r@xa)*exp(r@xa*(r@gamma-2)))

    ks <- r@k(r@x)
    # fft of offspring size distribution
    FqR <- fft(r@q(r@x))
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/r@L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)

    # Create vectors containing powers
    wsmxi <- r@ws^(-r@xi)
    womxi <- r@w^(1-r@xi)
    wmxi <- r@w^(-r@xi)

    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=r@Ns)
        Nu <- pN[length(pN)]
        pcp <- community_spectrum(p, r)

        # Calculate growth rate
        gp <- Re(fft(gkernel*pcp, inverse = TRUE))/r@Na  # from predation
        gr <- r@g(r@w, Nu) # from resource

        # Calculate death rate
        mp <- Re(fft(mkernel*pcp, inverse = TRUE))/r@Na # from predation

        # Calculate right-hand side of population balance equation
        f <- matrix(0, nrow = N, ncol = Ns)
        idx <- 1:N
        for (i in 1:Ns) {
            gs <- gr + womxi * gp[idx]  # growth rate
            ms <- r@m + wmxi * mp[idx]  # mortality rate
            f[,i] <- wsmxi[i] * (
                -(ks+ms)*p[,i] +  # linear part
                # birth part
                2*r@L/N*Re(fft(FqR*(fft(ks*p[,i])), inverse = TRUE)/N) +
                # growth part
                -Re(fft(fft(gs*p[,i])*k1, inverse=TRUE)/N)/r@w
            )
        }
        nutrientGrowth <- r@dNu(r@w, Nu, p, r@dxs)

        list(c(f, nutrientGrowth))
    }

    out <- ode(y=c(r@p0, r@Nu0), times=r@t, func=f)

    Nut <- out[ , ncol(out)]
    # The following is obsolete code that was necessary when I wanted to
    # include the right boundary in the return value
    # psit <- array(dim=c(length(t), N+1, Ns))
    # psit[ , -(N+1), ] <- out[ , 2:(ncol(out)-1)]
    # psit[ , N+1, ] <- psit[ , 1, ]
    psit <- out[ , 2:(N+1)]
    list(psit, Nut)
}

#' Determine community spectrum
#'
#' Adds together the population density of all species
#' wrapped around assuming periodicity in species size
#' Used in \code{\link{evolve_cel_pop}}
#'
#' @param p matrix of population densities (N x Ns)
#' @param r PlanktonParams object
#' @return A vector of community population densities at points \code{xa}
community_spectrum <- function(p, r) {
    pc <- vector("numeric", length=r@Nal)
    idx <- 1:r@N
    for (i in 1:r@Ns) {
        pc[idx] <- pc[idx] + r@wsmgamma[i]*p[,i]
        idx <- idx + r@ds
    }
    # Pull out a factor of w^{-\gamma} so that pc is constant in steady state
    pc <- r@walgamma*pc
    # Wrap around by moving the lowest Nal-Na entries to the top
    pcp <- pc[(r@Nal-r@Na+1L):r@Nal]
    top <- 2L*r@Na-r@Nal+1L:r@Na  # the top Nal-Na indices
    pcp[top] <- pcp[top] + pc[1:(r@Nal-r@Na)]
    pcp
}

#' Perform a Fourier interpolation
#'
#' Takes a discretisation at any number of equally-spaced points, performs
#' a Fast Fourier Transform on it and then does a slow inverse Fourier
#' transform evaluated at the desired number of equally-spaced points.
#'
#' @param p vector of values of function at equally spaced steps, excluding
#'          the right endpoint
#' @param n desired length of output vector
#' @return Vector of n values of Fourier interpolation at n equally
#'     spaced points, excluding the right endpoints.
fourier_interpolate <- function(p, n) {
    N <- length(p)
    x <- seq(0, 1-1/n, length.out=n)
    fp <- fft(p)
    f <- rep(Re(fp[1]), length(x))
    for (j in 2:(N/2+1)) {
        f <- f + 2*(Re(fp[j])*cos(2*pi*(j-1)*x) - Im(fp[j])*sin(2*pi*(j-1)*x))
    }
    f/N
}
