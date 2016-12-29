#' Evolve cell population density using the population balance equation
#'
#' see eq.(2.10)
#' @param p0 N x Ns matrix or vector with initial population for all Ns species
#' or a vector of length N with initial population for one species that is then
#' replicated over all species. If missing, the steady-state solution is used
#' and replicated over all species.
#' @param Nu0 Initial resource concentration. If missing, the steady-state value
#' is used
#' @param r Sim object
#' @return list of two elements:
#'     Nt x N x Ns  array of population densities
#'     vector of length Nt containing the nutrient densities
evolve_cell_pop <- function(p0, Nu0, r) {
    N <- r@N
    Ns <- r@Ns
    Na <- r@Na

    # Predation
    mkernel <- r@dx/Na*fft(r@s(-r@xa)*exp(r@xa*(-r@xi)))
    gkernel <- r@dx/Na*fft(r@epsilon*r@s(r@xa-r@xa[1])*exp(r@xa*(r@gamma-2)))

    ks <- r@k(r@w)
    # fft of offspring size distribution
    FqR <- 2*r@dx/N*fft(r@q(r@w))
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/r@L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)

    # Create vectors containing powers
    wsmxi <- r@ws^(-r@xi)
    womxi <- r@w^(1-r@xi)
    wmxi <- r@w^(-r@xi)

    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=r@Ns)
        Nu <- pN[length(pN)]
        pcp <- fft(community(p, r))

        # Calculate growth rate
        gp <- Re(fft(gkernel*pcp, inverse = TRUE))  # from predation
        gr <- r@g(r@w, Nu) # from resource

        # Calculate death rate
        mp <- Re(fft(mkernel*pcp, inverse = TRUE)) # from predation

        # Calculate right-hand side of population balance equation
        f <- matrix(0, nrow = N, ncol = Ns)
        idx <- 1:N
        for (i in 1:Ns) {
            gs <- gr + womxi * gp[idx]  # growth rate
            ms <- r@m + wmxi * mp[idx]  # mortality rate
            f[,i] <- wsmxi[i] * (
                -(ks+ms)*p[,i] +  # linear part
                # birth part
                Re(fft(FqR*(fft(ks*p[,i])), inverse = TRUE)) +
                # growth part
                -Re(fft(fft(gs*p[,i])*k1, inverse=TRUE)/N)/r@w
            )
        }
        nutrientGrowth <- r@dNu(r@w, Nu, p, r)

        list(c(f, nutrientGrowth))
    }

    out <- ode(y=c(p0, Nu0), times=r@t, func=f)

    Nut <- out[ , ncol(out)]
    psit <- array(dim=c(length(r@t), N, Ns))
    psit[ , , ] <- out[ , 2:(ncol(out)-1)]
    # The following is obsolete code that was necessary when I wanted to
    # include the right boundary in the return value
    # psit[ , N+1, ] <- psit[ , 1, ]
    list(psit, Nut)
}

#' Determine community spectrum
#'
#' Adds together the population density of all species
#' wrapped around assuming periodicity in species size.
#' Used in \code{\link{evolve_cel_pop}}.
#' This returns the $\tilde{p}_c$ from the vignette that
#' is constant in the steady-state.
#'
#' @param p matrix of population densities (N x Ns)
#' @param sim Sim object
#' @return A vector of community population densities at points \code{xa}
community <- function(p, sim) {
    pc <- vector("numeric", length=sim@Nal)
    idx <- 1:sim@N
    for (i in 1:sim@Ns) {
        pc[idx] <- pc[idx] + sim@wsmgamma[i]*p[,i]
        idx <- idx + sim@ds
    }
    # Pull out a factor of w^{-\gamma} so that pc is constant in steady state
    pc <- sim@walgamma*pc
    # Wrap around
    pcp <- pc[(sim@Nal-sim@Na+1L):sim@Nal]
    if (sim@Nal > sim@Na) {
    # Wrap around by adding the lowest Nal-Na entries to the top
        top <- (2L*sim@Na-sim@Nal+1L):sim@Na  # the top Nal-Na indices
        pcp[top] <- pcp[top] + pc[1:(sim@Nal-sim@Na)]
    }
    pcp
}

#' Get community spectrum
#'
#' This produces $\tilde{p}_c(t, w)$ as defined in the vignette.
#' In the steady state this should be constant in w.
#' @param sim Sim object
#' @return matrix Nt x Na
get_community <- function(sim) {
    aaply(sim@p, 1, "community", sim=sim)
}

#' Plot community spectrum against time and size
#'
#' This plots $\tilde{p}_c(t, w)$ as defined in the vignette.
#' @param sim Sim object
plot3d_community <- function(sim) {
    com <- get_community(sim)
    persp3d(sim@t, sim@xa, com, col = "lightblue",
            xlab="t", ylab="xa", zlab="p_c",
            main="Community size spectrum")
}

#' Plot community spectrum against size at one time
#'
#' This plots $\tilde{p}_c(t, w)$ as defined in the vignette.
#' @param sim Sim object
#' @param t Time at which to plot. If the value is not available
#' at that time, the last earlier time is used. Default: latest available time.
plot_community <- function(sim, t=NULL) {
    if (is.null(t)) {
        ti <- sim@Nt+1
    } else if (t >= 0) {
        ti <- which(sim@t >= t)[1] - 1
    } else {
        stop("The time can not be negative")
    }
    com <- community(sim@p[ti, , ], sim=sim)
    plot(sim@xa, com, type="l", xlab="xa", ylab="p_c",
         main=paste("Community spectrum at t=", sim@t[ti]))
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
