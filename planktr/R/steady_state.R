#' Determine steady state cell size distribution
#'
#' Using the analytic formulae in eqs.(4.16) to (4.21) of the paper
#' @param r Object of class PlanktonParams
#' @return List containing:
#    \enumerate{
#'     \item vector of steady-state population densities
#'     \item nutrient concentration at steady state
#'   }
steady_state <- function(r) {

    # Without loss of generality we can set w_* = 1 so that \tilde{w} = w/w_* = w
    # The population density for any other w_* can be obtained by scaling
    # see eqs.(6.1) and (6.2)

    wp <- min(r@wBar[r@wBar>=(0.5 + r@delta_q/2)])  # the point at which we glue
    kv <- r@k(r@wBar)
    qv <- r@q(r@wBar)

    # Calculate growth rate from predation
    # see eq.(5.9).
    # We assume that the coefficient of the community size spectrum is 1
    # Note change in exponent because we integrate with respect to log var
    gp <- r@wBar^(-r@xi)*s_moment(r, -r@xi)
    # Use constant mortality plus mortality from predation
    # see eq.(5.10)
    mv <- r@m + r@epsilon * r@wBar^(1-r@xi) * s_moment(r, r@gamma-2)

    p <- function(Nu0) {
        # Calculate the steady-state solution
        #
        # Args:
        #   Nu0: nutrient concentration
        #
        # Returns:
        #   List containing the solution and the integral over h/e
        #   For the correct value of Nu0 that integral will be equal to 1

        w <- r@wBar
        x <- log(r@wBar)

        # Add growth rates from nutrient for given Nu0
        gv <- gp + r@g(w, Nu0)

        # Calculate e(w) in eq.(4.16)
        ep <- (kv+mv)/gv
        # First calculate for w < wp
        es <- rev(intx(rev(ep[w<=wp]), rev(w[w<=wp])))
        # then for w >= wp
        el <- -intx(ep[w>=wp], w[w>=wp])
        # and put the results together
        e <- exp(c(es[-length(es)], el))

        # Calculate h(w) in eq.(4.17)
        hp <- kv*e/gv
        hp[w < r@w_th] <- 0
        h <- rep_len(0, length(w))
        for (i in 1:length(w)) {
            h[i] <- 2*intx(hp*r@q(w[i]/w)/w, w)[length(w)]
        }

        # Calculate Theta in eq.(4.21). Because h/e is zero beyond w=wp, we simply
        # integrate up to w=1. Then Theta is automatically 1 for w>wp when the
        # boundary condition (4.18) is satisfied.
        he <- h/e
        he[e==0] <- 0  # Removes the infinity from division by zero
        Theta <- intx(he, w)
        # The boundary condition (4.18) is satisfied iff b=1
        b <- Theta[length(w)]

        # Use eq.(4.20) to calculate the solution.
        Psi <- gv[w==wp]*e*Theta/gv

        return(list(Psi, b))
    }

    # Find the mortality rate constant that satisfies the boundary condition
    # in eq. (4.18)
    Nu <- uniroot(function(Nu0) p(Nu0)[[2]]-1, lower=1, upper=100)[["root"]]

    # Calculate and return the solution
    psi <- p(Nu)[[1]]
    list(psi, Nu)
}

intx <- function(f, x) {
    # Use trapezoidal rule for integration
    #
    # Args:
    #   f: value of function at points given by x
    #   x: the points at which the function values are given. Needs to be in
    #      either increasing or decreasing order.
    #
    # Returns:
    #   vector of values of integral over f at points given by x
    #
    # The integration starts at the first value in x, so first
    # value in the returned vector is zero.

    idx <- 1:(length(f)-1)
    dx <- abs(diff(x))
    c(0,cumsum((f[idx]+f[idx+1])*dx[idx]/2))
}
