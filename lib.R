library("deSolve")

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
    return(c(0,cumsum((f[idx]+f[idx+1])*dx[idx]/2)))
}

steady_state <- function(t, w, g, k, wa, q, delta) {
    # Determine steady state cell size distribution from analytic formulae
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   w: vector of equally spaced w-values (cell size)
    #   g: vector of growth rates
    #   k: vector of division rates
    #   wa: threshold for duplication (size of smallest cell that can divide)
    #   q: function giving offspring size distribution
    #
    # Value:
    #   List containing
    #     matrix of population densities (columns t, rows w)
    #     mortality rate
    
    Nw <- length(w)  # Number of w steps.
    wmin <- w[1]  # Smallest possible cell size
    wp <- min(w[w>=(0.5 + delta/2)])  # the point at which we glue
    
    p <- function(m0) {
        # Calculate the steady-state solution
        #
        # Args:
        #   m0: scalar giving the constant mortality rate
        #
        # Returns:
        #   List containing the solution and the integral over h/e 
        #   For the correct value of m0 that integral will be equal to 1
        
        # Use constant mortality 
        m <- rep_len(m0, length(w))
        
        # Calculate e(w)
        ep <- (k+m)/g
        # First calculate for w < wp
        es <- rev(intx(rev(ep[w<=wp]), rev(w[w<=wp])))
        # then for w >= wp
        el <- -intx(ep[w>=wp], w[w>=wp])
        # and put the results together
        e <- exp(c(es[-length(es)], el))
        
        # Calculate h(w)
        # TODO: Want to convert this to using spectral methods when the steps
        # are logarithmic
        hp <- k*e/g
        hp[w<wa] <- 0
        h <- rep_len(0, length(w))
        for (i in 1:length(w)) {
            h[i] <- 2*intx(hp*q(w[i]/w)/w, w)[length(w)-1]
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
        Psi <- g[w==wp]*e*Theta/g
        
        return(list(Psi, b))
    }
    
    # Find the mortality rate constant that satisfies the boundary condition
    # in eq. (4.18)
    m <- uniroot(function(m0) p(m0)[[2]]-1, lower=0.05, upper=10)[["root"]]
    
    # Calculate the solution
    psi <- p(m)[[1]]
    return(list(psi, m))
}

evolve_cell_pop <- function(t, x, p0, g, k, q, m) {
    # Evolve cell population density
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   x: vector of equally spaced steps in logarithmic cell size
    #      This should include both endpoints of the periodic interval.
    #      So length(x) is one larger than the number of steps
    #   p0 : vector on initial population density
    #   g: vector of growth rates
    #   k: vector of division rates
    #   q: vector giving offspring size distribution
    #   m: death rate
    #
    # Value:
    #   matrix of population densities (columns t, rows x)
    
    N <- length(x)-1  # Number of x steps. Also number of Fourier modes
    L <- max(x)-min(x)
    w <- exp(x)
    
    # We strip off the first value of everything because that is identical
    # to the last one by periodicity
    ks <- k[-1]
    gs <- g[-1]
    ws <- w[-1]
    # fft of offspring size distribution
    FqR <- fft(rev(q[-1]))
    
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)
    
    f <- function(t, p, parms) {
        linearPart <- -ks*p-m*p
        birthPart <- rev(2*L/N*Re(fft(
            FqR*(fft(rev(ks*p))), inverse = TRUE)/N))
        growthPart <- rev(Re(fft(fft(rev(gs*p))*k1, inverse=TRUE)/N))/ws
        return(list(linearPart + birthPart + growthPart))
    }
    
    out <- ode(y=p0[-1], times=t, func=f)
    # Replace the times contained in the first column
    # with the value at the boundary.
    out[, 1] <- out[, N+1]
    return(out)
}

fourier_interpolate <- function(p, n) {
    N <- length(p)-1
    x <- seq(0, 1, length.out=n)
    fp <- fft(p[1:N])
    f <- rep(Re(fp[1]), length(x))
    for (j in 2:(N/2+1)) {
        f <- f + 2*(Re(fp[j])*cos(2*pi*(j-1)*x) - Im(fp[j])*sin(2*pi*(j-1)*x))
    }
    return(f/N)
}
