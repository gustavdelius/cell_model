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
    c(0,cumsum((f[idx]+f[idx+1])*dx[idx]/2))
}

steady_state <- function(t, w, g, k, wa, q, delta) {
    # Determine steady state cell size distribution 
    # from the analytic formulae in eqs.(4.16) to (4.21)
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
    
    # Without loss of generality we can set w_* = 1 so that x=w/w_*=w
    # The population density for any other w_* can be obtained by scaling
    # see eqs.(6.1) and (6.2)
    
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
        
        # Calculate e(w) in eq.(4.16)
        ep <- (k+m)/g
        # First calculate for w < wp
        es <- rev(intx(rev(ep[w<=wp]), rev(w[w<=wp])))
        # then for w >= wp
        el <- -intx(ep[w>=wp], w[w>=wp])
        # and put the results together
        e <- exp(c(es[-length(es)], el))
        
        # Calculate h(w) in eq.(4.17)
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
    
    # Calculate and return the solution
    psi <- p(m)[[1]]
    list(psi, m)
}

evolve_cell_pop <- function(t, w, ws, p0, Nu0, g, k, q, m, dNu) {
    # Evolve cell population density using the population balance equation
    # see eq.(2.10)
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   w: vector of logarithmic spaced steps in cell size
    #      This should include both endpoints of the periodic interval.
    #      So length(W) is one larger than the number of steps
    #   ws: vector of characteristic cell sizes
    #   p0 : (N+1) x M matrix of initial population densities
    #   Nu0: initial nutrient concentration
    #   g: function giving growth rates g(w, Nu)
    #   k: vector of division rates
    #   q: vector giving offspring size distribution
    #   m: death rate
    #   dNu: function giving nutrient growth rate
    #
    # Value:
    #   list of two elements:
    #     Nt x (N+1) x M  array of population densities
    #     vector of length Nt containing the nutrient densities
    
    N <- length(w)-1  # Number of x steps. Also number of Fourier modes
    M <- length(ws)   # Number of species
    
    L <- log(max(w))-log(min(w))
    # create a matrix with N copies of column vector ws^(-xi). This is needed
    # later to implement the multiplication by ws^(-xi) in the fastest way 
    wsm <- rep(ws^(-xi), rep(N, M))
    # We strip off the first value of everything because that is identical
    # to the last one by periodicity
    ks <- k[-1]
    wsh <- w[-1]
    # fft of offspring size distribution
    FqR <- fft(rev(q[-1]))
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)
    
    ff <- function(p, gs) {
        -(ks+m)*p +  # linear part
            # birth part
            rev(2*L/N*Re(fft(FqR*(fft(rev(ks*p))), inverse = TRUE)/N)) +
            # growth part
            rev(Re(fft(fft(rev(gs*p))*k1, inverse=TRUE)/N))/wsh
    }
    
    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=M)
        Nu <- pN[length(pN)]
        gs <- g(wsh, Nu)
        f <- apply(p, 2, ff, gs) * wsm
        nutrientGrowth <- dNu(Nu, rbind(0, p))
        # above we added a zero at start of p to give it lenght N+1
        # Return
        list(c(f, nutrientGrowth))
    }
    
    out <- ode(y=c(p0[-1,], Nu0), times=t, func=f, parms=parms)
    
    Nut <- out[ , ncol(out)]
    psit <- array(dim=c(length(t), N+1, M))
    psit[ , -1, ] <- out[ , 2:(ncol(out)-1)]
    psit[ , 1, ] <- psit[ , N+1, ]
    list(psit, Nut)
}

fourier_interpolate <- function(p, n) {
    # Perform a Fourier interpolation
    # Args:
    #   p: vector of values of function at equally spaced steps,
    #      including both endpoints of the periodic interval
    #   n: desired length of output vector
    # Value:
    #   Vector of n values of Fourier interpolation at n equally
    #     spaced points, including both endpoints.
    if (p[length(p)] != p[1]) {
        error("The function is not periodic or you did not include both endpoints.")
    }
    N <- length(p)-1
    x <- seq(0, 1, length.out=n)
    fp <- fft(p[1:N])
    f <- rep(Re(fp[1]), length(x))
    for (j in 2:(N/2+1)) {
        f <- f + 2*(Re(fp[j])*cos(2*pi*(j-1)*x) - Im(fp[j])*sin(2*pi*(j-1)*x))
    }
    f/N
}
