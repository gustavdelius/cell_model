source("multi-species/lib.R")

steady_state <- function(w, g, k, wa, q, delta) {
    # Determine steady state cell size distribution 
    # from the analytic formulae in eqs.(4.16) to (4.21)
    #
    # Args:
    #   w: vector of equally spaced w-values (cell size)
    #   g: function giving growth rates
    #   k: vector of division rates
    #   wa: threshold for duplication (size of smallest cell that can divide)
    #   q: vector giving offspring size distribution
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
    
    p <- function(Nu0) {
        # Calculate the steady-state solution
        #
        # Args:
        #   m0: scalar giving the constant mortality rate
        #
        # Returns:
        #   List containing the solution and the integral over h/e 
        #   For the correct value of m0 that integral will be equal to 1
        
        # Calculate growth rates for given Nu0
        gv <- g(w, Nu0)
        
        # Calculate e(w) in eq.(4.16)
        ep <- (k+m)/gv
        # First calculate for w < wp
        es <- rev(intx(rev(ep[w<=wp]), rev(w[w<=wp])))
        # then for w >= wp
        el <- -intx(ep[w>=wp], w[w>=wp])
        # and put the results together
        e <- exp(c(es[-length(es)], el))
        
        # Calculate h(w) in eq.(4.17)
        hp <- k*e/gv
        hp[w<wa] <- 0
        h <- c(2*L/N*Re(fft(fft(q[1:N])*fft(hp[1:N]), inverse = TRUE))/N,0)
        
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