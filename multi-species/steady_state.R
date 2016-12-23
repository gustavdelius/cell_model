source("multi-species/lib.R")

steady_state <- function(w, g, k, wa, q, delta) {
    # Determine steady state cell size distribution 
    # from the analytic formulae in eqs.(4.16) to (4.21)
    #
    # Args:
    #   w: vector of equally spaced w-values (cell size)
    #   g: vector of growth rates
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
        hp <- k*e/g
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