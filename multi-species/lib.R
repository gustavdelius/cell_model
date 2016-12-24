library("deSolve")

evolve_cell_pop <- function(t, p0, Nu0, r) {
    # Evolve cell population density using the population balance equation
    # see eq.(2.10)
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   p0 : (N+1) x Ns matrix of initial population densities
    #   Nu0: initial nutrient concentration
    #   r: object with parameters
    #
    # Value:
    #   list of two elements:
    #     Nt x (N+1) x Ns  array of population densities
    #     vector of length Nt containing the nutrient densities
    N <- r@N
    Ns <- r@Ns
    
    # Predation
    mkernel <- fft(r@S(-r@xa)*exp(r@xa*(-r@xi)))
    gkernel <- fft(r@epsilon*r@S(r@xa)*exp(r@xa*(r@gamma-2)))
    
    # We strip off the last value of everything because that is identical
    # to the first one by periodicity
    ks <- r@k[-(N+1)]
    wsh <- r@w[-(N+1)]
    # fft of offspring size distribution
    FqR <- fft(r@q[-(N+1)])
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/r@L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)
    
    # Create vectors containing powers
    wsmxi <- r@ws^(-r@xi)
    womxi <- wsh^(1-r@xi)
    wmxi <- wsh^(-r@xi)
    
    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=r@Ns)
        Nu <- pN[length(pN)]
        pcp <- community_spectrum(p, r)
        
        # Calculate growth rate 
        gp <- Re(fft(gkernel*pcp, inverse = TRUE))/r@Na  # from predation
        gr <- r@g(Nu) # from resource
        
        # Calculate death rate
        mp <- Re(fft(mkernel*pcp, inverse = TRUE))/r@Na # from predation
        
        # Calculate right-hand side of population balance equation
        f <- matrix(nrow = N, ncol = Ns)
        idx <- 1:N
        for (i in 1:Ns) {
            gs <- gr + womxi * gp[idx]  # growth rate
            ms <- r@m + wmxi * mp[idx]  # mortality rate
            f[,i] <- wsmxi[i] * (
                -(ks+ms)*p[,i] +  # linear part
                # birth part
                2*r@L/N*Re(fft(FqR*(fft(ks*p[,i])), inverse = TRUE)/N) +
                # growth part
                -Re(fft(fft(gs*p[,i])*k1, inverse=TRUE)/N)/wsh
            )
        }
        nutrientGrowth <- r@dNu(Nu, rbind(0, p))
        # above we added a zero at start of p to give it lenght N+1
        # Return
        list(c(f, nutrientGrowth))
    }
    
    out <- ode(y=c(p0[-(N+1),], Nu0), times=t, func=f)
    
    Nut <- out[ , ncol(out)]
    psit <- array(dim=c(length(t), N+1, Ns))
    psit[ , -(N+1), ] <- out[ , 2:(ncol(out)-1)]
    psit[ , N+1, ] <- psit[ , 1, ]
    list(psit, Nut)
}

community_spectrum <- function(p, r) {
    # Determine community spectrum
    pc <- vector("numeric", length=r@Nal)
    idx <- 1:r@N
    for (i in 1:r@Ns) {
        pc[idx] <- pc[idx] + r@wsmgamma[i]*p[,i]
        idx <- idx + r@ds
    }
    # Pull out a factor of w^{-\gamma} so that pc is constant in steady state
    pc <- r@walgamma*pc
    # Wrap around everything below min(xs)-dxs
    pcp <- pc[(r@Nal-r@Na+1L):r@Nal]
    pcp[(2L*r@Na-r@Nal):r@Na] <- pcp[(2L*r@Na-r@Nal):r@Na] + pc[1:(r@Nal-r@Na+1L)]
    pcp
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
        stop("The function is not periodic or you did not include both endpoints.")
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