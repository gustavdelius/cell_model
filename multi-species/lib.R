library("deSolve")

evolve_cell_pop <- function(t, x, xs, p0, Nu0, g, k, q, m, dNu, S) {
    # Evolve cell population density using the population balance equation
    # see eq.(2.10)
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   x: vector of equally-spaced steps in log cell size
    #      This should include both endpoints of the periodic interval.
    #      So length(x) is one larger than the number of steps
    #   xs: vector of characteristic log cell sizes
    #   p0 : (N+1) x Ns matrix of initial population densities
    #   Nu0: initial nutrient concentration
    #   g: function giving growth rates g(w, Nu)
    #   k: vector of division rates
    #   q: vector giving offspring size distribution
    #   m: death rate
    #   dNu: function giving nutrient growth rate
    #   S: function giving predation kernel
    #
    # Value:
    #   list of two elements:
    #     Nt x (N+1) x Ns  array of population densities
    #     vector of length Nt containing the nutrient densities
    
    N <- length(x)-1  # Number of x steps. Also number of Fourier modes
    Ns <- length(xs)   # Number of species
    if (max(x) != 0) {
        stop("The largest x should be zero, by definition.")
    }
    if (min(x) >= 0) {
        stop("The smallest x must be negative.")
    }
    L <- -min(x)
    dx <- x[2]-x[1]  # Step size in x
    if (diff(range(diff(x))) > .Machine$double.eps ^ 0.5) {
        stop("Our code needs equally-spaced log-size steps")
    }
    if (max(xs) != 0) {
        stop("The largest xs should be zero, by convention.")
    }
    if (min(x) >= 0) {
        stop("The smallest xs must be negative.")
    }
    dxs <- xs[2]-xs[1]  # Step size in xs
    if (diff(range(diff(xs))) > .Machine$double.eps ^ 0.5) {
        stop("Our code needs equally-spaced log-size steps")
    }
    if (dxs %% dx  > .Machine$double.eps ^ 0.5) {
        stop("Our code needs the spacing in species to be a multiple of the size step size dx.")
    }
    
    # Create x steps for entire community spectrum
    # For now we impose a periodicity where the smallest species is also
    # assumed to sit above the largest by the same distance as between the
    # two smallest species.
    ds <- round(dxs/dx)  # number of steps separating species
    xa <- seq(xs[1]-dxs, 0, by=dx)
    Na <- length(xa)
    if (Na != Ns*ds+1) {
        stop("The length of xa is wrong.")
    }
    La <- min(xa)
    # Without wrapping around the size spectrum will be longer:
    Nal <- (Ns-1)*ds+N+1
    xal <- seq(xs[1]+x[1], 0, by=dx)
    # Create weight vectors from log-weight vectors
    w <- exp(x)
    ws <- exp(xs)
    wal <- exp(xal)
    
    # Predation
    mkernel <- fft(S(-xa)*wa^(-xi))
    gkernel <- fft(epsilon*S(xa)*wa^(gamma-2))
    
    # We strip off the last value of everything because that is identical
    # to the first one by periodicity
    ks <- k[-(N+1)]
    wsh <- w[-(N+1)]
    # fft of offspring size distribution
    FqR <- fft(q[-(N+1)])
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)
    
    # Create vectors containing powers
    wsmgamma <- ws^(-gamma)
    walgamma <- wal^(gamma)
    wmxi <- wsh^(-xi)
    womxi <- wsh^(1-xi)
    wsmxi <- ws^(-xi)
    
    ff <- function(p, gs, ms) {
        # Calculate right-hand side of population balance equation
        
    }
    
    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=Ns)
        Nu <- pN[length(pN)]
        pcp <- community_spectrum(p)
        
        # Calculate growth rate 
        gp <- Re(fft(gkernel*pcp, inverse = TRUE))/Na  # from predation
        gr <- g(wsh, Nu) # from resource
        
        # Calculate death rate
        mp <- Re(fft(mkernel*pcp, inverse = TRUE))/Na # from predation
        
        # Calculate right-hand side of population balance equation
        f <- matrix(nrow = N, ncol = Ns)
        idx <- 1:N
        for (i in 1:Ns) {
            gs <- gr + womxi * gp[idx]  # growth rate
            ms <- m + wmxi * mp[idx]  # mortality rate
            f[,i] <- wsmxi[i] * (
                -(ks+ms)*p[,i] +  # linear part
                # birth part
                2*L/N*Re(fft(FqR*(fft(ks*p[,i])), inverse = TRUE)/N) +
                # growth part
                -Re(fft(fft(gs*p[,i])*k1, inverse=TRUE)/N)/wsh
            )
        }
        nutrientGrowth <- dNu(Nu, rbind(0, p))
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

community_spectrum <- function(p) {
    # Determine community spectrum
    pc <- vector("numeric", length=Nal)
    idx <- 1:N
    for (i in 1:Ns) {
        pc[idx] <- pc[idx] + wsmgamma[i]*p[,i]
        idx <- idx + ds
    }
    # Pull out a factor of w^{-\gamma} so that pc is constant in steady state
    pc <- walgamma*pc
    # Wrap around everything below min(xs)-dxs
    pcp <- pc[(Nal-Na+1):Nal]
    pcp[(2*Na-Nal):Na] <- pcp[(2*Na-Nal):Na] + pc[1:(Nal-Na+1)]
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