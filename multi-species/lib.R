library("deSolve")

evolve_cell_pop <- function(t, x, xs, p0, Nu0, g, k, q, m, dNu) {
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
    dx <- x[2]-x[1]
    if (diff(range(diff(x))) > .Machine$double.eps ^ 0.5) {
        stop("Our code needs equally-spaced log-size steps")
    }
    if (max(xs) != 0) {
        stop("The largest xs should be zero, by convention.")
    }
    if (diff(range(diff(xs) %% dx))  > .Machine$double.eps ^ 0.5) {
        stop("Our code needs the spacing in species to be a multiple of the size step size dx.")
    }
    # Create x steps for entire community spectrum
    xa <- seq(min(xs)+min(x), 0, by=dx)
    Na <- length(xa)
    La <- min(xa)
    # Create index vector locating species maximal sizes in the xa vector
    indxs <- Na + xs / dx
    # Create weight vectors from log-weight vectors
    w <- exp(x)
    ws <- exp(xs)
    
    # Create vector containing ws^{-\gamma}
    wsgamma <- ws^{-gamma}
    # Create a matrix with N copies of column vector ws^(-xi). This is needed
    # later to implement the multiplication by ws^(-xi) in the fastest way 
    wsm <- rep(ws^(-xi), rep(N, Ns))
    # We strip off the last value of everything because that is identical
    # to the first one by periodicity
    ks <- k[-(N+1)]
    wsh <- w[-(N+1)]
    # fft of offspring size distribution
    FqR <- fft(q[-(N+1)])
    # For calculating first derivative by Fourier transform
    k1 <- (2*pi/L)*1i*c(0:(N/2-1),0,(-N/2+1):-1)
    
    ff <- function(p, gs, pc) {
        # Calculate right-hand side of population balance equation
        -(ks+m)*p +  # linear part
            # birth part
            2*L/N*Re(fft(FqR*(fft(ks*p)), inverse = TRUE)/N) +
            # growth part
            -Re(fft(fft(gs*p)*k1, inverse=TRUE)/N)/wsh
    }
    
    f <- function(t, pN, parms) {
        p <- matrix(pN[-length(pN)], ncol=Ns)
        Nu <- pN[length(pN)]
        # Determine community spectrum
        pc <- vector("numeric", length=Na)
        for (i in 1:Ns) {
            pc[(indxs[i]-N):(indxs[i]-1)] <- 
                pc[(indxs[i]-N):(indxs[i]-1)] + wsgamma[i]*p[,i]
        }
        # Calculate growth rate
        gs <- g(wsh, Nu)
        # Calculate right-hand side of population balance equation
        f <- apply(p, 2, ff, gs) * wsm
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