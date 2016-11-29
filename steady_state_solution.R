library("Matrix")

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

# Set parameters
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
xa <- 0.7;  # threshold for duplication
delta <- 0.2  # width of offspring size distribution
xmin <- xa*(1-delta)/2  # Smallest possible cell size

Nx <- 1024  # Choose number of steps
uselog <- TRUE
if (uselog) {
    y <- seq(log(xmin), 0, length.out = Nx+1)
    x <- exp(y)
} else {
    x <- seq(xmin, 1, length.out = Nx+1)
}

xp <- min(x[x>=(0.5 + delta/2)])  # the point at which we glue

# Growth rate
g <- a*x^alpha-b*x^beta

q <- function(x) {
    # Make q nonzero only between (1-delta)/2 and (1+delta)/2
    # Here we use a smooth bump function
    qr <- exp(-1/(1-(2/delta*(x-1/2))^2))/0.444*2/delta
    qr[abs(x-0.5)>=delta/2] <- 0  # Note that we need >= instead of just >
    # to avoid the singularity in the argument to the exponential
    
    # Use the following two lines if you want a step function
    # qr <- rep(0, length(x))
    # qr[abs(x-0.5)<delta/2] <- 1/delta
    
    return(qr)
}

# Use a k that stays finite but is large enough to ensure that
# almost all cells duplicate before reaching x=1
k <- 10000*(x-xa)^4
k[x<xa] <- 0

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
    m <- rep_len(m0, length(x))
    
    # Calculate e(x)
    ep <- (k+m)/g
    # First calculate for x < xp
    es <- rev(intx(rev(ep[x<=xp]), rev(x[x<=xp])))
    # then for x >= xp
    el <- -intx(ep[x>=xp], x[x>=xp])
    # and put the results together
    e <- exp(c(es[-length(es)], el))
    
    # Calculate h(x)
    # TODO: Want to convert this to using spectral methods when the steps
    # are logarithmic
    hp <- k*e/g
    hp[x<xa] <- 0
    h <- rep_len(0, length(x))
    for (i in 1:length(x)) {
        h[i] <- 2*intx(hp*q(x[i]/x)/x, x)[length(x)-1]
    }
    
    # Calculate Theta in eq.(4.21). Because h/e is zero beyond x=xp, we simply
    # integrate up to x=1. Then Theta is automatically 1 for x>xp when the
    # boundary condition (4.18) is satisfied.
    he <- h/e
    he[e==0] <- 0  # Removes the infinity from division by zero
    Theta <- intx(he, x) 
    # The boundary condition (4.18) is satisfied iff b=1
    b <- Theta[length(x)]
    
    # Use eq.(4.20) to calculate the solution.
    Psi <- g[x==xp]*e*Theta/g
    
    return(list(Psi, b))
}

# Find the mortality rate constant that satisfies the boundary condition
# in eq. (4.18)
m0 <- uniroot(function(m0) p(m0)[[2]]-1, lower=0.05, upper=10)[["root"]]

# Calculate the solution
psi <- p(m0)[[1]]
# Plot the solution
par(mar=c(5,5,1,1))
plot(x, psi, type="l", lwd=3,
     xlab="x", ylab=expression(Psi(x)))


###################### new code addition is below ##################

write(x,"solnsteplogx2.txt");
write(psi,"solnsteplogy2.txt");

#write(x,"solxRnewfun.txt");
#write(psi,"solyRnewfun.txt");
m0


#Get space step size (this only works when uselog==FALSE)

dx=x[2]-x[1]


#Use Riemann sum to work out birth term


birthterm <-sapply(x, function(w) 2*dx*sum(k*psi*sapply(w/x,q)/(x)))





#Use central differencing, assuming Z is zero at either end


#centraldiffwZeroes2 <- function(Z,dx){
#  outputtt <- rep(0,length(Z));
#  for (i in 1:length(Z)) {
#    outputtt[i] <- if (i==1) {
#      0
#    } else {
#      if (i==length(Z)) {
#        0
#      } else {
#        (Z[i+1]-Z[i-1])/(2*dx)
#      }
#    }
#  }
#return(outputtt)
#}

# compute linear term equal to minus psi(k(w)+m)

linterm <- -k*psi-m0*psi


#Make 4th order differentiation matrix D for computing growth term

NNN=length(x);
D <- sparseMatrix((1:NNN), c((2:NNN),1), x = 2*rep(1, NNN)/3)-sparseMatrix((1:NNN), c((3:NNN),1,2), x = rep(1, NNN)/12)
D <- (D-t(D))/dx;

#Multilply g*psi by D to get growth term

growthterm <- -D %*% (g*psi)
#Manually set growth term to zero for extreme w values

growthterm[1] <-0
growthterm[2] <-0
#growthterm[length(growthterm3)] <-0
#growthterm[length(growthterm3)-1] <-0

#growthterm[length(growthterm3)-2] <-0
#growthterm[length(growthterm3)-3] <-0


#note since it is a 4th order method, it may not be handelling the boundaries properly


#pde rhs is sum of birth, growth and linear terms

DpsiDt <- birthterm+growthterm+linterm

plot(x,DpsiDt)





######################computing birth term using Riemann sum################

xvals=log(x);
dxlog=xvals[2]-xvals[1]

# Really I should discard the last term when computing Riemann sum, but I think
# this term is zero anyway

birthlogRiemann=sapply(x, function(xchosen) 2*dxlog*sum(k*psi*q(xchosen/x)))

plot(x, birthlogRiemann, type="l")

############################################# computing birth term using fft #############


Fq=fft(sapply(xvals, function(xchosen) q(exp(xchosen))))
Fkpsi=fft(k*psi)

FqR=fft(rev(sapply(xvals, function(xchosen) q(exp(xchosen)))))
FkpsiR=fft(rev(k*psi))


### I should discard the last term again, also I should pad, but I have tried playing around with padding, and not seen much difference. 

birthlogF=2*dxlog*Re(fft(Fq*Fkpsi, inverse = TRUE)/(length(xvals)))

<<<<<<< HEAD
birthlogF2=rev(2*dxlog*Re(fft(FqR*FkpsiR, inverse = TRUE)/(length(xvals))))
=======
plot(x, birthlogRiemann, type="l")
plot(x, birthlogF-birthlogRiemann, type="l")
>>>>>>> a621e3e6752d35a2e978c600b46336f11bd4b0f9



plot(x, birthlogRiemann, type="l")
plot(x, birthlogF2-birthlogRiemann, type="l")


