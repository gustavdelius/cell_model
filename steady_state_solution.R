library("Matrix")
library("deSolve")
library("rgl")


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

Nx <- 2048  # Choose number of steps
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

#I asume uselog==False, so I switch off the non-spectral stuff

dx=x[2]-x[1]


#Use Riemann sum to work out birth term


#birthterm <-sapply(x, function(w) 2*dx*sum(k*psi*sapply(w/x,q)/(x)))


# compute linear term equal to minus psi(k(w)+m)

#linterm <- -k*psi-m0*psi


#Make 4th order differentiation matrix D for computing growth term

#NNN=length(x);
#D <- sparseMatrix((1:NNN), c((2:NNN),1), x = 2*rep(1, NNN)/3)-sparseMatrix((1:NNN), c((3:NNN),1,2), x = rep(1, NNN)/12)
#D <- (D-t(D))/dx;

#Multilply g*psi by D to get growth term

#growthterm <- -D %*% (g*psi)
#Manually set growth term to zero for extreme w values

#growthterm[1] <-0
#growthterm[2] <-0



#note since it is a 4th order method, it may not be handelling the boundaries properly


#pde rhs is sum of birth, growth and linear terms

#DpsiDt <- birthterm+growthterm+linterm

#plot(x,DpsiDt)

#plot(x,birthterm)




######################computing birth term using Riemann sum################

xvals=log(x);
dxlog=xvals[2]-xvals[1]

# Really I should discard the last term when computing Riemann sum, but I think
# this term is zero anyway

#I turn off this part for the sake of speed up, because it was mainly done for testing

#birthlogRiemann=sapply(x, function(xchosen) 2*dxlog*sum(k*psi*q(xchosen/x)))

#plot(x, birthlogRiemann, type="l")

############################################# computing birth term using fft #############


#computing fft s of reversals of components of convolution integral



#I comment out the code below that is mainly for testing

#birthlogF2=rev(2*dxlog*Re(fft(FqR*FkpsiR, inverse = TRUE)/(length(xvals))))
#plot(x, birthlogRiemann, type="l")
#plot(x, birthlogF2-birthlogRiemann, type="l")



############################# log growth term 4th order difference

# I comment out this code which computes the spatial derivatitives up to 4th order via matrix multiplication, although I still wonder if this is better than 
# the fft based approach I used in the code since there is a small discrepency between the two, although I'm hoping it is the fft approach which is more accurate.

#dxxx=xvals[2]-xvals[1]
#NNN=length(x);
#D <- sparseMatrix((1:NNN), c((2:NNN),1), x = 2*rep(1, NNN)/3)-sparseMatrix((1:NNN), c((3:NNN),1,2), x = rep(1, NNN)/12)
#D <- (D-t(D))/dxxx;
#growthtermlog <- (-D %*% (g*psi))/x
#growthtermlog[1] <-0
#growthtermlog[2] <-0


#function to determine growth term




# I dont use this anymore, but the code below uses matrix multiplication to get the derivatives, note the output object is a matrix, so it needs to be converted back 
# for some purposes.

#hpp=rep(1, length(x))
#dxxx=xvals[2]-xvals[1]
#NNN=length(x);
#D <- sparseMatrix((1:NNN), c((2:NNN),1), x = 2*rep(1, NNN)/3)-sparseMatrix((1:NNN), c((3:NNN),1,2), x = rep(1, NNN)/12)
#D <- (D-t(D))/dxxx;
#getgrowthterm2 <- function(psi) {
#    growthtermlog <- (-D %*% (g*psi))/x
#    for (i in 1:length(x)) {
#        hpp[i]=growthtermlog[i]
#    }
#    return(hpp)
#}


############################# log growth term, Fourier method (double reverse)


#I comment out this part, which is mainly for testing

#Fgrowthlog=(2*pi/(max(xvals)-min(xvals)))*rev(Re(fft(fft(rev(g*psi))*(1i*c(0:(length(g)/2-1),0,(-length(g)/2+1):-1)), inverse=TRUE)/length(g)))/x
#plot(x,Fgrowthlog-growthtermlog)

#############compute linear term and combine to get d psi/dt

#I comment out this part, which is mainly for testing


#lintermlog <- -k*psi-m0*psi
#derpsilog=lintermlog+growthtermlog+birthlogRiemann
#plot(x,derpsilog)


#I changed it so we throw away the first term, and take xval[2]... end as the Nx points over one period, note this change has not been done with commented out test code



FqR=fft(rev(sapply(xvals[-1], function(xchosen) q(exp(xchosen)))))
FkpsiR=fft(rev(k[-1]*psi[-1]))


derpsilog2=   (2*pi/(max(xvals)-min(xvals)))*rev(Re(fft(fft(rev(g[-1]*psi[-1]))*(1i*c(0:(length(g[-1])/2-1),0,(-length(g[-1])/2+1):-1)), inverse=TRUE)/length(g[-1])))/x[-1]+
    (-k[-1]*psi[-1]-m0*psi[-1])+
    rev(2*(max(xvals)-min(xvals))/(length(psi[-1]))*Re(fft(FqR*(fft(rev(k[-1]*psi[-1]))), inverse = TRUE)/(length(xvals[-1]))))


plot(x[-1],derpsilog2)



##################################pde solver

#this is the core of the PDE solver, it determines p psi/dt using spectral methods




fff <- function(t, psis, parms) {
    list(
        (2*pi/(max(xvals)-min(xvals)))*rev(Re(fft(fft(rev(g[-1]*psis))*(1i*c(0:(length(g[-1])/2-1),0,(-length(g[-1])/2+1):-1)), inverse=TRUE)/length(g[-1])))/x[-1]+
            (-k[-1]*psis-m0*psis)+
            rev(2*(max(xvals)-min(xvals))/(length(psis))*Re(fft(FqR*(fft(rev(k[-1]*psis))), inverse = TRUE)/(length(xvals[-1]))))
        
    )
}




#The rest of the code runs this and plots it

tmax <- 20  # final time
Nt <- 50    # number of time steps at which to store intermediate values
t <- seq(0, tmax, by=tmax/Nt)
parms <- list(dog=5)

#psis is psi with the first entry dropped, because of new period convention

spsi=psi[-1]

#Our input rpsi is a random pertubation of the solution generated by the code at the start

rpsi=spsi+0.05*runif(length(spsi))

#psigauss=exp(-100*(x - 0.6)^2)

out <- ode(y=rpsi, times=t, func=fff, parms=parms)

sselectx <- seq(2, length(spsi), by=1)
xshort=x[-1]
xss <- xshort[sselectx]
yss <- out[,sselectx]
persp3d(t, xss, yss, col = "lightblue")



#Testing spectral derivatives
#Ga=getgrowthterm2(psi);

#Gc=(2*pi/(max(xvals)-min(xvals)+xvals[2]-xvals[1]))*rev(Re(fft(fft(rev(g*psi))*(1i*c(0:(length(g)/2-1),0,(-length(g)/2+1):-1)), inverse=TRUE)/length(g)))/x


#plot(x,Ga-Gc)

