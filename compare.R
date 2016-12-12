# Here I explore how coarse the discretisation can be when working with
# spectral methods as opposed to finite-difference methods.

library("rgl")
source("lib.R")

evolve_cell_pop_fd <- function(t, x, p0, g, k, q, m) {
    # Evolve cell population density with second-order finite-difference scheme
    #
    # Args:
    #   t: vector of times at which to return the population density
    #   x: vector of equally spaced steps in logarithmic cell size
    #   p0 : vector on initial population density
    #   g: vector of growth rates
    #   k: vector of division rates
    #   q: function giving offspring size distribution
    #   m: death rate
    #
    # Value:
    #   matrix of population densities (columns t, rows x)
    N <- length(p)  # Number of x steps. Also number of Fourier modes
    dx <- x[2]-x[1]
    w <- exp(x)
    
    f <- function(t, p, parms) {
        N <- length(p)
        
        linearPart <- -k*p-m*p
        
        birthPart <- rep_len(0, N)
        for (i in 1:N) {
            birthPart[i] <- 2*intx(k*p*q(w[i]/w)/w, w)[N-1]
        }
        
        gp <- c(0, g*p, 0)
        growthPart <- (gp[1:N]-gp[3:(N+2)])/(2*dx)/w
        
        return(list(linearPart + birthPart + growthPart))
    }
    out <- ode(y=p0, times=t, func=f)
    return(out[, -1])
}

# Set parameters ----
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
wa <- 0.7;  # threshold for duplication
delta <- 0.2  # width of offspring size distribution
wmin <- wa*(1-delta)/2  # Smallest possible cell size

N <- 1024  # Choose number of steps
x <- seq(log(wmin), 0, length.out = N+1)
w <- exp(x)

# Growth rate
g <- a*w^alpha-b*w^beta

# Division rate
# Use a k that stays finite but is large enough to ensure that
# almost all cells duplicate before reaching w=1
k <- 10000*(w-wa)^4
k[w<wa] <- 0

# Offspring size distribution
q <- function(w) {
    # Make q nonzero only between (1-delta)/2 and (1+delta)/2
    # Here we use a smooth bump function
    qr <- exp(-1/(1-(2/delta*(w-1/2))^2))/0.444*2/delta
    qr[abs(w-0.5)>=delta/2] <- 0  # Note that we need >= instead of just >
    # to avoid the singularity in the argument to the exponential
    
    return(qr)
}

# Calculate the steady-state solution ----
sol <- steady_state(t, w, g, k, wa, q, delta)
psi <- sol[[1]]
m <- sol[[2]]

# Plot the solution
par(mar=c(5,5,1,1))
plot(w, psi, type="l", lwd=3, xlab="w", ylab=expression(Psi(w)))

# small N ----
select8 <- seq(1, N+1, by=N/8)
p8 <- approx(w[select8], psi[select8], n=1025)
lines(p8, col="blue", lwd=2)

fp8 <- fourier_interpolate(psi[select8], 1024)
lines(w, fp8, col="red", lwd=2)

plot(p8[[1]], approx(w, psi, n=1025)[[2]]-p8[[2]], col="blue", lwd=2, type="l")
lines(w, psi-fp8, col="red", lwd=2)

# medium N ----
plot(w, psi, type="l", lwd=3, xlab="w", ylab=expression(Psi(w)))
select32 <- seq(1, N+1, by=N/32)
p32 <- approx(w[select32], psi[select32], n=1025)
lines(p32, col="blue", lwd=2)

fp32 <- fourier_interpolate(psi[select32], 1024)
lines(w, fp32, col="red", lwd=2)

plot(p32[[1]], approx(w, psi, n=1025)[[2]]-p32[[2]], col="blue", lwd=2, type="l")
lines(w, psi-fp32, col="red", lwd=2)

# medium N ----
select128 <- seq(1, N+1, by=N/128)
w128 <- w[select128]
p128 <- psi[select128]
lines(w128, p128, col="blue")

fp128 <- fourier_interpolate(p128, 4096)
lines(w, fp128, col="red")

# Solve equation ----
tmax <- 10  # final time
Nt <- 100    # number of time steps at which to store intermediate values
t <- seq(0, tmax, by=tmax/Nt)

p0 <- exp(-100*(w - 0.6)^2)
plot(w, p0, type="l")

s <- seq(2, N+1, by=N/128)
p <- evolve_cell_pop(t, x[s], p0[s], g=g[s], k=k[s], q=q(w[s]), m)
persp3d(t, w[s], p, col = "lightblue")

pfd <- evolve_cell_pop_fd(t, x[s], p0[s], g=g[s], k=k[s], q, m)
persp3d(t, w[s], pfd, col = "lightblue")

s <- seq(2, N+1, by=N/32)
p <- evolve_cell_pop(t, x[s], p0[s], g=g[s], k=k[s], q=q(w[s]), m)
pp <- t(apply(p, 1, fourier_interpolate, 256))
sp <- seq(1, N+1, by=N/256)
persp3d(t, w[sp], pp, col = "lightblue")

pfd <- evolve_cell_pop_fd(t, x[s], p0[s], g=g[s], k=k[s], q, m)
persp3d(t, w[s], pfd, col = "lightblue")

select16 <- seq(1, N+1, by=N/16)
w16 <- w[select16]
p16 <- psi[select16]
lines(w16, p16, col="blue")

fp16 <- fourier_interpolate(p16, 4096)
lines(w, fp16, col="red")

select4 <- seq(1, N+1, by=N/4)
w4 <- w[select4]
p4 <- psi[select4]
lines(w4, p4, col="green")

fp4 <- fourier_interpolate(p4, 4096)
lines(w, fp4, col="red")