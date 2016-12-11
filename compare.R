# Here I explore how coarse the discretisation can be when working with
# spectral methods as opposed to finite-difference methods.

# Set parameters ----
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
wa <- 0.7;  # threshold for duplication
delta <- 0.2  # width of offspring size distribution
wmin <- wa*(1-delta)/2  # Smallest possible cell size

N <- 4096  # Choose number of steps
x <- seq(log(wmin), 0, length.out = N+1)
w <- exp(x)

# Growth rate
g <- a*w^alpha-b*w^beta

q <- function(w) {
    # Make q nonzero only between (1-delta)/2 and (1+delta)/2
    # Here we use a smooth bump function
    qr <- exp(-1/(1-(2/delta*(w-1/2))^2))/0.444*2/delta
    qr[abs(w-0.5)>=delta/2] <- 0  # Note that we need >= instead of just >
    # to avoid the singularity in the argument to the exponential
    
    return(qr)
}

# Use a k that stays finite but is large enough to ensure that
# almost all cells duplicate before reaching w=1
k <- 10000*(w-wa)^4
k[w<wa] <- 0

# Calculate the steady-state solution ----
sol <- steady_state(t, w, g, k, wa, q, delta)
psi <- sol[[1]]
m <- sol[[2]]

# Plot the solution
par(mar=c(5,5,1,1))
plot(w, psi, type="l", lwd=3,
     xlab="w", ylab=expression(Psi(w)))

# small N ----
select8 <- seq(1, N+1, by=N/8)
w8 <- w[select8]
p8 <- psi[select8]
lines(w8, p8, col="blue")

fp8 <- fourier_interpolate(p8, 4096)
lines(w, fp8, col="red")

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