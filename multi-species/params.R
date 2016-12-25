setParams <- setClass("PlanktonParams",
    slots = c(
        # Exponents
        xi    = "numeric",
        nu    = "numeric",
        gamma = "numeric",
        
        # Nutrient consumption
        a_inf = "numeric",
        rr    = "numeric",
        # Nutrient replenishment
        rho_0 = "numeric",
        Nu_0  = "numeric",
        
        # Duplication 
        w_th    = "numeric",  # threshold for duplication
        delta_q = "numeric",  # width of offspring size distribution
        k_0     = "numeric",  # scale
        ke      = "numeric",  # exponent
        
        # Cell growth rate (see eq.(2.5))
        alpha = "numeric",
        b     = "numeric",
        beta  = "numeric",
        
        # Death rate
        m = "numeric",
        
        # Predation
        epsilon = "numeric",  # Conversion efficiency
        s0      = "numeric",  # strength of predation
        beta_p  = "numeric" , # log of predator/prey mass ratio
        delta_p = "numeric",  # width of predation kernel
        
        # Rates
        a   = "function",
        dNu = "function",
        k   = "function",
        q   = "function",
        g   = "function",
        S   = "function",
        
        L = "numeric",  # length of single-cell size distribution
        
        # Steady-state solution
        wBar   = "numeric",
        psiBar = "numeric",
        NuBar  = "numeric"
    ),
    prototype = list(
        # Exponents
        xi = 0.15,
        nu = 0.85,
        
        N  = 32L,  # Number of steps within-species
        Ns = 40L,  # Number of species
        SpeciesSpacing = 8L,   # Number of steps between species
        
        # Nutrient consumption
        a_inf = 2,
        rr    = 100,
        abar  = 0.7,  # Nutrient consumption at steady state
        # Nutrient replenishment
        rho_0 = 100,
        Nu_0  = 100,
        
        # Duplication 
        w_th    = 0.7,    # threshold for duplication
        delta_q = 0.2,    # width of offspring size distribution
        k_0     = 10000,  # scale
        ke      = 4,      # exponent
        
        # Cell growth rate
        alpha = 0.85,
        b     = 0.5,
        beta  = 1,
        
        # Death rate
        m = 0.25,
        
        # Predation ----
        epsilon = 0.9,  # Conversion efficiency
        s0      = 0.2,  # strength of predation
        beta_p  = 2,    # log of predator/prey mass ratio
        delta_p = 1     # width of predation kernel
    )
)

setMethod("initialize", "PlanktonParams", 
    function(.Object, ...) {
        .Object <- callNextMethod()
        r <- .Object
        
        r@gamma <- 1 + r@nu + r@xi
        
        # Set rates ----
        # Nutrient dependent feeding coefficient in growth rate
        # See eq.(2.11)
        r@a <- function(Nu) {
            r@a_inf*Nu/(r@rr+Nu)
        }
        # Nutrient growth rate 
        # See eq.(2.12) and (2.13)
        r@dNu <- function(w, Nu, psi, dxs) {
            # Args:
            #   Nu: Nutrient concentration
            #   psi: N x Ns matrix with each column the scaled population of one 
            #        species
            integral <- colSums(w^(r@alpha+1)*psi)*dxs
            r@rho_0*(1-Nu/r@Nu_0) - (r@a(Nu)*sum(w^(2-r@xi-r@gamma)*integral))
        }
        
        # Duplication rate
        # Use a k that stays finite but is large enough to ensure that
        # almost all cells duplicate before reaching w=1
        r@k <- function(w) {
            k <- r@k_0*(w-r@w_th)^r@ke
            k[w<r@w_th] <- 0
            k
        }
        
        # Offspring size distribution 
        # See eq.(2.9) for the definition of q
        # Here we use a smooth bump function
        r@q <- function(w) {
            q <- exp(-1/(1-(2/r@delta_q*(w-1/2))^2))/0.444*2/r@delta_q
            # Make q nonzero only between (1-delta_q)/2 and (1+delta_q)/2
            q[abs(w-0.5)>=r@delta_q/2] <- 0
            q
        }
        
        # Cell growth rate from resource
        # See eq.(2.5)
        r@g <- function(w, Nu) {
            r@a(Nu)*w^r@alpha-r@b*w^r@beta
        }
        
        # Predation kernel
        r@S <- function(x) {
            # Here we use a smooth bump function
            s <- r@s0 * exp(-1/(1-(2*(x-r@beta_p)/r@delta_p)^2))
            # Make s nonzero only between betap-deltap/2 and betap+deltap/2
            s[abs(x-r@beta_p)>=r@delta_p/2] <- 0
            s
        }
        
        wmin <- r@w_th*(1-r@delta_q)/2  # Smallest possible cell size
        xmin <- log(wmin)
        r@L <- -xmin
        
        # Steady-state solution
        N <- 1024
        # equal step sizes in log size
        dx <- r@L/N
        x <- seq(log(wmin), 0, by=dx)
        r@wBar <- exp(x)  # vector of weights
        p <- steady_state(r)
        r@psiBar <- p[[1]]
        r@NuBar <- p[[2]]
        
        return(r)
    }
)

setMethod("plot", "PlanktonParams", 
    function(x, y, ...) {
        plot(x@wBar, x@psiBar, type="l",
             xlab="w/w*", ylab=expression(Psi))   
    }
)

setMethod("show", "PlanktonParams", 
    function(object) {
        cat("Parameters for the plankton model")
    }
)