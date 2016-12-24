setParams <- setClass("PlanktonParams",
         slots = c(
             # Exponents
             xi    = "numeric",
             nu    = "numeric",
             gamma = "numeric",
             
             N  = "integer",  # Number of steps within-species
             Ns = "integer",  # Number of species
             ds = "integer",  # Number of steps between species
             
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
             
             # Grids
             x  = "numeric",
             dx = "numeric",
             w  = "numeric",
             L  = "numeric",
             
             dxs = "numeric",
             xs =  "numeric",
             ws =  "numeric",
             
             xa =  "numeric",
             Na =  "integer",
             La =  "numeric",
             xal = "numeric",
             Nal = "integer",
             
             wsmgamma = "numeric",
             walgamma = "numeric",
             
             # Rates
             a   = "function",
             dNu = "function",
             k   = "numeric",
             q   = "numeric",
             g   = "function",
             S   = "function"
         ),
         prototype = list(
             # Exponents
             xi = 0.15,
             nu = 0.85,
             
             N  = 32L,  # Number of steps within-species
             Ns = 40L,  # Number of species
             ds = 4L,   # Number of steps between species
             
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
        
        # Create grids ----
        # Within-species cell size grid
        wmin <- r@w_th*(1-r@delta_q)/2  # Smallest possible cell size
        # equal step sizes in log size
        r@x <- seq(log(wmin), 0, length.out = r@N+1)
        r@dx <- r@x[2]-r@x[1]
        r@w <- exp(r@x)  # vector of weights
        r@L <- -r@x[1]
        
        # Characteristic cell size grid
        # ws denotes the species maximum cell size
        # xs denotes the log of the species maximum cell size
        # We space our species equidistant in log size
        r@dxs <- r@ds * r@dx  # Spacing of species in log size
        r@xs  <- seq(-(r@Ns-1)*r@dxs, 0, length.out=r@Ns)
        r@ws  <- exp(r@xs)
        
        # Create x steps for entire community spectrum
        # For now we impose a periodicity where the smallest species is also
        # assumed to sit above the largest by the same distance as between the
        # two smallest species.
        r@xa <- seq(r@xs[1]-r@dxs, 0, by=r@dx)
        r@Na <- length(r@xa)
        # Without wrapping around the size spectrum will be longer:
        r@Nal <- (r@Ns-1L)*r@ds+r@N+1L
        r@xal <- seq(r@xs[1]+r@x[1], 0, by=r@dx)
        
        # Create vectors containing powers
        r@wsmgamma <- exp(r@xs*(-r@gamma))
        r@walgamma <- exp(r@xal*(r@gamma))
        
        # Set rates ----
        # Nutrient dependent feeding coefficient in growth rate
        # See eq.(2.11)
        r@a <- function(Nu) {
            r@a_inf*Nu/(r@rr+Nu)
        }
        # Nutrient growth rate 
        # See eq.(2.12) and (2.13)
        r@dNu <- function(Nu, psi) {
            # Args:
            #   Nu: Nutrient concentration
            #   psi: (N+1) x M matrix with each column the scaled population of one 
            #        species
            integral <- colSums(r@w^(r@alpha+1)*psi)*r@dxs
            r@rho_0*(1-Nu/r@Nu_0) - (r@a(Nu)*sum(r@w^(2-r@xi-r@gamma)*integral))
        }
        
        # Duplication rate
        # Use a k that stays finite but is large enough to ensure that
        # almost all cells duplicate before reaching w=1
        r@k <- r@k_0*(r@w-r@w_th)^r@ke
        r@k[r@w<r@w_th] <- 0
        
        # Offspring size distribution 
        # See eq.(2.9) for the definition of q
        # Here we use a smooth bump function
        r@q <- exp(-1/(1-(2/r@delta_q*(r@w-1/2))^2))/0.444*2/r@delta_q
        # Make q nonzero only between (1-delta_q)/2 and (1+delta_q)/2
        r@q[abs(r@w-0.5)>=r@delta_q/2] <- 0  
        # Note that we needed >= instead of just >
        # to avoid the singularity in the argument to the exponential
        
        # Cell growth rate from resource
        # See eq.(2.5)
        r@g <- function(Nu) {
            r@a(Nu)*r@w^r@alpha-r@b*r@w^r@beta
        }
        
        # Predation kernel
        r@S <- function(x) {
            # Here we use a smooth bump function
            s <- r@s0 * exp(-1/(1-(2*(x-r@beta_p)/r@delta_p)^2))
            # Make s nonzero only between betap-deltap/2 and betap+deltap/2
            s[abs(x-r@beta_p)>=r@delta_p/2] <- 0
            s
        }
        
        return(r)
    }
)
