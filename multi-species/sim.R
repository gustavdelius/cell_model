#' An S4 class to represent the simulation of a plankton model
#' 
#' @slot N Number of steps within-species
#' @slot Ns Number of species
#' @slot SpeciesSpacing Number of species per length of species 
#' 
#' @slot x  Log sizes within a species.
#' @slot dx log step size within a species.
#' @slot w  Sizes within a species.
#' @slot ds Number of steps between species.
#' @slot dxs Log distance between species.
#' @slot xs Log maximum sizes of species.
#' @slot ws Maximum sizes of species.
#' @slot xa Log sizes within community.
#' @slot Na Number of steps in community spectrum.
#' @slot xal Log sizes in community spectrum before periodic wrapping.
#' @slot Nal Number of steps in community spectrum before periodic wrapping.
#' @slot wsmgamma \code{ws^(-gamma)}
#' @slot walgamma \code{wal^(gamma)}
#' @slot tmax Time until which to run simulation. If the vector \code{t} of
#' times is supplied than \code{tmax} is set from that.
#' @slot Nt Number of time steps at which to return population density.
#' If the vector \code{t} is supplied than \code{Nt} is set from that.
#' @slot t vector of \code{Nt} time steps between 0 and \code{tmax}. If this 
#' is not supplied then the default is a vector of \code{Nt} equally-spaced
#' time points between 0 and \code{tmax}.
#' @slot p0 Matrix (N x Ns) of initial population densities.
#' @slot Nu0 Initial resource concentration.
#' @slot p Array (Nt x N x Ns) of simulated population densities.
#' @slot Nu Vector of simulated nutrient concentrations.
simulate <- setClass("PlanktonSim",
    slots = c(
        N  = "integer",
        Ns = "integer",
        SpeciesSpacing = "integer",
        
        # Grids
        # within a species
        x  = "numeric",
        dx = "numeric",
        w  = "numeric",
        # maximal species sizes
        ds  = "integer",
        dxs = "numeric",
        xs =  "numeric",
        ws  =  "numeric",
        # entire spectrum
        xa =  "numeric",
        Na =  "integer",
        # entire spectrum before wrapping
        xal = "numeric",
        Nal = "integer",
        # some powers of weights
        wsmgamma = "numeric",
        walgamma = "numeric",
        
        # Times
        tmax = "numeric",
        Nt   = "integer",
        t    = "numeric",
        
        # Population
        p0  = "matrix",
        Nu0 = "numeric",
        p   = "array",
        Nu  = "numeric"
    ),
    contains = "PlanktonParams",
    prototype = list(
        N = 32L,
        Ns = 2L,
        SpeciesSpacing = 1L,
        tmax = 10,
        Nt = 100L
    )
)

setMethod("initialize", "PlanktonSim", 
    function(.Object, ...) {
        .Object <- callNextMethod()
        r <- .Object
        
        # Create grids ----
        # Within-species cell size grid
        wmin <- r@w_th*(1-r@delta_q)/2  # Smallest possible cell size
        xmin <- log(wmin)
        # equal step sizes in log size
        r@dx <- -xmin/r@N
        r@x <- seq(log(wmin), -r@dx, by=r@dx)
        r@w <- exp(r@x)  # vector of weights
        r@L <- -xmin
        
        # Characteristic cell size grid
        # ws denotes the species maximum cell size
        # xs denotes the log of the species maximum cell size
        # We space our species equidistant in log size
        if (r@N %% r@SpeciesSpacing) {
            stop("The number N of steps must be a multiple of the SpeciesSpacing")
        }
        r@ds = r@N%/%r@SpeciesSpacing
        r@dxs <- r@ds * r@dx  # Spacing of species in log size
        r@xs  <- seq(-(r@Ns-1)*r@dxs, 0, length.out=r@Ns)
        r@ws  <- exp(r@xs)
        
        # Create x steps for entire community spectrum
        # For now we impose a periodicity where the smallest species is also
        # assumed to sit above the largest by the same distance as between the
        # two smallest species.
        r@xa <- seq(r@xs[1]-r@dxs, -r@dx, by=r@dx)
        r@Na <- length(r@xa)
        # Without wrapping around the size spectrum will be longer:
        r@Nal <- (r@Ns-1L)*r@ds+r@N
        r@xal <- seq(r@xs[1]+r@x[1], -r@dx, by=r@dx)
        
        # Create vectors containing powers
        r@wsmgamma <- exp(r@xs*(-r@gamma))
        r@walgamma <- exp(r@xal*(r@gamma))
        
        # Times
        if (length(r@t) == 0) {
            r@t <- seq(0, r@tmax, by=r@tmax/r@Nt)
        } else {
            r@Nt <- length(r@t)
            r@tmax <- r@t[r@Nt]
        }
        
        # Population
        if (length(r@p0) == 0) {
            # if no initial population is provided use steady-state solution
            # First subsample at N points
            psi<-fourier_interpolate(r@psiBar[-1025], 32)
            # Then replicate this over all species
            p0 = matrix(psi, nrow=r@N, ncol=r@Ns)
            # Then normalise it so that the nutrient is at steady-state
            # For this we observe that in eq.(2.12) the sigma is proportional to psi
            # So we get \rho and \sigma to cancel by rescaling \psi -> psi * rho/sigma
            # Alternatively see eqs.(5.33)-(5.35)
            integral <- colSums(r@w^(r@alpha+1)*p0)*r@dx
            r@p0 <- p0 * r@rho_0*(1-r@NuBar/r@Nu_0) / 
                (r@a(r@NuBar)*sum(r@ws^(2-r@xi-r@gamma)*integral))
        } else if ((length(r@p0) != r@N) && (length(r@p0) != (r@N * r@Ns))) {
            # If only a single species is provided we will replicate this
            # but if a strange number of intial values is given we complain
            stop("p0 has the incorrect length")
        } else {
            r@p0 = matrix(p0, nrow=r@N, ncol=r@Ns)
        }
        
        if (length(r@Nu0) == 0) {
            # If no initial resource is given, use steady-state resource
            r@Nu0 <- r@NuBar
        }
        
        p <- evolve_cell_pop(r)
        r@p <- p[[1]]
        r@Nu <- p[[2]]
        
        return(r)
    }
)

