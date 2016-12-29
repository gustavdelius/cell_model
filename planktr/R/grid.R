#' An S4 class to represent a 3d grid on which to do a simulation
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
#' @include params.R
setClass("Grid",
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
             t    = "numeric"
         ),
         contains = "Params",
         prototype = list(
             N = 32L,
             Ns = 32L,
             SpeciesSpacing = 8L,
             tmax = 1,
             Nt = 99L
         )
)

#' Perform simulation of multi-species Plankton community
#'
Grid <- function(...) {
    r <- new("Grid", ...)

    # If necessary, call Params
    if (length(r@wBar) == 0) {
        params <- Params(...)
        return(Grid(params, ...))
    }
    # Create grids
    # Within-species cell size grid
    wmin <- r@w_th*(1-r@delta_q)/2  # Smallest possible cell size
    xmin <- log(wmin)
    # equal step sizes in log size
    r@dx <- -xmin/r@N
    r@x <- seq(xmin, -r@dx, by=r@dx)
    r@w <- exp(r@x)  # vector of weights
    r@L <- -xmin

    # Characteristic cell size grid
    # ws denotes the species maximum cell size
    # xs denotes the log of the species maximum cell size
    # We space our species equidistant in log size
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

    # Check that the feeding kernel has support within the community
    if ((r@s0) > 0 &&
        (-r@beta_p-r@delta_p) < (r@xs[1]+r@x[1])) {
        stop("The community spectrum is too short for the given feeding kernel.")
    }

    # Times
    if (length(r@t) == 0) {
        if (r@tmax > 0 && r@Nt > 0) {
            r@t <- seq(0, r@tmax, by=r@tmax/r@Nt)
        } else {
            stop("tmax and Nt must both be positive")
        }
    } else {
        if (r@t[1] != 0) {
            warning("The vector t did not start at 0. I added a zero to it.")
            r@t <- c(0, r@t)
        }
        r@Nt <- length(r@t)
        r@tmax <- r@t[r@Nt]
    }

    return(r)
}

setValidity("Grid", function(object) {
    err <- character()
    if (length(object@Ns) != 1) {
        err <- c(err, "Length of xi should be 1")
    }
    if (object@Ns < 1) {
        err <- c(err, "The number Ns of species can not be less than 1")
    }
    if (length(object@N) != 1) {
        err <- c(err, "Length of N should be 1")
    }
    if (object@N < 1) {
        err <- c(err, "The number of steps can not be less than 1")
    }
    if (length(object@SpeciesSpacing) != 1) {
        err <- c(err, "Length of SpeciesSpacing should be 1")
    }
    if (object@SpeciesSpacing < 1) {
        err <- c(err, "SpeciesSpacing can not be less than 1")
    }
    if (object@N %% object@SpeciesSpacing) {
        err <- c(err,
                 "The number N of steps must be a multiple of the SpeciesSpacing")
    }

    if (length(err) == 0) TRUE else err
})

#' @describeIn Grid Draws a picture of the grid
#' @export
setMethod("plot", "Grid",
          function(x, y, ...) {
              # TODO implement this
              "Not yet implemented"
          }
)

#' @describeIn Grid Show very short description of object
#' @export
setMethod("show", "Grid",
          function(object) {
              cat("Grid for the plankton model")
          }
)

#' @describeIn Params List grid and model parameter values
#' @export
setMethod("summary", "Grid",
    function(object) {
        cat("Grid:\nNs = ", object@Ns,
            ", N = ", object@N,
            ", SpeciesSpacing = ", object@SpeciesSpacing,
            ", tmax = ", object@tmax,
            ", Nt = ", object@Nt,
            "\n"
        )
        callNextMethod()
    }
)
