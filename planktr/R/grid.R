#' An S4 class to represent a 3d grid on which to do a simulation
#'
#' @slot Ns Number of species
#' @slot N Number of steps within-species
#' @slot ds Number of steps between species.
#'
#' @slot x  Log sizes within a species.
#' @slot dx log step size within a species.
#' @slot w  Sizes within a species.
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
             Ns = "integer",
             N  = "integer",
             ds  = "integer",

             # Grids
             # within a species
             x  = "numeric",
             dx = "numeric",
             w  = "numeric",
             # maximal species sizes
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
         contains = "Params"
)

#' Set up grid for multi-species Plankton simulation
#'
#' @param params Object of class Params. If missing, default params are used.
#' If an object of class Grid is supplied, then its Params subobject is used.
#' @param NS Number of species
#' @param N  Number of steps within-species
#' @param ds Number of steps between species
#' @param t  vector of time steps at which to return the solution. If this
#' is not supplied then the default is a vector of \code{Nt} equally-spaced
#' time points between 0 and \code{tmax}.
#' @param tmax Time until which to run simulation. If the vector \code{t} of
#' times is supplied than \code{tmax} is set from that.
#' @param Nt Number of time steps at which to return population density.
#' If the vector \code{t} is supplied than \code{Nt} is set from that.
#' @return Object of type Grid
Grid <- function(params=NULL, Ns = 32L, N = 32L, ds= 4L,
                 t=NULL, tmax = 1, Nt = 99L, ...) {
    if (is.null(params)) {
        params <- Params(...)
    }
    assert_that(is(params, "Params"))
    assert_that(is.count(Ns))
    assert_that(is.count(N))
    assert_that(is.count(ds))
    assert_that(N > 8)
    assert_that(ds <= N)

    # Times
    if (is.null(t)) {
        if (tmax > 0 && Nt > 0) {
            t <- seq(0, tmax, by=tmax/Nt)
        } else {
            stop("tmax and Nt must both be positive")
        }
    } else {
        if (t[1] != 0) {
            warning("The vector t did not start at 0. I added a zero to it.")
            t <- c(0, t)
        }
        Nt <- length(t)
        tmax <- t[Nt]
        assert_that(tmax > 0)
    }

    r <- new("Grid", params, N=as.integer(N), Ns=as.integer(Ns),
             ds=as.integer(ds), tmax=tmax, Nt=as.integer(Nt), t=t, ...)

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

    return(r)
}

setValidity("Grid", function(object) {
    err <- c(validate_that(is.count(object@Ns)),
             validate_that(is.count(object@N)),
             validate_that(object@N >= 8),
             validate_that(is.count(object@ds)),
             validate_that(object@ds <= object@N),
             validate_that(object@Nt > 0)
    )
    err <- err[err != "TRUE"]
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

#' @describeIn Grid List grid and model parameter values
#' @export
setMethod("summary", "Grid",
    function(object) {
        cat("Grid:\nNs = ", object@Ns,
            ", N = ", object@N,
            ", ds = ", object@ds,
            ", tmax = ", object@tmax,
            ", Nt = ", object@Nt,
            "\n"
        )
        callNextMethod()
    }
)
