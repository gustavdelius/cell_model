#' An S4 class to represent the simulation of a plankton model

#' @slot p Array (Nt x N x Ns) of simulated population densities.
#' @slot Nu Vector of simulated nutrient concentrations.
#' @slot rho_0 Nutrient replenishment rate when nutrient is low.
#' @slot Nu_0 Nutrient carrying capacity in absence of consumption.
#' @slot dNu Function giving nutrient replenishment rate.
#' @include params.R
#' @include grid.R
setClass("Sim",
    slots = c(
        p   = "array",
        Nu  = "numeric",

        # Nutrient replenishment
        rho_0 = "numeric",
        Nu_0  = "numeric",
        dNu = "function"
    ),
    contains = "Grid"
)

#' Perform simulation of multi-species Plankton community
#'
#' @param grid Object of class Grid. If missing or an object of class Params is
#' provided this is extended using the further arguments
#' @param params Object of class Params. Only used if grid is missing.
#' @param p0  N x Nx matrix with initial population density.
#' Defaults to steady-state solution
#' @param Nu0 Initial nutrient concentration. Defaults to steady state value.
#' @param ... Arguments that will be used to initialise a Grid object if non
#' is provided.
#' @export
Sim <- function(grid=NULL, params=NULL, p0=NULL, Nu0=NULL, ...) {
    if (is.null(grid)) {
        if(is.null(params)) {
            grid <- Grid(...)
        } else {
            assert_that(class(params) == "Params")
            grid <- Grid(params=params, ...)
        }
    }
    assert_that(is(grid, "Grid"))

    # Population
    if (is.null(p0)) {
        p0 <- make_p0(grid)
    } else if ((length(p0) != grid@N) && (length(p0) != (grid@N * grid@Ns))) {
        # If only a single species is provided we will replicate this
        # but if a strange number of intial values is given we complain
        stop("p0 has the incorrect length")
    } else {
        assert_that(is.numeric(p0))
        p0 <- matrix(p0, nrow=grid@N, ncol=grid@Ns)
    }

    if (is.null(Nu0)) {
        # If no initial resource is given, use steady-state resource
        Nu0 <- grid@NuBar
    }

    # Arbitrary choice for carrying capacity at twice the initial value
    Nu_0  <- 2 * Nu0
    # Set rho_0 so that Nu0 is the steady-state resource for the given p0
    integral <- colSums(grid@w^(grid@alpha+1)*p0)*grid@dx
    rho_0 <- (grid@a(Nu0)*sum(grid@ws^(2-grid@xi-grid@gamma)*integral)) /
        (1-Nu0/Nu_0)
    # Nutrient growth rate
    # See eq.(2.12) and (2.13)
    dNu <- function(w, Nu, psi, r) {
        # Args:
        #   Nu: Nutrient concentration
        #   psi: N x Ns matrix with each column the scaled population of one
        #        species
        integral <- colSums(w^(r@alpha+1)*psi)*r@dx
        r@rho_0*(1-Nu/r@Nu_0) - (r@a(Nu)*sum(r@ws^(2-r@xi-r@gamma)*integral))
    }

    sim <- new("Sim", grid, rho_0=rho_0, Nu_0=Nu_0, dNu=dNu)

    p <- evolve_cell_pop(p0, Nu0, sim)
    sim@p  <- p[[1]]
    sim@Nu <- p[[2]]

    sim
}

#' Extract grid from simulation
#'
#' This can be achieved just as easy with as(). I am providing this function
#' only lest I might forget how to use as()
#' @param sim Object of class Sim
#' @return Object of class Grid
#' @export
getGrid <- function(sim) {
    as(sim, "Grid")
}

#' Extract parameters from simulation
#' @param sim Object of class Sim
#' @return Object of class Grid
#' @export
getParams <- function(sim) {
    as(sim, "Params")
}

#' @describeIn Sim Plot the solution for one species at one time
#' @param sim Sim object
#' @param t Time at which to plot. If the value is not available
#' at that time, the next larger time is used. Default: latest available time.
#' @param s Index of species to plot. Default 1.
#' @export
setMethod("plot", "Sim",
    function(x, y, t=NULL, s=1, xlog=FALSE) {
        if (is.null(t)) {
            ti <- x@Nt+1
        } else {
            ti <- which(x@t >= t)[1]
            if (is.na(ti)) {
                ti <- length(r@t)
                message("Defaulting to t=", r@t[ti])
            }
        }
        xx <- if (xlog) x@x else x@w
        xlab <- if (xlog) "x" else "y"

        plot(xx, x@p[ti, , s], type="l",
             xlab=xlab, ylab=paste("p_1(t=", x@t[ti], ")"),
             main = paste("Size spectrum for species", s, "at t=", x@t[ti]))
    }
)

#' @describeIn Sim Show very short description of object
#' @export
setMethod("show", "Sim",
          function(object) {
              cat("A simulation of the plankton model")
          }
)

#' @describeIn Sim List grid and model parameter values
#' @export
setMethod("summary", "Sim",
          function(object) {
              callNextMethod()
              cat("\nNutrient replenishment:\n",
                  "  rho_0 = ", object@rho_0,
                  ", Nu_0 = ", object@Nu_0
              )
          }
)
