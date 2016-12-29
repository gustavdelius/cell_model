#' An S4 class to represent the simulation of a plankton model

#' @slot p Array (Nt x N x Ns) of simulated population densities.
#' @slot Nu Vector of simulated nutrient concentrations.
#' @include params.R
#' @include grid.R
setClass("Sim",
    slots = c(
        p   = "array",
        Nu  = "numeric"
    ),
    contains = "Grid"
)

#' Perform simulation of multi-species Plankton community
#'
Sim <- function(..., p0=NULL, Nu0=NULL) {
    r <- new("Sim", ...)

    # If necessary, call Grid
    if (length(r@x) == 0) {
        grid <- Grid(...)
        return(Sim(grid, ..., p0=p0, Nu0=Nu0))
    }
    # Population
    if (is.null(p0)) {
        p0 <- make_p0(r)
    } else if ((length(p0) != r@N) && (length(p0) != (r@N * r@Ns))) {
        # If only a single species is provided we will replicate this
        # but if a strange number of intial values is given we complain
        stop("p0 has the incorrect length")
    } else {
        p0 <- matrix(p0, nrow=r@N, ncol=r@Ns)
    }

    if (is.null(Nu0)) {
        # If no initial resource is given, use steady-state resource
        Nu0 <- r@NuBar
    }

    p <- evolve_cell_pop(p0, Nu0, r)
    r@p <- p[[1]]
    r@Nu <- p[[2]]

    return(r)
}

#' Extract grid from simulation
#'
#' This can be achieved just as easy with as(). I am providing this function
#' only lest I might forget how to use as()
#' @param sim Object of class Sim
#' @return Object of class Grid
getGrid <- function(sim) {
    as(sim, "Grid")
}

#' @describeIn Sim Plot the solution for one species at one time
#' @param sim Sim object
#' @param t Time at which to plot. If the value is not available
#' at that time, the last earlier time is used. Default: latest available time.
#' @param s Index of species to plot. Default 1.
#' @export
setMethod("plot", "Sim",
    function(x, y, t=NULL, s=1) {
        if (is.null(t)) {
            ti <- sim@Nt+1
        } else if (t >= 0) {
            ti <- which(sim@t >= t)[1] - 1
        } else {
            stop("The time can not be negative")
        }
        plot(sim@w, sim@p[ti, , s], type="l",
             xlab="w", ylab=paste("p_1(t=", sim@t[ti], ")"),
             main = paste("Size spectrum for species", s, "at t=", sim@t[ti]))
    }
)

#' @describeIn Sim Show very short description of object
#' @export
setMethod("show", "Sim",
          function(object) {
              cat("A simulation of the plankton model")
          }
)
