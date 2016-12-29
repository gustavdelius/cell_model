#' Produce 3 dimenstional plot of the simulation result
#'
#' @param r Object of type Sim
#' @param t Time at which to plot. If the time is not one of the values in r@t,
#' choose next larger time contained in r@t or, if there is none, choose
#' largest time in r@t.
#' @param xs Log size of species to plot. If there is no species at this size,
#' choose the next larger species, or if there is none, the largest species.
#' @include sim.R
setGeneric("plot3d", function(r, t, xs) {standardGeneric("plot3d")})

#' @describeIn plot3d 3d plot of population density against x and xa at given time
setMethod("plot3d", c(r="Sim", t="numeric", xs="missing"),
    function(r, t) {
        ti <- which(r@t >= t)[1]
        if (is.na(ti)) {
            ti <- length(r@t)
            warning("Defaulting to t=", r@t[ti])
        }
        persp3d(r@x, r@xs, r@p[ti, , ], col = "lightblue",
                xlab="x", ylab="xs", zlab="p",
                main=paste("Spectrum at t=", r@t[ti]))
    }
)

#' @describeIn plot3d 3d plot of population density of given species against t and x
setMethod("plot3d", c(r="Sim", t="missing", xs="numeric"),
    function(r, xs) {
        xsi <- which(r@xs >= xs)[1]
        if (is.na(xsi)) {
            xsi <- length(r@xs)
            warning("Defaulting to xs=", r@xs[xsi])
        }
        persp3d(r@t, r@x, r@p[ , , xsi], col = "lightblue",
                xlab="t", ylab="x", zlab="p",
                main=paste("Spectrum at xs=", r@xs[xsi]))
    }
)

#' @describeIn plot3d 3d plot of community density against x and xa
setMethod("plot3d", c(r="Sim", t="missing", xs="missing"),
          function(r) {
              com <- get_community(r)
              persp3d(r@t, r@xa, com, col = "lightblue",
                      xlab="t", ylab="xa", zlab="p_c",
                      main="Community size spectrum")
          }
)
