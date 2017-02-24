#' Produce 3 dimenstional plot of the simulation result
#'
#' @param r Object of type Sim
#' @param t Time at which to plot. If the time is not one of the values in r@t,
#' choose next larger time contained in r@t or, if there is none, choose
#' largest time in r@t.
#' @param xs Log size of species to plot. If there is no species at this size,
#' choose the next larger species, or if there is none, the largest species.
#' @include sim.R
setGeneric("plot3D", function(r, t, xs, ...) {standardGeneric("plot3D")})

#' @describeIn plot3D 3D plot of population density against x and xa at given time
setMethod("plot3D", c(r="Sim", t="numeric", xs="missing"),
    function(r, t, zlog=FALSE) {
        ti <- which(r@t >= t)[1]
        if (is.na(ti)) {
            ti <- length(r@t)
            message("Defaulting to t=", r@t[ti])
        }
        z <- if (zlog) log(pmax(r@p[ti, , ], 0.05)) else r@p[ti, , ]
        zlab <- if (zlog) "log(p)" else "p"
        persp3D(r@x, r@xs, z,
                xlab="x", ylab="x*", zlab=zlab,
                main=paste("Spectrum at t=", r@t[ti]),
                ticktype="detailed", nticks=4)
    }
)

#' @describeIn plot3D 3D plot of population density of given species against t and x
setMethod("plot3D", c(r="Sim", t="missing", xs="numeric"),
    function(r, xs, zlog=FALSE) {
        xsi <- which(r@xs >= xs)[1]
        if (is.na(xsi)) {
            xsi <- length(r@xs)
            message("Defaulting to xs=", r@xs[xsi])
        }
        z <- if (zlog) log(pmax(r@p[ , , xsi], 0.05)) else r@p[ , , xsi]
        zlab <- if (zlog) "log(p)" else "p"
        persp3D(r@t, r@x, z,
                xlab="t", ylab="x", zlab=zlab,
                main=paste("Spectrum at xs=", r@xs[xsi]),
                ticktype="detailed", nticks=4)

    }
)

#' @describeIn plot3D 3D plot of community density against x and xa
setMethod("plot3D", c(r="Sim", t="missing", xs="missing"),
          function(r, zlog=FALSE, smooth=TRUE) {
              com <- get_community(r, smooth = smooth)
              com <- if (zlog) log(pmax(com, 0.05)) else com
              zlab <- if (zlog) "log(p_c)" else "p_c"
              persp3D(r@t, r@xa, com,
                      xlab="t", ylab="x", zlab=zlab,
                      main="Community size spectrum",
                      ticktype="detailed", nticks=4)
          }
)
