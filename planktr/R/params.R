#' An S4 class to represent the parameters of a plankton model
#'
#' @slot xi Exponent for allometric scaling of mortality. Default 0.15
#' @slot nu Exponent in scaling of predation kernel. Default 0.85
#' @slot gamma Exponent in across-species steady-state power-law.
#' Equal to 1+nu+xi.
#' @slot a_inf Nutrient consumption rate when unlimited nutrient supply.
#' Default 2
#' @slot rr Half saturation level at which nutrient consumption is half of the
#' maximal value \code{a_inf}.
#' @slot w_th Threshold size for duplication. Cells smaller than \code{w_th}
#' can not duplicate.
#' @slot delta_q Width of offspring size distribution. Daughter cells can have
#' a weight that lies between \code{(1-delta_q)/2} and \code{(1+delta_q)/2}
#' of their parent's weight.
#' @slot k_0 Overall coefficient of duplication rate function.
#' @slot ke Exponent in duplication rate function.
#' @slot alpha Allometric exponent of nutrient intake rate.
#' @slot beta Allometric exponent of metabolic loss rate.
#' @slot b Coefficient of metabolic loss rate.
#' @slot m Coefficient of death rate.
#' @slot epsilon Conversion efficiency of prey mass into predator mass.
#' @slot s0 Coefficient of predation kernel.
#' @slot beta_p Center of predation kernel.
#' @slot delta_p Width of predation kernel.
#' @slot a Function \code{a(Nu)} giving nutrient intake rate coefficient as a
#' function of nutrient concentration.
#' @slot dNu Function \code{dNu(w, Nu, psi, r)} giving the nutrient growth
#' rate for a cell of weight \code{w} as a function of nutrient concentration
#' \code{Nu} and plankton population \code{psi}.
#' @slot k Function \code{k(w)} giving duplication rate for a cell of weight
#' \code{w}
#' @slot g Function \code{g(w, Nu)} giving the growth due to resource
#' consumption for a cell of weight \code{w} as a function of resource
#' concentration \code{Nu}.
#' @slot s Function \code{s(x)} giving the predation kernel as a function of
#' the logarithm \code{x} of the predator/prey mass ratio. The predation kernel
#' needs to be multiplied by the allometric scaling factor and the predator
#' and prey densities to give the actual predation rate.
#' @slot L Length of single-cell size distribution in log-size.
#' @slot wBar Weights at which the steady-state population is given.
#' @slot psiBar Population density at steady state.
#' @slot NuBar Nutrient concentration at steady state.
setClass("Params",
    slots = c(
        # Exponents
        xi    = "numeric",
        nu    = "numeric",
        gamma = "numeric",

        # Nutrient consumption
        a_inf = "numeric",
        rr    = "numeric",

        # Duplication
        w_th    = "numeric",
        delta_q = "numeric",
        k_0     = "numeric",
        ke      = "numeric",

        # Cell growth rate (see eq.(2.5))
        alpha = "numeric",
        b     = "numeric",
        beta  = "numeric",

        # Death rate
        m = "numeric",

        # Predation
        epsilon = "numeric",
        s0      = "numeric",
        beta_p  = "numeric" ,
        delta_p = "numeric",

        # Rates
        a   = "function",
        k   = "function",
        q   = "function",
        g   = "function",
        s   = "function",

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

        # Nutrient consumption
        a_inf = 2,
        rr    = 100,

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

        # Predation
        epsilon = 0.9,  # Conversion efficiency
        s0      = 0,  # strength of predation
        beta_p  = 2,    # log of predator/prey mass ratio
        delta_p = 1     # width of predation kernel
    )
)

# Initialise object ----
Params <- function(...) {
    r <- new("Params", ...)

    r@gamma <- 1 + r@nu + r@xi

    # Set rates ----
    # Nutrient dependent feeding coefficient in growth rate
    # See eq.(2.11)
    r@a <- function(Nu) {
        r@a_inf*Nu/(r@rr+Nu)
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
    r@s <- function(x) {
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

setValidity("Params", function(object) {
    err <- character()
    r <- object
    if (length(r@xi) != 1) {
        msg <- paste("Length of xi should be 1")
        err <- c(errors, msg)
    }
    if (r@beta_p-r@delta_p/2 < 0) {
        err <- c(err, "The feeding kernel is too wide.")
    }

    if (length(err) == 0) TRUE else err
})

#' @describeIn Params Plot the steady-state solution
#' @export
setMethod("plot", "Params",
    function(x, y) {
        plot(x@wBar, x@psiBar, type="l",
             xlab="w/w*", ylab=expression(Psi))
    }
)

#' @describeIn Params Show very short description of object
#' @export
setMethod("show", "Params",
    function(object) {
        cat("Parameters for the plankton model")
    }
)

#' @describeIn Params List model parameter values
#' @export
setMethod("summary", "Params",
    function(object) {
        cat("Exponents:\n",
            "  xi = ", object@xi,
            ", nu = ", object@nu,
            ", gamma = ", object@gamma,
            "\nNutrient consumption:\n",
            " a_inf = ", object@a_inf,
            ", rr = ", object@rr,
            "\nDuplication:\n",
            "  w_th = ", object@w_th,
            ", delta_q = ", object@delta_q,
            ", k_0  = ", object@k_0,
            ", ke = ", object@ke,
            "\nCell growth rate:\n",
            "  alpha = ", object@alpha,
            ", b = ", object@b,
            ", beta = ", object@beta,
            "\nDeath rate: \n",
            "  m = ", object@m,
            "\nPredation: \n",
            "  epsilon = ", object@epsilon,
            ", s0 = ", object@s0,
            ", beta_p = ", object@beta_p,
            ", delta_p = ", object@delta_p
        )
    }
)
