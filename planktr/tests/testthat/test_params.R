library(planktr)
context("Params")

test_that("steady state for default parameters is correct", {
    r <- Params()
    expect_equal(r@NuBar, 53.2794114834)
    expect_equal(mean(r@psiBar), 0.579918927586)
    expect_equal(max(r@psiBar), 1.31807481093)
})

test_that("steady state with predation is correct", {
    r <- Params(s0=5)
    expect_equal(r@NuBar, 36.6211602164)
    expect_equal(mean(r@psiBar), 0.532160558808)
    expect_equal(max(r@psiBar), 1.22087271895)
})

test_that("steady state with non-default parameters is correct", {
    r <- Params(xi = 0.2,
                nu = 0.7,
                # Nutrient consumption
                a_inf = 4,
                rr    = 50, # Duplication
                w_th    = 0.6,    # threshold for duplication
                delta_q = 0.1,    # width of offspring size distribution
                k_0     = 100000,  # scale
                ke      = 3.8,      # exponent

                # Cell growth rate
                alpha = 0.7,
                b     = 0.8,
                beta  = 1.1,

                # Death rate
                m = 0.5,

                # Predation
                epsilon = 0.5,
                s0      = 3,
                beta_p  = 1,
                delta_p = 1
                )
    expect_equal(r@NuBar, 22.6140431893)
    expect_equal(mean(r@psiBar), 0.636945584657)
    expect_equal(max(r@psiBar), 1.60510703523)
})
