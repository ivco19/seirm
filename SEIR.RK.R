ddo <-
  function(ou) {
    if (missing(ou)) ou <- tryRK()

    plot(ou[, "S"], type = "l", col = "blue", xlab = "Dias", ylab = "Poblacion")
    lines(ou[, "E"], col = "magenta", lty = "dotted", lwd = 2)
    lines(ou[, "I"], col = "red", lty = "dotted", lwd = 2)
    lines(ou[, "R"], col = "green", lty = "dotted", lwd = 2)
    lines(ou[, "M"], col = "black", lty = "dotted", lwd = 2)
    legend("right", c("Sanos", "Contagiados", "Infectados", "Recuperados", "Muertos"),
      col = c("blue", "magenta", "red", "green", "black"),
      pch = c(21, 21, 21, 21, 21)
    )
  }
#########################################################
# data:
# age : 0-9  / 10-19 / 20-29 / 30-39 / 40-49 / 50-59 / 60-69 / 70-79 / 80+
# mu  : 0.0  / 0.2   / 0.2   / 0.2   / 0.4   / 1.3   / 3.6   / 8.0   / 15.0
# %%% : 17.5 + 17    + 15.5  + 14.5  + 13    + 11.5  + 7.5   + 2     + 1.5
ddo2 <-
  function() {
    NN <- 3.5E6
    p0 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.002, N = NN)
    x0 <- c(S = NN - 11, E = 10, I = 1, R = 0, M = 0)
    out <- tryRK(parms = p0, xstart = x0)
    por <- 0.175
    p0 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.0, N = NN * por)
    x0 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    por <- 0.17 + 0.155 + 0.145
    p1 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.002, N = NN * por)
    x1 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    por <- 0.13
    p2 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.004, N = NN * por)
    x2 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    por <- 0.115
    p3 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.013, N = NN * por)
    x3 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    por <- 0.075
    p4 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.036, N = NN * por)
    x4 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    por <- 0.035
    p5 <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.1, N = NN * por)
    x5 <- c(S = NN * por - 11, E = 10, I = 1, R = 0, M = 0)
    out0 <- tryRK(parms = p0, xstart = x0)
    out1 <- tryRK(parms = p1, xstart = x1)
    out2 <- tryRK(parms = p2, xstart = x2)
    out3 <- tryRK(parms = p3, xstart = x3)
    out4 <- tryRK(parms = p4, xstart = x4)
    out5 <- tryRK(parms = p5, xstart = x5)

    pdf("salida.pdf")
    ddo(out)
    ddo(out0)
    ddo(out1)
    ddo(out2)
    ddo(out3)
    ddo(out4)
    ddo(out5)
    dev.off()
    return()
  }

SEIRMmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    Sdot <- -bet * I / N * S
    Edot <- bet * I / N * S - a * E - gam / 2 * E
    Idot <- a * E - gam / 2 * I - mu * I
    Rdot <- gam / 2 * I + gam / 2 * E
    Mdot <- mu * I
    res <- c(Sdot, Edot, Idot, Rdot, Mdot)
    list(res)
  })
}

tryRK <-
  function(parms, xstart, times) {
    ## The parameters
    library("deSolve")
    if (missing(parms)) parms <- c(bet = 10, a = 1.0 / 5.2, gam = 0.8, mu = 0.2, N = 3.5E6)

    ## vector of timesteps
    if (missing(xstart)) xstart <- c(S = 3.5E6 - 11, E = 10, I = 1, R = 0, M = 0)

    if (missing(times)) times <- seq(0, 180, length = 181)


    out <- rk(xstart, times, SEIRMmod, parms, hini = 1, method = "rk4")


    return(out)
  }
## Dormand-Prince method of order 5(4)
# out3 <- rk(xstart, times, SEIRMmod, parms, hmax = 1,
#            method = "rk45dp7")

# mf <- par("mfrow")
## deSolve plot method for comparing scenarios
# plot(out2, out3, which = c("S", "I", "R"),
#     main = c ("Substrate", "Producer", "Consumer"),
#     col =c("black", "red", "green"),
#     lty = c("solid", "dotted", "dotted"), lwd = c(1, 2, 1))
# }
## user-specified plot function
# plot (out1[,"P"], out1[,"C"], type = "l", xlab = "Producer", ylab = "Consumer")
# lines(out2[,"P"], out2[,"C"], col = "red",   lty = "dotted", lwd = 2)
# lines(out3[,"P"], out3[,"C"], col = "green", lty = "dotted")
#
# legend("center", legend = c("euler", "rk4", "rk45dp7"),
#  lty = c(1, 3, 3), lwd = c(1, 2, 1),
#  col = c("black", "red", "green"))
# par(mfrow = mf)
