###
## Funcion que llama a la intregracion RK4, del modelo de ecuaciones diferenciales
## luego plotea todas las poblaciones evolucionadas y las reales
##
ddo <-
  function(par, xstart, model) {
    source("src/models.R")
    source("src/load_data.R")
    if (missing(par)) {
      simu <- tryRK() # <- integro
    } else {
      simu <- tryRK(parms = par, xstart = xstart, model)
    } # <- integro
    # cargamos los datos reales
    Real <- load_data(download = FALSE)
    iday <- 20 # elijo dia 20 porque ya hay un buen nuemero de casos
    fday <- end(Real[, 1])[1]
    RInf <- Real[, 4][iday:fday]
    RMue <- Real[, 1][iday:fday]
    RRec <- Real[, 3][iday:fday]
    Rtime <- seq(iday, fday) #-8#+8

    # sumo los resultados para los 3 subgrupos
    # copio los datos homogeneos
    if (dim(simu)[2] == 6) {
      ST <- simu[, "S0"]
      ET <- simu[, "E0"]
      IT <- simu[, "I0"]
      RT <- simu[, "R0"]
      MT <- simu[, "M0"]
    } else {
      ST <- simu[, "S0"] + simu[, "S1"] + simu[, "S2"]
      ET <- simu[, "E0"] + simu[, "E1"] + simu[, "E2"]
      IT <- simu[, "I0"] + simu[, "I1"] + simu[, "I2"]
      RT <- simu[, "R0"] + simu[, "R1"] + simu[, "R2"]
      MT <- simu[, "M0"] + simu[, "M1"] + simu[, "M2"]
    }

    # ploteo los totales
    # plot (ST, type = "l", col='blue', xlab = "Dias", ylab = "Poblacion",ylim=c(0,3.5E6))
    plot(ST, type = "l", col = "blue", xlab = "Dias", ylab = "Poblacion", ylim = c(1, 3E7), log = "y")
    lines(ET, col = "magenta", lty = "dotted", lwd = 2)
    lines(IT, col = "red", lty = "dotted", lwd = 2)
    lines(RT, col = "green", lty = "dotted", lwd = 2)
    lines(MT, col = "black", lty = "dotted", lwd = 2)
    # overplot los reales
    points(Rtime, RInf, col = "red")
    points(Rtime, RMue, col = "black")
    points(Rtime, RRec, col = "green")
    legend("right", c("Sanos", "Contagiados", "Infectados", "Recuperados", "Muertos"),
      col = c("blue", "magenta", "red", "green", "black"),
      pch = c(21, 21, 21, 21, 21)
    )
    if (!missing(par)) {
      tlab01 <- sprintf("I0 %4.2f / E0 %4.2f", xstart[3], xstart[2])
      tlab02 <- sprintf("Beta:%2.4f / a:%2.4f / gamma:%2.4f / mu:%2.4f", par[1], par[2], par[3], par[4])
      text(60, 10, tlab01)
      text(60, 1, tlab02)
    }
  }

