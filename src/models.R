# SERIM
# es un mini modelo SEIR notas en el codigo
# Fuentes:
# https://web.stanford.edu/~jhj1/teachingdocs/Jones-on-R0.pdf
# http://indico.ictp.it/event/7960/session/3/contribution/19/material/slides/0.pdf
# https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
# Resumen:
# Hice un modelito SEIR , compartimente por simplicidad en 3 grupos etarios:
# [0] = 0-19
# [1] = 20-49
# [2] = 50+
# el modelo tiene mas modificaciones ver comments en el codigo

### Funciones generales
lseq <-
  function(from, to, len) {
    return(exp(seq(log(from), log(to), len = len)))
  }
##### Modelo SEIR Simple
SEIR <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    ####### Por simplicidad entran parametros generales y modifico los parametros para el vector de entrada X=5*3
    N0 <- N
    bet0 <- bet
    gam0 <- gam
    mu0 <- mu

    S0dot <- -bet0 * I0 * S0 / N0
    E0dot <- bet0 * S0 / N0 * I0 - a * E0
    I0dot <- a * E0 - (gam0 + mu0) * I0
    R0dot <- gam0 * I0
    M0dot <- mu0 * I0
    ###############################################################
    res <- c(S0dot, E0dot, I0dot, R0dot, M0dot)
    list(res)
  })
}

##### Modelo SEIR modificado
SEIRmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    ### parametros
    N0 <- N
    bet0 <- bet
    gam0 <- gam
    mu0 <- mu

    S0dot <- -bet0 * I0 * S0 / N0 - bet0 * E0 * S0 / N0
    E0dot <- bet0 * S0 / N0 * I0 - (a + gam0) * E0 + bet0 * E0 * S0 / N0
    I0dot <- a * E0 - (gam0 + mu0) * I0
    R0dot <- gam0 * I0 + gam0 * E0
    M0dot <- mu0 * I0
    ###############################################################
    res <- c(S0dot, E0dot, I0dot, R0dot, M0dot)
    list(res)
  })
}

SEIR_3 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    ####### Por simplicidad entran parametros generales y modifico los parametros para el vector de entrada X=5*3
    ### Nro de indiviuos por rangos etarios, mail dante de E Bologna
    N0 <- 0.304 * N
    N1 <- 0.430 * N
    N2 <- 0.266 * N
    ##### grupo [0-19], contagia igual que el resto, se reupera más rápido y no muere
    bet0 <- bet
    gam0 <- gam * 2
    mu0 <- 0.0
    a0 <- a * 2
    ##### grupo [20-49], contagia el doble (actividades todo el dia), se reupera la mitad de lo standard y muere como todos
    bet1 <- bet
    gam1 <- gam
    mu1 <- mu
    a1 <- a
    ##### grupo [50-inf], se contagia la mitad (asumo que estan guardados), se reupera la mitad de lo standard y el doble (o triple)
    bet2 <- bet / 2.
    gam2 <- gam / 2.
    mu2 <- mu * 2
    a2 <- a / 2
    #################################################

    #S0dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0
    #E0dot <- (bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0 - (a0 + gam0 * 0.0) * E0
    S0dot <- -(bet0 * I0 / N0)  * S0
    E0dot <- (bet0 * I0 / N0) * S0 - a0 * E0
    I0dot <- a0 * E0 - (gam0 + mu0) * I0
    R0dot <- gam0 * I0 + gam0 * E0
    M0dot <- mu0 * I0

    #S1dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1
    #E1dot <- (bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1 - (a1 + gam1 * 0.0) * E1
    S1dot <- -(bet1 * I1 / N1) * S1
    E1dot <- (bet1 * I1 / N1) * S1 - a1 * E1
    I1dot <- a1 * E1 - (gam1 + mu1) * I1
    R1dot <- gam1 * I1 + gam1 * 0.0 * E1
    M1dot <- mu1 * I1

    S2dot <- -(bet2 * I2 / N2) * S2
    E2dot <- (bet2 * I2 / N2) * S2 - a1 * E2
    I2dot <- a1 * E2 - (gam2 + mu2) * I2
    R2dot <- gam2 * I2
    M2dot <- mu2 * I2

    ###############################################################
    res <- c(
      S0dot, E0dot, I0dot, R0dot, M0dot,
      S1dot, E1dot, I1dot, R1dot, M1dot,
      S2dot, E2dot, I2dot, R2dot, M2dot
    )
    list(res)
  })
}

SEIR_3mod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {

    ####### Por simplicidad entran parametros generales y modifico los parametros para el vector de entrada X=5*3
    ### Nro de indiviuos por rangos etarios, mail dante de E Bologna
    N0 <- 0.304 * N
    N1 <- 0.430 * N
    N2 <- 0.266 * N
    ##### grupo [0-19], contagia igual que el resto, se reupera más rápido y no muere
    bet0 <- bet
    gam0 <- gam
    mu0 <- 0.0
    ##### grupo [20-49], contagia el doble (actividades todo el dia), se reupera la mitad de lo standard y muere como todos
    bet1 <- bet
    gam1 <- gam
    mu1 <- mu
    ##### grupo [50-inf], se contagia la mitad (asumo que estan guardados), se reupera la mitad de lo standard y el doble (o triple)
    bet2 <- bet / 2.
    gam2 <- gam / 2.
    mu2 <- mu * 2
    #################################################

    S0dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0 - bet0 * E0 * S0 / N0
    E0dot <- (bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0 - (a + gam0) * E0 + bet0 * E0 * S0 / N0
    I0dot <- a * E0 - (gam0 + mu0) * I0
    R0dot <- gam0 * I0 + gam0 * E0
    M0dot <- mu0 * I0

    S1dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1 - bet1 * E1 * S1 / N1
    E1dot <- (bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1 - (a + gam1) * E1 + bet1 * E1 * S1 / N1
    I1dot <- a * E1 - (gam1 + mu1) * I1
    R1dot <- gam1 * I1 + gam1 * E1
    M1dot <- mu1 * I1

    S2dot <- -(bet1 * I1 / N1 + bet2 * I2 / N2) / 2 * S2 - bet2 * E2 * S2 / N2
    E2dot <- (bet1 * I1 / N1 + bet2 * I2 / N2) / 2 * S2 - a * E2 + bet2 * E2 * S2 / N2
    I2dot <- a * E2 - (gam2 + mu2) * I2
    R2dot <- gam2 * I2
    M2dot <- mu2 * I2

    ###############################################################
    res <- c(
      S0dot, E0dot, I0dot, R0dot, M0dot,
      S1dot, E1dot, I1dot, R1dot, M1dot,
      S2dot, E2dot, I2dot, R2dot, M2dot
    )
    list(res)
  })
}

##### Integrador de modelos
tryRK <-
  function(model,parms, xstart, times) {
    ## The parameters
    library("deSolve")
    ####################### DO IT
    N <- 3.5E6
    if (missing(model)) model <- "SEIR"
    if (model == "SEIR" || model == "SEIRmod") {
      if (missing(parms)) parms <- c(bet = 0.44, a = 0.5 / 14.2, gam = 0.0355, mu = 0.02, N = N)
      if (missing(xstart)) xstart <- c(S0 = N - 13, E0 = 10, I0 = 1, R0 = 1, M0 = 1)
    } else if (model == "SEIR_3" || model == "SEIR_3mod") {
      if (missing(parms)) parms <- c(bet = 0.55, a = 0.5 / 14.2, gam = 0.0255, mu = 0.02, N = N)
      if (missing(xstart)) {
        xstart <- c(
          S0 = 0.3 * N - 1025, E0 = 1000, I0 = 20, R0 = 5, M0 = 0,
          S1 = 0.5 * N - 1196, E1 = 1000, I1 = 180, R1 = 15, M1 = 1,
          S2 = 0.2 * N - 125, E2 = 80, I2 = 35, R2 = 7, M2 = 3
        )
      }
    } else {
      print(c("Model ", model, " not implemented"))
      return(0)
    }
    # secuencia del tiempo de integracion
    if (missing(times)) times <- seq(0, 90, length = 91)
    model_in <- switch(model, "SEIRmod" = SEIRmod, "SEIR" = SEIR, "SEIR_3mod" = SEIR_3mod, "SEIR_3" = SEIR_3)

    out <- rk(xstart, times, model_in, parms, hini = 1, method = "rk4")
    # out es un tipo deSolve
    return(out)
  }
