###
## Funcion que llama a la intregracion RK4, del modelo de ecuaciones diferenciales
## luego plotea todas las poblaciones evolucionadas y las reales
##### Funciones generales
lseq <-
  function(from, to, len) {
    return(exp(seq(log(from), log(to), len = len)))
  }
#####
SEIR_2 <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    
    ####### Por simplicidad entran parametros generales y modifico los parametros para el vector de entrada X=5*3
    ### Nro de indiviuos por rangos etarios, mail dante de E Bologna
    N0 <- 0.82 * N
    N1 <- 0.18 * N
    tt = c(19,30,65)
    #bett = c(0.9  , 0.3 , 0.35 , 0.65  , 0.001 , 0.17)
    #gamm = c(0.0  , 0.04 , 0.06  , 0.    , 0.003  , 0.45)
    #muuu = c(0.0  , 0.0   , 0.0  , 0.0   , 0.05  , 0.19)
    bett = matrix(c(0.24   , 0.05  , 0.04  , 0.018  , 0.12  , 0.04   , 0.055 , 0.02 ),ncol=2)
    gamm = matrix(c(0.0    , 0.001 , 0.018 , 0.008 , 0.00  , 0.001  , 0.01  , 0.005),ncol=2)
    muuu = matrix(c(0.0    , 0.0   , 0.0   , 0.0   , 0.00  , 0.01   , 0.03  , 0.05 ),ncol=2)

    if(t <= tt[1]){
    ##### grupo [0-60], contagia igual que el resto, se reupera más rápido y no muere
       bet <- bett[1,]
       gam <- gamm[1,]
       mu  <- muuu[1,]
    ##### grupo [60-], contagia el doble (actividades todo el dia), se reupera la mitad de lo standard y muere como todos
    } else if( t > tt[1] && t <= tt[2]){
       bet <- bett[2,]
       gam <- gamm[2,]
       mu  <- muuu[2,]
    } else if( t > tt[2] && t <= tt[3]){
       bet <- bett[3,]
       gam <- gamm[3,]
       mu  <- muuu[3,]
    } else {
       bet <- bett[4,]
       gam <- gamm[4,]
       mu  <- muuu[4,]
    }
    #################################################

    S0dot <- -(bet[1] * I0 / N0)  * S0
    I0dot <- (bet[1] * I0 / N0) * S0 - (gam[1] + mu[1]) * I0
    R0dot <- gam[1] * I0 
    M0dot <- mu[1] * I0

    S1dot <- -(bet[2] * I1 / N1) * S1
    I1dot <- (bet[2] * I1 / N1) * S1 - (gam[2] + mu[2]) * I1
    R1dot <- gam[2] * I1
    M1dot <- mu[2] * I1


    ###############################################################
    res <- c(
      S0dot, I0dot, R0dot, M0dot,
      S1dot, I1dot, R1dot, M1dot
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
    if (missing(parms)) parms <- c(bet = 0.55, a = 0.5 / 14.2, gam = 0.0255, mu = 0.02, N = N)
    if (missing(xstart)) {
        xstart <- c(
          S0 = 0.82 * N, I0 = 1, R0 = 0, M0 = 0,
          S1 = 0.18 * N, I1 = 1, R1 = 0, M1 = 0
        )
      }
    # secuencia del tiempo de integracion
    if (missing(times)) times <- seq(0, 100, length = 101)

    out <- rk(xstart, times, SEIR_2, parms, hini = 1, method = "rk4")
    # out es un tipo deSolve
    return(out)
  }
###########################################
ddo <-
  function() {
    source('src/funciones_NG.R')
    simu <- tryRK() # <- integro
    # cargamos los datos reales
    dat <- prepara_datos('src/BM.csv','BM')
    Real <- todos(dat)

     iday <- 1 # elijo dia 20 porque ya hay un buen nuemero de casos
     fday <- length(Real$dates)
     RInf <- cumsum(Real$Sinf)
     RRec <- cumsum(Real$Srec)
     RMue <- cumsum(Real$Sfal)
     Rtime <- seq(iday, fday) #-8#+8

    # sumo los resultados para los 3 subgrupos
    # copio los datos homogeneos
      ST <- simu[, "S0"] + simu[, "S1"] 
   #   ET <- simu[, "E0"] + simu[, "E1"] 
      IT <- simu[, "I0"] + simu[, "I1"] 
      RT <- simu[, "R0"] + simu[, "R1"] 
      MT <- simu[, "M0"] + simu[, "M1"] 

    # ploteo los totales
    # plot (ST, type = "l", col='blue', xlab = "Dias", ylab = "Poblacion",ylim=c(0,3.5E6))
    plot(ST, type = "l", col = "blue", xlab = "Dias", ylab = "Poblacion", ylim = c(1, 3E7), log = "y")
  #  lines(ET, col = "magenta", lty = "dotted", lwd = 2)
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
   # if (!missing(par)) {
   #   tlab01 <- sprintf("I0 %4.2f / E0 %4.2f", xstart[3], xstart[2])
   #   tlab02 <- sprintf("Beta:%2.4f / a:%2.4f / gamma:%2.4f / mu:%2.4f", par[1], par[2], par[3], par[4])
   #   text(60, 10, tlab01)
   #   text(60, 1, tlab02)
   # }
  }

