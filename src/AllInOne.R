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
    tt = c(20,38,45,55)
    #bett = c(0.9  , 0.3 , 0.35 , 0.65  , 0.001 , 0.17)
    #gamm = c(0.0  , 0.04 , 0.06  , 0.    , 0.003  , 0.45)
    #muuu = c(0.0  , 0.0   , 0.0  , 0.0   , 0.05  , 0.19)
    #                t0    t1       t2      t3      t4 
    bett = matrix(c(0.3 , 0.045 , 0.045 , 0.035 , 0.04 ,
                    0.38, 0.045 , 0.11  , 0.035 , 0.03 ),ncol=2)
    gamm = matrix(c(0.0 , 0.007 , 0.026 , 0.03  , 0.015, 
                    0.0 , 0.005 , 0.015 , 0.025 , 0.01 ),ncol=2)
    muuu = matrix(c(0.0 , 0.0   , 0.0   , 0.0   , 0.00 , 
                    0.0 , 0.005 , 0.006 , 0.008 , 0.005),ncol=2)

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
    } else if( t > tt[3] && t <= tt[4]){
       bet <- bett[4,]
       gam <- gamm[4,]
       mu  <- muuu[4,]
    } else {
       bet <- bett[5,]
       gam <- gamm[5,]
       mu  <- muuu[5,]
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
          S0 = 0.82 * N, I0 = 13, R0 = 0, M0 = 0,
          S1 = 0.18 * N, I1 = 4, R1 = 0, M1 = 0
        )
      }
    # secuencia del tiempo de integracion
    bday=15
    if (missing(times)) times <- seq(bday, 100, length = 101-bday)

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
    dat  <- prepara_datos('src/BM.csv','BM')
    Real <- todos(dat)
    
    bday=15

     iday <- bday # elijo dia 20 porque ya hay un buen nuemero de casos
     fday <- length(Real$dates)

     Rtime <- seq(iday, fday) - bday
    
     RInf0 <- cumsum(Real$Sinf0[iday:fday])
     RInf1 <- cumsum(Real$Sinf1[iday:fday])
    
     RRec0 <- cumsum(Real$Srec0[iday:fday])
     RRec1 <- cumsum(Real$Srec1[iday:fday])
     
     RMue <- cumsum(Real$Sfal1[iday:fday])

      ST0 <- simu[, "S0"] 
      ST1 <- simu[, "S1"]
      IT0 <- simu[, "I0"]
      IT1 <- simu[, "I1"]
      RT0 <- simu[, "R0"] 
      RT1 <- simu[, "R1"] 
      MT0 <- simu[, "M0"] 
      MT1 <- simu[, "M1"] 

    # ploteo los totales
    # plot (ST, type = "l", col='blue', xlab = "Dias", ylab = "Poblacion",ylim=c(0,3.5E6))
    plot(ST0, type = "l", lty="dotted", lwd = 2, col = "blue", xlab = "Dias - 15 de FIS_0", ylab = "Poblacion", ylim = c(1,2E3), log = "y")
  #  lines(ET, col = "magenta", lty = "dotted", lwd = 2)
    lines(IT0, col = "red", lty = "dotted", lwd = 2)
    lines(RT0, col = "green", lty = "dotted", lwd = 2)
#    lines(MT0, col = "black", lty = "dotted", lwd = 2)
    lines(ST1, col = "blue" , lty = "dashed", lwd = 2)
    lines(IT1, col = "red"  , lty = "dashed", lwd = 2)
    lines(RT1, col = "green", lty = "dashed", lwd = 2)
    lines(MT1, col = "black", lty = "dashed", lwd = 2)
    # overplot los reales
    points(Rtime, RInf0, col = "red")
    points(Rtime, RInf1, col = "red",pch=5)
    points(Rtime, RMue,  col = "black",pch=5)
    points(Rtime, RRec0, col = "green")
    points(Rtime, RRec1, col = "green",pch=5)
    legend("topleft", c("Sanos", "Contagiados", "Infectados", "Recuperados", "Muertos"),
      col = c("blue", "magenta", "red", "green", "black"),
      pch = c(21, 21, 21, 21, 21)
    )
    legend("bottomright",c("Menores a 60","Mayores a 60"),pch=c(21,5))
   # if (!missing(par)) {
   #   tlab01 <- sprintf("I0 %4.2f / E0 %4.2f", xstart[3], xstart[2])
   #   tlab02 <- sprintf("Beta:%2.4f / a:%2.4f / gamma:%2.4f / mu:%2.4f", par[1], par[2], par[3], par[4])
   #   text(60, 10, tlab01)
   #   text(60, 1, tlab02)
   # }
  }

