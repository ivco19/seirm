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


############################################################
### data:
### age : 0-9  / 10-19 / 20-29 / 30-39 / 40-49 / 50-59 / 60-69 / 70-79 / 80+
### mu  : 0.0  / 0.2   / 0.2   / 0.2   / 0.4   / 1.3   / 3.6   / 8.0   / 15.0
### %%% : 17.5 + 17    + 15.5  + 14.5  + 13    + 11.5  + 7.5   + 2     + 1.5  % Porcentaje de la poblacion en Argentina
### faltan referencias aca pero son los primer hits del google en tasa de mortalidad y pirameida poblacional Argentina
############################################################
############################################################
## DEPRECATED
## Funcion que llama a la intregracion RK4, del modelo de ecuaciones diferenciales para una sola poblacion
## luego plotea todas las subpoblaciones evolucionadas
## Queda a modo de debug, por si es nesario
# ddo <-
# function(ou){
#
# if(missing(ou)) ou=tryRK()
#
# plot (ou[,"S"], type = "l", col='blue', xlab = "Dias", ylab = "Poblacion")
# lines(ou[,"E"], col = "magenta", lty = "dotted", lwd = 2)
# lines(ou[,"I"], col = "red", lty = "dotted", lwd = 2)
# lines(ou[,"R"], col = "green", lty = "dotted", lwd = 2)
# lines(ou[,"M"], col = "black", lty = "dotted", lwd = 2)
# legend('right', c('Sanos','Contagiados','Infectados','Recuperados','Muertos'),
#    col=c('blue','magenta','red','green','black'),
#    pch=c(21,21,21,21,21))
# }
##################################################################

###
## Funcion que llama a la intregracion RK4, del modelo de ecuaciones diferenciales para 3 grupos poblacionales
## luego plotea todas las subpoblaciones evolucionadas
##

ddo2 <-
  function() {
    a <- tryRK() # <- integro

    # sumo los resultados para los 3 subgrupos
    ST <- a[, "S0"] + a[, "S1"] + a[, "S2"]
    ET <- a[, "E0"] + a[, "E1"] + a[, "E2"]
    IT <- a[, "I0"] + a[, "I1"] + a[, "I2"]
    RT <- a[, "R0"] + a[, "R1"] + a[, "R2"]
    MT <- a[, "M0"] + a[, "M1"] + a[, "M2"]

    # ploteo los totales
    plot(ST, type = "l", col = "blue", xlab = "Dias", ylab = "Poblacion", ylim = c(0, 3.5E6))
    # plot (ST, type = "l", col='blue', xlab = "Dias", ylab = "Poblacion",ylim=c(0,10000),xlim=c(0,60))
    lines(ET, col = "magenta", lty = "dotted", lwd = 2)
    lines(IT, col = "red", lty = "dotted", lwd = 2)
    lines(RT, col = "green", lty = "dotted", lwd = 2)
    lines(MT, col = "black", lty = "dotted", lwd = 2)
    legend("right", c("Sanos", "Contagiados", "Infectados", "Recuperados", "Muertos"),
      col = c("blue", "magenta", "red", "green", "black"),
      pch = c(21, 21, 21, 21, 21)
    )
  }


## Funcion que define las ecuasiones diferenciales que quiero evolucionar
## t = entra la secuencia de tiempo en dias, en este caso no es usado, pero uno puede cambiar
##     cosas en funcion del tiempo en que se llama a la funcion.
##     Idea original: cambiar la porpagacion o las muertes segun las camas o aislamiento social en el tiempo
## x = vector de poblaciones,ver notacion mia:
## 	S: Sanos. es decir poblacion original. a t-> inf esto tiene que ser cero
## 	E: Contagiados. poblacion que tiene el virus asintomaticamente
## 	I: Infectados. Gente que se le detecto el virus
## 	R: recuperados. Gente que sobrevino a la epidemia
## 	M: Muertos. gente que NO sobrevino a la epidemia
## params = parametros generales del modelo:
## 	bet: 	beta en wikipedia, multiplicado por los contagiados es el rate de contagio, unidades 1/dias. Paso de S -> E
##      a: 	1/tiempo de incuvacion, esto es propio de la sepa de virus (y del cuerpo de la persona).
## 		yo lo fije 1.0/5.2 dias (Whutan), como propio del virus. Paso de E -> I
## 	gam:    1/tiempo de reuperacion (con esto tmb calculamos el R0), esto si varia en el tiempo y grupo poblacional.
## 		En mi "modelo", permito que tanto los contagiados como los infectados se recuperen. Pasa de E -> R y I-> R
## 	mu:	1/Tiempo de muertos. Paso de I-> M
## 	N: 	Numero total de la poblacion. 3.5E6 es como la pobacion de la provincia de cordoba
####################
SEIRMmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    # aca uno puede modificar los parametros segun el tiempo. ie: implementacion de politicas, cierre de fronteras, ailamiento, o cantidad de camas.
    #    if(t > 20)	bet=bet/10.0
    #    if(t > 40) mu=2*mu

    ####### Por simplicidad entran parametros generales y modifico los parametros para el vector de entrada X=5*3
    ### Nro de indiviuos por rangos etarios
    N0 <- 0.34 * N
    N1 <- 0.48 * N
    N2 <- 0.18 * N
    ##### grupo [0-19], contagia igual que el resto, se reupera más rápido y no muere
    bet0 <- bet
    gam0 <- gam * 2.
    mu0 <- mu * 0.0
    ##### grupo [20-49], contagia el doble (actividades todo el dia), se reupera la mitad de lo standard y muere como todos
    bet1 <- bet * 2.
    gam1 <- gam / 2
    mu1 <- mu
    ##### grupo [50-inf], se contagia la mitad (asumo que estan guardados), se reupera la mitad de lo standard y el doble (o triple)
    bet2 <- bet / 2.
    gam2 <- gam / 2.
    mu2 <- mu * 2
    #################################################
    ## Sistema de ecuaciones Ojo:
    ## Pienso que>
    ## 	la poblacion 0 y 1 estan en contacto y no con los viejos
    ## 	la poblacion 0 no muere
    ## 	la poblacion 1 esta en contacto con todos
    ## 	la poblacion 2 esta en contacto solo con la 1
    ## detalles: de los contagiados asintomaticos que pasan a recuperarse solos
    ## 	los contagiados, contagian 1/4 del rate que los sanos
    ## 	los contagiados viejos no se sanan, (si o si se infectan)


    S0dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0 - bet0 / 4 * E0 / N0
    E0dot <- (bet0 * I0 / N0 + bet1 * I1 / N1) / 2 * S0 - (a + gam0) * E0 + bet0 / 4 * E0 / N0
    I0dot <- a * E0 - (gam0 + mu0) * I0
    R0dot <- gam0 * I0 + gam0 * E0
    M0dot <- mu0 * I0

    S1dot <- -(bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1 - bet1 / 4 * E1 / N1
    E1dot <- (bet0 * I0 / N0 + bet1 * I1 / N1 + bet2 * I2 / N2) / 3 * S1 - a * E1 - gam1 * E1 + bet1 / 4 * E1 / N1
    I1dot <- a * E1 - (gam1 + mu1) * I1
    R1dot <- gam1 * I1 + gam1 * E1
    M1dot <- mu1 * I1

    S2dot <- -(bet1 * I1 / N1 + bet2 * I2 / N2) / 2 * S2 - bet2 / 4 * E2 / N2
    E2dot <- (bet1 * I1 / N1 + bet2 * I2 / N2) / 2 * S2 - a * E2 + bet2 / 4 * E2 / N2
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

tryRK <-
  function(parms, xstart, times) {
    ## The parameters
    library("deSolve")
    if (missing(parms)) parms <- c(bet = 2.5, a = 1.0 / 5.2, gam = 0.2, mu = 0.05, N = 3.5E6)
    N <- 3.5E6
    ## vector con las condiciones inciales
    if (missing(xstart)) {
      xstart <- c(
        S0 = 0.3 * N - 11, E0 = 10, I0 = 1, R0 = 0, M0 = 0,
        S1 = 0.5 * N - 11, E1 = 10, I1 = 1, R1 = 0, M1 = 0,
        S2 = 0.2 * N - 11, E2 = 10, I2 = 1, R2 = 0, M2 = 0
      )
    }

    # secuencia del tiempo de integracion
    if (missing(times)) times <- seq(0, 90, length = 91)

    out <- rk(xstart, times, SEIRMmod, parms, hini = 1, method = "rk4")
    # out es un tipo deSolve

    return(out)
  }
