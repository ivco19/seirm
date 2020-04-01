doall <-
   function() {
      s1=find_parms('SEIR',Nbet=14)
      s2=find_parms('SEIRmod',Nbet=14)
      s3=find_parms('SEIR_3',Nbet=14)
      s4=find_parms('SEIR_3mod',Nbet=14)
      return(list(seir=s1,seirm=s2,seir3=s3,seir3m=s4))
   }


find_parms <-
  function(model,Nbet) {
    library("deSolve")
    # source('getPar.R')
    source("src/models.R")
    source("src/plot.R")
    source("src/load_data.R")

    if (missing(model)) model <- "SEIR"
    if (missing(Nbet)) Nbet <- 10
    model_in <- switch(model, "SEIRmod" = SEIRmod, "SEIR" = SEIR, "SEIR_3mod" = SEIR_3mod, "SEIR_3" = SEIR_3)

    Real <- load_data(download = FALSE)
    iday <- 24
    fday <- end(Real[, 1])[1]
    RInf <- Real[, 4][iday:fday]
    RMue <- Real[, 1][iday:fday]
    RRec <- Real[, 3][iday:fday]
    Dday <- fday - iday + 1
    Rtime <- seq(iday, fday, len = Dday)

    N <- 3.5E6
    times <- seq(0, 60, length = 61)

    XXsig <- 1E9
    count <- 1
    for (BET in lseq(0.001, 2.5, Nbet)) {
      for (AAA in seq(1.0, 30, len = Nbet)) {
        for (GAM in seq(3.0, 30, len = Nbet)) {
          for (MUU in lseq(0.001, 1., Nbet)) {
            for (Eini in lseq(1, 1000, Nbet-1)) {
              for (Iini in lseq(1, 100, Nbet-1))
              {
                Eini <- as.integer(Eini)
                Iini <- as.integer(Iini)
                parms <- c(bet = BET, a = 1.0 / AAA, gam = 1.0 / GAM, mu = MUU, N = N)
    		if (model == "SEIR" || model == "SEIRmod") {
                  xstart <- c(S0 = N - Eini - Iini - 2, E0 = Eini, I0 = Iini, R0 = 1, M0 = 1)
                  out <- rk(xstart, times, model_in, parms, hini = 1, method = "rk4")
                  IT <- out[, "I0"]
                  RT <- out[, "R0"]
                  MT <- out[, "M0"]
    		} else if (model == "SEIR_3" || model == "SEIR_3mod") {
                  # porcent=c(0.304,0.430,0.266)
                  xstart <- c(
                    S0 = 0.304 * (N - Eini - Iini - 2), E0 = 0.304 * Eini, I0 = 0.304 * Iini, R0 = 1, M0 = 1,
                    S1 = 0.430 * (N - Eini - Iini - 2), E1 = 0.430 * Eini, I1 = 0.430 * Iini, R1 = 1, M1 = 1,
                    S2 = 0.266 * (N - Eini - Iini - 2), E2 = 0.266 * Eini, I2 = 0.266 * Iini, R2 = 1, M2 = 1
                  )
                  out <- rk(xstart, times, model_in, parms, hini = 1, method = "rk4")
                  IT <- out[, "I0"] + out[, "I1"] + out[, "I2"]
                  RT <- out[, "R0"] + out[, "R1"] + out[, "R2"]
                  MT <- out[, "M0"] + out[, "M1"] + out[, "M2"]
                }
                sig_inf <- 0
                sig_rec <- 0
                sig_mue <- 0
                mminf <- 0
                mmrec <- 0
                mmmue <- 0
                for (i in Rtime)
                {
                  # 	print(c(i,sig_inf, IT[i], RInf[i-iday+1]))
                  sig_inf <- sig_inf + (IT[i] - RInf[i - iday + 1])^2
                  sig_rec <- sig_rec + (RT[i] - RRec[i - iday + 1])^2
                  sig_mue <- sig_mue + (MT[i] - RMue[i - iday + 1])^2
                  mminf <- mminf + RInf[i - iday + 1]^2
                  mmrec <- mminf + RRec[i - iday + 1]^2
                  mmmue <- mminf + RMue[i - iday + 1]^2
                }
                mminf <- (mminf / Dday)
                mmrec <- (mmrec / Dday)
                mmmue <- (mmmue / Dday)
                sig_inf <- sqrt(sig_inf / mminf) ## por ahora le creemos a la estadistica de infectados
                sig_rec <- sqrt(sig_rec / mmrec)
                sig_mue <- sqrt(sig_mue / mmmue)

                XsigT <- sig_mue + sig_rec + sig_inf
                # print(c(XsigT,XXsig))
                if (!is.na(XsigT)) {
                  if (XXsig > XsigT) {
                    XXsig <- XsigT
                    Par_pref <- parms
                    Ini_pref <- xstart
                  }
                }
              }
            }
          }
        }
      }
      print("///////////////////////////")
      print(sprintf("Sigma: %8.4f / Porcent: %8.4f", XXsig, count / Nbet * 100))
      count <- count + 1
      if (!missing(Par_pref)) {
        print(c(Par_pref))
        print(c(Ini_pref))
      }
    }
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    if (missing(Par_pref)) {
      print("No se pudieron ajustar parametros")
      return(0)
    } else {
      print(c(XXsig, Par_pref, Ini_pref))

      ddo(model=model, par = Par_pref, xstart = Ini_pref)

      return(list(par = Par_pref, xstart = Ini_pref,Sigma=XXsig))
    }
  }
