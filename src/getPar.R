getPar <-
  function(download) {
    if (missing(download)) download <- FALSE

    source("src/load_data.R")

    a <- load_data(download = FALSE)

    iday <- 15
    fday <- end(a[, 1])[1]
    DelD <- fday - iday + 1
    Muer <- a[, 1][iday:fday]
    Recu <- a[, 3][iday:fday]
    Acti <- a[, 4][iday:fday]

    MMean <- (Acti[DelD] + Acti[1]) / 2
    Bet <- (Acti[DelD] - Acti[1]) / DelD # /Acti[DelD+1]
    Gam <- (Recu[DelD] - Recu[1]) / DelD / MMean
    Mu <- (Muer[DelD] - Muer[1]) / DelD / MMean
    print(c(Bet, Gam, Mu))
    aa <- 0
    for (i in c(2:DelD))
    {
      if (i < 4) {
        a <- (Acti[i] - Acti[i - 1])
      } else {
        a <- (Acti[i] - Acti[i - 3]) / 3
      }
      aa <- aa + a
      print(c("Ver", i, "sin", a, "tot", aa))
    }
    print(c("Final", aa / DelD))

    return(data.frame(Bet = Bet, Mu = Mu, Gam = Gam))
  }
