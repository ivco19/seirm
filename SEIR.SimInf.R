trySIR <-
function()
{
library('SimInf')
n <- 1E4
u0 <- data.frame(S = rep(n-1, n), I = rep(1, n), R = rep(0, n))
tspan <- seq(from = 1, to = 90, by = 1)
model <- SIR(u0 = u0, tspan = tspan, beta = 0.5, gamma = 0.077)
set.seed(345564)
set_num_threads(4)
res <- run(model = model)
plot(res)
return(res)
}
trySEIRM <-
function()
{
library('SimInf')
n <- 5.0E3
bet =2.5 #
gam =1.4 #0.75
mu =0.2
a0=1./14.2
print(c("R0:",bet*a0/(mu*(a0+mu)*(gam+mu))))
transitions <- c(
	"S -> bet*S*E/(S+E+I+R+M0+M1) -> E", 
	"S -> bet*a0*S*I/(S+E+I+R+M0+M1) -> I", 
	"E -> a0*E -> I", 
	"E -> gam*E -> R", 
	"I -> gam*I -> R",
	"I -> 0.05*I -> M0",
	"M0 -> mu*M0 -> M1",
	"M0 -> gam*M0 -> R"
	)
compartments <- c("S", "E","I", "R", "M0","M1")

u0 <- data.frame(S = rep(n-11, n), E=rep(10,n), I = rep(1, n), R = rep(0, n),M0=rep(0,n),M1=rep(0,n))
tspan <- seq(from = 1, to = 90, by = 1)
#model <- SIR(u0 = u0, tspan = tspan, beta = 0.5, gamma = 0.077)
model <- mparse(transitions=transitions, compartments=compartments, gdata=c(bet= bet, gam = gam, a0=a0,mu=mu), u0 = u0, tspan = tspan) 
set.seed(345564)
set_num_threads(4)
res <- run(model = model)
plot(res)
return(res)
}

