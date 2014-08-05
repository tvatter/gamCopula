#########
#########

# define 5-dimensional R-vine tree structure matrix
Matrix = c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
Matrix = matrix(Matrix,5,5)

# define R-vine pair-copula family matrix
family = c(0,1,3,4,4,0,0,3,4,1,0,0,0,4,1,0,0,0,0,3,0,0,0,0,0)
family = matrix(family,5,5)

# define R-vine pair-copula parameter matrix
par = c(0,0.2,0.9,1.5,3.9,0,0,1.1,1.6,0.9,0,0,0,1.9,0.5,
        0,0,0,0,4.8,0,0,0,0,0)
par = matrix(par,5,5)

# define second R-vine pair-copula parameter matrix
par2 = matrix(0,5,5)

# define RVineMatrix object
RVM = RVineMatrix(Matrix=Matrix,family=family,par=par,par2=par2,
                  names=c("V1","V2","V3","V4","V5"))

Nn <- rep(c(1,2,5), 4)*10^c(rep(1,3),rep(2,3), rep(3,3), rep(4,3))
t1 <- t2 <- rep(0, length(Nn))

for(i in 1:length(Nn)){
	N <- Nn[i]
	U <- matrix(runif(N*5), ncol = 5)
	t1[i] <- system.time(Sim1 <- RVineSim(N, RVM, U))[3]
	t2[i] <- system.time(Sim2 <- gamVineSim(N, RVM, U))[3]
	all.equal(Sim1,Sim2)
}

par(mfrow = c(1,2))
plot(Nn, t1, ylim = range(t1,t2), log = "xy")
points(Nn, t2, pch = 2)
plot(Nn, t2/t1, log = "xy")



