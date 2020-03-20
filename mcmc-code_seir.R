# -------------------------------------------------------------- #
#                MONTE CARLO MARKOV CHAIN PROJECT                #
# -------------------------------------------------------------- #


# Libraries ------------------------------
library(ggplot2)
library(reshape2)
library(pomp)
library(coda)
library(invgamma)
library(MASS)
library(mvtnorm)
library(foreach)
library(parallel)
library(doParallel)
library(gridExtra)
library(tmvtnorm)
library(matrixcalc)
library(plyr)
library(magrittr)

citation("pomp")
stopifnot(packageVersion("pomp")>="1.12")

# /!\ INSTALL RTools here : https://cran.r-project.org/bin/windows/Rtools/
      # Click on edith system path

# Check if your Rtools installation worked, ---------------------
# if yes then this should print hellow world!
cat("#include <R.h>
void hello (void) {
    Rprintf(\"hello world!\\n\");
    }",file="hello.c")
system("R CMD SHLIB hello.c")
dyn.load(paste0("hello",.Platform$dynlib.ext))
.C("hello",PACKAGE="hello")

# set path 

# 1. SIMULATION OF SEIR MODEL ------------------------------------

# /!\ READ ME ----------------------------------------------------

# The following code will make use of the pomp package (partially
# observed Markov process package). We build on this package to 
# develop the SEIR model with a contact rate driven by a Brownian
# motion. First of all, we simulate from this model data with fixed
# parameters. Then, we take those simulations and apply the particle
# filter and the pmcmc algorithm in order to estimate the parameters
# driving the model. 

# We observe Y, which is a noisy observation of dR, that is to say
# the people that enter the Removed compartment at time t. We have
# only a noisy measure of this latent variable. Moreover, we need to
# build an accumulator variable H that will record dR (built from
# the latent variable R).

# Note that in our states, we also put x_t, the Brownian motion
# driving beta_t. 

# ----------------------------------------------------------------

# There are four compartments:
#' @state S : susceptible
#' @state E : infective but not infectious
#' @state I : infectious
#' @state R : removed people 
#' @state H : accumulator variable, dR
#' @state x : Brownian motion driving beta


#' @param k parameter driving differential equations system
#' @param gamma parameter driving differential equations system
#' @param N size of population
#' @param S_0 initial susceptible population
#' @param E_0 initial infected population
#' @param I_0 initial infectious population
#' @param H_0 initial dR, will be zero
#' @param sigma standard deviation of Brownian motion
#' @param x_0 initial state of Brownian motion
#' @param tau standard deviation of normal in measurement equation

statenames <- c("S","E","I","R","H","x")
paramnames <- c("k","gamma","N","S_0","E_0","I_0","H_0","sigma","x_0")

# Differential equation specification:
closed.sir.ode <- Csnippet("
                           double dx = rnorm(0,sigma*dt);
                           double Beta = exp(x);
                           x += dx;
                           S -= Beta*S*I/N;
                           E += Beta*S*I/N-k*E;
                           I += k*E - gamma*I;
                           R += gamma*I;
                           H += gamma*I; 
                           ") 

# Initialization:
init1 <- Csnippet("
                  S = S_0;
                  E = E_0;
                  I = I_0;
                  R = N - S_0 - E_0 - I_0;
                  H = H_0;
                  x = x_0;
                  ")

# number of points for Euler approximation
m <- 6 

# Creation of partially observed markov process object ----------------
pomp(data= data.frame(time=1:50,data=NA),
              # 50 data points for epidemic observations, for now empty
     
     times = "time",t0=0,
     rprocess=euler.sim(closed.sir.ode,delta.t=1/(m+1)), 
              # euler.sim because of continuous time process, 
              # with euler approx, m
     
     initializer=init1,
     statenames=statenames,
     paramnames=paramnames) -> closed.sir

# Declaration of accumulator variable
pomp(closed.sir,zeronames="H") -> closed.sir 

# Declaration of density and generator of Y conditional on states
dmeas <- Csnippet("lik = dnorm(data,H,tau,give_log);") 
              # log-likelihood is given here
rmeas <- Csnippet("data = rnorm(H,tau);") 
              # process of observables conditional on states

              # vC)rifier dans C que le parametre est le sd et non la 
              # var


closed.sir <- pomp(closed.sir,rmeasure=rmeas,dmeasure=dmeas,
                   statenames="H",paramnames="tau")

# Simulations from the model ------------------------------------------------
set.seed(8)
# Fixed parameters
params1 <- c(k = 1/40, gamma=1/15, N = 800, S_0 = 730, E_0 = 50, 
             I_0=20, x_0 = -1, H_0 = 0, sigma = 0.1, tau = 0.1)
    # in the paper, they take tau = 0.1
sims1 <- simulate(closed.sir, params = params1,
                 nsim = 1, as.data.frame = TRUE, include.data = F)

# Plot of simulations -------------------------------------------------------
ggplot(data=sims1)+geom_line(mapping=aes(x = time,y = S, color = "S"))+
  geom_line(mapping=aes(x = time,y = E, color = "E"))+
  geom_line(mapping=aes(x = time,y = I, color = "I"))+
  geom_line(mapping=aes(x = time,y = R, color = "R"))
  
ggplot(data = sims1)+geom_line(mapping = aes(x=time,y=H, color = "H"))


# Simulation noisy observed data
data <- matrix(nrow = 50,ncol = 1)
H <- sims1$H
tau <- params1['tau']
for (i in 1:50){
  data[i] <- rnorm(1,mean = H[i],sd = tau)
}

sims1$data <- data 

# Noisy observations and latent variable:
ggplot(data=sims1)+geom_line(mapping=aes(x=time,y=data, color = "Y"))+
  geom_line(mapping=aes(x=time,y=H, color = "H"))
save(sims1, file = paste0(getwd(),"/sims1.RData"))

# PARTICLE FILTER ALGORITHM -------------------------------------------------

# We start by estimating only four parameters, sigma, tau, k and gamma.
# This leads to change a it the model and the initialization, we do as if
# we know the initialization states and the initialization of process xt
# as well as the size of the population.

dataset <- sims1[,1:2]
closed.sir.ode <- Csnippet("
                           double N = 800;
                           double dx = rnorm(0,sigma*dt);
                           double Beta = exp(x);
                           x += dx;
                           S -= Beta*S*I/N;
                           E += Beta*S*I/N-k*E;
                           I += k*E - gamma*I;
                           R += gamma*I;
                           H += gamma*I; 
                           ") # H accumulator var

# Initialisation des etats
init1 <- Csnippet("
                  S = 730;
                  E = 50;
                  I = 20;
                  R = 0;
                  H = 0;
                  x = -1;
                  ")
#specification of parameter transformation
fromEst <- Csnippet("
                    Tk = exp(k);
                    Tgamma = exp(gamma);
                    ")

toEst <- Csnippet("
                  Tk = log(k);
                  Tgamma = log(gamma);
                  ")

paramnames <- c("k","gamma", "sigma", "tau")
params2 <- c(k = 1/40,gamma=1/15, sigma = 0.1, tau = 0.1)
dmeas <- Csnippet("
                  lik = dnorm(data,H,tau,give_log);
                  ") # log-likelihood is given here
rmeas <- Csnippet("
                  data = rnorm(H,tau);
                  ") # process of observables conditional on states

pomp(
  data=dataset,
  times="time",t0=0,
  rprocess=euler.sim(closed.sir.ode,delta.t=1/(m+1)), # euler.sim because of continuous time process, with euler approx, m
  initializer=init1,
  zeronames="H",
  rmeasure=rmeas,dmeasure=dmeas,
  fromEstimationScale=fromEst,toEstimationScale=toEst,
  statenames=statenames,
  paramnames=paramnames) -> seir

set.seed(6)
# 10 simulations a partir des donnees
sims <- simulate(seir,params=params2,nsim=10,as.data.frame=TRUE, include.data = T)

# plot of noisy data:
ggplot(data=sims,mapping=aes(x=time,y=data,group=sim))+
  geom_line()

# Basic particle filter -----------------------------------------------
Nparticle <- 5000
pf <- pfilter(seir,params=params2,Np=Nparticle)
ll <- pf@loglik
ll
plot(pf)
# The plot shows the data, along with effective sample size of the particle filter (ess)
# and the log likelihood of each obs conditional on the preceding ones (cond.logLik)

# Effect of number of particles on estimated likelihood by pf -------
# parameters
#Nparticles <- c(1000,2000,5000)
Nparticles <- c(1000,2000,5000,10000,20000,50000)
Nrep <- 10

# particle filter replications /!\ Takes a minute to complete !

pf <- sapply(Nparticles, function(i) replicate(Nrep,
                pfilter(seir,params = params2, Np = i)))
ll <- sapply(pf,logLik)
ll <- data.frame(ll)
ll <- cbind(ll = ll, Particles = sort(rep(Nparticles,Nrep)))
ll$Particles <- as.factor(ll$Particles)
# we parallelize this calculation
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, c("Nrep","Nparticles","seir","params2"))
set.seed(998468235L,kind="L'Ecuyer")
foreach(i=1:6,
        .packages="pomp",
        .options.multicore=list(set.seed=TRUE)
) %dopar% {
  replicate(Nrep,pfilter(seir,Np = Nparticles[i],params=params2))
} -> pfs
stopCluster(cl)
pfs <- unlist(pfs)
ll2 <- sapply(pfs,logLik)
ll2 <- data.frame(ll2)
ll2 <- cbind(ll2 = ll2, Particles = sort(rep(Nparticles,Nrep)))
ll2$Particles <- as.factor(ll2$Particles)


# plot ------------------------------------------------------------------
ggplot(data = ll)+ 
  geom_boxplot(mapping = aes(x = Particles,y = ll, color = Particles))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Particle filt. est. log likelihood func. of Npart",
          caption = "MCMC project", 
          x = "Number of Particles",
          y = 'Estim. log lik')

# same thing but the computation have been parallelized
ggplot(data = ll2)+ 
  geom_boxplot(mapping = aes(x = Particles,y = ll2, color = Particles))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Particle filt. est. log likelihood func. of Npart",
          caption = "MCMC project", 
          x = "Number of Particles",
          y = 'Estim. log lik')

save(sims, file = "sims1.RData")
# Estimating the parameters ---------------------------------------------
# PMCMC -----------------------------------------------------------------
NItermcmc <- 2000
Nparticle <- 1000
# Prior information on parameters : 
seir <- pomp(seir,dprior=Csnippet("
                        lik = dnorm(k,1/40,0.1,1) + dnorm(gamma,1/15,0.1,1)
                          + dunif(sigma,0,1,1) + dunif(tau,0,1,1);
                        lik = (give_log) ? lik : exp(lik);"),
             paramnames=c("k","gamma","sigma", "tau"))
# The last 1 in each density is simply to indicate that we want the      log-likelihood

set.seed(73)
# Initial values of pmcmc
k0 <- rnorm(1,1/40,0.001)
gamma0 <- rnorm(1,1/15,0.001)
sigma0 <- runif(1,0.05,0.15)
tau0 <- runif(1,0.05,0.15)
Ntat <- 1000
t1 <- Sys.time()
pmcmc(seir,Nmcmc = Ntat,
      start = c(k = k0, gamma=gamma0, sigma = sigma0, tau = tau0),
      # Adaptative MCMC Proposal 
      proposal =  mvn.rw.adaptive(rw.sd=
                                    c(k = 10^(-5), gamma=10^(-5), sigma = 10^(-5), tau = 10^(-5)),scale.start=300,
                                  target = 0.2),
      # rw.sd is a named numeric vector; random-walk SDs for a 
      # multivariate normal randomwalk proposal with diagonal 
      # variance-covariance matrix.
      Np = Nparticle)-> chain



continue(chain,Nmcmc=NItermcmc ,proposal = mvn.rw(covmat(chain))) -> chain
t2 <- Sys.time()
pmcmc.t <- t2-t1
# During the first 1000 iterations, we allow no correlation between
# paramaters (variance covariance matrix is diagonal)
# Afterwards, we continue the chain by allowing
# correlation with a correlation estimated on the first 1000 iterations.

print(pmcmc.t)
trace1 <- window(conv.rec(chain,c("k","gamma")),start=1000)
trace2 <- window(conv.rec(chain,c("sigma","tau")),start=1500)
rejectionRate(trace1)
effectiveSize(trace1)
autocorr.diag(trace1)
summary(trace1)
print(paste0("True k0 is: ", 1/40, " and true gamma0 is: ", 1/15))
plot(trace1)
rejectionRate(trace2)
effectiveSize(trace2)
autocorr.diag(trace2)
summary(trace2)
plot(trace2)
# Effect of the number of particles on the path -----------------------
Npart <- c(500,1000,10000)
set.seed(73)
# Initial values of pmcmc
k0 <- rnorm(1,1/40,0.001)
gamma0 <- rnorm(1,1/15,0.001)
sigma0 <- runif(1,0.05,0.15)
tau0 <- runif(1,0.05,0.15)
Ntat <- 1000
Nitermcmc <- 2000

chain <- lapply(Npart, function(i) pmcmc(seir,Nmcmc = Ntat,
                                         start = c(k = k0, gamma=gamma0, sigma = sigma0, tau = tau0),
                                         # Adaptative MCMC Proposal 
                                         proposal =  mvn.rw.adaptive(rw.sd=
                                                                       c(k = 10^(-5), gamma=10^(-5), sigma = 10^(-5), tau = 10^(-5)),scale.start=300,
                                                                     target = 0.2),
                                         # rw.sd is a named numeric vector; random-walk SDs for a 
                                         # multivariate normal randomwalk proposal with diagonal 
                                         # variance-covariance matrix.
                                         Np = i))

chain.cont <- lapply(1:3, function(i) continue(chain[[i]],Nmcmc=NItermcmc ,proposal = mvn.rw(covmat(chain[[i]]))))
for (i in 1:3){
  assign(paste0("trace1.",i),window(conv.rec(chain.cont[[i]],c("k","gamma")),start=1000))
  assign(paste0("trace2.",i),window(conv.rec(chain.cont[[i]],c("sigma","tau")),start=1000))
}

summary(trace1.1)
plot(trace1.1) # 500 part
summary(trace1.2)
plot(trace1.2) # 1000 part
summary(trace1.3)
plot(trace1.3) # 10 000 part

summary(trace2.1)
plot(trace2.1) # 500 part
summary(trace2.2)
plot(trace2.2) # 1000 part
summary(trace2.3)
plot(trace2.3) # 10 000 part


# PROPER PARTICLE FLTER AND PMCMC NOT BASED ON POMP
#-----------------------------------------------------------------------------------------
sims1 <- load("sims1.RData")

# PARTICLE FILTER ALGORITHM -------------------------------------------------

# We start by estimating only two parameters, k and gamma. We code our
# own particle filter and pmcmc methods without pomp already written methods.


given.data <- sims1[,1:2]
# Initialization 
init1 <- Csnippet("
                  S = S_0;
                  E = E_0;
                  I = I_0;
                  R = R_0;
                  ")
closed.sir.ode <- Csnippet("
                           DS = -Beta*S*I/N;
                           DE = Beta*S*I/N-k*E;
                           DI = k*E - gamma*I;
                           DR = gamma*I;
                           ")

paramnames <- c("k","gamma","Beta" ,"S_0", "E_0", "I_0", "R_0","N")
statenames <- c("S","E","I","R")
tau <- 0.1# Resolution of ODE   



for (i in 1:49){
  assign(paste0("seir",i),pomp(
    data= data.frame(time = 1:2,data = given.data[i:(i+1),2]),
    times="time",t0=0,
    initializer=init1,
    skeleton=vectorfield(closed.sir.ode),
    statenames=statenames,
    paramnames=paramnames))
}
model <- list(seir1,seir2,seir3,seir4,seir5,seir6,seir7,seir8,
              seir9,seir10,seir11,seir12,seir13,seir14,seir15,
              seir16,seir17,seir18,seir19,seir20,seir21,seir22,
              seir23,seir24,seir25,seir26,seir27,seir28,seir29,
              seir30,seir31,seir32,seir33,seir34,seir35,seir36,
              seir37,seir38,seir39,seir40,seir41,seir42,seir43,
              seir44,seir45,seir46,seir47,seir48,seir49)

sigma <- 0.1
m <- 6

# Particle Filter algorithm ---------------------------------------------
particle_filter <- function(Npart, T_, k, gamma, m = 6, sigma = 0.1, tau = 0.1){
  # Initialization
  S0 <- array(730, Npart)
  E0 <- array(50, Npart)
  I0 <- array(20, Npart)
  R0 <- array(0, Npart)
  H0 <- array(0, Npart)
  V <- data.frame(S0, E0, I0, R0, H0)
  W <- array(1/Npart, Npart)
  X0 <- -1
  Xdat <- array(dim = c((T_),Npart))
  alpha <- array(0, Npart)
  L <- 1
  
  for (ki in 1:(T_)){
    dx <- replicate(Npart,sum(replicate((m+1),rnorm(1,0,sd = sqrt(sigma*(1/(m+1)))))))
    if (ki ==1){
      Xdat[ki,] <- X0 + dx
    }else{
      Xdat[ki,] <- Xdat[(ki-1),]+dx
    }
    for (j in 1:Npart){
      params2 <- c(k=k*(m+1) ,gamma=gamma*(m+1), Beta = exp(Xdat[ki,j])*((m+1)),S_0 = V$S0[j] ,
                   E_0 = V$E0[j], I_0=V$I0[j], 
                   R_0 = V$R0[j], N = 800)
      traj <- trajectory(model[[ki]],params=params2,as.data.frame=TRUE)
      # varaccumulator
      H <- traj$R[1] - V$R0[j]
      V[j,] <- cbind(traj$S[1],traj$E[1],traj$I[1],traj$R[1],H)
      alpha[j] <- dnorm(given.data[(ki),2], H, sqrt(tau)) 
    }
    alpha_tot <- sum(alpha)
    if((alpha_tot!=0 )&((!is.na(alpha_tot)))){
      W <- alpha/alpha_tot
      L <- L*(alpha_tot/Npart)
    } else{
      W <- rep(1/Npart, Npart)
      # print("ERROR sum of weight is null")
    }
    resampl <- as.numeric(rmultinom(1, Npart, W))
    l <- 1
    q <- 1
    Vcopy <- V
    Xcopy <- Xdat
    while (l <= Npart){ 
      jn <- as.integer(resampl[l]) # number of children of l
      # q position in new vector
      # l position in old vector
      if (jn > 0){ # if the number of children is positive
        
        VV <- do.call("rbind", replicate(jn, Vcopy[l, ], simplify = FALSE ))
        V[q:(q+jn-1), ] <- VV 
        Xdat[ ,q:(q+jn-1)] <- matrix(rep(Xcopy[, l], jn), ncol = jn)
        l <- l+1 
        q <- q+jn
      }
      else{
        l <- l+1 
      }
    }
  }
  #select_part <- sample.int(Npart, size = 1)
  #X_res <- Xdat[, select_part]
  #return(L = L, Xres = Xres)
  return(L = L)
}

# parameters
Npart <- 500
T_ <- 49

# output ----------------------------------------------------
print("Test of the particle filter algorithm")
set.seed(34)
# We start in a neighborhood of k and gamma
k0 <- runif(1,0.021,0.027)
gamma0 <- runif(1,0.062,0.07)
t1  <- Sys.time()
pfilter_output <- particle_filter(Npart = Npart, T_ = T_,k = k0, gamma = gamma0) 
t2 <- Sys.time()
pfilter.time <- t2 - t1
print(pfilter.time)
pfilter_output

# Effect of number of particles on est. likelihood
Nparticles <- c(50,100,200,500,1000)
Nrep <- 50
T_ <- 15

set.seed(34)
t1 <- Sys.time()
pf <- lapply(Nparticles, function(i) data.frame(lik =  replicate(Nrep,
                                                                 particle_filter(Npart = i , T_ = T_,k = k0, gamma = gamma0)),
                                                Particles = i))
t2 <- Sys.time()
ll <- do.call(rbind,pf)
ll$Particles <- as.factor(ll$Particles)
#save(ll, file = "llv2.Rdata")

#load("llv2.Rdata")
# plot ------------------------------------------------------------------
ggplot(data = ll)+ 
  geom_boxplot(mapping = aes(x = Particles,y = lik, color = Particles))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Particle filt. est. likelihood func. of Npart",
          caption = "MCMC project", 
          x = "Number of Particles",
          y = 'Estim. lik')

# Pmcmc algorithm -------------------------------------------
#' @param Nit number of iterations of Metropolis Hastings algorithm
#' @param sd variance covariance matrix of the proposal Q distribution (normal random walk)
#' @param theta list of parameters to estimate
#' @param y data we observe
#' @param Npart number of particles in particle filter



PMCMC <- function(Nit, sd, Npart, T_, k, gamma, Start, End = (Nit+1)){
  theta <- c(k, gamma)
  Theta <- theta
  #X <- list()
  print("Initialization")
  pf <- particle_filter(Npart, T_, theta[1], theta[2])
  Like <- pf
  #x <- pf[2:length(pf)]
  r <- 0
  for (i in 1:Nit){
    
    print(paste0('Iteration of Metropolis Hastings: ',i))
    up_theta <- rtmvnorm(1,as.numeric(theta),  sd, lower = rep(0,2))
    print(paste0("up_theta: ",up_theta))
    pf <- particle_filter(Npart, T_, k = up_theta[1], gamma = up_theta[2])
    up_Like <- pf
    #up_x <- pf[2:length(pf)]
    
    up_Q <- dtmvnorm.marginal2(up_theta[1],up_theta[2],1,2, as.numeric(theta), sd,lower = rep(0,2))
    Q <- dtmvnorm.marginal2(theta[1],theta[2],1,2, as.numeric(up_theta), sd,lower = rep(0,2))
    up_prob <- (up_Like*Q)/(Like*up_Q)
    print(paste0("Like: ", Like, " Up_like : ", up_Like))
    up_prob <- min(1, up_prob)
    u <- runif(1, 0, 1)
    if ((u < up_prob)&(up_Like!=0)){
      theta <- up_theta
      # x <- up_x
      Like <- up_Like
    }else{
      print("rejet")
      r <- r+1
      
    }
    Theta <- rbind(Theta, theta)
    Theta <- matrix(Theta, ncol = 2)
    colnames(Theta) <- c("k","gamma")
    #X <- rbind(X, x)
    if ((i > Start)&(i<End)){
      v1 <- var(Theta[,1])
      v2 <- var(Theta[,2])
      sd <- diag(c(v1,v2))
    }
  }
  return(list(Theta = Theta,  accept =1 - r/Nit))
}


# parameters --------------------------------------------------------
Nit <- 1000
Npart <- 500
T_ <- 49
Start <- 10
set.seed(34) # then for second chain set.seed(92)
# We start in a neighborhood of k and gamma
k0 <- runif(1,0.021,0.027)
gamma0 <- runif(1,0.062,0.07)

# Output ---------------------------------------------------------------------
res <- PMCMC(Nit = Nit, sd = diag(x = c(10^(-5),10^(-5))) , 
             Npart = Npart, T_ = T_, k = k0, gamma = gamma0, Start = Start)

# We save it in order to use it later to estimate a variance covariance matrix
# for the two coefficients
# save(res, "Result1000It500Part49TPMCMCv1.RData")
res.seed34 <- res
# res.seed92 <- res

res.seed34$accept
res.seed92$accept

# plot -----------------------------------------------------------------------

theta.c <- data.frame(res.seed34$Theta,res.seed92$Theta,0:Nit)
colnames(theta.c)<- c("k1","gamma1","k2","gamma2","It")

# summary stat chains
mean(theta.c$k1)
mean(theta.c$k2)
sd(theta.c$k1)
sd(theta.c$k2)

mean(theta.c$gamma1)
mean(theta.c$gamma2)
sd(theta.c$gamma1)
sd(theta.c$gamma2)

burnin <- 50
theta.c <- theta.c[-(1:burnin),]



g1 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = k1, color = "k1"))+
  geom_line(mapping = aes(x = It,y = k2, color = "k2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot k ",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'k')

g2 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = k1, color = "k1"))+
  geom_histogram(mapping = aes(x = k2, color = "k2"))+
  geom_line(aes(x = k1, y = ..density.., colour = 'Empirical density k1'), stat = 'density',size=1) + 
  geom_line(aes(x = k2, y = ..density.., colour = 'Empirical density k2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram k",
          caption = "MCMC project", 
          x = "k")

grid.arrange(g1,g2, nrow = 2)

g3 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = gamma1, color ="gamma1"))+
  geom_line(mapping = aes(x = It,y = gamma2, color ="gamma2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot gamma ",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'gamma')

g4 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = gamma1, color ="gamma1"))+
  geom_histogram(mapping = aes(x = gamma2, color ="gamma2"))+
  geom_line(aes(x = gamma1, y = ..density.., colour = 'Empirical density gamma1'), stat = 'density',size=1) + 
  geom_line(aes(x = gamma2, y = ..density.., colour = 'Empirical density gamma2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram gamma ",
          caption = "MCMC project", 
          x = "gamma")
grid.arrange(g3,g4, nrow = 2)

# Evaluation of the variation in results ---------------------------------------------
# parameters
Nit <- 500
Npart <- 300
T_1 <- 5
T_2 <- 10
Start <- 10
End <- 100
set.seed(34)
# We start in a neighborhood of k and gamma
k0 <- runif(1,0.021,0.027)
gamma0 <- runif(1,0.062,0.07)
Nrep <- 2

# Output ---------------------------------------------------------------------
res.chains <- replicate(Nrep, PMCMC(Nit = Nit, sd = diag(x = c(10^(-5),10^(-5))) , 
                                    Npart = Npart, T_ = T_1, k = k0, gamma = gamma0, Start = Start, End = End),
                        simplify = FALSE)

# Output ---------------------------------------------------------------------
res.chains <- replicate(Nrep, PMCMC(Nit = Nit, sd = diag(x = c(10^(-5),10^(-5))) , 
                                    Npart = Npart, T_ = T_2, k = k0, gamma = gamma0, Start = Start, End = End),
                        simplify = FALSE)

# plot -----------------------------------------------------------------------
chains <- cbind(res.chains[[1]]$Theta,res.chains[[2]]$Theta)
theta.c <- data.frame(chains,0:Nit)
colnames(theta.c)<- c("k1","gamma1","k2","gamma2","It")
# summary stat chains
mean(theta.c$k1)
mean(theta.c$k2)
sd(theta.c$k1)
sd(theta.c$k2)

mean(theta.c$gamma1)
mean(theta.c$gamma2)
sd(theta.c$gamma1)
sd(theta.c$gamma2)


burnin <- 11
theta.c <- theta.c[-(1:burnin),]




g1 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = k1, color = "k1"))+
  geom_line(mapping = aes(x = It,y = k2, color = "k2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot k (T=5)",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'k')

g2 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = k1, color = "k1"))+
  geom_histogram(mapping = aes(x = k2, color = "k2"))+
  geom_line(aes(x = k1, y = ..density.., colour = 'Empirical density k1'), stat = 'density',size=1) + 
  geom_line(aes(x = k2, y = ..density.., colour = 'Empirical density k2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram k (T=5)",
          caption = "MCMC project", 
          x = "k")

grid.arrange(g1,g2, nrow = 2)

g3 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = gamma1, color ="gamma1"))+
  geom_line(mapping = aes(x = It,y = gamma2, color ="gamma2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot gamma (T=5)",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'gamma')

g4 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = gamma1, color ="gamma1"))+
  geom_histogram(mapping = aes(x = gamma2, color ="gamma2"))+
  geom_line(aes(x = gamma1, y = ..density.., colour = 'Empirical density gamma1'), stat = 'density',size=1) + 
  geom_line(aes(x = gamma2, y = ..density.., colour = 'Empirical density gamma2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram gamma (T=5)",
          caption = "MCMC project", 
          x = "gamma")
grid.arrange(g3,g4, nrow = 2)
res.chains2 <- res.chains

# plot -----------------------------------------------------------------------
chains <- cbind(res.chains2[[1]]$Theta,res.chains2[[2]]$Theta)
theta.c <- data.frame(chains,0:Nit)
colnames(theta.c)<- c("k1","gamma1","k2","gamma2","It")
# summary stat chains
mean(theta.c$k1)
mean(theta.c$k2)
sd(theta.c$k1)
sd(theta.c$k2)

mean(theta.c$gamma1)
mean(theta.c$gamma2)
sd(theta.c$gamma1)
sd(theta.c$gamma2)

burnin <- 11
theta.c <- theta.c[-(1:burnin),]



g1 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = k1, color = "k1"))+
  geom_line(mapping = aes(x = It,y = k2, color = "k2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot k (T=10)",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'k')

g2 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = k1, color = "k1"))+
  geom_histogram(mapping = aes(x = k2, color = "k2"))+
  geom_line(aes(x = k1, y = ..density.., colour = 'Empirical density k1'), stat = 'density',size=1) + 
  geom_line(aes(x = k2, y = ..density.., colour = 'Empirical density k2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram k (T=10)",
          caption = "MCMC project", 
          x = "k")

grid.arrange(g1,g2, nrow = 2)

g3 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = gamma1, color ="gamma1"))+
  geom_line(mapping = aes(x = It,y = gamma2, color ="gamma2"))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot gamma (T=10)",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'gamma')

g4 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = gamma1, color ="gamma1"))+
  geom_histogram(mapping = aes(x = gamma2, color ="gamma2"))+
  geom_line(aes(x = gamma1, y = ..density.., colour = 'Empirical density gamma1'), stat = 'density',size=1) + 
  geom_line(aes(x = gamma2, y = ..density.., colour = 'Empirical density gamma2'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram gamma (T=10)",
          caption = "MCMC project", 
          x = "gamma")
grid.arrange(g3,g4, nrow = 2)



# Refined PMCMC ------------------------------------------------------------------

PMCMC.adaptative <- function(Nit, sd, Npart, T_, k, gamma, Start){
  theta <- c(k, gamma)
  Theta <- theta
  #X <- list()
  print("Initialization")
  pf <- particle_filter(Npart, T_, theta[1], theta[2])
  Like <- pf
  #x <- pf[2:length(pf)]
  r <- 0
  for (i in 1:Nit){
    
    print(paste0('Iteration of Metropolis Hastings: ',i))
    up_theta <- rtmvnorm(1,as.numeric(theta),  sd, lower = rep(0,2))
    print(paste0("up_theta: ",up_theta))
    pf <- particle_filter(Npart, T_, k = up_theta[1], gamma = up_theta[2])
    up_Like <- pf
    #up_x <- pf[2:length(pf)]
    
    up_Q <- dtmvnorm.marginal2(up_theta[1],up_theta[2],1,2, as.numeric(theta), sd,lower = rep(0,2))
    Q <- dtmvnorm.marginal2(theta[1],theta[2],1,2, as.numeric(up_theta), sd,lower = rep(0,2))
    up_prob <- (up_Like*Q)/(Like*up_Q)
    print(paste0("Like: ", Like, " Up_like : ", up_Like))
    up_prob <- min(1, up_prob)
    u <- runif(1, 0, 1)
    if ((u < up_prob)&(up_Like!=0)){
      theta <- up_theta
      # x <- up_x
      Like <- up_Like
    }else{
      print("rejet")
      r <- r+1
      
    }
    Theta <- rbind(Theta, theta)
    Theta <- matrix(Theta, ncol = 2)
    colnames(Theta) <- c("k","gamma")
    #X <- rbind(X, x)
    if (i > Start){
      v1 <- var(Theta[,1])
      v2 <- var(Theta[,2])
      cov12 <-cov(Theta[,1],Theta[,2])
      Sigma.hat.new <- matrix(c(v1,cov12,cov12,v2),ncol = 2, nrow = 2)
      if (is.positive.semi.definite(Sigma.hat.new)){
        sd <- Sigma.hat.new
      }
    }
  }
  return(list(Theta = Theta,  accept =1 - r/Nit))
}

## parameters -----------------------------------------------------------

# Estimation of Sigma.hat on previous chain 
load("Result1000It500Part49TPMCMCv1.RData")
Nit <- 1000
chains <- cbind(res$Theta)
theta.c <- data.frame(chains,0:Nit)
colnames(theta.c)<- c("k","gamma","It")

v1 <-var(theta.c$k)
v2 <-var(theta.c$gamma)
cov12 <- cov(theta.c$k,theta.c$gamma)
Sigma.hat <- matrix(c(v1,cov12,cov12,v2),ncol = 2, nrow = 2)
# We check if it is positive semi definite in order to use it as a 
# variance covariance matrix
is.positive.semi.definite(Sigma.hat)


# parameters
Nit <- 500
Npart <- 300
T_ <- 25
Start <- 50
set.seed(34)
# We start in a neighborhood of k and gamma
k0 <- runif(1,0.021,0.027)
gamma0 <- runif(1,0.062,0.07)
sd <- Sigma.hat

# Output ---------------------------------------------------------------------
chain.adapt <- PMCMC.adaptative(Nit = Nit, sd = sd, 
                                Npart = Npart, T_ = T_, k = k0, gamma = gamma0, 
                                Start = Start)
chain.adapt$accept

# plot -----------------------------------------------------------------------

theta.c <- data.frame(chain.adapt$Theta,0:Nit)
colnames(theta.c)<- c("k","gamma","It")
# summary stat chains
mean(theta.c$k)
sd(theta.c$k)

mean(theta.c$gamma)
sd(theta.c$gamma)

burnin <- 20
theta.c <- theta.c[-(1:burnin),]



g1 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = k))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot k",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'k')

g2 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = k, color = "k"))+
  geom_line(aes(x = k, y = ..density.., colour = 'Empirical density k'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram k",
          caption = "MCMC project", 
          x = "k")

grid.arrange(g1,g2, nrow = 2)

g3 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = gamma))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot gamma",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'gamma')

g4 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = gamma, color ="gamma1"))+
  geom_line(aes(x = gamma, y = ..density.., colour = 'Empirical density gamma'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram gamma",
          caption = "MCMC project", 
          x = "gamma")

grid.arrange(g3,g4, nrow = 2)
# Starting values effect --------------------------------------------------------
# parameters
Nit <- 500
Npart <- 300
T_ <- 25
Start <- 50
set.seed(34)

# We start in a neighborhood of k and gamma (wider than before)
k0 <- runif(1,0.02,0.04)
gamma0 <- runif(1,0.03,0.1)
sd <- Sigma.hat

# Output ---------------------------------------------------------------------
chain.adapt <- PMCMC.adaptative(Nit = Nit, sd = sd, 
                                Npart = Npart, T_ = T_, k = k0, gamma = gamma0, 
                                Start = Start)
chain.adapt$accept
# plot -----------------------------------------------------------------------

theta.c <- data.frame(chain.adapt$Theta,0:Nit)
colnames(theta.c)<- c("k","gamma","It")

# summary stat chains
mean(theta.c$k)
sd(theta.c$k)

mean(theta.c$gamma)
sd(theta.c$gamma)

g1 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = k))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot k",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'k')

g2 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = k, color = "k"))+
  geom_line(aes(x = k, y = ..density.., colour = 'Empirical density k'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram k",
          caption = "MCMC project", 
          x = "k")

grid.arrange(g1,g2, nrow = 2)

g3 <- ggplot(data = theta.c)+
  geom_line(mapping = aes(x = It,y = gamma))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Traceplot gamma",
          caption = "MCMC project", 
          x = "Niter (metropolis hastings)",
          y = 'gamma')

g4 <- ggplot(data = theta.c)+
  geom_histogram(mapping = aes(x = gamma, color ="gamma1"))+
  geom_line(aes(x = gamma, y = ..density.., colour = 'Empirical density gamma'), stat = 'density',size=1) + 
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Histogram gamma",
          caption = "MCMC project", 
          x = "gamma")

grid.arrange(g3,g4, nrow = 2)

# APPLICATION TO CASE OF MEASLES IN LONDON ---------------------------------------
set.seed(594709947L)
# load data
datfile <- load("twentycities.rda")
measles %>% 
  mutate(year=as.integer(format(date,"%Y"))) %>%
  subset(town=="London" & year>=1950 & year<1964) %>%
  mutate(time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950) %>%
  subset(time>1950 & time<1964, select=c(time,cases)) -> dat
demog %>% subset(town=="London",select=-town) -> demogLondon
## plot -------------------------------------------------
# plot of cases
dat %>% ggplot(aes(x=time,y=cases))+geom_line()
# plot of the demography (population and birth)
demogLondon %>% melt(id="year") %>%
  ggplot(aes(x=year,y=value))+geom_point()+
  facet_wrap(~variable,ncol=1,scales="free_y")
## ----prep-covariates-----------------------------------------------------
demogLondon %>% 
  summarize(
    time=seq(from=min(year),to=max(year),by=1/12),
    pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
    birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
  ) -> covar

## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
                  double beta, br, seas, foi, dw, births;
                  double rate[6], trans[6];
                  
                  // cohort effect
                  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
                  br = cohort*birthrate/dt + (1-cohort)*birthrate;
                  else 
                  br = (1.0-cohort)*birthrate;
                  
                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;
                  
                  // transmission rate
                  beta = R0*(gamma+mu)*seas;
                  // expected force of infection
                  foi = beta*pow(I+iota,alpha)/pop;
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  rate[0] = foi*dw/dt;  // stochastic force of infection
                  rate[1] = mu;			    // natural S death
                  rate[2] = sigma;		  // rate of ending of latent stage
                  rate[3] = mu;			    // natural E death
                  rate[4] = gamma;		  // recovery
                  rate[5] = mu;			    // natural I death
                  
                  // Poisson births
                  births = rpois(br*dt);
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  
                  S += births   - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  ")

## ----initializer---------------------------------------------------------
initlz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")

## ----dmeasure------------------------------------------------------------
dmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")


## ----rmeasure------------------------------------------------------------
rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

## ----pomp-construction---------------------------------------------------
dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
       time="time",
       rprocess=euler.sim(rproc,delta.t=1/365.25),
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covar,
       tcovar="time",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> m1


## ----mles,include=FALSE--------------------------------------------------
read.csv(text="
         town,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
         Bedwellty,-1125.1,0.14,0.02,4,57.9,146,0.311,24.7,0.16,0.937,0.0396,0.351,0.951,0.0396,2.64e-05,2.45e-05,0.96,0.0611
         Birmingham,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
         Bradford,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
         Bristol,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
         Cardiff,-2364.9,0.73,0.02,4,39,143,0.602,34.4,0.223,0.996,0.141,0.267,0.27,0.0317,1.01e-05,9.21e-06,0.968,0.0539
         Consett,-1362.9,0.73,0.02,4,42.6,172,0.65,35.9,0.2,1.01,0.0731,0.31,0.406,0.0322,1.83e-05,1.97e-05,0.968,0.0712
         Dalton.in.Furness,-726.1,0.3,0.02,4,73.6,257,0.455,28.3,0.203,0.989,0.0386,0.421,0.818,0.0387,2.23e-05,2.36e-05,0.961,0.0779
         Halesworth,-318.6,0.51,0.02,4,49.6,210,0.754,33.1,0.381,0.948,0.00912,0.547,0.641,0.0526,1.99e-05,2.82e-05,0.947,0.0748
         Hastings,-1583.7,0.21,0.02,4,56.3,74.1,0.695,34.2,0.299,1,0.186,0.329,0.396,0.0233,5.61e-06,3.4e-06,0.977,0.0955
         Hull,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
         Leeds,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
         Lees,-548.1,1.1,0.02,4,45.6,244,0.612,29.7,0.153,0.968,0.0311,0.648,0.681,0.0477,2.66e-05,2.08e-05,0.952,0.0802
         Liverpool,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
         London,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
         Manchester,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
         Mold,-296.5,0.25,0.02,4,67.4,301,0.131,21.4,0.271,1.04,0.0145,0.436,2.87,0.064,2.61e-05,2.27e-05,0.936,0.0544
         Northwich,-1195.1,2.25,0.02,4,45.6,147,0.795,30.1,0.423,0.948,0.0602,0.236,0.402,0.0213,1.32e-05,1.58e-05,0.979,0.0857
         Nottingham,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
         Oswestry,-696.1,0.49,0.02,4,37.3,168,0.631,52.9,0.339,1.04,0.0298,0.263,0.476,0.0218,1.56e-05,1.61e-05,0.978,0.0699
         Sheffield,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
         ",stringsAsFactors=FALSE) -> mles

## ----mle-----------------------------------------------------------------
mles %>% subset(as.character(mles$town)=="         London") -> mle
paramnames <- c("R0","mu","sigma","gamma","alpha","iota",
                "rho","sigmaSE","psi","cohort","amplitude",
                "S_0","E_0","I_0","R_0")
mle %>% extract(paramnames) %>% unlist() -> theta
mle %>% subset(select=-c(S_0,E_0,I_0,R_0)) %>%
  knitr::kable(row.names=FALSE)


## ----transforms----------------------------------------------------------
toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tcohort = logit(cohort);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  TR0 = log(R0);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")
fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    Tcohort = expit(cohort);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    TR0 = exp(R0);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")


pomp(m1,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     statenames=c("S","E","I","R","C","W"),
     paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                  "rho","sigmaSE","psi","cohort","amplitude",
                  "S_0","E_0","I_0","R_0")) -> m1


# PMCMC -----------------------------------------------------------------
Ntat<- 400
Nparticle <- 500

sd0 <- c(0,0,rep(0.1,length(theta)-13),rep(0,11)) # only two parameters move
theta0 <- theta
names(sd0) <- names(theta)
set.seed(11)
theta0[3] <- runif(1,28,29)
theta0[4] <- runif(1,30,30.5)

pmcmc(m1,Nmcmc = Ntat,
      start = theta0,
      # Adaptative MCMC Proposal 
      proposal =  mvn.rw.adaptive(rw.sd=
                                    sd0,scale.start=300),
      # rw.sd is a named numeric vector; random-walk SDs for a 
      # multivariate normal randomwalk proposal with diagonal 
      # variance-covariance matrix.
      Np = Nparticle)-> chain
trace1 <- window(conv.rec(chain,c("sigma","gamma")),start=20)
summary(trace1)
plot(trace1)