################ Tarea 10 ##################################
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(fda)
##### base de datos ####
set.seed(114041)
mu <- 5
n <- 100
x <- rnorm(n,mu,1)
base <- data.frame(x=x)
real<-exp(mu)
real
# el error estandar via método delta es : 1/sqrt(n) y por lo tanto
# es posible obtener intervalosde confianza para mu

se_delta <- 1/sqrt(n)
IC_delta <- c(mean(base$x)-1.96*(se_delta), mean(base$x)+1.96*(se_delta))
matrix(IC_delta,nrow=1,ncol=2,byrow=TRUE,dimnames=list(c("delta"),c("lim inf","lim sup")))

##### bootstrap parametrico #####

thetaboot<-function(){
  return(mean(rnorm(100,mean(base$x),1)))
}
sims_boot <- rdply(1000,thetaboot)
IC_boot_param <- c(mean(sims_boot$V1)-1.96*(sd(sims_boot$V1)),
                   mean(sims_boot$V1)+1.96*(sd(sims_boot$V1)))
IC_boot_param
IC<-c(IC_delta , IC_boot_param)
matrix(IC,nrow=2,ncol=2,byrow=TRUE,dimnames=list(c("delta","boot_param"),c("lim inf","lim sup")))

##### bootstrap no parametrico #####
thetaboot_np<-function(){
  return(mean(sample(base$x,size=100,replace=TRUE)))
}

sims_boot_np <- rdply(1000,thetaboot_np)
IC_boot_np <- c(mean(sims_boot_np$V1)-1.96*(sd(sims_boot_np$V1)),
                   mean(sims_boot_np$V1)+1.96*(sd(sims_boot_np$V1)))
IC_boot_np
IC<-c(IC_delta , IC_boot_param, IC_boot_np)
comp<-as.data.frame(matrix(IC,nrow=3,ncol=2,byrow=TRUE,
                           dimnames=list(c("delta","boot param","boot no param"),c("lim_inf","lim_sup"))))
comp$dif<-comp$lim_sup-comp$lim_inf
comp

# el metodo con menor intervalo de confianza es el del bootstrap parametrico

sims_delta<-rnorm(1000,mean(base$x),se_delta)
delta<-as.data.frame(sims_delta)
ggplot(delta,aes(x =sims_delta)) +  geom_histogram(fill = "blue")


boot <- data.frame(parametrico = sims_boot$V1, no_parametrico = sims_boot_np$V1)
boot$indice <- 1:1000
boots <- gather(boot,bootstrap, value, -indice)

ggplot(boots, aes(x=value,color=bootstrap)) +
  facet_wrap(~bootstrap) +
  geom_histogram() +
  geom_vline(xintercept = c(IC_boot_np[1],IC_boot_np[2],5,IC_boot_param[1],IC_boot_param[2]),
             size = 0.5,
             color = c("cyan3","cyan3","black","indianred2","indianred2")) +
  labs(title = "Bootstrap parametrico y no parametrico") +
  theme(plot.title = element_text(size = rel(1.8)))


####################### Vancouver ##############################



vancouver <- read.csv("C:/Users/Administrador/Desktop/tarea10/vancouver.csv")
vancouver <- arrange(vancouver,id)

knots <- quantile(vancouver$id)
base <- create.bspline.basis(
  norder = 4, # polinomios cubicos
  breaks = knots # nodos en los cuartiles de x
)
plot(base, lty = "solid")
H <- eval.basis(vancouver$id, base)
head(H)
beta_hat <- as.vector(solve(t(H) %*% H) %*% t(H) %*% vancouver$prec)
beta_hat
mu <- function(x, betas){
  as.numeric(betas %*% t(eval.basis(x, base)))
}
ggplot(vancouver, aes(x = id, y = prec)) +
  geom_point(color = "red") +
  stat_function(fun = mu, arg = list(betas = beta_hat), color = "blue") +
  labs(title = "B-splines")
n <- 365
mu_hat <- mu(vancouver$id, beta_hat)
sigma_hat <- sqrt(1 / n * sum((vancouver$prec - mu_hat) ^ 2))
# creamos las muestras bootstrap (parametrico)
splinesBootP <- function(){
  vancouver_boot <- data.frame(id = vancouver$id, prec = mu_hat + rnorm(n, 0,
                                                                        sigma_hat))
  H <- eval.basis(vancouver_boot$id, base)
  as.vector(solve(t(H) %*% H) %*% t(H) %*% vancouver_boot$prec)
}
sims_splines <- rdply(2000, splinesBootP())
betas <- as.matrix(sims_splines[, -1])
splines_boot <- ggplot(vancouver, aes(x = id, y = prec))
for(i in 1:100){
  splines_boot <- splines_boot +
    stat_function(fun = mu, arg = list(betas = betas[i, ]), alpha = 0.1)
}
splines_boot + geom_point(color = "red")






