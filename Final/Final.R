library(nullabor)
library(plyr)
library(ggplot2)
library(MCMCpack)
library(dplyr)
library(reshape2)
library(expm)
library(tidyr)
library(R2jags)
library(rjags)
library(gridExtra)

######################################### Pregunta 1 ############################################## 
places <- read.table("places.csv", header=TRUE, quote="\"")
places_ch <- places[,c("Climate","HousingCost")]
head(places_ch)
# valores bajos de clima implican climas inconvenientes
# valores altos en housing implican costos altos

# la hipotesis nula considera que no hay relacion entre los datos

# los datos originales se ven como siguen
ggplot(places_ch, aes(x = HousingCost, y = Climate)) +
  geom_point()

# inferencia por graficas
nrow(places_ch)

clim_null <- lineup(null_permute('HousingCost'), n = 10, places_ch)

ggplot(clim_null, aes(x = HousingCost, y = Climate)) +
  facet_wrap(~ .sample) +
  geom_jitter(position = position_jitter(width = 0.1, height = 1), 
              size = 1, alpha = 0.5)

decrypt("xrgp bViV NE njJNiNjE 0")
decrypt("xrgp bViV NE njJNiNjE G")


# la grafica real esta en la ventana 1 por lo tanto rechazamos la hipotesis nula
# de que las variables no estan relacionadas entre si.

######################################### Pregunta 2 #######################################

# preparacion de la informacion
wages <- read.csv("wages.csv")
wages <- subset(wages, hgc >= 9 & hgc <= 11 & lnw < 3.5)
wages$raza <- ifelse((wages$hispanic == 1), "hisp",ifelse((wages$black==1), "black","white"))

# usamos el mismo numpero de individuos por raza
u.hisp  <- unique(subset(wages, raza=="hisp",select=c(id,raza)))
u.black <- unique(subset(wages, raza=="black",select=c(id,raza)))
u.white <- unique(subset(wages, raza=="white",select=c(id,raza)))

ind <- c(sample(u.hisp$id,size=110,replace=FALSE),
         sample(u.black$id,size=110,replace=FALSE),
         sample(u.white$id,size=110,replace=FALSE))

wages_r <- wages[wages$id %in% ind, ]


# generamos datos nulos 
wages_null <- lineup(null_permute('raza'), n = 20, wages_r)

ggplot(wages_null, aes(x = exper, y = lnw)) +
  facet_wrap(~ .sample) +
  geom_point(alpha = 0.02, size = 2) + 
  geom_smooth(aes(group = raza, color = raza), method = "loess", se = FALSE) 

decrypt("xrgp bViV NE njJNiNjE 01")

######################################### Pregunta 3 #######################################


###################################a
########## bootstrap
###################################a

load("x.Rdata")

n<-length(x)

var <- (1 / n * sum((x - 0) ^ 2)) 
# 131.291

# Paso 1: calculamos la varianza
sigma2 <- 1 / n * sum((x - 0) ^ 2)

# Pasos 2 y 3
thetaBoot <- function(){
  # Simular X_1*,...X_N* con distribución N(mu_hat, sigma_hat^2) 
  x_boot <- rnorm(n, mean = 0, sd = sqrt(sigma2)) 
  # Calcular mu*, sigma* y theta*
  sigma2_boot <- 1 / n * sum((x_boot - 0) ^ 2)
  sigma2_boot
}

# Paso 4: Repetimos B = 3000 veces estimamos el error estándar
sims_boot <- rdply(3000, thetaBoot())

# sigma2 a partir de bootstrap
sigma2_b <- mean(sims_boot$V1)
sigma2_b

# erro estandar bootstrap
es_b <- sqrt(1 / 2999 * sum((sims_boot$V1 - sigma2) ^ 2))

# histograma de simulaciones bootstrap
ggplot(sims_boot, aes(x = V1)) +
  geom_histogram(aes(y = ..density..), binwidth = 2, fill = "gray") +
  geom_vline(xintercept = c(sigma2_b), 
             color = c("red"))



###################################a
########## analisis bayesiano
###################################a

# definimos la seleccion de los parametros 

shape0 <- 1
scale0 <- 1

# 1. simulamos valores porvenientes de una distribución gamma
x_gamma <- rgamma(100000, shape = shape0, rate = scale0)
# 2. invertimos los valores simulados
x_igamma <- 1 / x_gamma
x_igamma <- data.frame(x_igamma)

# grafica de la funcion de densidad de la gamma inversa
ggplot(x_igamma, aes(x = x_igamma)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "gray") + 
  stat_function(fun = dinvgamma, args = list(shape = shape0, scale = scale0), 
                color = "blue") +
  scale_x_continuous(limits = c(0, 7))



# calculo de la distribucion posterior
alpha1 <- shape0 + n/2
beta1 <- scale0 + sum((x - 0)^2)/2

posterior <- rinvgamma(10000, alpha1, beta1)

# valor de sigma2 estimado por analisis bayesianos
sigma2_bayes<- mean(posterior)
sigma2_bayes


# histograma posterior
ggplot(as.data.frame(posterior), aes(x = posterior)) +
  geom_histogram(aes(y = ..density..), binwidth = 2, fill = "gray") + 
  stat_function(fun = dinvgamma, args = list(shape = alpha1, scale = beta1), 
                color = "blue") +
  geom_vline(xintercept = c(sigma2_bayes), 
             color = c("red"))

# error estandar
es_bayes <- sqrt(var(posterior)) 
es_bayes

###################################a
########## enfoque JAGS
###################################a

N <- n
modelo_normal_sigma.txt <-
  '
model{
  for(i in 1:N){
    x[i] ~ dnorm(0, nu)
  }
  # iniciales
  nu <- 1/sigma2
  sigma2 ~ dunif(0.1, 300)
}
'

cat(modelo_normal_sigma.txt, file = 'modelo_normal_sigma.bugs')

# valores iniciales para los parámetros, si no se especifican la función jags
# generará valores iniciales
jags.inits <- function(){
  list("nu" = runif(1, 0.1, 300))
}

jags_fit_sigma <- jags(
  model.file = "modelo_normal_sigma.bugs",    # modelo de JAGS
  #inits = jags.inits,   # valores iniciales
  data = list(x = x, N = N),    # lista con los datos
  parameters.to.save = c("nu", "sigma2"),  # parámetros por guardar
  n.chains = 3,   # número de cadenas
  n.iter = 10000,    # número de pasos
  n.burnin = 1000,   # calentamiento de la cadena
  n.thin = 1
)

jags_fit_sigma

# podemos ver las cadenas
traceplot(jags_fit_sigma, varname = c("nu", "sigma2"))

head(jags_fit_sigma$BUGSoutput$summary)

sigmas <- (jags_fit_sigma$BUGSoutput$sims.matrix[, 3])
sigma2_jags <- mean(sigmas)
sigma2_jags

# histograma posterior
ggplot(as.data.frame(sigmas), aes(x = sigmas)) +
  geom_histogram(aes(y = ..density..), binwidth = 2, fill = "gray")  +
  geom_vline(xintercept = c(sigma2_jags), 
             color = c("red"))


# error estandar
es_jags <- sqrt(var(sigmas)) 
es_jags


################a
## comparacion
################a

length(sims_boot$V1)
length(posterior)
length(sigmas)

histogramas_b <- as.data.frame(sample(sims_boot$V1,size=3000,replace=FALSE))
colnames(histogramas_b) <- c("Hist")
histogramas_b$tipo <- "bootstrap"
head(histogramas_b)

histogramas_bayes <- as.data.frame(sample(posterior,size=3000,replace=FALSE))
colnames(histogramas_bayes) <- c("Hist")
histogramas_bayes$tipo <- "bayes"
head(histogramas_bayes)

histogramas_jags <- as.data.frame(sample(sigmas,size=3000,replace=FALSE))
colnames(histogramas_jags) <- c("Hist")
histogramas_jags$tipo <- "jags"
head(histogramas_jags)

histogramas <- rbind(histogramas_b, histogramas_bayes, histogramas_jags)



hist_cut <- ggplot(histogramas, aes(x=Hist, fill=tipo)) + 
  geom_bar(alpha=0.4, binwidth=3) +
  geom_vline(xintercept = c(sigma2_jags,sigma2_bayes,sigma2_b), 
             color = c("blue","red","green"))

hist_cut

trial <- matrix(c(sigma2_b,es_b,sigma2_bayes,es_bayes,sigma2_jags,es_jags), ncol=3)
colnames(trial) <- c('Bootstrap', 'Bayes', 'JAGS')
rownames(trial) <- c('Sigma2', 'Error Estandar')
trial.table <- as.table(trial)
trial.table



####################################a
######## analisis de log(sigma)
#####################################a

################a
## bootstrap
################a

# estimador de max verosimilitud es
lnsigma <- log(sqrt(1 / n * sum((x - 0) ^ 2)))

# utilizamos bootstrap parametrico 
lnsigmaBoot <- function(){
  # Simular X_1*,...X_N* con distribución N(mu_hat, sigma_hat^2) 
  x_boot <- rnorm(n, mean = 0, sd = sqrt(sigma2)) 
  # Calcular mu*, sigma* y theta*
  lnsigma_boot <- log(sqrt(1 / n * sum((x_boot - 0) ^ 2)))
  lnsigma_boot
}

# Paso 4: Repetimos B = 3000 veces estimamos el error estándar
sims_lnsigma_boot <- rdply(3000, lnsigmaBoot())
lnsigma_boot <- mean(sims_lnsigma_boot$V1)
lnsigma_boot

# histograma de las replicaciones bootstrap
ggplot(sims_lnsigma_boot,aes(x=V1))+
  geom_histogram(binwidth = .05)


# como los cuantiles de la distribucion de la estimacion de lnsigma son muy parecidos 
# a los cuantiles de una normal, usaremos los cuantiles del histograma bootstrap para
# definir los intervalos de confianza

round(quantile(sims_lnsigma_boot$V1, probs = c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975)), 2)
round(qnorm(p = c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975), mean = lnsigma, 
            sd = sd(sims_lnsigma_boot$V1)), 2)

# intervalo de confianza por percentiles
ls_per <- round(quantile(sims_lnsigma_boot$V1, prob = 0.975), 2)
li_per <- round(quantile(sims_lnsigma_boot$V1, prob = 0.025), 2)


# hisograma con intervalos de confianza
ggplot(sims_lnsigma_boot, aes(x = V1)) + 
  geom_histogram(binwidth = 0.02, alpha=.7,fill = "gray30") +
  geom_vline(xintercept = c(li_per, ls_per, lnsigma_boot), 
             color = c("black", "black", "red"))


################a
## jags
################a


N <- n
modelo_normal_tau.txt <-
  '
model{
for(i in 1:N){
x[i] ~ dnorm(0, nu)
}
# iniciales
nu <- 1/sigma2
sigma2 ~ dunif(0.1, 300)
tau <- log(sqrt(sigma2))
}
'

cat(modelo_normal_tau.txt, file = 'modelo_normal_tau.bugs')

# valores iniciales para los parámetros, si no se especifican la función jags
# generará valores iniciales
#jags.inits <- function(){
#  list("nu" = runif(1, 0.1, 300))
#}

jags_fit_tau <- jags(
  model.file = "modelo_normal_tau.bugs",    # modelo de JAGS
  #inits = jags.inits,   # valores iniciales
  data = list(x = x, N = N),    # lista con los datos
  parameters.to.save = c("nu", "sigma2", "tau"),  # parámetros por guardar
  n.chains = 3,   # número de cadenas
  n.iter = 10000,    # número de pasos
  n.burnin = 1000,   # calentamiento de la cadena
  n.thin = 1
)

jags_fit_tau

# podemos ver las cadenas
traceplot(jags_fit_tau, varname = c("nu", "sigma2"))

head(jags_fit_tau$BUGSoutput$summary)

lnsigmas_jags <- (jags_fit_tau$BUGSoutput$sims.matrix[, 4])
lnsigma_jags <- mean(lnsigmas_jags)
lnsigma_jags

# intervalo de confianza por percentiles
ls_per_jags <- round(quantile(lnsigmas_jags, prob = 0.975), 2)
li_per_jags <- round(quantile(lnsigmas_jags, prob = 0.025), 2)



# histograma posterior
ggplot(as.data.frame(lnsigmas_jags), aes(x = lnsigmas_jags)) +
  geom_histogram(aes(y = ..density..), binwidth = .01, fill = "gray")  +
  geom_vline(xintercept = c(li_per_jags, ls_per_jags, lnsigma_jags), 
             color = c("black", "black", "red"))

######################################### Pregunta 4 #######################################


prior <- function(mu, tau){
  function(theta) {
    1/sqrt(2*pi*tau^2)*exp(-1/(2*tau^2)*((theta-mu)^2))
  }
}

mi_prior<-prior(150,15)

likeNorm <- function(S, S2, N){
  function(theta){
    (1/(800*pi))^(N/2)*exp(-(1/800)*(S2-2*theta*S+N*(theta^2)))
  }
}

mi_like<-likeNorm(13000,1700000,100)

postRelProb <- function(theta){
  mi_like(theta) * mi_prior(theta)
}


### a) rechazados para N(0,0.2)
#Metrópolis

set.seed(114041)
# para cada paso decidimos el movimiento de acuerdo a la siguiente función
caminaAleat <- function(theta){ # theta: valor actual
  salto_prop <- rnorm(1, 0, sd = sqrt(0.2)) # salto propuesto
  theta_prop <- theta + salto_prop # theta propuesta
  #if(theta_prop < 0 | theta_prop > 1){ # si el salto implica salir del dominio
  #  return(theta)
  #}
  u <- runif(1) 
  p_move = min(postRelProb(theta_prop) / postRelProb(theta), 1) # prob mover
  if(p_move  > u){
    return(theta_prop) # aceptar valor propuesto
  }
  else{
    return(theta) # rechazar
  }
}

pasos <- 6000
camino_A <- numeric(pasos) # vector que guardará las simulaciones
camino_A[1] <- 150 # valor inicial

rechazo = 0

# Generamos la caminata aleatoria
for (j in 2:pasos){
  camino_A[j] <- caminaAleat(camino_A[j-1])
  rechazo <- rechazo + 1 * (camino_A[j] == camino_A[j - 1])
}


caminata <- data.frame(pasos = 1:pasos, theta = camino_A)

pasosN02 <- ggplot(caminata[1:2000, ], aes(x = pasos, y = theta)) +
  geom_point(size = 0.8) +
  geom_path(alpha = 0.5) +
  scale_y_continuous(expression(theta), limits = c(120, 160)) +
  scale_x_continuous("Tiempo") +
  geom_hline(yintercept = mean(caminata$theta), color = "red", alpha = 0.5)


rechazo / pasos
### 0.0785



### b) rechazados para N(0,5)
#Metrópolis

set.seed(114041)
# para cada paso decidimos el movimiento de acuerdo a la siguiente función
caminaAleat <- function(theta){ # theta: valor actual
  salto_prop <- rnorm(1, 0, sd = sqrt(5)) # salto propuesto
  theta_prop <- theta + salto_prop # theta propuesta
  #if(theta_prop < 0 | theta_prop > 1){ # si el salto implica salir del dominio
  #  return(theta)
  #}
  u <- runif(1) 
  p_move = min(postRelProb(theta_prop) / postRelProb(theta), 1) # prob mover
  if(p_move  > u){
    return(theta_prop) # aceptar valor propuesto
  }
  else{
    return(theta) # rechazar
  }
}

pasos <- 6000
camino_B <- numeric(pasos) # vector que guardará las simulaciones
camino_B[1] <- 150 # valor inicial

rechazo = 0

# Generamos la caminata aleatoria
for (j in 2:pasos){
  camino_B[j] <- caminaAleat(camino_B[j-1])
  rechazo <- rechazo + 1 * (camino_B[j] == camino_B[j - 1])
}


caminata <- data.frame(pasos = 1:pasos, theta = camino_B)

pasosN05 <- ggplot(caminata[1:2000, ], aes(x = pasos, y = theta)) +
  geom_point(size = 0.8) +
  geom_path(alpha = 0.5) +
  scale_y_continuous(expression(theta), limits = c(120, 160)) +
  scale_x_continuous("Tiempo") +
  geom_hline(yintercept = mean(caminata$theta), color = "red", alpha = 0.5)
pasosN05

rechazo / pasos
### 0.3216667

### c) rechazados para N(0,20)
#Metrópolis

set.seed(114041)
# para cada paso decidimos el movimiento de acuerdo a la siguiente función
caminaAleat <- function(theta){ # theta: valor actual
  salto_prop <- rnorm(1, 0, sd = sqrt(20)) # salto propuesto
  theta_prop <- theta + salto_prop # theta propuesta
  #if(theta_prop < 0 | theta_prop > 1){ # si el salto implica salir del dominio
  #  return(theta)
  #}
  u <- runif(1) 
  p_move = min(postRelProb(theta_prop) / postRelProb(theta), 1) # prob mover
  if(p_move  > u){
    return(theta_prop) # aceptar valor propuesto
  }
  else{
    return(theta) # rechazar
  }
}

pasos <- 6000
camino_C <- numeric(pasos) # vector que guardará las simulaciones
camino_C[1] <- 150 # valor inicial

rechazo = 0

# Generamos la caminata aleatoria
for (j in 2:pasos){
  camino_C[j] <- caminaAleat(camino_C[j-1])
  rechazo <- rechazo + 1 * (camino_C[j] == camino_C[j - 1])
}

caminata <- data.frame(pasos = 1:pasos, theta = camino_C)

pasosN020 <- ggplot(caminata[1:2000, ], aes(x = pasos, y = theta)) +
  geom_point(size = 0.8) +
  geom_path(alpha = 0.5) +
  scale_y_continuous(expression(theta), limits = c(120, 160)) +
  scale_x_continuous("Tiempo") +
  geom_hline(yintercept = mean(caminata$theta), color = "red", alpha = 0.5)
pasosN020

rechazo / pasos
### 0.5348333


### a) valores de la distribucion posterior y grafica para N(0,0.2)

caminata_A <- data.frame(pasos = 1:pasos, theta = camino_A)
grafica_A <- ggplot(caminata_A[1:2000, ], aes(x = pasos, y = theta)) +
                  geom_point(size = 0.8) +
                  geom_path(alpha = 0.3) 
grafica_A

### b) valores de la distribucion posterior y grafica para N(0,5)

caminata_B <- data.frame(pasos = 1:pasos, theta = camino_B)
grafica_B <- ggplot(caminata_B[1:2000, ], aes(x = pasos, y = theta)) +
  geom_point(size = 0.8) +
  geom_path(alpha = 0.3) 
grafica_B


### c) valores de la distribucion posterior y grafica para N(0,20)

caminata_C <- data.frame(pasos = 1:pasos, theta = camino_C)
grafica_C <- ggplot(caminata_C[1:2000, ], aes(x = pasos, y = theta)) +
  geom_point(size = 0.8) +
  geom_path(alpha = 0.3) 
grafica_C


grid.arrange(grafica_A, grafica_B, grafica_C, ncol = 1)
grid.arrange(pasosN02, pasosN05, pasosN020, ncol = 1)








### a) histograma de la distribucion posterior para N(0,0.2)
hist_A <- ggplot(filter(caminata_A, pasos > 1000), aes(x = theta)) +
                geom_histogram(aes(y = ..density..),binwidth=1,alpha=.4) 
hist_A

### b) histograma de la distribucion posterior para N(0,5)
hist_B <- ggplot(filter(caminata_B, pasos > 1000), aes(x = theta)) +
                 geom_histogram(aes(y = ..density..),binwidth=1,alpha=.4) 
hist_B

### c) histograma de la distribucion posterior para N(0,20)
hist_C <- ggplot(filter(caminata_C, pasos > 1000), aes(x = theta)) +
                  geom_histogram(aes(y = ..density..),binwidth=1,alpha=.4) 
hist_C

grid.arrange(hist_A, hist_B, hist_C, ncol = 3)



### histograma de la probabilidad predictiva posterior
caminata_f <- filter(caminata_B, pasos > 1000)

sims_y <- rnorm(nrow(caminata_f), mean = mean(caminata_f$theta), sd=sqrt(5) )
mean(sims_y) # p(y = 1 | x) probabilidad predictiva

ls_per_pred <- round(quantile(as.data.frame(sims_y)$sims_y, prob = 0.975), 2)
li_per_pred <- round(quantile(as.data.frame(sims_y)$sims_y, prob = 0.025), 2)

ggplot(as.data.frame(sims_y), aes(x = sims_y)) +
  geom_histogram(aes(y = ..density..),binwidth=1,alpha=.4) +
  geom_vline(xintercept = c(mean(sims_y),li_per_pred,ls_per_pred), 
             color = c("red","blue","blue"))

######################################### Pregunta 5 #######################################


modelo_regresion.txt <-
  '
model{
  for(i in 1 : n) {
    y[i] ~ dnorm(y.hat[i], tau.y) 
    y.hat[i] <- a + b * x[i]
  }
  a ~ dnorm(0, 0.001)
  b ~ dnorm(0, 0.001)
  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 100)
}
'
cat(modelo_regresion.txt, file = 'modelo_regresion.bugs')

### Radon
load("radon.Rdata")


# Iniciamos preparando los datos para el análisis, trabajaremos en
# escala logarítmica, hay algunos casos con medición cero, para éstos
# hacemos una pequeña correción redondeándolos a 0.1.
y <- log(ifelse (radon.2$activity == 0, 0.1, radon.2$activity))
n <- nrow(radon.2)
x <- radon.2$floor

# jags
radon1.data <- list("n", "y", "x")
radon1.inits <- function(){
  list(a = rnorm(1), 
       b = rnorm(1), 
       sigma.y = runif(1))}

radon1.parameters <- c("a", "b", "sigma.y")

radon1.jags <- jags(
  model.file = "modelo_regresion.bugs",
  inits= radon1.inits, 
  data = radon1.data, 
  parameters.to.save = radon1.parameters, 
  n.chains = 3,   # número de cadenas
  n.iter = 10000,    # número de pasos
  n.burnin = 1000,   # calentamiento de la cadena
  n.thin = 1)

traceplot(radon1.jags, varname = c("a", "b", "sigma.y"))

######################################### Pregunta 6 #######################################



################a
## 6.1
################a

# la variable y_ij la observación en el i-ésimo conejo perteneciente al j-ésimo experimento
# 1 indicando que el conejo desarrolló tumor y 0 en el caso contrario. Corresponde a 
# la realizacion de lanzar la moneda

# la variable theta_j que es la probabilidad de desarrollar un tumor en el jésimo grupo, 
# corresponde a las diferentes monedas que participan en el experimento

################a
## 6.2
################a

load("rabbits.RData")
head(rabbits)
# hay 21 realizaciones en cada experimento

modelo_jer.txt <-
  '
model{
for(t in 1:N){
x[t] ~ dbern(theta[grupo[t]])
}
for(j in 1:nGrupos){
theta[j] ~ dbeta(a, b)
}
a <- mu * kappa
b <- (1 - mu) * kappa
mu ~ dbeta(1, 1)
kappa ~ dgamma(1, 0.1)
}
'

cat(modelo_jer.txt, file = 'modelo_jer.bugs')

jags_fit_rab <- jags(
  model.file = "modelo_jer.bugs",    # modelo de JAGS
  #inits = jags.inits,   # valores iniciales
  data = list(x = rabbits$tumor, grupo = rabbits$experiment, nGrupos = 71,  N = nrow(rabbits)),    # lista con los datos
  parameters.to.save = c("mu", "kappa", "theta"),  # parámetros por guardar
  n.chains = 3,   # número de cadenas
  n.iter = 10000,    # número de pasos
  n.burnin = 1000   # calentamiento de la cadena
)

traceplot(jags_fit_rab, varname = c("kappa", "mu"))

# validacion de convergencia
jags_fit_rab


# histograma de la distribucion posterior para mu y para kappa

sims_df <- data.frame(n_sim = 1:jags_fit_rab$BUGSoutput$n.sims,
                      jags_fit_rab$BUGSoutput$sims.matrix) %>% 
  dplyr::select(-deviance) %>%
  gather(parametro, value, -n_sim)

mu_rab_62 <- ggplot(filter(sims_df, parametro == "mu"), aes(x = value)) +
  geom_histogram(alpha = 0.8, binwidth=.003)+ 
  xlab("Mu")

kappa_rab_62 <- ggplot(filter(sims_df, parametro == "kappa"), aes(x = value)) +
  geom_histogram(alpha = 0.8, binwidth=4)+ 
  xlab("kappa")

grid.arrange(mu_rab_62, kappa_rab_62, ncol = 2)

################a
## 6.3
################a

ggplot(subset(sims_df, parametro != "kappa" & parametro != "mu"), aes(x = parametro, y = value)) +
  geom_boxplot()

################a
## 6.4
################a

modelo_jer_64.txt <-
  '
model{
for(t in 1:N){
x[t] ~ dbern(theta[grupo[t]])
}
for(j in 1:nGrupos){
theta[j] ~ dbeta(a, b)
}
a <- mu * kappa
b <- (1 - mu) * kappa
mu ~ dbeta(10, 10)
kappa ~ dgamma(.51, 0.1)
}
'

cat(modelo_jer_64.txt, file = 'modelo_jer_64.bugs')

jags_fit_rab_64 <- jags(
  model.file = "modelo_jer_64.bugs",    # modelo de JAGS
  #inits = jags.inits,   # valores iniciales
  data = list(x = rabbits$tumor, grupo = rabbits$experiment, nGrupos = 71,  N = nrow(rabbits)),    # lista con los datos
  parameters.to.save = c("mu", "kappa", "theta"),  # parámetros por guardar
  n.chains = 3,   # número de cadenas
  n.iter = 10000,    # número de pasos
  n.burnin = 1000   # calentamiento de la cadena
)

traceplot(jags_fit_rab_64, varname = c("kappa", "mu"))

# validacion de convergencia
head(jags_fit_rab_64)


#grafica comoparacion medias posteriores
mod1 <- as.data.frame(jags_fit_rab$BUGSoutput$summary)
mod2 <- as.data.frame(jags_fit_rab_64$BUGSoutput$summary)
mod1 <- mod1[-c(1,2,3),1]
mod2 <- mod2[-c(1,2,3),1]
thetas_modelos <- as.data.frame(cbind(mod1,mod2))


ggplot(thetas_modelos, aes(x=mod1, y=mod2)) +
  geom_point(alpha=.5) + 
  geom_abline(intercept = 0, colour="red")








