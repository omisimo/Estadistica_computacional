library(ggplot2)
library(plyr)
library(dplyr)
library(arm)
library(tidyr)
library(Hmisc)
library(fda) 
set.seed(918739837)

############################ SIMULACION DISCRETA PREDICTIVA ##########################

n_ninas <- rbinom(1, 400, 0.488)
sims_ninas <- rdply(1000, rbinom(1, 400, 0.488))
mean(sims_ninas$V1)
ggplot(sims_ninas, aes(x = V1)) + geom_histogram(binwidth = 2)
tipo_nacimiento <- sample(c("unico", "fraternal", "identicos"), 
                          size = 400, replace = TRUE, prob = c(1 - 1 / 125 - 1 / 300, 1 / 125, 1 / 300))
n_unico <- sum(tipo_nacimiento == "unico")  # número de nacimientos únicos
n_fraternal <- sum(tipo_nacimiento == "fraternal")
n_identicos <- 400 - n_unico - n_fraternal
n_ninas <- rbinom(1, n_unico, 0.488) +
  rbinom(1, 2 * n_fraternal, 0.495) + # en cada nacimiento hay 2 bebés
  2 * rbinom(1, n_identicos, 0.495)
n_ninas

modelo2 <- function(){
  tipo_nacimiento <- sample(c("unico", "fraternal", "identicos"), 
                            size = 400, replace = TRUE, 
                            prob = c(1 - 1 / 125 - 1 / 1300, 1 / 125, 1 / 300))
  # número de nacimientos de cada tipo
  n_unico <- sum(tipo_nacimiento == "unico")  # número de nacimientos únicos
  n_fraternal <- sum(tipo_nacimiento == "fraternal")
  n_identicos <- 400 - n_unico - n_fraternal
  # simulamos para cada tipo de nacimiento
  n_ninas <- rbinom(1, n_unico, 0.488) +
    rbinom(1, 2 * n_fraternal, 0.495) + # en cada nacimiento hay 2 bebés
    2 * rbinom(1, n_identicos, 0.495)
  n_ninas
}

sims_ninas2 <- rdply(1000, modelo2)
mean(sims_ninas2$V1)
ggplot(sims_ninas2, aes(x = V1)) + geom_histogram(binwidth = 2)


############################ SIMULACION CONTINUA PREDICTIVA ##########################

sexo <- rbinom(10, 1, 0.52)
altura <- rnorm(sexo, mean = 161.8 * (sexo == 1) + 175 * (sexo == 0), 
                sd = 6.86 * (sexo == 1) + 7.37 * (sexo == 0))
mean(altura)
mediaAltura <- function(){
  sexo <- rbinom(10, 1, 0.52)
  altura <- rnorm(sexo, mean = 161.8 * (sexo == 1) + 175 * (sexo == 0), 
                  sd = 6.86 * (sexo == 1) + 7.37 * (sexo == 0))
  mean(altura)
}
sims_alturas <- rdply(1000, mediaAltura)
mean(sims_alturas$V1)
ggplot(sims_alturas, aes(x = V1)) + geom_histogram(binwidth = 1)






# Verosimilitud X_1,...,X_n ~ Bernoulli(theta)
L_bernoulli <- function(n, S){
  function(theta){
    theta ^ S * (1 - theta) ^ (n - S)
  }  
}
# log-verosimilitud
l_bernoulli <- function(n, S){
  function(theta){
    S * log(theta) + (n - S) * log(1 - theta)
  }  
}
xy <- data.frame(x = 0:1, y = 0:1)
ggplot(xy, aes(x = x, y = y)) +
  stat_function(fun = L_bernoulli(n = 20, S = 12)) +
  xlab(expression(theta)) +
  ylab(expression(L(theta))) +
  ggtitle("Verosimilitud (n=20, S = 12)")

### BOOTSTRAP PARAMETRICO

n <- 150
x <- rnorm(n, mean = 10, sd = 5)  # observaciones normales

# Paso 1: calcular mu_hat y sigma_hat
mu_hat <- mean(x)  
sigma_hat <- sqrt(1 / n * sum((x - mu_hat) ^ 2)) 

# Pasos 2 y 3
thetaBoot <- function(){
  # Simular X_1*,...X_N* con distribución N(mu_hat, sigma_hat^2) 
  x_boot <- rnorm(n, mean = mu_hat, sd = sigma_hat) 
  # Calcular mu*, sigma* y theta*
  mu_boot <- mean(x_boot)
  sigma_boot <- sqrt(1 / n * sum((x_boot - mu_boot) ^ 2)) 
  sigma_boot / mu_boot # theta*
}

# Paso 4: Repetimos B = 2000 veces estimamos el error estándar
sims_boot <- rdply(3000, thetaBoot())
sqrt(1 / 2999 * sum((sims_boot$V1 - sigma_hat/mu_hat) ^ 2))
1 / sqrt(n) * (1 / mu_hat ^ 4 + sigma_hat ^ 2 / (2 * mu_hat ^ 2)) ^ (1 / 2)


#### BOOTSTRAP PARAMETRICO,NO PARAMETRICO Y VEROSIMILITUD

set.seed(90984)
shm <- function(t, A = 1.5, omega = 4){ # Esta es una función sinusoidal
  t * A * sin(omega * t)
}
n <- 90
x <- sample(seq(0, 3, 0.02), n) # creamos una base con n observaciones
y <- shm(x) + rnorm(length(x), sd = 1)
outliers <- sample(1:length(y), 4)  # elijo 4 puntos al azar a los que agrego ruido
y[outliers] <- y[outliers] + rnorm(4, sd = 2)
toy <- data.frame(x, y)

ggplot(toy, aes(x, y)) +
  geom_point()


toy_1 <- toy
toy_1 <- toy %>%
  mutate(int = cut2(x, g = 4)) %>%
  group_by(int) %>%
  mutate(media = mean(y))

ggplot(toy_1, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, y = media, group = int), color = "red")

ggplot(toy_1, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", aes(x, y = y, group = int), color = "red", se = FALSE)

knots <- quantile(x)
# usamos la función create.bspline.basis para crear la base
base <- create.bspline.basis(
  norder = 4, # polinomios cúbicos
  breaks = knots # nodos en los cuartiles de x
)
plot(base, lty = "solid")
H <- eval.basis(x, base)
head(H)
beta_hat <- as.vector(solve(t(H) %*% H) %*% t(H) %*% toy$y)
beta_hat

# creamos una función que calcula mu(x)
mu <- function(x, betas){
  as.numeric(betas %*% t(eval.basis(x, base)))
}

ggplot(toy, aes(x = x, y = y)) +
  geom_point(color = "red") + 
  stat_function(fun = mu, arg = list(betas = beta_hat), color = "blue") +
  labs(title = "B-splines")


splinesBoot <- function(){
  toy_boot <- sample_n(toy, size = n, replace = TRUE)
  H <- eval.basis(toy_boot$x, base)
  as.vector(solve(t(H) %*% H) %*% t(H) %*% toy_boot$y)
}
sims_splines <- rdply(2000, splinesBoot())
betas <- as.matrix(sims_splines[, -1]) # quitamos la columna de ínidices

splines_boot <- ggplot(toy, aes(x = x, y = y)) 
for(i in 1:100){
  splines_boot <- splines_boot + 
    stat_function(fun = mu, arg = list(betas = betas[i, ]), alpha = 0.1) 
}

splines_boot + geom_point(color = "red")

amex<-2162
bbva<-1934
deb<-3817
