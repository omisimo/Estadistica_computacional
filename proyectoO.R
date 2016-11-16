library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(boot)
library(bootstrap)
library(xlsx)

###### problema 1 parte 1 #####
datosaretes <- read.csv("C:/Users/Administrador/Desktop/parcial est comp/datosaretes.csv")
arete <- read.table("C:/Users/Administrador/Desktop/parcial est comp/arete.csv", header=T, quote="\"")
arete = as.data.frame(arete)
nopar<-function(){
  muestra <- sample_n(as.data.frame(arete),size=8,replace=FALSE)
  indices = as.numeric(rownames(muestra))%%10
  par=length(unique(indices))<length(indices)
  return(par)
}


E=50000
nop<-rdply(E,nopar)
1-(sum(nop$V1)/E)
#.0919

##### problema 1 parte 2 #### 

expar<-function(){
  muestra <- sample_n(arete,size=8,replace=FALSE)
  indices = as.numeric(rownames(muestra))%%10
  par=length(unique(indices))==(length(indices)-1)
  return(par)
}

E=50000
nop<-rdply(E,expar)
(sum(nop$V1)/E)
#.43004




###### problema 2 parte 1 ####

#### p(t4nA)=0.1152
#### p(t4nB)=0.2592
#### .3744


juego<-function(){
  A<-0
  B<-0
  k<-0
  
  acabo=FALSE
  while(!acabo){
  u=runif(1)
  if(u<.4){
    A=A+1
  }else{
    B=B+1
  }
  
  if(A==B+2 | B==A+2){
    acabo=TRUE
  }
  k<-k+1
}
return(k)  
}

E=50000
juada<-rdply(E,juego)
length(which(juada$V1==4))/E

#### problema 2 inciso 2 ####

juego_a<-function(){
  A<-0
  B<-0
  k<-0
  
  acabo=FALSE
  while(!acabo){
    u=runif(1)
    if(u<.4){
      A=A+1
    }else{
      B=B+1
    }
    
    if(A==B+2){
      k=1
      acabo=TRUE
    }else{
      if(B==A+2){
        acabo=TRUE
        K=0
      }
    }
  }
  return(k)  
}

E=1000000
juada<-rdply(E,juego_a)
sum(juada$V1)/E



#### problema 3 #####
amis <- read.csv("C:/Users/Administrador/Desktop/parcial est comp/amis.csv")

head(amis)

##### a)
# La muestra no es aleatoria, debido a que la base de datos
# se obtuvo a través de un diseño de experimentos.
# La variable period indica si se hizo la medición antes y
# después, para la variable warning se tienen valores 1 y 2 
# para referir a los valores de carreteras y la variable period
# corresponde a las 11 carreteras. Por cada combinacion de medición 
# (antes y después) para cada warning (1 y 2) se tomaron 100 datos
# por lo que se tienen 400 mediciones de velocidad por cada 
# carretera. Hay 11 carreteras, por lo cual el total de datos es
# de 11*2*2*100 = 4400.


##### b)

aux11 <- subset(amis, amis$warning == 1 & amis$period ==2)
n11 <- quantile(aux11$speed,0.85)

aux21 <- subset(amis, amis$warning == 2 & amis$period ==2)
n21 <- quantile(aux21$speed,0.85)

aux12 <- subset(amis, amis$warning == 1 & amis$period ==1)
n12 <- quantile(aux12$speed,0.85)

aux22 <- subset(amis, amis$warning == 2 & amis$period ==1)
n22 <- quantile(aux22$speed,0.85)

eta = (1/11)*((as.numeric(n11) - as.numeric(n21)) - (as.numeric(n12) - as.numeric(n22)))

# eta es el estimador plug-in


##### c)

boot <- function(aux11,aux21,aux12,aux22){ 
  n <- nrow(aux11)
  muestra_boot_aux11 <- sample_n(aux11, size = n, replace = TRUE)
  muestra_boot_aux21 <- sample_n(aux21, size = n, replace = TRUE)
  muestra_boot_aux12 <- sample_n(aux12, size = n, replace = TRUE)
  muestra_boot_aux22 <- sample_n(aux22, size = n, replace = TRUE)
  
  n11 <- quantile(muestra_boot_aux11$speed,0.85)
  n21 <- quantile(muestra_boot_aux21$speed,0.85)
  n12 <- quantile(muestra_boot_aux12$speed,0.85)
  n22 <- quantile(muestra_boot_aux22$speed,0.85)
  
  eta = (1/11)*((as.numeric(n11) - as.numeric(n21)) - (as.numeric(n12) - as.numeric(n22)))
  
  return(eta)
  
}

B=3000
thetas_boot <- rdply(B, boot(aux11,aux21,aux12,aux22))
hist(thetas_boot[,2],20)

ggplot(thetas_boot, aes(x = V1)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.015, fill = "darkgray") + 
  geom_vline(aes(xintercept = mean(V1), color = "red"))

##### d)

li_normal <- round(eta - 1.96 * sd(thetas_boot$V1), 2)
ls_normal <- round(eta + 1.96 * sd(thetas_boot$V1), 2)
c(li_normal,ls_normal) #intervalos normales


# percentiles
ls_per <- round(quantile((thetas_boot$V1), prob = 0.975), 2)
li_per <- round(quantile((thetas_boot$V1), prob = 0.025), 2)
c(li_per,ls_per)


eta
c(li_normal,ls_normal) #intervalos normales
c(li_per,ls_per)


##### problema 4 #####

#### a) la suma de geometricas se distribuye binomial negativo

#### b) 
set.seed(221285)

geom<-function(p){
  g<-floor((log(runif(1))/log(1-p)))+1
  return(g)
}
geometricas<-rdply(100000,geom(.4))
ggplot(geometricas,aes(x=V1))+
  geom_histogram(aes(y=..density..),binwidth=1,fill="darkgray")



binomialN<-function(p,r){
  u<-runif(r)
  g<-floor(log(u)/log(1-p))+1
  return(sum(g))
}

# primeros 10 valores de la binomial negativa
set.seed(221285)
r=20
p=0.4
sim1<-rdply(10,binomialN(p,r))$V1
sim1





##### c)

set.seed(221285)
binoneg=function(p,r) {
u<-runif(1)
j<-20
x<-j
prob<-p**r
acum<-prob
while(u>acum){
  prob<-j/(j-r+1)*(1-p)*prob
  acum<-acum+prob
  x<-j
  j<-j+1
}
return(x)
}

set.seed(221285)
sim2<-rdply(10,binoneg(p=0.4,r=20))$V1
sim2

set.seed(221285)
rnbinom(10,r,p)





##### d)

set.seed(221285)
m1<-system.time({
  rdply(10000,binomialN(p,r))$V1
})
set.seed(221285)
m2<-system.time({
  rdply(10000,binoneg(0.4,20))$V1
})


m1 #tiempo con el metodo 1
m2 #tiempo con el metodo 2


#es mas lento el segundo metodo

##### e)
set.seed(221285)
m1<-rdply(1000,binomialN(0.4,20))$V1
set.seed(221285)
m2<- rdply(1000,binoneg(0.4,20))$V1
set.seed(221285)
rbin<-rnbinom(1000,size=20,prob=0.4)+20

#devtools::install_github("hadley/tidyr")
datos<-data.frame(Metodo1=m1,R=rbin)
datos1<-gather(datos,metodo,var)
ggplot(datos1, aes(var, fill = metodo)) + geom_histogram(alpha = 0.7, aes(y = ..density..), position = 'identity')

datos<-data.frame(Metodo2=m2,R=rbin)
datos1<-gather(datos,metodo,var)
ggplot(datos1, aes(var, fill = metodo)) + geom_histogram(alpha = 0.7, aes(y = ..density..), position = 'identity')



##### problema 5 #####

set.seed(221285)
muestra<-rnorm(10)

theta<-exp(mean(muestra))

thetaboot<-function(x){
  n<-length(x)
  x_boot <- sample(x, n, replace = TRUE)
  return(exp(mean(x_boot)))
}

B <- 6000
theta_b <- rdply(B, thetaboot(muestra))
li_normal <- round(theta - 1.96 * sd(theta_b$V1), 2)
ls_normal <- round(theta + 1.96 * sd(theta_b$V1), 2)
c(li_normal,ls_normal) #intervalos normales


ggplot(theta_b, aes(x = V1)) + 
  geom_histogram(binwidth = 0.05, fill = "gray30") +
  geom_vline(xintercept = c(li_normal, ls_normal, theta), 
             color = c("black", "black", "red"))


# percentiles
ls_per <- round(quantile((theta_b$V1), prob = 0.975), 2)
li_per <- round(quantile((theta_b$V1), prob = 0.025), 2)
c(li_per,ls_per)

#bca
# definimos la funcion necesaria para utilizar boot
set.seed(221285)
theta_est<-function(data,indices){
  d<-data[indices,]
  theta<-exp(mean(d))
return(theta)
}

#la funcion boot usa la funcion anterior en la que se definio el estadistico de interes
# para el cual se desea conocer el intervalo de confianza
theta_boot<-boot(data=as.data.frame(muestra),statistic=theta_est,R=6000)
bca<-boot.ci(theta_boot,statistic=theta_est,type="bca",index=1)
bca_ci<-round(bca$bca[4:5],2)

c(li_normal,ls_normal)
c(li_per,ls_per)
bca_ci
# en los tres casos se encuentra el valor real theta=1
# en la normal el interavlo no es bueno pues la distribuicon es asimetrica
# mientras que los de percentiles y bca son similares debido a que el primera
# considera la distribucion acumulada unicamente y el segundo realiza una 
# correccion por sesgo

# los intervalos de confianza normales podrian mejorarse si se aplicara la funcion inversa
# que normalice la distribucion del estadistico



# proceso para comparar los tres intervalos de confianza estimados, en el cual se pretende observar
# si incluyen al valor real, si hay fallo por la izquierda o fallo por la derecha


n<-500
fi_n<-0
fi_p<-0
fi_bca<-0
fd_n<-0
fd_p<-0
fd_bca<-0
cob_n<-0
cob_p<-0
cob_bca<-0
li_normal<-rep(0,n)
ls_normal<-rep(0,n)
li_per<-rep(0,n)
ls_per<-rep(0,n)
li_bca<-rep(0,n)
ls_bca<-rep(0,n)
theta<-rep(0,n)

set.seed(221285)

for(i in 1:n){
  # para medias

  muestra<-rnorm(10)
  theta[i]<-exp(mean(muestra))
  theta_b <- rdply(6000, thetaboot(muestra))
  
  # se calcula el intervalo de confiaza normal
  li_normal[i] <- theta[i] - 1.96 * sd(theta_b$V1)
  ls_normal[i] <- theta[i] + 1.96 * sd(theta_b$V1)
  
  # se calcula el intervalo de confianza por percentiles
  ls_per[i] <- as.numeric(quantile((theta_b$V1), prob = 0.975))
  li_per[i] <- as.numeric(quantile((theta_b$V1), prob = 0.025))
  
  # se calcula el intervalo de confianza por bca
  theta_boot<-boot(data=as.data.frame(muestra),statistic=theta_est,R=6000)
  bca<-boot.ci(theta_boot,statistic=theta_est,type="bca",index=1)
  bca_ci<-bca$bca[4:5]
  li_bca[i]<-bca_ci[1]
  ls_bca[i]<-bca_ci[2]  
  
}


datos_ic<-data.frame(theta,li_normal,ls_normal,li_per,ls_per,li_bca,ls_bca)

cob_normal<-sum(datos_ic$li_normal<=1 & datos_ic$ls_normal>=1)
cob_per<-sum(datos_ic$li_per<=1 & datos_ic$ls_per>=1)
cob_bca<-sum(datos_ic$li_bca<=1 & datos_ic$ls_bca>=1)

fi_normal<-sum(datos_ic$li_normal>1)
fi_per<-sum(datos_ic$li_per>1)
fi_bca<-sum(datos_ic$li_bca>1)

fd_normal<-sum(datos_ic$ls_normal<1)
fd_per<-sum(datos_ic$ls_per<1)
fd_bca<-sum(datos_ic$ls_bca<1)


tabla<-data.frame(fizq=c(fi_normal,fi_per,fi_bca),fder=c(fd_normal,fd_per,fd_bca),cob=c(cob_normal,cob_per,cob_bca))
tabla_p<-(tabla/n)*100

#write.csv(tabla,"tabla.csv")
#write.csv(tabla_p,"tabla_p.csv")

##### b)


coef_df<- datos_ic %>%
  gather(metodo,limite,-theta) %>%
  separate(metodo,c("int","metodo"),3) %>%
  mutate(lim=substr(int,1,2))%>%
  select(-int) %>%
  spread(lim,limite) %>%
  mutate(ind=rep(1:n, each=3)) 



ggplot(coef_df, aes(x = ind, y = theta, ymin = li, 
                    ymax= ls)) +
  geom_pointrange(size = 0.4) + 
  facet_wrap(~ metodo, scales = "free_y", nrow = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))+
  geom_hline(yintercept=1,color="red")





##### c)
##### C) parte a)

set.seed(221285)
muestra<-rpois(20,2.5)
theta<-exp(-2*mean(muestra))

thetaboot<-function(x){
  n<-length(x)
  x_boot <- sample(x, n, replace = TRUE)
  return(exp(-2*mean(x_boot)))
}

B <- 6000
theta_b <- rdply(B, thetaboot(muestra))
li_normal <-theta - 1.96 * sd(theta_b$V1)
ls_normal <-theta + 1.96 * sd(theta_b$V1)
c(li_normal,ls_normal) #intervalos normales


ggplot(theta_b, aes(x = V1)) + 
  geom_histogram(binwidth = 0.005, fill = "gray30") +
  geom_vline(xintercept = c(li_normal, ls_normal, theta), 
             color = c("black", "black", "red"))


# percentiles
ls_per <- quantile((theta_b$V1), prob = 0.975)
li_per <- quantile((theta_b$V1), prob = 0.025)
c(li_per,ls_per)

#bca
# definimos la funcion necesaria para utilizar boot
set.seed(221285)
theta_est<-function(data,indices){
  d<-data[indices,]
  theta<-exp(-2*mean(d))
  return(theta)
}

#la funcion boot usa la funcion anterior en la que se definio el estadistico de interes
# para el cual se desea conocer el intervalo de confianza
theta_boot<-boot(data=as.data.frame(muestra),statistic=theta_est,R=6000)
bca<-boot.ci(theta_boot,statistic=theta_est,type="bca",index=1)
bca_ci<-bca$bca[4:5]

c(li_normal,ls_normal)
c(li_per,ls_per)
bca_ci



n<-500
fi_n<-0
fi_p<-0
fi_bca<-0
fd_n<-0
fd_p<-0
fd_bca<-0
cob_n<-0
cob_p<-0
cob_bca<-0
li_normal<-rep(0,n)
ls_normal<-rep(0,n)
li_per<-rep(0,n)
ls_per<-rep(0,n)
li_bca<-rep(0,n)
ls_bca<-rep(0,n)
theta<-rep(0,n)
set.seed(221285)
for(i in 1:n){
  # para medias
  muestra<-rpois(20,2.5)
  theta[i]<-exp(-2*mean(muestra))
  theta_b <- rdply(6000, thetaboot(muestra))
  
  # se calcula el intervalo de confiaza normal
  li_normal[i] <- theta[i] - 1.96 * sd(theta_b$V1)
  ls_normal[i] <- theta[i] + 1.96 * sd(theta_b$V1)
  
  # se calcula el intervalo de confianza por percentiles
  ls_per[i] <- as.numeric(quantile((theta_b$V1), prob = 0.975))
  li_per[i] <- as.numeric(quantile((theta_b$V1), prob = 0.025))
  
  # se calcula el intervalo de confianza por bca
  theta_boot<-boot(data=as.data.frame(muestra),statistic=theta_est,R=6000)
  bca<-boot.ci(theta_boot,statistic=theta_est,type="bca",index=1)
  bca_ci<-bca$bca[4:5]
  li_bca[i]<-bca_ci[1]
  ls_bca[i]<-bca_ci[2]  
  
}


datos_ic<-data.frame(theta,li_normal,ls_normal,li_per,ls_per,li_bca,ls_bca)

real<-exp(-2*2.5)

cob_normal<-sum(datos_ic$li_normal<=real & datos_ic$ls_normal>=real)
cob_per<-sum(datos_ic$li_per<=real & datos_ic$ls_per>=real)
cob_bca<-sum(datos_ic$li_bca<=real & datos_ic$ls_bca>=real)

fi_normal<-sum(datos_ic$li_normal>real)
fi_per<-sum(datos_ic$li_per>real)
fi_bca<-sum(datos_ic$li_bca>real)

fd_normal<-sum(datos_ic$ls_normal<real)
fd_per<-sum(datos_ic$ls_per<real)
fd_bca<-sum(datos_ic$ls_bca<real)


tabla<-data.frame(fizq=c(fi_normal,fi_per,fi_bca),fder=c(fd_normal,fd_per,fd_bca),cob=c(cob_normal,cob_per,cob_bca))
tabla_p<-(tabla/n)*100

#write.csv(tabla,"tabla_pois.csv")
#write.csv(tabla_p,"tabla_p_pois.csv")

coef_df<- datos_ic %>%
  gather(metodo,limite,-theta) %>%
  separate(metodo,c("int","metodo"),3) %>%
  mutate(lim=substr(int,1,2))%>%
  select(-int) %>%
  mutate(ind=rep(1:n,times=6)) %>%
  spread(lim,limite) %>%
  arrange(metodo) %>%
  select(-ind) %>%
  mutate(ind=rep(1:n, times=3))
  
  
ggplot(coef_df, aes(x = ind, y = theta, ymin = li, 
                    ymax= ls)) +
  geom_pointrange(size = 0.4) + 
  facet_wrap(~ metodo, scales = "free_y", nrow = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  geom_hline(yintercept=real,color="red")

##### problema 6 #####

datos <- read.xlsx("Cuadros_2010.xlsx",sheetName="PEA",header=TRUE,colClasses=c("character",rep("double",14)))
datos<-datos[4:36,1:11]
datos<-datos[complete.cases(datos),]
rownames(datos) <- 1:33

datos.pea <- data.frame(datos[1], datos[3], datos[5], datos[8], datos[10])
colnames(datos.pea) <- c('Estado','act_mig','ina_mig','act_nomig',
                         'ina_nomig')



pea_df<- datos.pea %>%
  gather(tipo,valor,-Estado) %>%
  separate(tipo,c("migrante","Estado_migra"),4) %>%
  mutate(Estado_migracion=substr(migrante,1,3)) %>%
  select(-migrante) %>%
  spread(Estado_migra,valor)




ggplot(pea_df, aes(x=valor, fill=Estado_migra)) + 
  geom_bar() +
  facet_wrap(~Estado_migracion ) + 
  xlab('Porcentaje relativo de poblaci�n econ�micamente activa')

qplot(x=Estado, y=valor, fill=Estado_migracion,data=pea_df, 
      geom="bar", stat="identity",position="identity",alpha=0.3) + 
  facet_wrap(~Estado_migra) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(pea_df,aes(x=Estado,y=valor,fill=Estado_migracion)) + 
  geom_bar(stat="identity",position = "identity", alpha=.3) +
  facet_wrap(~Estado_migra) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





