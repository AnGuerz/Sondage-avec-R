# III. Les types de sondage 

## 3.1. Sondage aléatoire simple (SAS) 

library(sampling)
srswr(4,8) #Plan simple avec remise
srswor(4,8) #Plan simple sans remise
srswor1(4,8) #Plan simple sans remise

#Replace permet de dire s'il y aura remise ou pas dans le tirage. C'est un
#logical qui prend TRUE quand on veut un tirage avec remise et FALSE pour un
#tirage sans remise.

sample(8,4) #Par défaut, replace est défini comme étant FALSE
sample(8,4,replace = TRUE)

## 3.2. Sondage à probabilités inégales

x = 1:9
x
n = 4
N=length(x)
pik=inclusionprobabilities(x,n)
pik
UPmultinomial(pik = pik) #plan à probabilités inégales de taille fixe avec remise
UPpoisson(pik = pik) #plan  à probabilités inégales de taille aléatoire sans remise
UPbrewer(pik = pik) #plan à probabilités inégales, de taille fixe
UPsystematic(pik = pik) #plan à probabilités inégales, de taille fixe
s=UPsystematic(pik = pik)
(1:N)[s==1]
getdata(x,s)

## 3.3. Sondage aléatoire stratifié

data=rbind(matrix(rep("nc",165),165,1,byrow=TRUE),matrix(rep("sc",70),70,1,byrow=TRUE))
data=cbind.data.frame(data,c(rep(1,100), rep(2,50), rep(3,15), rep(1,30),rep(2,40)),
                      1000*runif(235))
names(data)=c("state","region","income")
View(data)
table(data$state,data$region)
s = strata(data, stratanames = c("region","state"), 
           size = c(10,5,10,4,6), method = "srswor")
s1 = strata(data, stratanames = c("region","state"), 
            size = c(5,5,5,2,3), method = "srswr")
s2 = strata(data, stratanames = c("region","state"), 
            size = c(10,5,10,4,6), method = "poisson", pik = data$income)
s3 = strata(data, stratanames = c("region","state"), 
            size = c(10,5,10,4,6), method = "systematic", pik = data$income)
print(s)
print(s1)
print(s2)
print(s3)
sample_sastr = getdata(data,s)
print(sample_sastr)

## 3.4. Le sondage par grappe

data=rbind(matrix(rep("nc",165),165,1,byrow=TRUE),
           matrix(rep("sc",70),70,1,byrow=TRUE))
data=cbind.data.frame(data,c(rep(1,100), rep(2,50), rep(3,15), rep(1,30),rep(2,40)),
                      1000*runif(235))
names(data)=c("state","region","income")
View(data)
cl = cluster(data, c("state"), size = 2, method = "srswor")
cl1 = cluster(data, c("state"), size = 2, method = "srswr")
cl2 = cluster(data, c("state"), size = 2, 
              method = "poisson", pik = data$income)
cl3 = cluster(data, c("state"), size = 2, 
              method = "systematic", pik = data$income)
print(cl)
getdata(data,cl)

## 3.5. Le Sondage à plusieurs degrés

data=rbind(matrix(rep("nc",165),165,1,byrow=TRUE),
           matrix(rep("sc",70),70,1,byrow=TRUE))
data=cbind.data.frame(data,c(rep(1,100), rep(2,50), rep(3,15), rep(1,30),rep(2,40)),
                      1000*runif(235))
names(data)=c("state","region","income")
data1=data[order(data$state,data$region),]
table(data1$state,data1$region)
View(data1)
m = mstage(data1, size =list(25, 10), method = list("srswor","srswor"))
View(m)
getdata(data1, m)

## 3.6. Sondage équilibré

X=cbind(c(1,1,1,1,1,1,1,1,1,1), 
        c(1.1,2.2,3.1,4.2,5.1,6.3,7.1,8.1,9.1,10), 
        c(2,3,4,6,1,2,4,5,6,4))
# probabilités d'inclusion
# taille de l'échantillon n=5
pik=c(1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2)
# sélection d'un échantillon
s=samplecube(X,pik,order=1,comment=TRUE)
print(s)

# IV. Les estimateurs et les méthodes de calcul de précision
## 4.1. L'estimateur de Horvitz-Thompson

x = 1:9
n = 4
N=length(x)
pik=inclusionprobabilities(x,n)
pik
s=UPsystematic(pik = pik)
(1:N)[s==1]
y=c(2,4,3,2)  # variable d'intérêt connue sur l'échantillon s
HTestimator(y,pik[s==1])  # l'estimateur HT du total
#ou utiliser
sum(y/pik[s==1])

## 4.2. L'estimateur par calage

# on suppose que s a été tiré
# variables de calage au niveau de s
Xs=cbind(c(1,1,1,1,1,0,0,0,0,0), c(0,0,0,0,0,1,1,1,1,1), c(1,2,3,4,5,6,7,8,9,10))
# probabilités d'inclusion au niveau de s
pik=rep(0.2,times=10)
# totaux de la population pour X
total=c(24,26,290)
# les poids g en utilisant la méthode linéaire tronquée
g=calib(Xs,d=1/pik,total,method="truncated",
        bounds=c(0.75,1.2))
# les poids g sont entre 0.75 et 1.2
g
# l'estimateur de Horvitz-Thompson de X
colSums(Xs/pik)
# l'estimateur calé de X
colSums(Xs*g/pik)
#check the calibration
checkcalibration(Xs, d=1/pik, total, g)
# variable d'intérêt connue sur l'échantillon s
ys=c(3,4,5,6,1,2,4,5,2,1)
# l'estimateur calé de Y
sum(ys*g/pik)

############# Pratique avec la base swissmunicipalities ###############

data("swissmunicipalities")
View(swissmunicipalities)
data = data.frame(swissmunicipalities)

#Echantillon:
##On détermine la taille optimale de l'échantillon à partir des théories du sondage.
#Supposons que cette taille est n = 500
n = 500
N = 2896 #Nombre d'observations de la base de données

# Echantillon tiré à partir du SAS

##Plan sans remise
sampling_SAS_wor = srswor(n, N)
print(sampling_SAS_wor)
sampling_SAS_wor = getdata(data, sampling_SAS_wor)
View(sampling_SAS_wor)
## Plan avec remise
sampling_SAS_wr = srswr(n, N)
print(sampling_SAS_wr)
sampling_SAS_wr = getdata(data, sampling_SAS_wr)
View(sampling_SAS_wr)

# Echantillon tiré à partir du Sondage à probabilités inégales

#On choisit en premier la variable de stratification. Ici prenons la variable REG
REG = data$REG
pik_REG = inclusionprobabilities(REG, n)
# Tirage systématique
sampling_SPI_sys = UPsystematic(pik = pik_REG)
print(sampling_SPI_sys)
sampling_SPI_sys = getdata(data, sampling_SPI_sys)
View(sampling_SPI_sys)

# Echantillon tiré à partir du Sondage aléatoire stratifié

## On veut former des strates à partir des variables CT et POPTOT
POPTOT = data$POPTOT
POPTOT_Class = list()
for (i in 1:length(POPTOT)){
  if (POPTOT[i]<=1000){
    POPTOT_Class = append(POPTOT_Class,1)
  }else if (POPTOT[i]<=3000){
    POPTOT_Class = append(POPTOT_Class,2)
  }else {
    POPTOT_Class = append(POPTOT_Class,3)
  }
}
POPTOT_Class = as.data.frame(POPTOT_Class)
POPTOT_Class = t(POPTOT_Class)
View(POPTOT_Class)
data_strata = cbind(data,POPTOT_Class)
View(data_strata)

table(data_strata$REG,data_strata$POPTOT_Class)

sampling_SAStr = strata(data_strata, stratanames = c("REG","POPTOT_Class"), 
           size = rep(15, times=21), method = "srswor")
sampling_SAStr = getdata(data_strata, sampling_SAStr)

sampling_SAStr1 = strata(data_strata, stratanames = c("REG","POPTOT_Class"), 
                         size = rep(15, times=21), method = "srswr")
sampling_SAStr1 = getdata(data_strata, sampling_SAStr1)

sampling_SAStr2 = strata(data_strata, stratanames = c("REG","POPTOT_Class"), 
                         size = rep(15, times=21), method = "poisson", pik = data_strata$CT)
sampling_SAStr2 = getdata(data_strata, sampling_SAStr2)

sampling_SAStr3 = strata(data_strata, stratanames = c("REG","POPTOT_Class"), 
                         size = rep(15, times=21), method = "systematic", pik = data_strata$CT)
sampling_SAStr3 = getdata(data_strata, sampling_SAStr3)

## Echantillon à  partir du sondage par grappe
# On veut effectuer un sondage par grappe sur la variable REG

cl = cluster(data, c("REG"), size = 7, method = "srswor")
cl1 = cluster(data, c("REG"), size = 7, method = "srswr")
cl2 = cluster(data, c("REG"), size = 7, 
              method = "poisson", pik = data$income)
cl3 = cluster(data, c("REG"), size = 7, 
              method = "systematic", pik = data$income)
getdata(data,cl)
getdata(data,cl1)
getdata(data,cl2)
getdata(data,cl3)

## Echantillon à partir du sondage à plusieurs degrés

## On veut effectuer un sondage à plusieurs dégrés sur les variables REG et CT avec des tailles de 500 et 250 suivant les degrés:

m = mstage(data, size =list(500, 250), method = list("srswor","srswor"))
getdata(data, m)

## Estimateurs de HT

# On considère l'échantillon tiré à partir sondage à probabilités inégales méthode systématique
# On souhaite estimer la population totale
y = rep(seq(from = 1100, to = 2000, by = 100),50)
length(y)
pik=inclusionprobabilities(data$POPTOT,n)
s = UPsystematic(pik = pik)
HTestimator(y,pik[s==1])

## Estimateurs par le calage
data_cal = cbind(data$CT,data$REG,data$COM,data$HApoly,data$Surfacesbois,data$Surfacescult,data$POPTOT)
colnames(data_cal) = c("CT", "REG","COM","HApoly","Surfacesbois","Surfacescult","POPTOT")

View(data_cal)
Xs = getdata(data_cal,s)
Xs = Xs[,-1]
Xs = as.matrix(Xs)
total = c(41202, 9248, 9934534, 3998831, 1270996, 987317, 7288010)
View(Xs)
g=calib(Xs,d=1/pik,total,method="truncated",
        bounds=c(0.75,1.2))





