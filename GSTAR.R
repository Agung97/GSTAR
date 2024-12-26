## package ##
library(dplyr)
library(timeSeries)
library(gstar)
library(tseries)
library(MASS)
library(pracma)
library(forecast)
library(Metrics)

## persiapan data ##
data = read.delim("clipboard")
data


## statistika deskriptif ##
summary(data)


## membagi data training (80) dan testing (20)
s <- round(nrow(data) * 0.8) ## split into training and testing (80:20)
data_train <- data[1:s, ]
data_train
data_test <- data[-c(1:s), ]
data_test


## korelasi ##
cor.test(data_train$Jabar, data_train$Jateng, alternative = "two.sided", method = "pearson", conf.level = 0.95)
cor.test(data_train$Jabar, data_train$Jatim, alternative = "two.sided", method = "pearson", conf.level = 0.95)
cor.test(data_train$Jateng, data_train$Jatim, alternative = "two.sided", method = "pearson", conf.level = 0.95)


## cek kestasioneran data dalam rataan ##
adf.test(data_train$Jabar)
adf.test(data_train$Jateng)
adf.test(data_train$Jatim)

## membuat plot acf ##
acf(data_train$Jabar)
acf(data_train$Jateng)
acf(data_train$Jatim)


## proses differencing ##
Data1=diff(data_train$Jabar)
Data2=diff(data_train$Jateng)
Data3=diff(data_train$Jatim)


## cek stasioneritas data setelah diff ##
adf.test(Data1)
adf.test(Data2)
adf.test(Data3)


## membuat plot acf setelah diff ##
acf(Data1)
acf(Data2)
acf(Data3)

data.diff1 = data.frame(Data1, Data2, Data3)
data.diff1
write.table(data.diff1, "D:/data diff.csv", sep = ",",row.names = FALSE)


## cek kestasioneran data dalam ragam ##

index <- seq(1:197)
bc1 = boxcox(data_train$Jabar~index, lambda = seq(0,4,by=0.01))
lambda1 <- bc1$x[which.max(bc1$y)]
lambda1
bc1$x[bc1$y > max(bc1$y) - 1/2 * qchisq(.95,1)]


bc2 = boxcox(data_train$Jateng~index, lambda = seq(0,4,by=0.01))
lambda2 <- bc2$x[which.max(bc2$y)]
lambda2
bc2$x[bc2$y > max(bc2$y) - 1/2 * qchisq(.95,1)]


bc3 = boxcox(data_train$Jatim~index, lambda = seq(0,4,by=0.01))
lambda3 <- bc3$x[which.max(bc3$y)]
lambda3
bc3$x[bc3$y > max(bc3$y) - 1/2 * qchisq(.95,1)]

## transformasi logaritma ##
new1 <- log(data_train$Jabar)
bc4 = boxcox(new1~index, lambda = seq(0,4,by=0.01))
lambda4 <- bc4$x[which.max(bc4$y)]
lambda4
bc4$x[bc4$y > max(bc4$y) - 1/2 * qchisq(.95,1)]

new2 <- log(new1)
bc5 = boxcox(new2~index, lambda = seq(0,4,by=0.01))
lambda5 <- bc5$x[which.max(bc5$y)]
lambda5
bc5$x[bc5$y > max(bc5$y) - 1/2 * qchisq(.95,1)]

new3 <- log(data_train$Jateng)
bc6 = boxcox(new3~index, lambda = seq(0,4,by=0.01))
lambda5 <- bc6$x[which.max(bc6$y)]
lambda5
bc6$x[bc6$y > max(bc6$y) - 1/2 * qchisq(.95,1)]


## penentuan bobot lokasi seragam ##
bobot = matrix(c(0,1,1,
                 1,0,1,
                 1,1,0), ncol=3, nrow=3)
bobot.srgm = bobot/(ncol(data_train)-1)
bobot.srgm


## penentuan bobot lokasi invers jarak ##
Jabar = '6 51 36S, 107 36 00E'
Jateng = '6 59 24S, 110 25 21E'
Jatim = '7 14 45S, 112 44 16E'

Jrk12 = haversine(Jabar, Jateng)
Jrk12

Jrk13 = haversine(Jabar, Jatim)
Jrk13

Jrk23 = haversine(Jateng, Jatim)
Jrk23

bobot.inv = matrix(c(0,0.451794873, 0.311327646,
                     0.645770338, 0, 0.688672354,
                     0.354229662, 0.548205127,0),3,3)
bobot.inv


## penentuan bobot lokasi normalisasi korelasi silang ##
z1=data_train$Jabar
z2=data_train$Jateng
z3=data_train$Jatim

sz1z2 = sum((z1[2:197]-mean(z1))*(z2[1:196]-mean(z2)))
sz1z3 = sum((z1[2:197]-mean(z1))*(z3[1:196]-mean(z3)))
sz2z3 = sum((z2[2:197]-mean(z2))*(z3[1:196]-mean(z3)))

z1kuadrat = sum((z1-mean(z1))^2)
z2kuadrat = sum((z2-mean(z2))^2)
z3kuadrat = sum((z3-mean(z3))^2)

uz1z2 = sqrt(z1kuadrat*z2kuadrat)
uz1z3 = sqrt(z1kuadrat*z3kuadrat)
uz2z3 = sqrt(z2kuadrat*z3kuadrat)

rz1z2 = sz1z2/uz1z2
rz1z3 = sz1z3/uz1z3
rz2z3 = sz2z3/uz2z3

r <-matrix(c(1,rz1z2,rz1z3,rz1z2,1,rz2z3,rz1z3,rz2z3,1),3,3)
bobot.kor <- matrix(c(0,abs(rz1z2)/(abs(rz1z2)+abs(rz2z3)),abs(rz1z3)/(abs(rz1z3)+abs(rz2z3)),
               abs(rz1z2)/(abs(rz1z2)+abs(rz1z3)),0,abs(rz2z3)/(abs(rz1z3)+abs(rz2z3)),
               abs(rz1z3)/(abs(rz1z3)+abs(rz1z2)),abs(rz2z3)/(abs(rz1z2)+abs(rz2z3)),0),3,3)
bobot.kor



## estimasi parameter bobot seragam ##
fit.srgm = gstar(data_train, weight = bobot.srgm, p=2, est = "OLS")
summary(fit.srgm)

seragam = matrix(c(-0.005142, 0,0, 
                   0, 0.051604,0,
                   0,0,-0.006065),3,3)
seragam%*%bobot.srgm

seragam1 = matrix(c(0.011539, 0,0, 
                   0,0.090587,0,
                   0,0, 0.061642),3,3)
seragam1%*%bobot.srgm

write.table(fit.srgm$fitted_values, "D:/fitted values seragam.csv", sep = ",",row.names = FALSE)

## estimasi parameter bobot invers jarak ##
fit.inv = gstar(data_train, weight = bobot.inv, p=2, est = "OLS")
summary(fit.inv)
invers = matrix(c(-0.007924, 0,0, 
                   0,0.034025,0,
                   0,0,0.005198),3,3)
invers%*%bobot.inv

invers1 = matrix(c( 0.014453, 0,0, 
                  0,0.072698 ,0,
                  0,0,0.047299),3,3)
invers1%*%bobot.inv

write.table(fit.inv$fitted_values, "D:/fitted values invers.csv", sep = ",",row.names = FALSE)

## estimasi parameter bobot invers jarak ##
fit.kor = gstar(data_train, weight = bobot.kor, p=2, est = "OLS")
summary(fit.kor)
korelasi = matrix(c(-0.005022,0,0, 
                  0,0.051213,0,
                  0,0,-0.006547),3,3)
korelasi%*%bobot.kor

korelasi1 = matrix(c(0.011556, 0,0, 
                    0,0.094950,0,
                    0,0,0.061015),3,3)
korelasi1%*%bobot.kor

write.table(fit.kor$fitted_values, "D:/fitted values korelasi2.csv", sep = ",",row.names = FALSE)


## uji kelayakan model ##
res.srgm <- read.delim("clipboard")
Box.test(res.srgm$Jabar, type = "Ljung-Box")
Box.test(res.srgm$Jateng, type = "Ljung-Box")
Box.test(res.srgm$Jatim, type = "Ljung-Box")

res.inv <- read.delim("clipboard")
Box.test(res.inv$Jabar, type = "Ljung-Box")
Box.test(res.inv$Jateng, type = "Ljung-Box")
Box.test(res.inv$Jatim, type = "Ljung-Box")

res.kor <- read.delim("clipboard")
Box.test(res.kor$Jabar, type = "Ljung-Box")
Box.test(res.kor$Jateng, type = "Ljung-Box")
Box.test(res.kor$Jatim, type = "Ljung-Box")


## Evaluasi model terbaik ##
performance(fit.srgm)
performance(fit.inv)
performance(fit.kor)


## prediksi ##
prediksi <- predict(fit.kor, n=5)
prediksi
write.table(prediksi, "D:/prediksi2.csv", sep = ",",row.names = FALSE)

dataakhir <- read.delim("clipboard")
dataakhir$pred_jabar <- as.numeric(dataakhir$pred_jabar)
dataakhir$pred_jatim <- as.numeric(dataakhir$pred_jatim)
dataakhir$test_jabar <- as.numeric(dataakhir$test_jabar)
dataakhir$test_jatim <- as.numeric(dataakhir$test_jatim)


##  akurasi peramalan ##
n=5
MAPEjabar <- 1/n*sum(abs((dataakhir$test_jabar - dataakhir$pred_jabar)/dataakhir$test_jabar))*100
MAPEjabar
MAPEjateng <- 1/n*sum(abs((dataakhir$test_jateng - dataakhir$pred_jateng)/dataakhir$test_jateng))*100
MAPEjateng
MAPEjatim <- 1/n*sum(abs((dataakhir$test_jatim - dataakhir$pred_jatim)/dataakhir$test_jatim))*100
MAPEjatim


