
# data provided by the Argentina goverment in https://datos.gob.ar/
# Data under analysis is the production of natural gas since 1996 to 2023
# it is a time series

# read data
gn = read.csv('C:/Users/marie/Documents/Mariela/data_science/produccion_gas_natural.csv', header = TRUE)
attach(gn)
gasnatural <- ts(gn, start = c(1996, 1), freq = 12)
plot(gasnatural, col = rainbow(5), ylab="Millones de metros cúbicos", xlab="Año", main="Producción de gas natural (1996-2023)")

head(gasnatural)
tail(gasnatural)

# anual behavior
annual.ts <- aggregate(gasnatural)/12
plot(annual.ts, col = rainbow(5), ylab="Mm³", xlab="Año", main="Comportamiento anual")

# mensual behavior
gn.month.ts <- gasnatural
plot(gn.month.ts, col = rainbow(5), ylab="Mm³", xlab="Año", main="Comportamiento mensual")

#histograms
hist(gasnatural, col=gray.colors(9), main="Histograma produccion de gas natural", xlab="Mm³")
hist(gasnatural, prob = T, col=gray.colors(9), main="Histograma produccion de gas natural", xlab="Mm³")
plot(density(gasnatural), col=gray.colors(9), main="Histograma de densidad", xlab="Mm³")

# anual mean 
annual.mean <- aggregate(gasnatural, FUN = mean)
plot(annual.mean, col=rainbow(5), ylab="Mm³", xlab="Año", main="Promedio anual de producción de gas natural")

# boxplots
boxplot(gasnatural~cycle(gasnatural), col=rainbow(12), ylab="Mm³", xlab="Mes", main="Boxplot de la Producción de gas natural mensual")

# time series decompose
plot(decompose(gasnatural, type=c("additive")), col=rainbow(5), xlab="Año")
plot(decompose(gasnatural, type=c("multiplicative")), col=rainbow(5), xlab="Año")
plot(decompose(gasnatural)$trend, col=rainbow(5), xlab="Año", main="Componente trend")
plot(decompose(gasnatural)$random, col=rainbow(5), xlab="Año", main="Componente random")

plot(decompose(diff(gasnatural), type=c("additive")), col=rainbow(5), xlab="Año")
plot(decompose(diff(gasnatural))$trend, col=rainbow(5), xlab="Año", main="Componente trend aplicado a diff(gasnatural)")

# Stocastic models
# Dickey-Fuller test to verify if the series is stationary or not
library(tseries)
adf.test(gasnatural)

acf(gasnatural)
# the correlations are not significant
# the range of correlations is between -0.2 to 0.2

pacf(gasnatural)

# Serie should not be stationary, so I gonna transformate the serie
# this is because the varianza should be constance accross the time

library(tseries)
adf.test(diff(gasnatural))

acf(diff(gasnatural))
pacf(diff(gasnatural))

install.packages('forecast')
library(forecast)
auto.arima(gasnatural)
auto.arima(diff(gasnatural))

gasnatural.ar <- ar(gasnatural)
mean(gasnatural)
gasnatural.ar$order
gasnatural.ar$ar


best.order <- c(0,0,0)
best.aic <- Inf
for (i in 0:2) for (j in 0:2){
  fit.aic <- AIC(arima(gasnatural, order=c(i,0,j)))
  if (fit.aic < best.aic){
    best.order <- c(i,0,j)
    best.arma <- arima(gasnatural, order = best.order)
    best.aic <- fit.aic
  }
}

# residuals analysis
Time <- 1:length(gasnatural)
Imth <- cycle(gasnatural)
gasnatural.lm <- lm(log(gasnatural) ~ Time + I(Time^2) + factor(Imth))

pacf(residuals(gasnatural.lm)) # una cota razonable es 2

best.order <- c(0,0,0)
best.aic <- Inf
for (i in 0:2) for (j in 0:2){
  fit.aic <- AIC(arima(resid(gasnatural.lm), order=c(i,0,j)))
  if (fit.aic < best.aic){
    best.order <- c(i,0,j)
    best.arma <- arima(resid(gasnatural.lm), order = best.order)
    best.aic <- fit.aic
  }
}

acf(resid(best.arma))

new.time <- seq(length(gasnatural), length = 36)
new.data <- data.frame(Time = new.time, Imth = rep(1:12, 3))

predict.lm <- predict(gasnatural.lm, new.data)
predict.arma <- predict(best.arma, n.ahead = 36)
gasnatural.pred <- ts(exp(predict.lm + predict.arma$pred), start=2024, freq=12)

# prediction of the production of gas natural for next 3 years
ts.plot(cbind(gasnatural, gasnatural.pred), lty=1:2)

layout(1:3)
plot(gasnatural,main="Producción de gas natural.")
plot(diff(gasnatural),main="Diferencia de la producción de gas natural.")
plot(diff(log(gasnatural)),main="Diferencia del logaritmo de la producción de gas natural.")

# SARIMA
get.best.arima <- function(gasnatural, maxord = c(1,1,1,1,1,1))
{
  best.aic <- 1e8
  n <-length(gasnatural)
  for (p in 0:maxord[1]) for (d in 0:maxord[2]) for (q in 0:maxord[3])
    for (P in 0:maxord[4]) for (D in 0:maxord[5]) for (Q in 0:maxord[6])
    {
      fit <- arima(gasnatural, order = c(p,d,q),
                   seas = list(order = c(P,D,Q),
                               frequency(gasnatural)), method = "CSS")
      fit.aic <- -2 * fit$loglik + (log(n) + 1) * length(fit$coef)
      if (fit.aic < best.aic){
        best.aic <- fit.aic
        best.fit <- fit
        best.model <- c(p,d,q,P,D,Q)
      }}
  list(best.aic, best.fit, best.model)
}

best.aic

best.fit
best.model

# fitting SARIMA
best.arima.gasnatural <- get.best.arima(log(gasnatural), maxord = c(2,2,2,2,2,2))
best.fit.gasnatural <- best.arima.gasnatural[[2]]
acf(resid(best.fit.gasnatural))
best.arima.gasnatural[[3]]
ts.plot(cbind(window(gasnatural,start=2005), exp(predict(best.fit.gasnatural,12)$pred)),lty=1:2)

# possibles models to analyze and compare
#SARIMA(1,1,0)(0,1,1)12
#SARIMA(0,1,1)(0,1,1)12
#SARIMA(1,1,1)(0,1,1)12
#SARIMA(1,0,1)(0,1,1)12
auto.arima(diff(gasnatural))

# I chose the following 3 models
mod1 <- arima(diff(gasnatural),c(1,1,0),c(0,1,1))
mod1 #aic = 1458.2

mod2 <- arima(diff(gasnatural),c(0,1,1),c(0,1,1))
mod2 #aic = 1406.2

mod3 <- arima(diff(gasnatural),c(1,1,1),c(0,1,1))
mod3 #aic = 1403.7

mod4 <- arima(diff(gasnatural),c(1,0,1),c(0,1,1))
mod4 #aic = 1400.94

# the best model is the one with the lowest AIC, this is mod4
summary(mod1$residuals)
summary(mod2$residuals)
summary(mod3$residuals)
summary(mod4$residuals)

par(mfcol = c(1,4))
hist(mod1$residuals)
hist(mod2$residuals)
hist(mod3$residuals)
hist(mod4$residuals)

jarque.bera.test(mod1$residuals)
jarque.bera.test(mod2$residuals)
jarque.bera.test(mod3$residuals)
jarque.bera.test(mod4$residuals)

forecastDiff <- forecast(mod4, level= c(95), h=24)
autoplot(forecastDiff)


# *********************************************************************

acf(log(gasnatural))
acf(diff(log(gasnatural)))
pacf(diff(log(gasnatural)))

#Series: gasnatural 
#ARIMA(1,0,0)(0,1,1)[12] with drift 

#Coefficients:
#  ar1     sma1    drift
#0.7958  -0.7190  -7.0375
#s.e.  0.0602   0.1006   1.1438

#sigma^2 = 6737:  log likelihood = -702.47
#AIC=1412.94   AICc=1413.29   BIC=1424.09

(fit <- arima(log(gasnatural), c(1, 0, 0),seasonal = list(order = c(1, 0, 0), period = 12)))
#Call:
#  arima(x = log(gasnatural), order = c(1, 0, 0), seasonal = list(order = c(1, 0, 0), period = 12))
#
#Coefficients:
#  ar1    sar1  intercept
#0.7787  0.8741     8.2646
#s.e.  0.0527  0.0359     0.0498
#
#sigma^2 estimated as 0.0005889:  log likelihood = 294.39,  aic = -580.77

GNforecast <- predict(fit, n.ahead=10*12)
#We now visualize the prediction along with the training data.
#
ts.plot(gasnatural,2.718^GNforecast$pred, log = "y", lty = c(1,3)) 

forecast = e ^ GNforecast$pred

# pruebo tb para las transformadas y diff
# estos 3 mod los tenemos que probar y estimar
mod1T <- arima(gasnaturalT,c(1,1,0),c(0,1,1))
mod1T #aic = -571.41
mod2T <- arima(gasnaturalT,c(0,1,1),c(0,1,1))
mod2T #aic = -573.03
mod3T <- arima(gasnaturalT,c(1,1,1),c(0,1,1))
mod3T #aic = -574.35
mod4T <- arima(gasnaturalT,c(1,0,0),c(0,1,1))
mod4T #aic = -572.6

mod1T1d <- arima(gasnaturalT1d,c(1,1,0),c(0,1,1))
mod1T1d #aic = -501.85
mod2T1d <- arima(gasnaturalT1d,c(0,1,1),c(0,1,1))
mod2T1d #aic = -552.95
mod3T1d <- arima(gasnaturalT1d,c(1,1,1),c(0,1,1))
mod3T1d #aic = -556.46
mod4T1d <- arima(gasnaturalT1d,c(1,0,0),c(0,1,1))
mod4T1d #aic = -571.42


##########  *************
# elijo este
# mod3T <- arima(gasnaturalT,c(1,1,1),c(0,1,1))
# mod3T #aic = -574.35
## xq es el de menor aic

# opcion de ajuste automatico
mod_1T <- auto.arima (gasnaturalT1d)
mod_1T
#Series: gasnaturalT1d 
#ARIMA(1,0,1)(0,1,1)[12] 
#
#Coefficients:
#  ar1      ma1     sma1
#0.4975  -0.7766  -0.6942
#s.e.  0.1653   0.1168   0.0978
#
#sigma^2 = 0.0004201:  log likelihood = 291.18
#AIC=-574.36   AICc=-574.01   BIC=-563.25

forecastGNT1d <- forecast(mod_1T, level= c(95), h=10)
autoplot(forecastGNT1d)

# una vez q tenemos la prediccion tenemos q ver si hemos conseguido estos modelos c/la serie
# transformada o sin transformar, en este caso lo consegui con la diferencia del log
# x lo tanto para evaluar el error y obtener la prediccion correcta debemos hacer la antitransformada
# al valor predicho
forecastGN <- forecast(e^(mod_1T), level= c(95), h=10)
autoplot(forecastGN)