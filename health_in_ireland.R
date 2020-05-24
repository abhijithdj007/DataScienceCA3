library(dplyr)
library(naniar)
library(VIM)
library(FactoMineR)
library(missMDA)
library("Hmisc")
library(corrplot)
library("PerformanceAnalytics")
library("factoextra")
library(pwr)
library(tidyverse)
library(lubridate)
library(DataCombine)
library(forecast)
library("fUnitRoots")
library(tseries)
library(lmtest)
library(FitAR)
library(DAAG)

.covidprofile_file_path <- file.path('data/CovidStatisticsProfileHPSCIrelandOpenData.csv')

covidprofile_df <- read.csv(.covidprofile_file_path, header=TRUE, stringsAsFactors=FALSE)

header_names <- c("Date", "ProfileDate", "ConfirmedCases", "TotalConfirmedCases", "ConfirmedDeaths", 
                  "TotalConfirmedDeaths", "ConfirmedRecovered", "TotalConfirmedRecovered", "CasesConfirmed", 
                  "HospitalisedCases", "RequiringICUCases", "HealthcareWorkersCases", "ClustersNotified", 
                  "HospitalisedAged5", "HospitalisedAged5to14", "HospitalisedAged15to24", "HospitalisedAged25to34", 
                  "HospitalisedAged35to44", "HospitalisedAged45to54", "HospitalisedAged55to64", "HospitalisedAged65up", 
                  "Male", "Female", "Unknown", "Aged1", "Aged1to4", "Aged5to14", "Aged15to24", "Aged25to34", "Aged35to44", 
                  "Aged45to54", "Aged55to64", "Aged65up", "Median_Age", "CommunityTransmission", "CloseContact", 
                  "TravelAbroad", "UnderInvestigation", "FID")

colnames(covidprofile_df) <- header_names 

covidprofile_df$Date <- as.Date(covidprofile_df$Date, format = '%Y/%m/%d')

# view(covidprofile_df)

summary(covidprofile_df)

dim(covidprofile_df)

dim(na.omit(covidprofile_df))

gg_miss_var(covidprofile_df)

subset_covid_prof_df <- covidprofile_df[, c("Date", "ConfirmedCases", "TotalConfirmedCases", "ConfirmedDeaths", 
                                            "TotalConfirmedDeaths", "ConfirmedRecovered", "TotalConfirmedRecovered")]

missing_data_res <- summary(aggr(subset_covid_prof_df, sortVar=TRUE))$combinations

head(missing_data_res[rev(order(missing_data_res[,2])),])

# view(subset_covid_prof_df)

matrixplot(subset_covid_prof_df, sortby = 2)

str(subset_covid_prof_df)

subset_covid_prof_df[is.na(subset_covid_prof_df)] <- 0

subset_covid_prof_df$Date <- as.numeric(subset_covid_prof_df$Date)

str(subset_covid_prof_df)

subset_covid_prof_df.cor = cor(subset_covid_prof_df)

subset_covid_prof_df.cor

subset_covid_prof_df.rcorr = rcorr(as.matrix(subset_covid_prof_df))

subset_covid_prof_df.coeff = subset_covid_prof_df.rcorr$r
subset_covid_prof_df.p = subset_covid_prof_df.rcorr$P

corrplot(subset_covid_prof_df.cor)

palette = colorRampPalette(c("blue", "white", "red")) (20)
heatmap(x = subset_covid_prof_df.cor, col = palette, symm = TRUE)

chart.Correlation(subset_covid_prof_df, histogram=TRUE, pch=19)

##########################################################################################################
pca <- prcomp(subset_covid_prof_df, center = TRUE, scale. = TRUE)
summary(pca)

eig_values <- get_eigenvalue(pca)
eig_values

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))

pca_for_variables <- get_pca_var(pca)
pca_for_variables

fviz_pca_var(pca, col.var = "black")

head(pca_for_variables$cos2, 10)

fviz_cos2(pca, choice = "var", axes = 1:5)

fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("red", "Blue", "Green"), 
             repel = TRUE # Avoid text overlapping
)

head(pca_for_variables$contrib, 10)

fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("red", "Blue", "Green"),
)

###########################################################################################################

# group by month(using floor_date to extract month from data). 
# since counts already combined we are taking MAX count from grouped value.
# Arranging data via month
grp_month_cases_df <- covidprofile_df %>%
  group_by(month=floor_date(Date, "month")) %>%
  summarise(Confirmed = max(TotalConfirmedCases), Recovered = max(TotalConfirmedRecovered), 
            Deaths = max(TotalConfirmedDeaths)) %>%
  arrange(factor(month, levels = month.name))

# gather column name inside table and dividing it by category as column name
grp_subset_cases <- grp_month_cases_df %>%
  gather(key = case, value = Value, Confirmed:Deaths)

# view(grp_subset_cases)

# Ploting the graph to see confirmed, Deaths and recovered cases month wise 
ggplot(grp_subset_cases, aes(month, Value, fill = case)) + geom_col(position = "dodge")


subset_gender_cases_df <- covidprofile_df[, c("Date","Male", "Female", "Unknown")]

subset_gender_cases_df[is.na(subset_gender_cases_df)] <- 0


# group by month(using floor_date to extract month from data). 
# since counts already combined we are taking MAX count from grouped value.
# Arranging data via month
grp_gender_cases_df <- subset_gender_cases_df %>%
  group_by(month=floor_date(Date, "month")) %>%
  summarise(Male = max(Male), Female = max(Female), Unknown = max(Unknown)) %>%
  arrange(factor(month, levels = month.name))

colnames(grp_gender_cases_df)

# gather column name inside table and dividing it by category as column name
grp_gender_subset_cases <- grp_gender_cases_df %>%
  gather(key = case, value = Value, Male:Unknown)

# Ploting the graph to see confirmed, Deaths and recovered cases month wise 
ggplot(grp_gender_subset_cases, aes(month, Value, fill = case)) + geom_col(position = "dodge")


###########################################################################################################

t.test(grp_gender_cases_df$Male,grp_gender_cases_df$Female, paired = TRUE, alt = "greater")

t.test(grp_month_cases_df$Confirmed,grp_month_cases_df$Recovered, paired = TRUE, alt = "greater")

test <- cor.test(grp_gender_cases_df$Male, grp_gender_cases_df$Female,
                 method = 'spearman', exact = FALSE) 


##########################################################################################################

# subsetting columns Date and confirmed cases
mydata <- covidprofile_df[, c("Date","ConfirmedCases")]
# calculating growth of infected cases and adding new column in dataframe
mydata <- change(mydata,"ConfirmedCases",NewVar = "Growth",slideBy = -1, type="percent")

# Plotting graph to visualise growth of confirmed cases
a<-ggplot(data=mydata,aes(x=Date))+ggtitle('Number of Covid Infected Case: Ireland')
b<-a+geom_line(aes(y=ConfirmedCases,color="Total Infected"),size=1)+geom_line(aes(y=ConfirmedCases*100,color="% change"),size=1)
b + labs(y='Total Infected',x='Date',color='Legend')+scale_y_continuous(name="Total_Infected",sec.axis=sec_axis(~./100,name="% change"))


ConfirmedCasesdata <- covidprofile_df[, c("ConfirmedCases")]

# transforming data frame to time series
tsData = ts(ConfirmedCasesdata, start = decimal_date(as.Date("2020-02-29")), frequency = 365.25)
# plot of ts graph
plot(tsData)


boxplot(tsData ~ cycle(tsData),
        xlab="Year",
        ylab = "Confirmed Cases (1000's)" ,
        main ="Dates of Confirmed Cases Boxplot from Feb to May")

#decompose data frame to understand data 
components.ts = decompose(tsData)
plot(components.ts)

#Test stationarity of the time series
adf.test(tsData, alternative = "stationary")

#Kwiatkowski-Phillips-Schmidt-Shin unit test
urkpssTest(tsData, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
# remove non-stationary values
tsstationary = diff(tsData, differences=1)

plot(tsstationary)

#autoregressive, integrated, and moving average parts of the model
acf(tsstationary)
pacf(tsstationary)

#fitting of p, q and d value into ARIMA model
fitARIMA <- arima(tsData, order=c(1,0,1),seasonal = list(order = c(1,0,0), period = 7),method="ML")
fitARIMA

coeftest(fitARIMA)

confint(fitARIMA)

acf(fitARIMA$residuals)



#predicting 7 days infected count
prediction <- predict(fitARIMA, n.ahead = 7)
prediction

#forecast future value of time series
forecast_infected_cases <- forecast(fitARIMA, level = c(95), h = 36)
forecast_infected_cases

autoplot(forecast_infected_cases)

plot(forecast(forecast_infected_cases, 6), xlab = "Days", ylab = "Annual Infected Cases")

#find optimal and better combination of order parameter
auto_arima_model <- auto.arima(tsData)
auto_arima_model

accuracy(auto_arima_model)

accuracy(fitARIMA)

plot(forecast(auto_arima_model, 6), xlab = "Days", ylab = "Annual Infected Cases")

qqnorm(fitARIMA$residuals)
qqline(fitARIMA$residuals)

qqnorm(auto_arima_model$residuals)
qqline(auto_arima_model$residuals)

Box.test(fitARIMA$residuals, type = "Ljung-Box")

Box.test(auto_arima_model$residuals, type = "Ljung-Box")

View(ConfirmedCasesdata)

# For training and test of ARIMA model
confirmed_cases_train <- window(x = tsData, start=decimal_date(as.Date("2020-02-29")), end=decimal_date(as.Date("2020-04-28")))
confirmed_cases_test <- window(x = tsData, start=decimal_date(as.Date("2020-04-28")))

# training the model with 70% data points
fit <- arima(confirmed_cases_train, c(1,1,1), seasonal = list(order = c(1,1,1), period = 7))
fit

#checking auto arima values
auto_arima_model <- auto.arima(confirmed_cases_train)
auto_arima_model

#forecasting from auto arima the data points
predict_auto_ARIMA <- forecast(auto_arima_model, 24)
predict_auto_ARIMA

#forecasting from manual arima the data points
precict_manual_ARIMA <- forecast(fit, 7)
precict_manual_ARIMA

ts.union(confirmed_cases_test, predict_auto_ARIMA)

#combining and test the trained model
actuals_predictions <- data.frame(cbind(actuals = confirmed_cases_test, predicted = predict_auto_ARIMA))
head(actuals_predictions)

#checking the accuracy of the model
correlation_accuracy <- cor(actuals_predictions)
correlation_accuracy


############################################################################################################

mydata <- covidprofile_df[, c("Date","ConfirmedCases")]
mydata <- change(mydata,"ConfirmedCases",NewVar = "Growth",slideBy = -1, type="percent")


myts=ts(log(mydata$Growth)[3:n],start = as.Date(mydata$Date[3]),end = as.Date(mydata$Date[n]), frequency = 1)

fit<-arima(myts,order = c(0,1,1))

predi=predict(fit,n.ahead=14)

ts.plot(diff(myts),ylab="diff(log(Growth))")

Predicted_case<-c(0,rep(NA,14)) 
Predicted_case[1]<-mydata$ConfirmedCases[length(mydata$ConfirmedCases)]
for (i in 2:15) {
  Predicted_case[i] = Predicted_case[i-1]*(1+exp(predi$pred[i-1])/100)
}

Fitted_case<-c(rep(NA,n))
for (i in 7:n) {
  Fitted_case[i] = round(mydata$ConfirmedCases[i-1]*(1+exp(myts[i-7+1]-fit$residuals[i-7+1])/100),digits=0)
}

#plot the results
report_data <- data.frame("Date"=c(time(predi$pred)),"ConfirmedCases"=c(Predicted_case[2:15]),"Growth"= c(exp(predi$pred)))
full<-data.frame("Date"=c(as.Date(mydata$Date),as.Date(as.numeric(time(predi$pred)),origin="1970-01-01")),"Infected_case_fit"=c(Fitted_case,round(report_data$ConfirmedCases,digits=0)),"Growth"=c(mydata$Growth,report_data$Growth),"raw_dt"=c(mydata$ConfirmedCases,rep(NA,14)))

a<-ggplot() + geom_jitter(data=full[1:42,], aes(x = Date, y = raw_dt), color = "blue")
b<-a + ylab('Total Infected case') + xlab('Date') + ggtitle('Ireland Covid-19 Infected Forecast')
c<-b + geom_jitter(data=full[(n+1):(n+7),], aes(x = Date, y = Infected_case_fit), color = "red")
d<-c + geom_text(data=full[(n+1):(n+7),], aes(x = Date, y = Infected_case_fit, label=Infected_case_fit), size=3, vjust = 0.25)
d

##########################################################################################################

# # subset_cases_df <- covidprofile_df[, c("Date","ConfirmedCases","ConfirmedDeaths", "ConfirmedRecovered")]
# # 
# # plot(subset_cases_df$month, subset_cases_df$ConfirmedCases)
# # 
# # plot(subset_cases_df$Date, subset_cases_df$ConfirmedRecovered)
# # 
# # plot(subset_cases_df$Date, subset_cases_df$ConfirmedDeaths)
# 
# subset_Confirmedcases_df <- grp_month_cases_df[, c("Confirmed")]
# 
# ts_subset_cases_df <- ts(subset_Confirmedcases_df)
# 
# cat("Start of air passengers : ", start(ts_subset_cases_df), "\n")
# 
# cat("End of air passengers : ", end(ts_subset_cases_df), "\n")
# 
# cat("Frequency of air passengers : ", frequency(ts_subset_cases_df), "\n")
# 
# print(summary(ts_subset_cases_df))
# 
# frequency(ts_subset_cases_df)
# 
# cycle(ts_subset_cases_df)
# 
# View(subset_Confirmedcases_df)
# 
# air_passengers <- ts(subset_Confirmedcases_df, start=c(2019, 2), end=c(2020, 6), frequency=12)
# 
# cycle(air_passengers)
# 
# # View(air_passengers)
# 
# frequency(air_passengers)
# 
# na_records <- air_passengers[!complete.cases(air_passengers)]
# sum(na_records)
# 
# # View(air_passengers)
# 
# plot(air_passengers)
# 
# abline(reg=lm(air_passengers~time(air_passengers)))
# 
# 
# # y_subset_cases_df <- ts(subset_cases_df, start=c(2020, 1), frequency=12)
# # 
# # cycle(y_subset_cases_df)
# # 
# # frequency(y_subset_cases_df)
# # 
# # na_records <- y_subset_cases_df[!complete.cases(y_subset_cases_df)]
# # 
# # sum(na_records)
# 
# plot(aggregate(ts_subset_cases_df,FUN=mean))
# 
# # ggplot(covidprofile_df, aes(Date, ConfirmedCases)) + geom_bar(stat = "identity")
# # 
# # plot(covidprofile_df$ConfirmedCases)
# # 
# # d.PosCases <- diff(covidprofile_df$ConfirmedCases)
# # plot(d.PosCases)
# # 
# # summary(covidprofile_df$ConfirmedCases)
# # 
# # adf.test(covidprofile_df$ConfirmedCases,"stationary", k=0)
# 
# # tscovid <- ts(covidprofile_df$ConfirmedCases, frequency = 2, start = c(29/02/2020,1))
# # 
# # plot(tscovid)
# # 
# # cycle(tscovid)
# 
# # tscovid <- ts(covidprofile_df$ConfirmedRecovered, frequency = 2, start = c(15/02/2020,1))
# # 
# # plot(tscovid)
#  
# # tscovid <- ts(covidprofile_df$ConfirmedDeaths, frequency = 1, start = c(15/02/2020,1))
# # 
# # plot(tscovid)








