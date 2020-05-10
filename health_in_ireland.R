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

.covidprofile_file_path <- file.path('data/CovidStatisticsProfileHPSCIrelandOpenData.csv')

covidprofile_df <- read.csv(.test_file_path, header=TRUE, stringsAsFactors=FALSE)

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

grp_month_cases_df <- covidprofile_df %>%
  group_by(month=floor_date(Date, "month")) %>%
  summarise(Confirmed = max(TotalConfirmedCases), Recovered = max(TotalConfirmedRecovered), Deaths = max(TotalConfirmedDeaths)) %>%
  arrange(factor(month, levels = month.name))

grp_subset_cases <- grp_month_cases_df %>%
  gather(key = case, value = Value, Confirmed:Deaths)

# view(grp_subset_cases)

ggplot(grp_subset_cases, aes(month, Value, fill = case)) + geom_col(position = "dodge")


subset_gender_cases_df <- covidprofile_df[, c("Date","Male", "Female", "Unknown")]

subset_gender_cases_df[is.na(subset_gender_cases_df)] <- 0


grp_gender_cases_df <- subset_gender_cases_df %>%
  group_by(month=floor_date(Date, "month")) %>%
  summarise(Male = max(Male), Female = max(Female), Unknown = max(Unknown)) %>%
  arrange(factor(month, levels = month.name))

colnames(grp_gender_cases_df)

grp_gender_subset_cases <- grp_gender_cases_df %>%
  gather(key = case, value = Value, Male:Unknown)

ggplot(grp_gender_subset_cases, aes(month, Value, fill = case)) + geom_col(position = "dodge")


###########################################################################################################

t.test(grp_gender_cases_df$Male,grp_gender_cases_df$Female, paired = TRUE, alt = "greater")

t.test(grp_month_cases_df$Confirmed,grp_month_cases_df$Recovered, paired = TRUE, alt = "greater")

test <- cor.test(grp_gender_cases_df$Male, grp_gender_cases_df$Female,
                 method = 'spearman', exact = FALSE) 

test











