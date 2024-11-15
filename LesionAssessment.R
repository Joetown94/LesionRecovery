#############################################################
#Thesis- Final Cleaned Script                               #
#Joseph Townsend                                            #
#Last Modified 11/28/2018                                   #
#############################################################

setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis/Hobo Data/2018 Data")
library(stringr)
library(plyr)
library(data.table)
library(reshape2)
library(ncdf4)
library(naniar)
library(ggplot2)
library(vegan)
library(grid)
library(ggrepel)
library(psycho)
library(Hmisc)
library(lme4)
library(lmerTest)
library(gridExtra)
############################Manual Temperature Data Import For Sample Period#################################
setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis/Hobo Data/2018 Data")


list.filenames = list.files(pattern=".csv$")
Exp_temp = list()

for(i in 1:length(list.filenames)){
  Exp_temp[[i]] = read.csv(list.filenames[i])
}

Exp_temp[[1]]$Site = "SCPD"
Exp_temp[[2]]$Site = "BUC"
Exp_temp[[3]]$Site = "CSE"
Exp_temp[[4]]$Site = "FLC"
Exp_temp[[5]]$Site = "GKT"
Exp_temp[[6]]$Site = "HND"
Exp_temp[[7]]$Site = "SCPS"
Exp_temp[[8]]$Site = "SHR"

names(Exp_temp[[1]]) = names(Exp_temp[[2]]) = names(Exp_temp[[3]]) = names(Exp_temp[[4]]) = names(Exp_temp[[5]]) = names(Exp_temp[[6]]) = names(Exp_temp[[7]]) = names(Exp_temp[[8]])

Exp_T = rbindlist(Exp_temp)

names(Exp_T) = c("Num","Time","Temp","Site")

Exp_split = as.data.frame(str_split_fixed(Exp_T$Time, " ", 2))

Exp_sp1 = cbind(Exp_T, Exp_split)

Exp_split2 = as.data.frame(str_split_fixed(Exp_sp1$V1, "/", 3))

Exp_sep = cbind(Exp_sp1, Exp_split2)

names(Exp_sep) = c("number", "date_combined", "temperature", "site", "date", "time", "month", "day", "year")
Exp_sep$year = revalue(Exp_sep$year, c("18" = "2018", "17" = "2017", "16" = "2016", "15" = "2015" , "14" = "2014", "13" = "2013", "12" = "2012", "11" = "2011", "10" = "2010", "09" = "2009", "08" = "2008", "07" = "2007", "06" = "2006", "05" = "2005", "04" = "2004", "03" = "2003"))

Temp_samp2 = data.frame(Exp_sep)
Temp_samp2$month = as.numeric(levels(Temp_samp2$month))[Temp_samp2$month]
Temp_samp2$day = as.numeric(levels(Temp_samp2$day))[Temp_samp2$day]
Temp_samp2$year = as.numeric(levels(Temp_samp2$year))[Temp_samp2$year]
Temp_samp2$date_f = as.Date(with(Temp_samp2, paste(year, month, day, sep = '-')), "%Y-%m-%d")
Temp_samp2$date_f = as.POSIXlt(Temp_samp2$date_f, format = "%y-%m-%d")
Temp_samp2$julian = Temp_samp2$date_f$yday + 18000
Temp_samp2$date_f = as.Date(with(Temp_samp2, paste(year, month, day, sep = '-')), "%Y-%m-%d")
Temp_samp2$site = as.factor(Temp_samp2$site)
#Average it together by daily
#Daily Averages
Temp_samp2_daily = ddply(Temp_samp2, .(site,julian), summarise,
                         mean(temperature),
                         max(temperature),
                         min(temperature),
                         sd(temperature))

names(Temp_samp2_daily) = c("site","julian","avgtemp","maxtemp","mintemp","sdtemp")
Temp_samp3_daily = data.frame(Temp_samp2_daily[Temp_samp2_daily$julian > 18131 & Temp_samp2_daily$julian < 18252,])


#Now to get analytics on the actual sample windows#

#Buck Island STT#
Temp_samp3_daily.buc = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "BUC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "BUC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "BUC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "BUC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "BUC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "BUC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "BUC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "BUC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "BUC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "BUC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "BUC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "BUC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "BUC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "BUC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "BUC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "BUC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "BUC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "BUC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "BUC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "BUC" & julian > 18134 & julian < 18252]))))

Temp_samp3_daily.buc$site = "BUC"

names(Temp_samp3_daily.buc) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

#College Shoal East#
Temp_samp3_daily.cse = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "CSE" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "CSE" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "CSE" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "CSE" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "CSE" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "CSE" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "CSE" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "CSE" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "CSE" & julian > 18202 & julian < 18251]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "CSE" & julian > 18202 & julian < 18251]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "CSE" & julian > 18202 & julian < 18251]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "CSE" & julian > 18202 & julian < 18251]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "CSE" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "CSE" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "CSE" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "CSE" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "CSE" & julian > 18142 & julian < 18251]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "CSE" & julian > 18142 & julian < 18251]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "CSE" & julian > 18142 & julian < 18251]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "CSE" & julian > 18142 & julian < 18251]))))

Temp_samp3_daily.cse$site = "CSE"

names(Temp_samp3_daily.cse) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")


#Flat Cay#
Temp_samp3_daily.flc = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "FLC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "FLC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "FLC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "FLC" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "FLC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "FLC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "FLC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "FLC" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "FLC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "FLC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "FLC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "FLC" & julian > 18194 & julian < 18252]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "FLC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "FLC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "FLC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "FLC" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "FLC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "FLC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "FLC" & julian > 18134 & julian < 18252]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "FLC" & julian > 18134 & julian < 18252]))))

Temp_samp3_daily.flc$site = "FLC"

names(Temp_samp3_daily.flc) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

#Grammanik Bank#
Temp_samp3_daily.gkt = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "GKT" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "GKT" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "GKT" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "GKT" & julian > 18142 & julian < 18172]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "GKT" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "GKT" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "GKT" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "GKT" & julian > 18172 & julian < 18202]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "GKT" & julian > 18202 & julian < 18248]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "GKT" & julian > 18202 & julian < 18248]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "GKT" & julian > 18202 & julian < 18248]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "GKT" & julian > 18202 & julian < 18248]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "GKT" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "GKT" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "GKT" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "GKT" & julian > 18142 & julian < 18202]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "GKT" & julian > 18142 & julian < 18248]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "GKT" & julian > 18142 & julian < 18248]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "GKT" & julian > 18142 & julian < 18248]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "GKT" & julian > 18142 & julian < 18248]))))

Temp_samp3_daily.gkt$site = "GKT"

names(Temp_samp3_daily.gkt) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

#Hind Bank#
Temp_samp3_daily.hnd = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "HND" & julian > 18131 & julian < 18161]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "HND" & julian > 18131 & julian < 18161]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "HND" & julian > 18131 & julian < 18161]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "HND" & julian > 18131 & julian < 18161]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "HND" & julian > 18161 & julian < 18191]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "HND" & julian > 18161 & julian < 18191]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "HND" & julian > 18161 & julian < 18191]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "HND" & julian > 18161 & julian < 18191]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "HND" & julian > 18191 & julian < 18248]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "HND" & julian > 18191 & julian < 18248]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "HND" & julian > 18191 & julian < 18248]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "HND" & julian > 18191 & julian < 18248]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "HND" & julian > 18131 & julian < 18191]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "HND" & julian > 18131 & julian < 18191]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "HND" & julian > 18131 & julian < 18191]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "HND" & julian > 18131 & julian < 18191]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "HND" & julian > 18131 & julian < 18248]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "HND" & julian > 18131 & julian < 18248]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "HND" & julian > 18131 & julian < 18248]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "HND" & julian > 18131 & julian < 18248]))))

Temp_samp3_daily.hnd$site = "HND"

names(Temp_samp3_daily.hnd) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")


#South Capella Deep#
Temp_samp3_daily.scpd = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "SCPD" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPD" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPD" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPD" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPD" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPD" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPD" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPD" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPD" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPD" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPD" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPD" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPD" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPD" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPD" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPD" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPD" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPD" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPD" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPD" & julian > 18134 & julian < 18247]))))

Temp_samp3_daily.scpd$site = "SCPD"

names(Temp_samp3_daily.scpd) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

#South Capella Shallow#
Temp_samp3_daily.scps = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "SCPS" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPS" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPS" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPS" & julian > 18134 & julian < 18164]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPS" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPS" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPS" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPS" & julian > 18164 & julian < 18194]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPS" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPS" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPS" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPS" & julian > 18194 & julian < 18247]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPS" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPS" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPS" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPS" & julian > 18134 & julian < 18194]))),
                             (with(Temp_samp3_daily, mean(avgtemp[site == "SCPS" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, max(avgtemp[site == "SCPS" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, min(avgtemp[site == "SCPS" & julian > 18134 & julian < 18247]))),
                             (with(Temp_samp3_daily, sd(avgtemp[site == "SCPS" & julian > 18134 & julian < 18247]))))

Temp_samp3_daily.scps$site = "SCPS"

names(Temp_samp3_daily.scps) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

#Seahorse Cottage Shoal#
Temp_samp3_daily.shr = data.frame((with(Temp_samp3_daily, mean(avgtemp[site == "SHR" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "SHR" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "SHR" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "SHR" & julian > 18134 & julian < 18164]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "SHR" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "SHR" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "SHR" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "SHR" & julian > 18164 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "SHR" & julian > 18194 & julian < 18247]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "SHR" & julian > 18194 & julian < 18247]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "SHR" & julian > 18194 & julian < 18247]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "SHR" & julian > 18194 & julian < 18247]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "SHR" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "SHR" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "SHR" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "SHR" & julian > 18134 & julian < 18194]))),
                            (with(Temp_samp3_daily, mean(avgtemp[site == "SHR" & julian > 18134 & julian < 18247]))),
                            (with(Temp_samp3_daily, max(avgtemp[site == "SHR" & julian > 18134 & julian < 18247]))),
                            (with(Temp_samp3_daily, min(avgtemp[site == "SHR" & julian > 18134 & julian < 18247]))),
                            (with(Temp_samp3_daily, sd(avgtemp[site == "SHR" & julian > 18134 & julian < 18247]))))

Temp_samp3_daily.shr$site = "SHR"

names(Temp_samp3_daily.shr) = c("T1temp", "T1max", "T1min","T1std","T2temp","T2max","T2min","T2std","T3temp","T3max","T3min","T3std", "T12temp","T12max","T12min","T12std", "Totaltemp","Totalmax","Totalmin","Totalstd","Site")

samp_temp = rbind(Temp_samp3_daily.buc, Temp_samp3_daily.cse, Temp_samp3_daily.flc, Temp_samp3_daily.gkt, Temp_samp3_daily.hnd, Temp_samp3_daily.scpd, Temp_samp3_daily.scps, Temp_samp3_daily.shr)

samp_temp$Site = as.factor(samp_temp$Site)

samp_temp.T1 = data.frame(samp_temp$T1temp, samp_temp$T1max, samp_temp$T1min, samp_temp$T1std, samp_temp$Site)
samp_temp.T1$timepoint = "T1"
names(samp_temp.T1) = c("avgtemp","maxtemp","mintemp","stdtemp","site","timepoint")

samp_temp.T2 = data.frame(samp_temp$T2temp, samp_temp$T2max, samp_temp$T2min, samp_temp$T2std, samp_temp$Site)
samp_temp.T2$timepoint = "T2"
names(samp_temp.T2) = c("avgtemp","maxtemp","mintemp","stdtemp","site","timepoint")

samp_temp.T3 = data.frame(samp_temp$T3temp, samp_temp$T3max, samp_temp$T3min, samp_temp$T3std, samp_temp$Site)
samp_temp.T3$timepoint = "T3"
names(samp_temp.T3) = c("avgtemp","maxtemp","mintemp","stdtemp","site","timepoint")

samp_temp.T12 = data.frame(samp_temp$T12temp, samp_temp$T12max, samp_temp$T12min, samp_temp$T12std, samp_temp$Site)
samp_temp.T12$timepoint = "T12"
names(samp_temp.T12) = c("avgtemp","maxtemp","mintemp","stdtemp","site","timepoint")

samp_temp.Tot = data.frame(samp_temp$Totaltemp, samp_temp$Totalmax, samp_temp$Totalmin, samp_temp$Totalstd, samp_temp$Site)
samp_temp.Tot$timepoint = "Total"
names(samp_temp.Tot) = c("avgtemp","maxtemp","mintemp","stdtemp","site","timepoint")

samp_temp.tall = rbind(samp_temp.T1, samp_temp.T2, samp_temp.T3, samp_temp.T12, samp_temp.Tot)
samp_temp.tall$timepoint = as.factor(samp_temp.tall$timepoint)

#######################Benthic Orbital Velocity Import#############################
setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis/BOV data")
library(dplyr)

flist = list.files(pattern = "^.*\\.(nc|NC|Nc|Nc)$")
nc = nc_open(flist[[1]])

#fun with .netcdf files!
attributes(nc)
print(paste("The file has",nc$nvars, "variables", nc$ndims, "dimensions, and", nc$natts, "NetCDF attributes"))

attributes(nc$var)$names
ncatt_get(nc, attributes(nc$var)$names[1])

bov_may = nc_open("ubenth_new_may.nc")
bov_may.data = data.frame(ncvar_get(bov_may))
bov.may = data.frame(t(bov_may.data))

bov_june = nc_open("ubenth_new_june.nc")
bov_june.data = data.frame(ncvar_get(bov_june))
bov.june = data.frame(t(bov_june.data))

bov_july = nc_open("ubenth_new_july.nc")
bov_july.data = data.frame(ncvar_get(bov_july))
bov.july = data.frame(t(bov_july.data))

bov_aug = nc_open("ubenth_new_aug.nc")
bov_aug.data = data.frame(ncvar_get(bov_aug))
bov.aug = data.frame(t(bov_aug.data))

bov_sept = nc_open("ubenth_new_sept.nc")
bov_sept.data = data.frame(ncvar_get(bov_sept))
bov.sept = data.frame(t(bov_sept.data))

bov_all = rbind(bov.may, bov.june, bov.july, bov.aug, bov.sept)

names(bov_all) = c("BUC","CSE","FLC","GKT","HND","SHR","SCPS","SCPD")

bov_all = bov_all %>% replace_with_na_all(condition = ~.x == 0)
bov_all$julian = seq(18120, 18274, 1)

all_bov.tall = data.frame(as.vector(rbind(bov_all$BUC,bov_all$CSE,bov_all$FLC,bov_all$GKT,bov_all$HND,bov_all$SHR,bov_all$SCPS,bov_all$SCPD)))
all_bov.tall$julian = rep(seq(18120,18274,1),8)
all_bov.tall$site = "NA"
names(all_bov.tall) = c("bov","julian","site")
all_bov.tall$site[1:153] = "BUC"
all_bov.tall$site[154:306] = "CSE"
all_bov.tall$site[307:459] = "FLC"
all_bov.tall$site[460:612] = "GKT"
all_bov.tall$site[613:765] = "HND"
all_bov.tall$site[766:918] = "SHR"
all_bov.tall$site[919:1071] = "SCPS"
all_bov.tall$site[1072:1224] = "SCPD"

all_bov.tall = data.frame(all_bov.tall[complete.cases(all_bov.tall),])
bov_all = data.frame(bov_all[complete.cases(bov_all),]) #Removes all rows with NA values

###################Process BOV data like temperature data############################
#Buck Island STT#
bov_all.buc = data.frame(with(bov_all, mean(BUC[julian > 18134 & julian < 18164])),
                         (with(bov_all, max(BUC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, min(BUC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, sd(BUC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, mean(BUC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, max(BUC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, min(BUC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, sd(BUC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, mean(BUC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, max(BUC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, min(BUC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, sd(BUC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, mean(BUC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, max(BUC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, min(BUC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, sd(BUC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, mean(BUC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, max(BUC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, min(BUC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, sd(BUC[julian > 18134 & julian < 18252]))))

bov_all.buc$site = "BUC"

names(bov_all.buc) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

#College Shoal East#
bov_all.cse = data.frame((with(bov_all, mean(CSE[julian > 18142 & julian < 18172]))),
                         (with(bov_all, max(CSE[julian > 18142 & julian < 18172]))),
                         (with(bov_all, min(CSE[julian > 18142 & julian < 18172]))),
                         (with(bov_all, sd(CSE[julian > 18142 & julian < 18172]))),
                         (with(bov_all, mean(CSE[julian > 18172 & julian < 18202]))),
                         (with(bov_all, max(CSE[julian > 18172 & julian < 18202]))),
                         (with(bov_all, min(CSE[julian > 18172 & julian < 18202]))),
                         (with(bov_all, sd(CSE[julian > 18172 & julian < 18202]))),
                         (with(bov_all, mean(CSE[julian > 18202 & julian < 18251]))),
                         (with(bov_all, max(CSE[julian > 18202 & julian < 18251]))),
                         (with(bov_all, min(CSE[julian > 18202 & julian < 18251]))),
                         (with(bov_all, sd(CSE[julian > 18202 & julian < 18251]))),
                         (with(bov_all, mean(CSE[julian > 18142 & julian < 18202]))),
                         (with(bov_all, max(CSE[julian > 18142 & julian < 18202]))),
                         (with(bov_all, min(CSE[julian > 18142 & julian < 18202]))),
                         (with(bov_all, sd(CSE[julian > 18142 & julian < 18202]))),
                         (with(bov_all, mean(CSE[julian > 18142 & julian < 18251]))),
                         (with(bov_all, max(CSE[julian > 18142 & julian < 18251]))),
                         (with(bov_all, min(CSE[julian > 18142 & julian < 18251]))),
                         (with(bov_all, sd(CSE[julian > 18142 & julian < 18251]))))

bov_all.cse$site = "CSE"

names(bov_all.cse) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")


#Flat Cay#
bov_all.flc = data.frame((with(bov_all, mean(FLC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, max(FLC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, min(FLC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, sd(FLC[julian > 18134 & julian < 18164]))),
                         (with(bov_all, mean(FLC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, max(FLC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, min(FLC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, sd(FLC[julian > 18164 & julian < 18194]))),
                         (with(bov_all, mean(FLC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, max(FLC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, min(FLC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, sd(FLC[julian > 18194 & julian < 18252]))),
                         (with(bov_all, mean(FLC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, max(FLC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, min(FLC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, sd(FLC[julian > 18134 & julian < 18194]))),
                         (with(bov_all, mean(FLC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, max(FLC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, min(FLC[julian > 18134 & julian < 18252]))),
                         (with(bov_all, sd(FLC[julian > 18134 & julian < 18252]))))

bov_all.flc$site = "FLC"

names(bov_all.flc) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

#Grammanik Bank#
bov_all.gkt = data.frame((with(bov_all, mean(GKT[julian > 18142 & julian < 18172]))),
                         (with(bov_all, max(GKT[julian > 18142 & julian < 18172]))),
                         (with(bov_all, min(GKT[julian > 18142 & julian < 18172]))),
                         (with(bov_all, sd(GKT[julian > 18142 & julian < 18172]))),
                         (with(bov_all, mean(GKT[julian > 18172 & julian < 18202]))),
                         (with(bov_all, max(GKT[julian > 18172 & julian < 18202]))),
                         (with(bov_all, min(GKT[julian > 18172 & julian < 18202]))),
                         (with(bov_all, sd(GKT[julian > 18172 & julian < 18202]))),
                         (with(bov_all, mean(GKT[julian > 18202 & julian < 18248]))),
                         (with(bov_all, max(GKT[julian > 18202 & julian < 18248]))),
                         (with(bov_all, min(GKT[julian > 18202 & julian < 18248]))),
                         (with(bov_all, sd(GKT[julian > 18202 & julian < 18248]))),
                         (with(bov_all, mean(GKT[julian > 18142 & julian < 18202]))),
                         (with(bov_all, max(GKT[julian > 18142 & julian < 18202]))),
                         (with(bov_all, min(GKT[julian > 18142 & julian < 18202]))),
                         (with(bov_all, sd(GKT[julian > 18142 & julian < 18202]))),
                         (with(bov_all, mean(GKT[julian > 18142 & julian < 18248]))),
                         (with(bov_all, max(GKT[julian > 18142 & julian < 18248]))),
                         (with(bov_all, min(GKT[julian > 18142 & julian < 18248]))),
                         (with(bov_all, sd(GKT[julian > 18142 & julian < 18248]))))

bov_all.gkt$site = "GKT"

names(bov_all.gkt) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

#Hind Bank#
bov_all.hnd = data.frame((with(bov_all, mean(HND[julian > 18131 & julian < 18161]))),
                         (with(bov_all, max(HND[julian > 18131 & julian < 18161]))),
                         (with(bov_all, min(HND[julian > 18131 & julian < 18161]))),
                         (with(bov_all, sd(HND[julian > 18131 & julian < 18161]))),
                         (with(bov_all, mean(HND[julian > 18161 & julian < 18191]))),
                         (with(bov_all, max(HND[julian > 18161 & julian < 18191]))),
                         (with(bov_all, min(HND[julian > 18161 & julian < 18191]))),
                         (with(bov_all, sd(HND[julian > 18161 & julian < 18191]))),
                         (with(bov_all, mean(HND[julian > 18191 & julian < 18248]))),
                         (with(bov_all, max(HND[julian > 18191 & julian < 18248]))),
                         (with(bov_all, min(HND[julian > 18191 & julian < 18248]))),
                         (with(bov_all, sd(HND[julian > 18191 & julian < 18248]))),
                         (with(bov_all, mean(HND[julian > 18131 & julian < 18191]))),
                         (with(bov_all, max(HND[julian > 18131 & julian < 18191]))),
                         (with(bov_all, min(HND[julian > 18131 & julian < 18191]))),
                         (with(bov_all, sd(HND[julian > 18131 & julian < 18191]))),
                         (with(bov_all, mean(HND[julian > 18131 & julian < 18248]))),
                         (with(bov_all, max(HND[julian > 18131 & julian < 18248]))),
                         (with(bov_all, min(HND[julian > 18131 & julian < 18248]))),
                         (with(bov_all, sd(HND[julian > 18131 & julian < 18248]))))

bov_all.hnd$site = "HND"

names(bov_all.hnd) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")


#South Capella Deep#
bov_all.scpd = data.frame((with(bov_all, mean(SCPD[julian > 18134 & julian < 18164]))),
                          (with(bov_all, max(SCPD[julian > 18134 & julian < 18164]))),
                          (with(bov_all, min(SCPD[julian > 18134 & julian < 18164]))),
                          (with(bov_all, sd(SCPD[julian > 18134 & julian < 18164]))),
                          (with(bov_all, mean(SCPD[julian > 18164 & julian < 18194]))),
                          (with(bov_all, max(SCPD[julian > 18164 & julian < 18194]))),
                          (with(bov_all, min(SCPD[julian > 18164 & julian < 18194]))),
                          (with(bov_all, sd(SCPD[julian > 18164 & julian < 18194]))),
                          (with(bov_all, mean(SCPD[julian > 18194 & julian < 18247]))),
                          (with(bov_all, max(SCPD[julian > 18194 & julian < 18247]))),
                          (with(bov_all, min(SCPD[julian > 18194 & julian < 18247]))),
                          (with(bov_all, sd(SCPD[julian > 18194 & julian < 18247]))),
                          (with(bov_all, mean(SCPD[julian > 18134 & julian < 18194]))),
                          (with(bov_all, max(SCPD[julian > 18134 & julian < 18194]))),
                          (with(bov_all, min(SCPD[julian > 18134 & julian < 18194]))),
                          (with(bov_all, sd(SCPD[julian > 18134 & julian < 18194]))),
                          (with(bov_all, mean(SCPD[julian > 18134 & julian < 18247]))),
                          (with(bov_all, max(SCPD[julian > 18134 & julian < 18247]))),
                          (with(bov_all, min(SCPD[julian > 18134 & julian < 18247]))),
                          (with(bov_all, sd(SCPD[julian > 18134 & julian < 18247]))))

bov_all.scpd$site = "SCPD"

names(bov_all.scpd) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

#South Capella Shallow#
bov_all.scps = data.frame((with(bov_all, mean(SCPS[julian > 18134 & julian < 18164]))),
                          (with(bov_all, max(SCPS[julian > 18134 & julian < 18164]))),
                          (with(bov_all, min(SCPS[julian > 18134 & julian < 18164]))),
                          (with(bov_all, sd(SCPS[julian > 18134 & julian < 18164]))),
                          (with(bov_all, mean(SCPS[julian > 18164 & julian < 18194]))),
                          (with(bov_all, max(SCPS[julian > 18164 & julian < 18194]))),
                          (with(bov_all, min(SCPS[julian > 18164 & julian < 18194]))),
                          (with(bov_all, sd(SCPS[julian > 18164 & julian < 18194]))),
                          (with(bov_all, mean(SCPS[julian > 18194 & julian < 18247]))),
                          (with(bov_all, max(SCPS[julian > 18194 & julian < 18247]))),
                          (with(bov_all, min(SCPS[julian > 18194 & julian < 18247]))),
                          (with(bov_all, sd(SCPS[julian > 18194 & julian < 18247]))),
                          (with(bov_all, mean(SCPS[julian > 18134 & julian < 18194]))),
                          (with(bov_all, max(SCPS[julian > 18134 & julian < 18194]))),
                          (with(bov_all, min(SCPS[julian > 18134 & julian < 18194]))),
                          (with(bov_all, sd(SCPS[julian > 18134 & julian < 18194]))),
                          (with(bov_all, mean(SCPS[julian > 18134 & julian < 18247]))),
                          (with(bov_all, max(SCPS[julian > 18134 & julian < 18247]))),
                          (with(bov_all, min(SCPS[julian > 18134 & julian < 18247]))),
                          (with(bov_all, sd(SCPS[julian > 18134 & julian < 18247]))))

bov_all.scps$site = "SCPS"

names(bov_all.scps) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

#Seahorse Cottage Shoal#
bov_all.shr = data.frame((with(bov_all, mean(SHR[julian > 18134 & julian < 18164]))),
                         (with(bov_all, max(SHR[julian > 18134 & julian < 18164]))),
                         (with(bov_all, min(SHR[julian > 18134 & julian < 18164]))),
                         (with(bov_all, sd(SHR[julian > 18134 & julian < 18164]))),
                         (with(bov_all, mean(SHR[julian > 18164 & julian < 18194]))),
                         (with(bov_all, max(SHR[julian > 18164 & julian < 18194]))),
                         (with(bov_all, min(SHR[julian > 18164 & julian < 18194]))),
                         (with(bov_all, sd(SHR[julian > 18164 & julian < 18194]))),
                         (with(bov_all, mean(SHR[julian > 18194 & julian < 18247]))),
                         (with(bov_all, max(SHR[julian > 18194 & julian < 18247]))),
                         (with(bov_all, min(SHR[julian > 18194 & julian < 18247]))),
                         (with(bov_all, sd(SHR[julian > 18194 & julian < 18247]))),
                         (with(bov_all, mean(SHR[julian > 18134 & julian < 18194]))),
                         (with(bov_all, max(SHR[julian > 18134 & julian < 18194]))),
                         (with(bov_all, min(SHR[julian > 18134 & julian < 18194]))),
                         (with(bov_all, sd(SHR[julian > 18134 & julian < 18194]))),
                         (with(bov_all, mean(SHR[julian > 18134 & julian < 18247]))),
                         (with(bov_all, max(SHR[julian > 18134 & julian < 18247]))),
                         (with(bov_all, min(SHR[julian > 18134 & julian < 18247]))),
                         (with(bov_all, sd(SHR[julian > 18134 & julian < 18247]))))

bov_all.shr$site = "SHR"

names(bov_all.shr) = c("T1bov", "T1bovmax", "T1min","T1std","T2bov","T2max","T2min","T2std","T3bov","T3max","T3min","T3std", "T12bov","T12max","T12min","T12std", "Totalbov","Totalmax","Totalmin","Totalstd","Site")

samp_bov = rbind(bov_all.buc, bov_all.cse, bov_all.flc, bov_all.gkt, bov_all.hnd, bov_all.scpd, bov_all.scps, bov_all.shr)

samp_bov$Site = as.factor(samp_bov$Site)

samp_bov.T1 = data.frame(samp_bov$T1bov, samp_bov$T1bovmax, samp_bov$T1min, samp_bov$T1std, samp_bov$Site)
samp_bov.T1$timepoint = "T1"
names(samp_bov.T1) = c("avgbov","maxbov","minbov","stdbov","site","timepoint")

samp_bov.T2 = data.frame(samp_bov$T2bov, samp_bov$T2max, samp_bov$T2min, samp_bov$T2std, samp_bov$Site)
samp_bov.T2$timepoint = "T2"
names(samp_bov.T2) = c("avgbov","maxbov","minbov","stdbov","site","timepoint")

samp_bov.T3 = data.frame(samp_bov$T3bov, samp_bov$T3max, samp_bov$T3min, samp_bov$T3std, samp_bov$Site)
samp_bov.T3$timepoint = "T3"
names(samp_bov.T3) = c("avgbov","maxbov","minbov","stdbov","site","timepoint")

samp_bov.T12 = data.frame(samp_bov$T12bov, samp_bov$T12max, samp_bov$T12min, samp_bov$T12std, samp_bov$Site)
samp_bov.T12$timepoint = "T12"
names(samp_bov.T12) = c("avgbov","maxbov","minbov","stdbov","site","timepoint")

samp_bov.Tot = data.frame(samp_bov$Totalbov, samp_bov$Totalmax, samp_bov$Totalmin, samp_bov$Totalstd, samp_bov$Site)
samp_bov.Tot$timepoint = "Total"
names(samp_bov.Tot) = c("avgbov","maxbov","minbov","stdbov","site","timepoint")

samp_bov.tall = rbind(samp_bov.T1, samp_bov.T2, samp_bov.T3, samp_bov.T12, samp_bov.Tot)
samp_bov.tall$timepoint = as.factor(samp_bov.tall$timepoint)


################################PART TWO: PHYSICAL MODELING AND ANALYSIS###############################
setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis")
library(plyr)

#####DATA IMPORT AND RECOVERY METRIC CALCULATION
rel_LLR = read.csv("Lesion_Area_full.csv", header = T, na.strings = c('NA'," NA", "NA ", "#VALUE!","#DIV/0!", "#REF!"))
rel_LLR2 = data.frame(rel_LLR[1:48,])
rel_LLR2 = data.frame(rel_LLR[-28,])

#removing unnecessary columns
rel_LLR3 = rel_LLR2[-c(1,6,8,10,12)]

names(rel_LLR3)[1] = "site"
names(rel_LLR3)[2] = "type"
names(rel_LLR3)[3] = "sitedepth"
names(rel_LLR3)[4] = "depth"
names(rel_LLR3)[9] = "species"

rel_LLR3$colony_size = rel_LLR3$Length*rel_LLR3$Width*rel_LLR3$Height

Rel_rec_t1 = data.frame(
  rel_LLR3$site, rel_LLR3$type, rel_LLR3$sitedepth, rel_LLR3$depth, 
  rel_LLR3$T0.Julian, rel_LLR3$T1.Julian, 
  rel_LLR3$species, rel_LLR3$Colony, rel_LLR3$colony_size, rel_LLR3$Longest.RLRR, 
  rel_LLR3$Initial.lesion.size, rel_LLR3$T0.Perimeter, rel_LLR3$SA.P.0, rel_LLR3$T1.per.recov, 
  rel_LLR3$T1.RLRR, rel_LLR3$T1.PRR, rel_LLR3$SA.P.R1)

Rel_rec_t1$timepoint = as.factor("T1")
names(Rel_rec_t1) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony", "colony_size", "longest LRR", "initial_size","initial_perimeter","SA_ratio_initial","percent_recov","RLRR","PRR","SAP_recov","timepoint")

Rel_rec_t2 = data.frame(
  rel_LLR3$site, rel_LLR3$type, rel_LLR3$sitedepth, rel_LLR3$depth, 
  rel_LLR3$T1.Julian, rel_LLR3$T2.Julian, 
  rel_LLR3$species, rel_LLR3$Colony, rel_LLR3$colony_size, rel_LLR3$Longest.RLRR,
  rel_LLR3$Initial.lesion.size, rel_LLR3$T0.Perimeter, rel_LLR3$SA.P.0, rel_LLR3$T2.per.recov, 
  rel_LLR3$T2.RLRR, rel_LLR3$T2.PRR, rel_LLR3$SA.P.R2)

Rel_rec_t2$timepoint = as.factor("T2")
names(Rel_rec_t2) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony", "colony_size", "longest LRR", "initial_size","initial_perimeter","SA_ratio_initial","percent_recov","RLRR","PRR","SAP_recov", "timepoint")

Rel_rec_t3 = data.frame(
  rel_LLR3$site, rel_LLR3$type, rel_LLR3$sitedepth, rel_LLR3$depth, 
  rel_LLR3$T2.Julian, rel_LLR3$T3.Julian, 
  rel_LLR3$species, rel_LLR3$Colony, rel_LLR3$colony_size, rel_LLR3$Longest.RLRR,
  rel_LLR3$Initial.lesion.size, rel_LLR3$T0.Perimeter, rel_LLR3$SA.P.0, rel_LLR3$T3.per.recov, 
  rel_LLR3$T3.RLRR, rel_LLR3$T3.PRR, rel_LLR3$SA.P.R3)

Rel_rec_t3$timepoint = as.factor("T3")
names(Rel_rec_t3) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony", "colony_size", "longest LRR", "initial_size","initial_perimeter","SA_ratio_initial","percent_recov","RLRR","PRR","SAP_recov","timepoint")

Rel_rec_t12 = data.frame(
  rel_LLR3$site, rel_LLR3$type, rel_LLR3$sitedepth, rel_LLR3$depth,
  rel_LLR3$T0.Julian, rel_LLR3$T1.Julian,
  rel_LLR3$species, rel_LLR3$Colony, rel_LLR3$colony_size, rel_LLR3$Longest.RLRR,
  rel_LLR3$Initial.lesion.size, rel_LLR3$T0.Perimeter, rel_LLR3$SA.P.0, rel_LLR3$T12.per.recov,
  rel_LLR3$T12.RLRR, rel_LLR3$T12.PRR, rel_LLR3$SA.P.R12)

Rel_rec_t12$timepoint = as.factor("T12")
names(Rel_rec_t12) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony", "colony_size", "longest LRR", "initial_size","initial_perimeter","SA_ratio_initial","percent_recov","RLRR","PRR","SAP_recov","timepoint")

Rel_rec_tot = data.frame(
  rel_LLR3$site, rel_LLR3$type, rel_LLR3$sitedepth, rel_LLR3$depth, 
  rel_LLR3$T0.Julian, rel_LLR3$T3.Julian, 
  rel_LLR3$species, rel_LLR3$Colony,rel_LLR3$colony_size, rel_LLR3$Longest.RLRR,
  rel_LLR3$Initial.lesion.size, rel_LLR3$T0.Perimeter, rel_LLR3$SA.P.0, rel_LLR3$Total.per.recov, 
  rel_LLR3$Total.RLRR, rel_LLR3$Total.PRR, rel_LLR3$SA.P.Rtotal)

Rel_rec_tot$timepoint = as.factor("Total")
names(Rel_rec_tot) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony", "colony_size", "longest LRR", "initial_size","initial_perimeter","SA_ratio_initial","percent_recov","RLRR","PRR","SAP_recov","timepoint")

Rel_rec_tall = rbind(Rel_rec_t1,Rel_rec_t2,Rel_rec_t3,Rel_rec_t12,Rel_rec_tot)

###########ABSOLUTE RECOVERY (ALRR)
setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis")
les.ab = read.csv("Lesion_Area_ab.csv", header = T, na.strings = c('NA'," NA", "NA ", "#VALUE!","#DIV/0!",""))
les.ab2 = data.frame(les.ab[1:48,])
les.ab2 = data.frame(les.ab2[-28,])


#removing unnecessary columns
les.ab3 = les.ab2[,-c(1,5,7,9,11,36:42)]

names(les.ab3)[1] = "site"
names(les.ab3)[2] = "type"
names(les.ab3)[3] = "sitedepth"
# names(rel_LLR3)[4] = "depth"
names(les.ab3)[8] = "species"
les.ab3$depth = rel_LLR3$depth

#adding in columns not added in excel
les.ab3$T1.ALRR = ((les.ab3$Initial.lesion.size-les.ab3$T1.Lesion.Size)/(les.ab3$T1.Julian-les.ab3$T0.Julian))
les.ab3$T2.ALRR = ((les.ab3$T1.Lesion.Size-les.ab3$T2.Lesion.Size)/(les.ab3$T2.Julian-les.ab3$T1.Julian))
les.ab3$T3.ALRR = ((les.ab3$T2.Lesion.Size-les.ab3$T3.Lesion.Size)/(les.ab3$T3.Julian-les.ab3$T2.Julian))
les.ab3$T12.ALRR = ((les.ab3$Initial.lesion.size-les.ab3$T2.Lesion.Size)/(les.ab3$T2.Julian-les.ab3$T0.Julian))
les.ab3$Total.ALRR = ((les.ab3$Initial.lesion.size-les.ab3$T3.Lesion.Size)/(les.ab3$T3.Julian-les.ab3$T0.Julian))

les.ab3$T1.APRR = ((les.ab3$T0.Perimeter-les.ab3$T1.Perimeter)/(les.ab3$T1.Julian-les.ab3$T0.Julian))
les.ab3$T2.APRR = ((les.ab3$T1.Perimeter-les.ab3$T2.Perimeter)/(les.ab3$T2.Julian-les.ab3$T1.Julian))
les.ab3$T3.APRR = ((les.ab3$T2.Perimeter-les.ab3$T3.Perimeter)/(les.ab3$T3.Julian-les.ab3$T2.Julian))
les.ab3$T12.APRR = ((les.ab3$T0.Perimeter-les.ab3$T2.Perimeter)/(les.ab3$T2.Julian-les.ab3$T0.Julian))
les.ab3$Total.APRR = ((les.ab3$T0.Perimeter-les.ab3$T3.Perimeter)/(les.ab3$T3.Julian-les.ab3$T0.Julian))

les.ab3$ASA.P.R.T1 = ((les.ab3$SA.P.1-les.ab3$SA.P.0)/(les.ab3$T1.Julian-les.ab3$T0.Julian))
les.ab3$ASA.P.R.T2 = ((les.ab3$SA.P.2-les.ab3$SA.P.1)/(les.ab3$T2.Julian-les.ab3$T1.Julian))
les.ab3$ASA.P.R.T3 = ((les.ab3$SA.P.3-les.ab3$SA.P.2)/(les.ab3$T3.Julian-les.ab3$T2.Julian))
les.ab3$ASA.P.R.T12 = ((les.ab3$SA.P.2-les.ab3$SA.P.0)/(les.ab3$T2.Julian-les.ab3$T0.Julian))
les.ab3$ASA.P.R.Total = ((les.ab3$SA.P.3-les.ab3$SA.P.0)/(les.ab3$T3.Julian-les.ab3$T0.Julian))

##Stacking Columns to match All_Data
Ab_rec_t1 = data.frame(
  les.ab3$site, les.ab3$type, les.ab3$sitedepth, les.ab3$depth,
  les.ab3$T0.Julian, les.ab3$T1.Julian, 
  les.ab3$species, les.ab3$Colony, 
  les.ab3$Initial.lesion.size, les.ab3$T0.Perimeter, les.ab3$SA.P.0, les.ab3$T1.Lesion.Size, 
  les.ab3$T1.ALRR, les.ab3$T1.APRR, les.ab3$ASA.P.R.T1)

Ab_rec_t1$timepoint = as.factor("T1")
names(Ab_rec_t1) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony","initial_size","initial_perimeter","ASA_ratio_initial","lesion_size","ALRR","APRR","ASAP_rate","timepoint")

Ab_rec_t2 = data.frame(
  les.ab3$site, les.ab3$type, les.ab3$sitedepth, les.ab3$depth,
  les.ab3$T1.Julian, les.ab3$T2.Julian, 
  les.ab3$species, les.ab3$Colony, 
  les.ab3$Initial.lesion.size, les.ab3$T0.Perimeter, les.ab3$SA.P.0, les.ab3$T2.Lesion.Size, 
  les.ab3$T2.ALRR, les.ab3$T2.APRR, les.ab3$ASA.P.R.T2)

Ab_rec_t2$timepoint = as.factor("T2")
names(Ab_rec_t2) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony","initial_size","initial_perimeter","ASA_ratio_initial","lesion_size","ALRR","APRR","ASAP_rate","timepoint")

Ab_rec_t3 = data.frame(
  les.ab3$site, les.ab3$type, les.ab3$sitedepth, les.ab3$depth,
  les.ab3$T2.Julian, les.ab3$T3.Julian, 
  les.ab3$species, les.ab3$Colony, 
  les.ab3$Initial.lesion.size, les.ab3$T0.Perimeter, les.ab3$SA.P.0, les.ab3$T3.Lesion.Size, 
  les.ab3$T3.ALRR, les.ab3$T3.APRR, les.ab3$ASA.P.R.T3)

Ab_rec_t3$timepoint = as.factor("T3")
names(Ab_rec_t3) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony","initial_size","initial_perimeter","ASA_ratio_initial","lesion_size","ALRR","APRR","ASAP_rate","timepoint")

Ab_rec_t12 = data.frame(
  les.ab3$site, les.ab3$type, les.ab3$sitedepth, les.ab3$depth,
  les.ab3$T0.Julian, les.ab3$T2.Julian, 
  les.ab3$species, les.ab3$Colony, 
  les.ab3$Initial.lesion.size, les.ab3$T0.Perimeter, les.ab3$SA.P.0, les.ab3$T2.Lesion.Size, 
  les.ab3$T12.ALRR, les.ab3$T12.APRR, les.ab3$ASA.P.R.T12)

Ab_rec_t12$timepoint = as.factor("T12")
names(Ab_rec_t12) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony","initial_size","initial_perimeter","ASA_ratio_initial","lesion_size","ALRR","APRR","ASAP_rate","timepoint")

Ab_rec_tot = data.frame(
  les.ab3$site, les.ab3$type, les.ab3$sitedepth, les.ab3$depth,
  les.ab3$T0.Julian, les.ab3$T3.Julian, 
  les.ab3$species, les.ab3$Colony, 
  les.ab3$Initial.lesion.size, les.ab3$T0.Perimeter, les.ab3$SA.P.0, les.ab3$T3.Lesion.Size, 
  les.ab3$Total.ALRR, les.ab3$Total.APRR, les.ab3$ASA.P.R.Total)

Ab_rec_tot$timepoint = as.factor("Total")
names(Ab_rec_tot) = c("site","sitetype","sitedepth","depth","Julianstart","Julianend","species","colony","initial_size","initial_perimeter","ASA_ratio_initial","lesion_size","ALRR","APRR","ASAP_rate","timepoint")

Ab_rec_tall = rbind(Ab_rec_t1, Ab_rec_t2, Ab_rec_t3, Ab_rec_t12, Ab_rec_tot)

rec_all = cbind(Ab_rec_tall,Rel_rec_tall)
rec_all2 = data.frame(rec_all[,-c(17:24)])
rec_all3 = data.frame(rec_all2[,-c(19,20,26)])


######PIGMENTATION RECOVERY RATE (IRR AND/OR PRR)
setwd("H:/Backups/UVI/Smith Lab/Thesis/R Analysis")
les.col = read.csv("Lesion_color.csv", header = T, na.strings = c('NA'," NA", "NA ", "#VALUE!","#DIV/0!",""))

les.col2 = data.frame(les.col[1:48,])
les.col2 = data.frame(les.col2[-23,])


#removing unnecessary columns
les.col3 = les.col2[,-c(1,6,8,10,12,39:44)]

names(les.col3)[1] = "site"
names(les.col3)[2] = "type"
names(les.col3)[3] = "depth"
names(les.col3)[4] = "sitedepth"
names(les.col3)[9] = "species"

les.col3$T1.IRR = ((les.col3$T1.Intensity/les.col3$T1.Reference)-(les.col3$T0.Intensity/les.col3$T0.Reference))/(les.col3$T1.Julian-les.col3$T0.Julian)
les.col3$T2.IRR = ((les.col3$T2.Intensity/les.col3$T2.Reference)-(les.col3$T1.Intensity/les.col3$T1.Reference))/(les.col3$T2.Julian-les.col3$T1.Julian)
les.col3$T3.IRR = ((les.col3$T3.Intensity/les.col3$T3.Reference)-(les.col3$T2.Intensity/les.col3$T2.Reference))/(les.col3$T3.Julian-les.col3$T2.Julian)
les.col3$T12.IRR = ((les.col3$T2.Intensity/les.col3$T2.Reference)-(les.col3$T0.Intensity/les.col3$T0.Reference))/(les.col3$T2.Julian-les.col3$T0.Julian)
les.col3$Tot.IRR = ((les.col3$T3.Intensity/les.col3$T3.Reference)-(les.col3$T0.Intensity/les.col3$T0.Reference))/(les.col3$T3.Julian-les.col3$T0.Julian)

#Stacking to tall
Col_rec_T1 = data.frame(
  les.col3$site, les.col3$type, les.col3$sitedepth, les.col3$depth,
  les.col3$T0.Julian, les.col3$T1.Julian,
  les.col3$species, les.col3$Colony,
  les.col3$T0.Reference, les.col3$T0.Value, les.col3$T1.Value, les.col3$T1.IRR)

Col_rec_T1$timepoint = as.factor("T1")
names(Col_rec_T1) = c("site","type","sitedepth","depth","Julianstart","Julianend","species","colony","Reference","initial_intensity","timepoint_intensity","IRR","timepoint")

Col_rec_T2 = data.frame(
  les.col3$site, les.col3$type, les.col3$sitedepth, les.col3$depth,
  les.col3$T1.Julian, les.col3$T2.Julian,
  les.col3$species, les.col3$Colony,
  les.col3$T0.Reference, les.col3$T0.Value, les.col3$T2.Value, les.col3$T2.IRR)

Col_rec_T2$timepoint = as.factor("T2")
names(Col_rec_T2) = c("site","type","sitedepth","depth","Julianstart","Julianend","species","colony","Reference","initial_intensity","timepoint_intensity","IRR","timepoint")


Col_rec_T3 = data.frame(
  les.col3$site, les.col3$type, les.col3$sitedepth, les.col3$depth,
  les.col3$T2.Julian, les.col3$T3.Julian,
  les.col3$species, les.col3$Colony,
  les.col3$T0.Reference, les.col3$T0.Value, les.col3$T3.Value, les.col3$T3.IRR)

Col_rec_T3$timepoint = as.factor("T3")
names(Col_rec_T3) = c("site","type","sitedepth","depth","Julianstart","Julianend","species","colony","Reference","initial_intensity","timepoint_intensity","IRR","timepoint")

Col_rec_T12 = data.frame(
  les.col3$site, les.col3$type, les.col3$sitedepth, les.col3$depth,
  les.col3$T0.Julian, les.col3$T2.Julian,
  les.col3$species, les.col3$Colony,
  les.col3$T0.Reference, les.col3$T0.Value, les.col3$T2.Value, les.col3$T12.IRR)

Col_rec_T12$timepoint = as.factor("T12")
names(Col_rec_T12) = c("site","type","sitedepth","depth","Julianstart","Julianend","species","colony","Reference","initial_intensity","timepoint_intensity","IRR","timepoint")

Col_rec_Tot = data.frame(
  les.col3$site, les.col3$type, les.col3$sitedepth, les.col3$depth,
  les.col3$T0.Julian, les.col3$T3.Julian,
  les.col3$species, les.col3$Colony,
  les.col3$T0.Reference, les.col3$T0.Value, les.col3$T3.Value, les.col3$Tot.IRR)

Col_rec_Tot$timepoint = as.factor("Total")
names(Col_rec_Tot) = c("site","type","sitedepth","depth","Julianstart","Julianend","species","colony","Reference","initial_intensity","timepoint_intensity","IRR","timepoint")

Col_rec_tall = rbind(Col_rec_T1,Col_rec_T2,Col_rec_T3,Col_rec_T12,Col_rec_Tot)

Les.all = cbind(Ab_rec_tall,Rel_rec_tall,Col_rec_tall)
Les.all2 = data.frame(Les.all[,c(1:16,25:33,43:46)])



#################FULL DATA MERGE- PRE-EMPTIVE MANIPULATION TO MULTIPLE SHAPES/FORMATIONS
All_dat = merge(Les.all2, samp_temp.tall, by = c("site","timepoint"))
All_dat = merge(All_dat, samp_bov.tall, by = c("site","timepoint"))
All_dat$RLRR = All_dat$RLRR*100
All_dat$species = revalue(All_dat$species, c("OFRA" = "Orbicella franksi", "AL" = "Agaricia lamarcki"))
All_dat$IRR = All_dat$IRR*(-100)
All_dat.spec = split(All_dat, All_dat$species)
All_dat.spec$`Agaricia lamarcki`$Spec_value = as.numeric(0)
All_dat.spec$`Orbicella franksi`$Spec_value = as.numeric(1)
All_dat = rbind(All_dat.spec$`Agaricia lamarcki`,All_dat.spec$`Orbicella franksi`)

#########Timepoints only- no long periods##########
All_dat.tp1 = data.frame(All_dat[All_dat$timepoint == "T1",])
All_dat.tp2 = data.frame(All_dat[All_dat$timepoint == "T2",])
All_dat.tp3 = data.frame(All_dat[All_dat$timepoint == "T3",])
All_dat.total = data.frame(All_dat[All_dat$timepoint == "Total",])
All_dat.total.spec = split(All_dat.total, All_dat.total$species)
All_dat.tp = data.frame(rbind(All_dat.tp1,All_dat.tp2,All_dat.tp3))
All_dat.tp.spec = split(All_dat.tp, All_dat.tp$species)
All_dat.meso = data.frame(All_dat[All_dat$sitetype == "Mesophotic Site",])
All_dat.shal = data.frame(All_dat[All_dat$sitetype == "Shallow Site",])

#########REMOVE NA ROWS FOR RESPECTIVE ANALYSES
All_dat.total.NA.ALRR = data.frame(All_dat.total[!is.na( All_dat.total$ALRR),])
All_dat.total.NA.RLRR = data.frame(All_dat.total[!is.na( All_dat.total$RLRR),])
All_dat.total.NA.IRR = data.frame(All_dat.total[!is.na( All_dat.total$IRR),])

All_dat.tot.NA.ALRR.AL = data.frame(All_dat.total.spec$`Agaricia lamarcki`[!is.na( All_dat.total.spec$`Agaricia lamarcki`$ALRR),])
All_dat.tot.NA.ALRR.OFRA = data.frame(All_dat.total.spec$`Orbicella franksi`[!is.na( All_dat.total.spec$`Orbicella franksi`$ALRR),])
All_dat.tot.NA.RLRR.AL = data.frame(All_dat.total.spec$`Agaricia lamarcki`[!is.na( All_dat.total.spec$`Agaricia lamarcki`$RLRR),])
All_dat.tot.NA.RLRR.OFRA = data.frame(All_dat.total.spec$`Orbicella franksi`[!is.na( All_dat.total.spec$`Orbicella franksi`$RLRR),])
All_dat.tot.NA.IRR.AL = data.frame(All_dat.total.spec$`Agaricia lamarcki`[!is.na( All_dat.total.spec$`Agaricia lamarcki`$IRR),])
All_dat.tot.NA.IRR.OFRA = data.frame(All_dat.total.spec$`Orbicella franksi`[!is.na( All_dat.total.spec$`Orbicella franksi`$IRR),])

All_dat.tp.NA.ALRR = data.frame(All_dat.tp[!is.na( All_dat.tp$ALRR),])
All_dat.tp.NA.RLRR = data.frame(All_dat.tp[!is.na( All_dat.tp$RLRR),])
All_dat.tp.NA.IRR = data.frame(All_dat.tp[!is.na( All_dat.tp$IRR),])

All_dat.tp.NA.ALRR.AL = data.frame(All_dat.tp.spec$`Agaricia lamarcki`[!is.na( All_dat.tp.spec$`Agaricia lamarcki`$ALRR),])
All_dat.tp.NA.ALRR.OFRA = data.frame(All_dat.tp.spec$`Orbicella franksi`[!is.na( All_dat.tp.spec$`Orbicella franksi`$ALRR),])
All_dat.tp.NA.RLRR.AL = data.frame(All_dat.tp.spec$`Agaricia lamarcki`[!is.na( All_dat.tp.spec$`Agaricia lamarcki`$RLRR),])
All_dat.tp.NA.RLRR.OFRA = data.frame(All_dat.tp.spec$`Orbicella franksi`[!is.na( All_dat.tp.spec$`Orbicella franksi`$RLRR),])
All_dat.tp.NA.IRR.AL = data.frame(All_dat.tp.spec$`Agaricia lamarcki`[!is.na( All_dat.tp.spec$`Agaricia lamarcki`$IRR),])
All_dat.tp.NA.IRR.OFRA = data.frame(All_dat.tp.spec$`Orbicella franksi`[!is.na( All_dat.tp.spec$`Orbicella franksi`$IRR),])


# #################################PCA OF PHYSICAL ENVIRONMENTAL FACTORS
# 
# # PCA WITH TOTAL SAMPLE PERIOD
# sites.phys = All_dat.total[,c(5,30:37)]
# 
# corr = rcorr(as.matrix(sites.phys))
# corr
# 
# sites.phys.pca = rda(sites.phys, scale = TRUE)
# summary(sites.phys.pca)
# 
# biplot(sites.phys.pca)
# plot(sites.phys.pca, scaling = 2)
# 
# smry = summary(sites.phys.pca)
# df1 = data.frame(smry$sites[,1:2])
# df2 = data.frame(smry$species[,1:2])
# 
# df1$site = sapply(strsplit(as.character(row.names(df1)), "_"), "[[",1)
# 
# df1$location = sapply(strsplit(as.character(All_dat.total$site), "_"),"[[",1)
# 
# df1$type = sapply(strsplit(as.character(All_dat.total$sitetype), "_"),"[[",1)
# 
# df2$variable = sapply(strsplit(as.character(row.names(df2)), "_"), "[[",1)
# 
# ggplot(df1, aes(x = PC1, y = PC2))+ 
#   geom_point(cex = 5.5, aes(shape = type, color = location))+ #cex is size
#   theme_bw()+
#   geom_text_repel(data = df2, aes(x = PC1, y = PC2, label = variable), size = 8)+
#   geom_segment(data = df2, aes(x=0, xend = PC1, y=0, yend= PC2), color = "blue", arrow = arrow(length=unit(0.01, "npc")))+
#   scale_colour_discrete(
#     name = "Site", labels = c("Buck Island", "College Shoal","Flat Cay","Grammanik Bank","Hind Bank","South Capella Deep", "South Capella Shallow","Seahorse Cottage Shoal"))+
#   scale_shape_discrete(
#     name = "Site Designation")+
#   xlab(label = "PC1 (81%)") + ylab(label = "PC2 (13%)")+
#   ggtitle(label = "PCA of Environmental Variables, Type II scaling")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=18), legend.text=element_text(size=18))
# 
# 
# #DEPTH REMOVED TOTAL SAMPLING PERIOD PCA
# sites.phys = All_dat.total[,c(30:37)]
# 
# corr = rcorr(as.matrix(sites.phys))
# corr
# 
# sites.phys.pca = rda(sites.phys, scale = TRUE)
# summary(sites.phys.pca)
# 
# biplot(sites.phys.pca)
# plot(sites.phys.pca, scaling = 2)
# 
# smry = summary(sites.phys.pca)
# df1 = data.frame(smry$sites[,1:2])
# df2 = data.frame(smry$species[,1:2])
# 
# df1$site = sapply(strsplit(as.character(row.names(df1)), "_"), "[[",1)
# 
# df1$location = sapply(strsplit(as.character(All_dat.total$site), "_"),"[[",1)
# 
# df1$type = sapply(strsplit(as.character(All_dat.total$sitetype), "_"),"[[",1)
# 
# df2$variable = sapply(strsplit(as.character(row.names(df2)), "_"), "[[",1)
# 
# ggplot(df1, aes(x = PC1, y = PC2))+ 
#   geom_point(cex = 5.5, aes(shape = type, color = location))+ #cex is size
#   theme_bw()+
#   geom_text_repel(data = df2, aes(x = PC1, y = PC2, label = variable), size = 8)+
#   geom_segment(data = df2, aes(x=0, xend = PC1, y=0, yend= PC2), color = "blue", arrow = arrow(length=unit(0.01, "npc")))+
#   scale_colour_discrete(
#     name = "Site", labels = c("Buck Island", "College Shoal","Flat Cay","Grammanik Bank","Hind Bank","South Capella Deep", "South Capella Shallow","Seahorse Cottage Shoal"))+
#   scale_shape_discrete(
#     name = "Site Designation")+
#   xlab(label = "PC1 (81%)") + ylab(label = "PC2 (13%)")+
#   ggtitle(label = "PCA of Environmental Variables Without Depth, Type II scaling")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=18), legend.text=element_text(size=18))
# 
# 
# ############SEPARATE TIMEPOINT ANALYSIS
# sites.phys = All_dat.tp[,c(30:37)]
# 
# corr = rcorr(as.matrix(sites.phys))
# corr
# 
# sites.phys.pca = rda(sites.phys, scale = TRUE)
# summary(sites.phys.pca)
# 
# biplot(sites.phys.pca)
# plot(sites.phys.pca, scaling = 2)
# 
# smry = summary(sites.phys.pca)
# df1 = data.frame(smry$sites[,1:2])
# df2 = data.frame(smry$species[,1:2])
# 
# df1$site = sapply(strsplit(as.character(row.names(df1)), "_"), "[[",1)
# 
# df1$timepoint = sapply(strsplit(as.character(All_dat.tp$timepoint), "_"), "[[",1)
# 
# df1$location = sapply(strsplit(as.character(All_dat.tp$site), "_"),"[[",1)
# 
# df1$type = sapply(strsplit(as.character(All_dat.tp$sitetype), "_"),"[[",1)
# 
# df2$variable = sapply(strsplit(as.character(row.names(df2)), "_"), "[[",1)
# 
# ggplot(df1, aes(x = PC1, y = PC2))+ 
#   geom_point(cex = 5.5, aes(shape = type, color = location))+ #cex is size
#   facet_grid(~timepoint)+
#   theme_bw()+
#   geom_text_repel(data = df2, aes(x = PC1, y = PC2, label = variable), size = 6)+
#   geom_segment(data = df2, aes(x=0, xend = PC1, y=0, yend= PC2), color = "blue", arrow = arrow(length=unit(0.01, "npc")))+
#   scale_colour_discrete(
#     name = "Site", labels = c("Buck Island", "College Shoal","Flat Cay","Grammanik Bank","Hind Bank","South Capella Deep", "South Capella Shallow","Seahorse Cottage Shoal"))+
#   scale_shape_discrete(
#     name = "Site Designation")+
#   xlab(label = "PC1 (81%)") + ylab(label = "PC2 (13%)")+
#   ggtitle(label = "PCA of Environmental Variables through time, Type II scaling")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=18), legend.text=element_text(size=18))






##########################################################ALRR#################################################
# ALRR.lm.tot.desat = lm(ALRR ~ depth*Spec_value, data = All_dat.total)
# summary(ALRR.lm.tot.desat)
# confint(ALRR.lm.tot.desat)
# 
# ALRR.lm.tot = lm(ALRR ~ depth*Spec_value + avgbov + 
#                    avgtemp +
#                    initial_size + 
#                    colony_size +
#                    maxtemp + 
#                    mintemp + 
#                    stdtemp + 
#                    maxbov + 
#                    minbov + 
#                    stdbov,
#                  data = All_dat.total)
# 
# summary(ALRR.lm.tot)
# plot(ALRR.lm.tot)
# anova(ALRR.lm.tot)
# confint(ALRR.lm.tot)
# step(ALRR.lm.tot, direction = 'both')
# 
# ALRR.lm.tot.fitted = lm(ALRR ~ depth*Spec_value + 
#                            avgbov+
#                            avgtemp+
#                            initial_size+
#                            stdtemp+
#                            maxbov,
#                          data = All_dat.total)
# 
# summary(ALRR.lm.tot.fitted)
# plot(ALRR.lm.tot.fitted)
# anova(ALRR.lm.tot.fitted,ALRR.lm.tot)
# confint(ALRR.lm.tot.fitted)
# 
# ggplot(All_dat.total.NA.ALRR, aes(x = depth, y = ALRR, color = sitetype))+
#   facet_wrap(~species)+
#   theme(strip.text.x = element_text(size = 20))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(ALRR.lm.tot.fitted), group = species), cex = 2)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
#   xlab(label = "Depth (m)") + ylab(label = "Absolute Lesion Recovery Rate (cm^2/day)")+
#   ggtitle(label = "Recovery Rate of Corals Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

#AL
AL.ALRR.lm = lm(ALRR~depth, data = All_dat.total.spec$`Agaricia lamarcki`)
summary(AL.ALRR.lm)

ALRR.lm.tot.AL = lm(ALRR ~ depth + avgbov + 
                  avgtemp +
                  colony_size +
                  initial_size +
                  maxtemp + 
                  mintemp + 
                  stdtemp + 
                  maxbov + 
                  minbov + 
                  stdbov,
                 data = All_dat.total.spec$`Agaricia lamarcki`)

summary(ALRR.lm.tot.AL)
#Initial Size and colony size
anova(AL.ALRR.lm, ALRR.lm.tot.AL)

ALRR.lm.tot.AL.fitted = lm(ALRR ~ depth + initial_size + colony_size, data = All_dat.total.spec$`Agaricia lamarcki`)
summary(ALRR.lm.tot.AL.fitted)
anova(AL.ALRR.lm,ALRR.lm.tot.AL.fitted)

ALRR.AL.graph = ggplot(All_dat.tot.NA.ALRR.AL, aes(x = depth, y = ALRR, color = sitetype))+
  theme(strip.text.x = element_text(size = 16))+
  geom_hline(yintercept = 0, color = 'black')+
  geom_line(aes(x = depth, y = predict(AL.ALRR.lm), group = 'sitetype'), cex = 2, alpha = 0.4)+
  geom_line(aes(x = depth, y = predict(ALRR.lm.tot.AL.fitted), group = 'sitetype'), cex = 2)+
  geom_point(stat = "identity")+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = "Depth (m)") + ylab(label = "Regeneration (cm^2/day)")+ xlim(12,43)+ ylim(-0.01,0.06)+
  ggtitle(label = "Tissue Regeneration Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), text = element_text(size=16), legend.text=element_text(size=16), legend.position = 'none')

ALRR.AL.graph
#OFRA
OFRA_no_outlier = data.frame(All_dat.total.spec$`Orbicella franksi`)
OFRA_no_outlier = data.frame(OFRA_no_outlier[!is.na(OFRA_no_outlier$ALRR),])
OFRA_no_outlier = data.frame(OFRA_no_outlier[!(OFRA_no_outlier$site=="GKT"),])

OFRA.ALRR.lm = lm(ALRR~depth, data = OFRA_no_outlier)
summary(OFRA.ALRR.lm)
confint(OFRA.ALRR.lm)

ALRR.lm.tot.OFRA = lm(ALRR ~ depth + avgbov + 
                      avgtemp +
                      initial_size + 
                      colony_size +
                      maxtemp + 
                      mintemp + 
                      stdtemp + 
                      maxbov + 
                      minbov + 
                      stdbov,
                    data = OFRA_no_outlier)

summary(ALRR.lm.tot.OFRA)
confint(ALRR.lm.tot.OFRA)
#no significant figures
anova(OFRA.ALRR.lm, ALRR.lm.tot.OFRA)
#no significant additional terms

ALRR.OFRA.graph = ggplot(OFRA_no_outlier, aes(x = depth, y = ALRR, color = sitetype))+
  theme(strip.text.x = element_text(size = 16))+
  geom_hline(yintercept = 0, color = 'black')+
  geom_line(aes(x = depth, y = predict(OFRA.ALRR.lm), group = 'sitetype'), cex = 2)+
  geom_point(stat = "identity")+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = " ") + ylab(label = "Regeneration (cm^2/day)")+ xlim(12,43)+ ylim(-0.01,0.06)+
  ggtitle(label = "Tissue Regeneration Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), text = element_text(size=16), legend.text=element_text(size=16), legend.position = 'none')

ALRR.OFRA.graph
# ##########################Mixed Effects Analysis#################################
# ##Together##
# #random intercept model
# ALRR.lmm.random = lmer(ALRR ~ 1 +  (1 | site) + (1 |timepoint) + (1|sitetype), REML = FALSE, data = All_dat.tp)
# summary(ALRR.lmm.random)
# 
# #Level 1 predictors
# ALRR.lmm.tp = lmer(ALRR ~ 1+ depth*Spec_value + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp)
# summary(ALRR.lmm.tp)
# anova(ALRR.lmm.random, ALRR.lmm.tp)
# rand(ALRR.lmm.tp)
# 
# ALRR.lmm.tp.sat = lmer(ALRR ~ 1+ depth*Spec_value + (1 | site) + (1 |timepoint) + 
#                        avgbov + 
#                        avgtemp +
#                        initial_size + 
#                        colony_size +
#                        maxtemp + 
#                        mintemp + 
#                        stdtemp + 
#                        maxbov + 
#                        minbov + 
#                        stdbov, 
#                        REML = FALSE, data = All_dat.tp)
# 
# step(ALRR.lmm.tp.sat, direction = 'both')
# summary(ALRR.lmm.tp.sat)
# anova(ALRR.lmm.tp.sat)
# anova(ALRR.lmm.tp, ALRR.lmm.tp.sat)
# rand(ALRR.lmm.tp.sat)
# Anova(ALRR.lmm.tp.sat)
# analyze(ALRR.lmm.tp.sat)
# 
# ALRR.lmm.tp.fitted = lmer(ALRR ~ 1 + depth*Spec_value + (1 |timepoint) +
#                           avgbov + 
#                           maxbov+
#                           stdbov, 
#                           REML = FALSE, data = All_dat.tp)
# 
# summary(ALRR.lmm.tp.fitted)
# anova(ALRR.lmm.tp.sat, ALRR.lmm.tp.fitted)
# confint(ALRR.lmm.tp.fitted)
# analyze(ALRR.lmm.tp.fitted)
# rand(ALRR.lmm.tp.fitted)
# 
# ggplot(All_dat.tp.NA.ALRR, aes(x = depth, y = ALRR, color = sitetype))+
#   facet_grid(species~timepoint)+
#   theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(ALRR.lmm.tp.fitted), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
#   xlab(label = "Depth (m)") + ylab(label = "Absolute Lesion Recovery Rate (cm^2/day)")+
#   ggtitle(label = "Recovery Rate of Corals Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

##########################################################RLRR#################################################
# RLRR.lm.tot.desat = lm(RLRR ~ depth*Spec_value, data = All_dat.total)
# summary(RLRR.lm.tot.desat)
# 
# RLRR.lm.tot = lm(RLRR ~ depth*Spec_value + avgbov + 
#                    avgtemp +
#                    initial_size + 
#                    colony_size +
#                    maxtemp + 
#                    mintemp + 
#                    stdtemp + 
#                    maxbov + 
#                    minbov + 
#                    stdbov,
#                  data = All_dat.total)
# 
# summary(RLRR.lm.tot)
# plot(RLRR.lm.tot)
# Anova(RLRR.lm.tot)
# anova(RLRR.lm.tot)
# anova(RLRR.lm.tot.desat,RLRR.lm.tot)
# step(RLRR.lm.tot, direction = 'backward')
# 
# RLRR.lm.tot.fitted2 = lm(RLRR ~ depth*Spec_value +
#                            initial_size,
#                          data = All_dat.total)
# 
# summary(RLRR.lm.tot.fitted2)
# anova(RLRR.lm.tot.fitted2)
# anova(RLRR.lm.tot,RLRR.lm.tot.fitted2)
# confint(RLRR.lm.tot.fitted2)
# 
# 
# All_dat.total.NA.RLRR = data.frame(All_dat.total[!is.na( All_dat.total$RLRR),])
# 
# ggplot(All_dat.total.NA.RLRR, aes(x = depth, y = RLRR, color = sitetype))+
#   facet_wrap(~species)+
#   theme(strip.text.x = element_text(size = 20))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(RLRR.lm.tot.fitted2), group = species), cex = 2)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
#   xlab(label = "Depth (m)") + ylab(label = "Relative Lesion Recovery Rate (%/day)")+
#   ggtitle(label = "Recovery Rate of Corals Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))


#AL-specific model
RLRR.lm.tot.AL.desat = lm(RLRR~depth, data = All_dat.total.spec$`Agaricia lamarcki`)
summary(RLRR.lm.tot.AL.desat)

RLRR.lm.tot.AL = lm(RLRR ~ depth + avgbov + 
                      avgtemp +
                      initial_size + 
                      colony_size +
                      maxtemp + 
                      mintemp + 
                      stdtemp + 
                      maxbov + 
                      minbov + 
                      stdbov,
                    data = All_dat.total.spec$`Agaricia lamarcki`)

summary(RLRR.lm.tot.AL)
#Colony size
anova(RLRR.lm.tot.AL, RLRR.lm.tot.AL.desat)
#significant term

RLRR.lm.tot.AL.fitted = lm(RLRR ~ depth +  
                      colony_size,
                    data = All_dat.total.spec$`Agaricia lamarcki`)
summary(RLRR.lm.tot.AL.fitted)
anova(RLRR.lm.tot.AL.fitted, RLRR.lm.tot.AL.desat)

All_dat.total.spec.NA.AL.RLRR = data.frame(All_dat.total.spec$`Agaricia lamarcki`[!is.na( All_dat.total.spec$`Agaricia lamarcki`$RLRR),])

RLRR.AL.graph = ggplot(All_dat.total.spec.NA.AL.RLRR, aes(x = depth, y = RLRR, color = sitetype))+
  geom_point(stat = "identity")+
  # geom_line(aes(x = depth, y = predict(RLRR.lm.tot.AL)), cex = 2)+
  # geom_line(aes(x = depth, y = predict(RLRR.lm.tot.AL.desat), group = species), cex = 2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = 'black')+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = "Depth (m)") + ylab(label = "Closure (% Recovered/day)")+ xlim(12,43)+ ylim(-0.5,1)+
  ggtitle(label = "Lesion Closure Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), text = element_text(size=16), legend.text=element_text(size=16), legend.position = 'none')

RLRR.AL.graph
#OFRA
OFRA_no_outlier = data.frame(All_dat.total.spec$`Orbicella franksi`)
OFRA_no_outlier = data.frame(OFRA_no_outlier[!is.na(OFRA_no_outlier$RLRR),])

RLRR.lm.tot.OFRA.desat = lm(RLRR~depth, data = OFRA_no_outlier)
summary(RLRR.lm.tot.OFRA.desat)
#no significant effect
confint(RLRR.lm.tot.OFRA.desat)

RLRR.lm.tot.OFRA = lm(RLRR ~ depth + avgbov + 
                        avgtemp +
                        colony_size +
                        maxtemp + 
                        mintemp + 
                        stdtemp + 
                        maxbov + 
                        minbov + 
                        stdbov,
                      data = OFRA_no_outlier)

summary(RLRR.lm.tot.OFRA)
#no significant terms
anova(RLRR.lm.tot.OFRA.desat,RLRR.lm.tot.OFRA)
#no siginificant additional effect

RLRR.OFRA.graph = ggplot(OFRA_no_outlier, aes(x = depth, y = RLRR, color = sitetype))+
  geom_hline(yintercept = 0, color = 'black')+
  # geom_line(aes(x = depth, y = predict(RLRR.lm.tot.OFRA.desat), group = 'sitetype'), alpha = 0.4, cex = 2)+
  geom_point(stat = "identity")+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = " ") + ylab(label = "Closure (% Recovered/day)")+ xlim(12,43)+ ylim(-0.5,1)+
  ggtitle(label = "Lesion Closure Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), text = element_text(size=16), legend.text=element_text(size=16), legend.position = 'none')

RLRR.OFRA.graph
#Combine to one figure
# AL.lm.plot = 
#   ggplot(All_dat.tot.NA.RLRR.AL, aes(x = depth, y = RLRR, color = sitetype))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(RLRR.lm.tot.AL), group = species), cex = 1.5, alpha = 0.3)+
#   geom_line(aes(x = depth, y = predict(RLRR.lm.tot.AL.desat), group = species), cex = 1.5)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Depth Category") +
#   theme(panel.border = element_rect(colour = "black", fill=NA),
#         text = element_text(size = 20),
#         legend.justification=c(1,1), 
#         legend.position=c(0.558, 0.147),
#         legend.background = element_rect(size = 0.5, colour = 'black'),
#         legend.title = element_text(size = 20), 
#         legend.text = element_text(size = 20))+
#   ylim(-0.5, 1.0)+
#   xlab(label = "Depth (m)") + ylab(label = "Relative Lesion Recovery Rate (%/day)")+
#   ggtitle(label = "Recovery Rate of Agaricia lamarcki")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))
# 
# OFRA.lm.plot = ggplot(OFRA_no_outlier, aes(x = depth, y = RLRR, color = sitetype))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(RLRR.lm.tot.OFRA.fitted), group = species), cex = 1.5, alpha = 0.3)+
#   geom_line(aes(x = depth, y = predict(RLRR.lm.tot.OFRA.desat), group = species), cex = 1.5)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), guide = FALSE)+
#   theme(panel.border = element_rect(colour = "black", fill=NA),
#         text = element_text(size = 20))+
#   ylim(-0.5, 1.0)+
#   xlab(label = "Depth (m)") + ylab(label = " ")+
#   ggtitle(label = "Recovery Rate of Orbicella franksi")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = 'none')
# 
# grid.arrange(AL.lm.plot,OFRA.lm.plot, ncol = 2)
# ##########################Mixed Effects Analysis#################################
# ##Together##
# #random intercept model
# RLRR.lmm.random = lmer(RLRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp)
# summary(RLRR.lmm.random)
# rand(RLRR.lmm.random)
# 
# #Level 1 predictors
# RLRR.lmm.tp = lmer(RLRR ~ 1+ depth*Spec_value + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp)
# summary(RLRR.lmm.tp)
# 
# anova(RLRR.lmm.random, RLRR.lmm.tp)
# 
# RLRR.lmm.tp.sat = lmer(RLRR ~ 1+ depth*Spec_value + (1 | site) + (1 |timepoint) + 
#                          avgbov + 
#                          avgtemp +
#                          initial_size + 
#                          colony_size +
#                          maxtemp + 
#                          mintemp + 
#                          stdtemp + 
#                          maxbov + 
#                          minbov + 
#                          stdbov, 
#                        REML = FALSE, data = All_dat.tp)
# 
# step(RLRR.lmm.tp.sat, direction = 'both')
# summary(RLRR.lmm.tp.sat)
# anova(RLRR.lmm.tp.sat)
# anova(RLRR.lmm.tp, RLRR.lmm.tp.sat)
# rand(RLRR.lmm.tp.sat)
# analyze(RLRR.lmm.tp.sat)
# 
# ##AL-Specific Mixed Effects
# RLRR.lmm.random.AL = lmer(RLRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(RLRR.lmm.random.AL)
# rand(RLRR.lmm.random.AL)
# 
# #Level 1 predictors
# RLRR.lmm.tp.AL = lmer(RLRR ~ 1+ depth + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(RLRR.lmm.tp.AL)
# 
# anova(RLRR.lmm.random.AL, RLRR.lmm.tp.AL)
# 
# RLRR.lmm.tp.sat.AL = lmer(RLRR ~ 1 + avgbov + 
#                             avgtemp +
#                             maxtemp + 
#                             mintemp + 
#                             stdtemp + 
#                             maxbov + 
#                             minbov + 
#                             stdbov + (1 | site) + (1 |timepoint) + 
#                             + depth +
#                             initial_size + 
#                             colony_size, 
#                           REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# 
# step(RLRR.lmm.tp.sat.AL, direction = 'backward')
# summary(RLRR.lmm.tp.sat.AL)
# anova(RLRR.lmm.tp.AL, RLRR.lmm.tp.sat.AL)
# anova(RLRR.lmm.tp.sat.AL)
# rand(RLRR.lmm.tp.sat.AL)
# analyze(RLRR.lmm.tp.sat.AL)
# 
# RLRR.lmm.tp.fitted.AL = lmer(RLRR ~ 1 + depth + (1 |timepoint) + 
#                                avgbov, 
#                                # avgtemp +
#                                # initial_size,
#                              # maxtemp + 
#                              # mintemp + 
#                              # stdtemp + 
#                              # maxbov + 
#                              # minbov + 
#                              # stdbov, 
#                              REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(RLRR.lmm.tp.fitted.AL)
# anova(RLRR.lmm.tp.sat.AL,RLRR.lmm.tp.fitted.AL)
# analyze(RLRR.lmm.tp.fitted.AL)
# Anova(RLRR.lmm.tp.fitted.AL)
# rand(RLRR.lmm.tp.fitted.AL)
# 
# All_dat.tp.NA.RLRR.AL = data.frame(All_dat.tp.spec$`Agaricia lamarcki`[!is.na( All_dat.tp.spec$`Agaricia lamarcki`$RLRR),])
# 
# ggplot(All_dat.tp.NA.RLRR.AL, aes(x = depth, y = RLRR, color = sitetype))+
#   facet_grid(~timepoint)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(RLRR.lmm.tp.fitted.AL), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_discrete(name = "Depth Category", labels = c("Mesophotic","Shallow"))+
#   xlab(label = "Depth (m)") + ylab(label = "Relative Lesion Recovery Rate (%/day)")+
#   ggtitle(label = "Recovery Rate of Agaricia lamarcki Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# ##OFRA-Specific Mixed Effects
# RLRR.lmm.random.OFRA = lmer(RLRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(RLRR.lmm.random.OFRA)
# 
# #Level 1 predictors
# RLRR.lmm.tp.OFRA = lmer(RLRR ~ 1+ depth + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(RLRR.lmm.tp.OFRA)
# 
# anova(RLRR.lmm.random.OFRA, RLRR.lmm.tp.OFRA)
# 
# RLRR.lmm.tp.sat.OFRA = lmer(RLRR ~ 1+ depth + (1 | site) + (1 |timepoint) +
#                               avgbov + 
#                               avgtemp +
#                               initial_size + 
#                               colony_size +
#                               maxtemp + 
#                               mintemp + 
#                               stdtemp + 
#                               maxbov + 
#                               minbov + 
#                               stdbov +
#                               initial_size,
#                             REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# 
# step(RLRR.lmm.tp.sat.OFRA, direction = 'backward')
# summary(RLRR.lmm.tp.sat.OFRA)
# anova(RLRR.lmm.tp.OFRA, RLRR.lmm.tp.sat.OFRA)
# anova(RLRR.lmm.tp.sat.OFRA)
# rand(RLRR.lmm.tp.sat.OFRA)
# analyze(RLRR.lmm.tp.sat.OFRA)
# 
# RLRR.lmm.tp.fitted.OFRA = lmer(RLRR ~ 1 + depth + (1 |timepoint) + 
#                                  # avgbov + 
#                                  avgtemp,
#                                  # initial_size,
#                                  # maxtemp + 
#                                  # mintemp + 
#                                  # stdtemp + 
#                                  # maxbov + 
#                                  # minbov, 
#                                # stdbov, 
#                                REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(RLRR.lmm.tp.fitted.OFRA)
# anova(RLRR.lmm.tp.sat.OFRA,RLRR.lmm.tp.fitted.OFRA)
# analyze(RLRR.lmm.tp.fitted.OFRA)
# anova(RLRR.lmm.tp.fitted.OFRA)
# rand(RLRR.lmm.tp.fitted.OFRA)
# 
# All_dat.tp.NA.RLRR.OFRA = data.frame(All_dat.tp.spec$`Orbicella franksi`[!is.na( All_dat.tp.spec$`Orbicella franksi`$RLRR),])
# 
# ggplot(All_dat.tp.NA.RLRR.OFRA, aes(x = depth, y = RLRR, color = sitetype))+
#   facet_grid(~timepoint)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(RLRR.lmm.tp.fitted.OFRA), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_discrete(name = "Depth Category", labels = c("Mesophotic","Shallow"))+
#   xlab(label = "Depth (m)") + ylab(label = "Relative Lesion Recovery Rate (%/day)")+
#   ggtitle(label = "Recovery Rate of Orbicella franksi Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# 
# 
# 
# 
# 
# 
##########################################################IRR#################################################
# IRR.lm.tot.desat = lm(IRR ~ depth*Spec_value, data = All_dat.total)
# summary(IRR.lm.tot.desat)
# confint(IRR.lm.tot.desat)
# 
# IRR.lm.tot = lm(IRR ~ depth*Spec_value + avgbov + 
#                   avgtemp +
#                   initial_size + 
#                   colony_size +
#                   maxtemp + 
#                   mintemp + 
#                   stdtemp + 
#                   maxbov + 
#                   minbov + 
#                   stdbov,
#                 data = All_dat.total)
# 
# summary(IRR.lm.tot)
# plot(IRR.lm.tot)
# anova(IRR.lm.tot, IRR.lm.tot.desat)
# anova(IRR.lm.tot)
# step(IRR.lm.tot, direction = 'backward')
# 
# IRR.lm.tot.fitted = lm(IRR ~ depth*Spec_value + 
#                          avgtemp, 
#                        data = All_dat.total)
# 
# summary(IRR.lm.tot.fitted)
# anova(IRR.lm.tot.fitted)
# plot(IRR.lm.tot.fitted)
# anova(IRR.lm.tot.fitted,IRR.lm.tot)
# confint(IRR.lm.tot.fitted)
# 
# ggplot(All_dat.total.NA.IRR, aes(x = depth, y = IRR, color = sitetype))+
#   facet_wrap(~species)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lm.tot.fitted), group = species), cex = 2)+
#   scale_color_discrete(name = "Site Type")+
#   xlab(label = "Depth (m)") + ylab(label = "Pigmentation Lesion Recovery Rate (color value/day)")+
#   ggtitle(label = "Recovery Rate of Corals Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
#AL-specific model
IRR.lm.tot.AL.desat = lm(IRR~depth, data = All_dat.total.spec$`Agaricia lamarcki`)
summary(IRR.lm.tot.AL.desat)
#positive effect of depth
confint(IRR.lm.tot.AL.desat)

IRR.lm.tot.AL = lm(IRR ~ depth + avgbov + 
                     avgtemp +
                     initial_size + 
                     colony_size +
                     maxtemp + 
                     mintemp + 
                     stdtemp + 
                     maxbov + 
                     minbov + 
                     stdbov,
                   data = All_dat.total.spec$`Agaricia lamarcki`)

summary(IRR.lm.tot.AL)
#no significant terms
anova(IRR.lm.tot.AL.desat, IRR.lm.tot.AL)
#no significant additions
step(IRR.lm.tot.AL, direction = 'backward')

All_dat.total.spec.NA.AL.IRR = data.frame(All_dat.total.spec$`Agaricia lamarcki`[!is.na( All_dat.total.spec$`Agaricia lamarcki`$IRR),])

IRR.AL.graph = ggplot(All_dat.total.spec.NA.AL.IRR, aes(x = depth, y = IRR, color = sitetype))+
  geom_point(stat = "identity")+
  geom_hline(yintercept = 0, color = 'black')+
  geom_line(aes(x = depth, y = predict(IRR.lm.tot.AL.desat), group = species), cex = 2)+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = "Depth (m)") + ylab(label = "Pigmentation Recovery (% healthy color/day)")+ xlim(12,43) + ylim(-1.0,0.7)+
  ggtitle(label = "Pigmentation Recovery Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y = element_text(size = 12), legend.text=element_text(size=16), legend.position = 'none', axis.text.x = element_text(size = 13))

IRR.AL.graph
#OFRA
# OFRA_no_outlier = data.frame(All_dat.total.spec$`Orbicella franksi`)
# OFRA_no_outlier = data.frame(OFRA_no_outlier[!is.na(OFRA_no_outlier$IRR),])
# OFRA_no_outlier$IRR = OFRA_no_outlier$IRR*100


IRR.lm.tot.OFRA.desat = lm(IRR~depth, data = All_dat.tot.NA.IRR.OFRA)
summary(IRR.lm.tot.OFRA.desat)
#no effect of depth

IRR.lm.tot.OFRA = lm(IRR ~ depth + avgbov + 
                       avgtemp +
                       initial_size + 
                       colony_size +
                       maxtemp + 
                       mintemp + 
                       stdtemp + 
                       maxbov + 
                       minbov + 
                       stdbov,
                     data = All_dat.tot.NA.IRR.OFRA)

summary(IRR.lm.tot.OFRA)
#no effects
anova(IRR.lm.tot.OFRA.desat, IRR.lm.tot.OFRA)
#no significant additions
step(IRR.lm.tot.OFRA, direction = 'both')

IRR.lm.tot.OFRA.fitted = lm(IRR ~ depth + avgtemp, 
                            data = All_dat.tot.NA.IRR.OFRA)

summary(IRR.lm.tot.OFRA.fitted)
#significant effect of both
anova(IRR.lm.tot.OFRA.desat,IRR.lm.tot.OFRA.fitted)
#significantly improved model fit

IRR.lm.tot.OFRA.fitted2 = lm(IRR ~ avgbov + avgtemp + stdtemp, data = 
                               All_dat.tot.NA.IRR.OFRA)

summary(IRR.lm.tot.OFRA.fitted2)
#no significance
confint(IRR.lm.tot.OFRA.fitted2)
anova(IRR.lm.tot.OFRA.fitted, IRR.lm.tot.OFRA.fitted2)
anova(IRR.lm.tot.OFRA.desat, IRR.lm.tot.OFRA.fitted2)
#not significant

IRR.OFRA.graph = ggplot(All_dat.tot.NA.IRR.OFRA, aes(x = depth, y = IRR, color = sitetype))+
  geom_point(stat = "identity")+
  geom_hline(yintercept = 0, color = 'black')+
  # geom_line(aes(x = depth, y = predict(IRR.lm.tot.OFRA.desat), alpha = 0.3, group = species), cex = 2)+
  scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Site Type") +
  xlab(label = " ") + ylab(label = "Pigmentation Recovery (% Healthy color/day)")+ xlim(12,43) + ylim(-1.0, 0.7)+
  ggtitle(label = "Pigmentation Recovery Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.text.y = element_text(size=12), legend.text=element_text(size=16), legend.position = 'none', axis.text.x = element_text(size = 13))
IRR.OFRA.graph

Fig1 = grid.arrange(ALRR.OFRA.graph, RLRR.OFRA.graph, IRR.OFRA.graph, ALRR.AL.graph, RLRR.AL.graph, IRR.AL.graph, ncol = 3)
Fig1
# #MULTIPLOT OF INDEPENDENT MODELS
# AL.lm.plot = ggplot(All_dat.tot.NA.IRR.AL, aes(x = depth, y = IRR, color = sitetype))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lm.tot.AL.desat), group = species), cex = 1.5)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), name = "Depth Category") +
#   ylim(-0.010,0.007)+
#   theme(panel.border = element_rect(colour = "black", fill=NA),
#         text = element_text(size = 20),
#         legend.justification=c(1,1), 
#         legend.position=c(0.588, 0.147),
#         legend.background = element_rect(size = 0.5, colour = 'black'),
#         legend.title = element_text(size = 20), 
#         legend.text = element_text(size = 20))+
#   theme(panel.border = element_rect(colour = "black", fill=NA),
#         text = element_text(size = 20))+
#   xlab(label = "Depth (m)") + ylab(label = "Pigmentation Lesion Recovery Rate (% healthy pixel color/day)")+
#   ggtitle(label = "Pigmentation Recovery- A. lamarcki")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20))
# 
# OFRA.lm.plot = ggplot(OFRA_no_outlier, aes(x = depth, y = IRR, color = sitetype))+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lm.tot.OFRA.fitted), group = species), cex = 1.5)+
#   scale_color_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ), guide = FALSE) +
#   ylim(-0.010,0.007)+
#   theme(panel.border = element_rect(colour = "black", fill=NA),
#         text = element_text(size = 20))+
#   xlab(label = "Depth (m)") + ylab(label = " ")+
#   ggtitle(label = "Pigmentation Recovery- O. franksi")+
#   theme(plot.title = element_text(hjust = 0.5, size = 20))
# 
# grid.arrange(AL.lm.plot,OFRA.lm.plot, ncol = 2)
# 
# ##########################Mixed Effects Analysis#################################
# ##Together##
# #random intercept model
# IRR.lmm.random = lmer(IRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp)
# summary(IRR.lmm.random)
# 
# #Level 1 predictors
# IRR.lmm.tp = lmer(IRR ~ 1+ depth*Spec_value + (1 |timepoint), REML = FALSE, data = All_dat.tp)
# analyze(IRR.lmm.tp)
# summary(IRR.lmm.tp)
# 
# anova(IRR.lmm.random, IRR.lmm.tp)
# 
# IRR.lmm.tp.sat = lmer(IRR ~ 1+ depth*Spec_value + (1 | site) + (1 |timepoint) + 
#                         avgbov + 
#                         avgtemp +
#                         initial_size + 
#                         colony_size +
#                         maxtemp + 
#                         mintemp + 
#                         stdtemp + 
#                         maxbov + 
#                         minbov + 
#                         stdbov, 
#                       REML = FALSE, data = All_dat.tp)
# 
# step(IRR.lmm.tp.sat, direction = 'both')
# summary(IRR.lmm.tp.sat)
# anova(IRR.lmm.tp, IRR.lmm.tp.sat)
# analyze(IRR.lmm.tp.sat)
# rand(IRR.lmm.tp.sat)
# 
# IRR.lmm.tp.fitted = lmer(IRR ~ 1 + depth + (1 |timepoint) + 
#                            # avgbov + 
#                            # avgtemp +
#                            # initial_size,
#                            # maxtemp + 
#                            mintemp, 
#                          # stdtemp + 
#                          # maxbov + 
#                          # minbov + 
#                          # stdbov, 
#                          REML = FALSE, data = All_dat.tp)
# summary(IRR.lmm.tp.fitted)
# anova(IRR.lmm.tp,IRR.lmm.tp.fitted)
# analyze(IRR.lmm.tp.fitted)
# Anova(IRR.lmm.tp.fitted)
# 
# All_dat.tp.NA.IRR = data.frame(All_dat.tp[!is.na(All_dat.tp$IRR),])
# 
# ggplot(All_dat.tp.NA.IRR, aes(x = depth, y = IRR, color = sitetype))+
#   facet_grid(species~timepoint)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lmm.tp.fitted), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_discrete(name = "Depth Category", labels = c("Mesophotic","Shallow"))+
#   xlab(label = "Depth (m)") + ylab(label = "Pigmentation Lesion Recovery Rate (color value/day)")+
#   ggtitle(label = "Recovery Rate of Corals Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# ##AL-Specific Mixed Effects
# IRR.lmm.random.AL = lmer(IRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(IRR.lmm.random.AL)
# 
# #Level 1 predictors
# IRR.lmm.tp.AL = lmer(IRR ~ 1+ depth + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(IRR.lmm.tp.AL)
# 
# anova(IRR.lmm.random.AL, IRR.lmm.tp.AL)
# 
# IRR.lmm.tp.sat.AL = lmer(IRR ~ 1+ depth + (1 | site) + (1 |timepoint) + 
#                            avgbov + 
#                            avgtemp +
#                            initial_size + 
#                            colony_size +
#                            maxtemp + 
#                            mintemp + 
#                            stdtemp + 
#                            maxbov + 
#                            minbov + 
#                            stdbov, 
#                          REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# 
# step(IRR.lmm.tp.sat.AL, direction = 'both')
# summary(IRR.lmm.tp.sat.AL)
# anova(IRR.lmm.tp.sat.AL)
# analyze(IRR.lmm.tp.sat.AL)
# 
# IRR.lmm.tp.fitted.AL = lmer(IRR ~ 1 + depth + (1 |timepoint) + 
#                               # avgbov + 
#                               # avgtemp +
#                               # initial_size,
#                               # maxtemp + 
#                               mintemp, 
#                             # stdtemp + 
#                             # maxbov + 
#                             # minbov + 
#                             # stdbov, 
#                             REML = FALSE, data = All_dat.tp.spec$`Agaricia lamarcki`)
# summary(IRR.lmm.tp.fitted.AL)
# anova(IRR.lmm.tp.sat.AL,IRR.lmm.tp.fitted.AL)
# analyze(IRR.lmm.tp.fitted.AL)
# Anova(IRR.lmm.tp.fitted.AL)
# rand(IRR.lmm.tp.fitted.AL)
# 
# All_dat.tp.NA.IRR.AL = data.frame(All_dat.tp.spec$`Agaricia lamarcki`[!is.na( All_dat.tp.spec$`Agaricia lamarcki`$IRR),])
# 
# ggplot(All_dat.tp.NA.IRR.AL, aes(x = depth, y = IRR, color = sitetype))+
#   facet_grid(~timepoint)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lmm.tp.fitted.AL), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_discrete(name = "Depth Category", labels = c("Mesophotic","Shallow"))+
#   xlab(label = "Depth (m)") + ylab(label = "Pigmentation Lesion Recovery Rate (color value/day)")+
#   ggtitle(label = "Recovery Rate of Agaricia lamarcki Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# ##OFRA-Specific Mixed Effects
# IRR.lmm.random.OFRA = lmer(IRR ~ 1 +  (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(IRR.lmm.random.OFRA)
# 
# #Level 1 predictors
# IRR.lmm.tp.OFRA = lmer(IRR ~ 1+ depth + (1 | site) + (1 |timepoint), REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(IRR.lmm.tp.OFRA)
# rand(IRR.lmm.tp.OFRA)
# anova(IRR.lmm.random.OFRA, IRR.lmm.tp.OFRA)
# 
# IRR.lmm.tp.sat.OFRA = lmer(IRR ~ 1+ depth + (1 | site) + (1 |timepoint) +
#                              avgbov + 
#                              avgtemp +
#                              initial_size + 
#                              colony_size +
#                              maxtemp + 
#                              mintemp + 
#                              stdtemp + 
#                              maxbov + 
#                              minbov + 
#                              stdbov +
#                              initial_size,
#                            REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# 
# step(IRR.lmm.tp.sat.OFRA, direction = 'both')
# summary(IRR.lmm.tp.sat.OFRA)
# anova(IRR.lmm.tp.OFRA, IRR.lmm.tp.sat.OFRA)
# rand(IRR.lmm.tp.sat.OFRA)
# analyze(IRR.lmm.tp.sat.OFRA)
# 
# IRR.lmm.tp.fitted.OFRA = lmer(IRR ~ 1 + depth + (1 |timepoint) + 
#                                 # avgbov + 
#                                 avgtemp+
#                                 # initial_size,
#                                 maxtemp + 
#                                 mintemp + 
#                               # stdtemp + 
#                               maxbov + 
#                               # minbov, 
#                               stdbov, 
#                               REML = FALSE, data = All_dat.tp.spec$`Orbicella franksi`)
# summary(IRR.lmm.tp.fitted.OFRA)
# anova(IRR.lmm.tp.sat.OFRA,IRR.lmm.tp.fitted.OFRA)
# analyze(IRR.lmm.tp.fitted.OFRA)
# Anova(IRR.lmm.tp.fitted.OFRA)
# rand(IRR.lmm.tp.fitted.OFRA)
# 
# ggplot(All_dat.tp.NA.IRR.OFRA, aes(x = depth, y = IRR, color = sitetype))+
#   facet_grid(~timepoint)+
#   geom_point(stat = "identity")+
#   geom_line(aes(x = depth, y = predict(IRR.lmm.tp.fitted.OFRA), group = species), cex = 2)+
#   scale_shape_discrete(name = "Timepoint", labels = c("T1","T2","T3"))+
#   scale_color_discrete(name = "Depth Category", labels = c("Mesophotic","Shallow"))+
#   xlab(label = "Depth (m)") + ylab(label = "Pigmentation Lesion Recovery Rate (color value/day)")+
#   ggtitle(label = "Recovery Rate of Orbicella franksi Across Depth")+
#   theme(plot.title = element_text(hjust = 0.5))



#####################################################################################################################################################################################################################################################################################CHAPTER 2: COMPARISON OF LESION CHARACTERISTICS############################################

######NMDS#########
#Remove NAs
All_dat.NA = na.omit(All_dat)
All_dat.tot.NA = na.omit(All_dat.total)
All_dat.tp.NA = na.omit(All_dat.tp)

###Total Only
rec_matrix = as.matrix(cbind(All_dat.tot.NA$ALRR, 
                             All_dat.tot.NA$ASAP_rate, 
                             All_dat.tot.NA$RLRR, 
                             All_dat.tot.NA$SAP_recov, 
                             All_dat.tot.NA$IRR))
rec.mds = metaMDS(rec_matrix, k = 2, distance = 'euclidean')


#Method2
data<- as.data.frame(rec_matrix)
head(data)
ncol(data)
head(data)
d_matrix=data.matrix(data)
head(d_matrix)
d_veg.tot=vegdist(d_matrix, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE)

#Running each possible stress of NMDS- these don't all have to be run if stress selection is predetermined
mds1 <- metaMDS(d_veg.tot, k=1, autotransform=FALSE)
mds2 <- metaMDS(d_veg.tot, k=2, autotransform=FALSE)
mds3 <- metaMDS(d_veg.tot, k=3, autotransform=FALSE)
mds4 <- metaMDS(d_veg.tot, k=4, autotransform=FALSE)
mds5 <- metaMDS(d_veg.tot, k=5, autotransform=FALSE)
mds6 <- metaMDS(d_veg.tot, k=6, autotransform=FALSE)
mds7 <- metaMDS(d_veg.tot, k=7, autotransform=FALSE)
mds8 <- metaMDS(d_veg.tot, k=8, autotransform=FALSE)
mds9 <- metaMDS(d_veg.tot, k=9, autotransform=FALSE)
mds10 <- metaMDS(d_veg.tot, k=10, autotransform=FALSE)
mds11 <- metaMDS(d_veg.tot, k=11, autotransform=FALSE)
mds12 <- metaMDS(d_veg.tot, k=12, autotransform=FALSE)
mds13 <- metaMDS(d_veg.tot, k=13, autotransform=FALSE)
mds14 <- metaMDS(d_veg.tot, k=14, autotransform=FALSE)
mds15 <- metaMDS(d_veg.tot, k=15, autotransform=FALSE)
mds16 <- metaMDS(d_veg.tot, k=16, autotransform=FALSE)
mds17 <- metaMDS(d_veg.tot, k=17, autotransform=FALSE)
mds18 <- metaMDS(d_veg.tot, k=18, autotransform=FALSE)
mds19 <- metaMDS(d_veg.tot, k=19, autotransform=FALSE)
mds20 <- metaMDS(d_veg.tot, k=20, autotransform=FALSE)

#Graph Stress and Determine Desired Stress Level
stress=c(mds1$stress, mds2$stress, mds3$stress, mds4$stress, mds5$stress, mds6$stress, mds7$stress, mds8$stress, mds9$stress, mds10$stress, mds11$stress, mds12$stress, mds13$stress, mds14$stress, mds15$stress, mds16$stress, mds17$stress, mds18$stress, mds19$stress, mds20$stress)
dim=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
plot(dim, stress, xlab="Number of Dimensions (k)",ylab="Stress Value")

#Making the NMDS#
summary(mds2)
names(mds2)
#Isolate Points- use the desired stress
sampleScores <- mds2$points
sampleScores

#Input Bray Curtis points with transect information
dat=data.frame(NMDS1=sampleScores[,1], NMDS2=sampleScores[,2], site=All_dat.tot.NA$site, type = All_dat.tot.NA$sitetype, species=All_dat.tot.NA$species, timepoint=All_dat.tot.NA$timepoint)


#Plot the ordination
ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  geom_point(stat = 'identity', aes(colour = type, shape = species), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ))+
  scale_shape_manual("Species", values = c("Agaricia lamarcki" = 15, "Orbicella franksi" = 17))+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  facet_grid(~species)+
  theme(strip.text.x = element_text(size = 20))+
  geom_point(stat = 'identity', aes(colour = type), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ))+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

###Timpoints without Summary
rec_matrix = as.matrix(cbind(All_dat.tp.NA$ALRR, 
                             All_dat.tp.NA$ASAP_rate, 
                             All_dat.tp.NA$RLRR, 
                             All_dat.tp.NA$SAP_recov, 
                             All_dat.tp.NA$IRR))
rec.mds = metaMDS(rec_matrix, k = 2, distance = 'euclidean')

data<- as.data.frame(rec_matrix)
head(data)
ncol(data)
head(data)
d_matrix=data.matrix(data)
d_veg.tp=vegdist(d_matrix, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

#Running each possible stress of NMDS- these don't all have to be run if stress selection is predetermined
mds1 <- metaMDS(d_veg.tp, k=1, autotransform=FALSE)
mds2 <- metaMDS(d_veg.tp, k=2, autotransform=FALSE)
mds3 <- metaMDS(d_veg.tp, k=3, autotransform=FALSE)
mds4 <- metaMDS(d_veg.tp, k=4, autotransform=FALSE)
mds5 <- metaMDS(d_veg.tp, k=5, autotransform=FALSE)
mds6 <- metaMDS(d_veg.tp, k=6, autotransform=FALSE)
mds7 <- metaMDS(d_veg.tp, k=7, autotransform=FALSE)
mds8 <- metaMDS(d_veg.tp, k=8, autotransform=FALSE)
mds9 <- metaMDS(d_veg.tp, k=9, autotransform=FALSE)
mds10 <- metaMDS(d_veg.tp, k=10, autotransform=FALSE)
mds11 <- metaMDS(d_veg.tp, k=11, autotransform=FALSE)
mds12 <- metaMDS(d_veg.tp, k=12, autotransform=FALSE)
mds13 <- metaMDS(d_veg.tp, k=13, autotransform=FALSE)
mds14 <- metaMDS(d_veg.tp, k=14, autotransform=FALSE)
mds15 <- metaMDS(d_veg.tp, k=15, autotransform=FALSE)
mds16 <- metaMDS(d_veg.tp, k=16, autotransform=FALSE)
mds17 <- metaMDS(d_veg.tp, k=17, autotransform=FALSE)
mds18 <- metaMDS(d_veg.tp, k=18, autotransform=FALSE)
mds19 <- metaMDS(d_veg.tp, k=19, autotransform=FALSE)
mds20 <- metaMDS(d_veg.tp, k=20, autotransform=FALSE)

#Graph Stress and Determine Desired Stress Level
stress=c(mds1$stress, mds2$stress, mds3$stress, mds4$stress, mds5$stress, mds6$stress, mds7$stress, mds8$stress, mds9$stress, mds10$stress, mds11$stress, mds12$stress, mds13$stress, mds14$stress, mds15$stress, mds16$stress, mds17$stress, mds18$stress, mds19$stress, mds20$stress)
dim=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
plot(dim, stress, xlab="Number of Dimensions (k)",ylab="Stress Value")

#Making the NMDS#
summary(mds2)
names(mds2)
#Isolate Points- use the desired stress
sampleScores <- mds2$points
sampleScores

#Input Bray Curtis points with transect information
dat=data.frame(NMDS1=sampleScores[,1], NMDS2=sampleScores[,2], site=All_dat.tp.NA$site, type = All_dat.tp.NA$sitetype, species=All_dat.tp.NA$species, timepoint=All_dat.tp.NA$timepoint)


#Plot the ordination
ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  facet_wrap(~timepoint)+
  theme(strip.text.x = element_text(size = 20))+
  geom_point(stat = 'identity', aes(colour = type, shape = species), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ))+
  scale_shape_manual("Species", values = c("Agaricia lamarcki" = 15, "Orbicella franksi" = 17))+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  facet_grid(species~timepoint)+
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20))+
  geom_point(stat = 'identity', aes(colour = type), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1" ))+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size=20), legend.text=element_text(size=20))

###All Data

rec_matrix = as.matrix(cbind(All_dat.NA$ALRR, 
                             All_dat.NA$ASAP_rate, 
                             All_dat.NA$RLRR, 
                             All_dat.NA$SAP_recov, 
                             All_dat.NA$IRR))
rec.mds = metaMDS(rec_matrix, k = 2, distance = 'euclidean')

data<- as.data.frame(rec_matrix)
head(data)
ncol(data)
head(data)
d_matrix=data.matrix(data)
head(d_matrix)
d_veg=vegdist(d_matrix, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

#Running each possible stress of NMDS- these don't all have to be run if stress selection is predetermined
mds1 <- metaMDS(d_veg, k=1, autotransform=FALSE)
mds2 <- metaMDS(d_veg, k=2, autotransform=FALSE)
mds3 <- metaMDS(d_veg, k=3, autotransform=FALSE)
mds4 <- metaMDS(d_veg, k=4, autotransform=FALSE)
mds5 <- metaMDS(d_veg, k=5, autotransform=FALSE)
mds6 <- metaMDS(d_veg, k=6, autotransform=FALSE)
mds7 <- metaMDS(d_veg, k=7, autotransform=FALSE)
mds8 <- metaMDS(d_veg, k=8, autotransform=FALSE)
mds9 <- metaMDS(d_veg, k=9, autotransform=FALSE)
mds10 <- metaMDS(d_veg, k=10, autotransform=FALSE)
mds11 <- metaMDS(d_veg, k=11, autotransform=FALSE)
mds12 <- metaMDS(d_veg, k=12, autotransform=FALSE)
mds13 <- metaMDS(d_veg, k=13, autotransform=FALSE)
mds14 <- metaMDS(d_veg, k=14, autotransform=FALSE)
mds15 <- metaMDS(d_veg, k=15, autotransform=FALSE)
mds16 <- metaMDS(d_veg, k=16, autotransform=FALSE)
mds17 <- metaMDS(d_veg, k=17, autotransform=FALSE)
mds18 <- metaMDS(d_veg, k=18, autotransform=FALSE)
mds19 <- metaMDS(d_veg, k=19, autotransform=FALSE)
mds20 <- metaMDS(d_veg, k=20, autotransform=FALSE)

#Graph Stress and Determine Desired Stress Level
stress=c(mds1$stress, mds2$stress, mds3$stress, mds4$stress, mds5$stress, mds6$stress, mds7$stress, mds8$stress, mds9$stress, mds10$stress, mds11$stress, mds12$stress, mds13$stress, mds14$stress, mds15$stress, mds16$stress, mds17$stress, mds18$stress, mds19$stress, mds20$stress)
dim=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
plot(dim, stress, xlab="Number of Dimensions (k)",ylab="Stress Value")

#Making the NMDS#
summary(mds2)
names(mds2)
#Isolate Points- use the desired stress
sampleScores <- mds2$points
sampleScores

#Input Bray Curtis points with transect information
dat=data.frame(NMDS1=sampleScores[,1], NMDS2=sampleScores[,2], site=All_dat.NA$site, type = All_dat.NA$sitetype, species=All_dat.NA$species, timepoint=All_dat.NA$timepoint)

#Plot the ordination
ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  facet_wrap(~timepoint)+
  geom_point(stat = 'identity', aes(colour = type, shape = species), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "dodgerblue1", "Shallow Site" = "goldenrod1" ))+
  scale_shape_manual("Species", values = c("Agaricia lamarcki" = 15, "Orbicella franksi" = 17))

ggplot(dat, aes(x = NMDS1, y = NMDS2))+
  facet_grid(species~timepoint)+
  geom_point(stat = 'identity', aes(colour = type), cex = 4.0)+
  scale_colour_manual("Depth Category", values = c("Mesophotic Site" = "dodgerblue1", "Shallow Site" = "goldenrod1" ))

#########################PERMANOVA#######################
Rec_corr = cor(All_dat.tp.NA[,c('ALRR','RLRR','IRR')])
Rec_corr #Correlation of recovery variables
#####Print#######
#         ALRR      RLRR       IRR
# ALRR 1.0000000 
# RLRR 0.7377088 1.0000000 
# IRR  0.3445599 0.2570211 1.0000000
#ALRR and RLRR are highly correlated, IRR is not correlated with either#

rec_perm = adonis(d_veg.tp ~ species*sitetype*timepoint, data = All_dat.tp.NA) #looking at centroids of the data
rec_perm #species, no interaction

rec_dist.spec = betadisper(d_veg.tp, group = All_dat.tp.NA$species) #looking at variance of the data
anova(rec_dist.spec) #nothing

rec_dist.depth = betadisper(d_veg.tp, group =  All_dat.tp.NA$sitetype) #looking at variance of the data among site type
anova(rec_dist.depth) #not much going on
rec_dist.timepoint = betadisper(d_veg.tp, group =  All_dat.tp.NA$timepoint) #looking at variance of the data
anova(rec_dist.timepoint) #no changes in variance

##In order to assess interactions, I can look at subsamples of the data

#Agaricia Lamarcki
All_dat.tp.spec$`Agaricia lamarcki` = na.omit(All_dat.tp.spec$`Agaricia lamarcki`)
rec_matrix.AL = as.matrix(cbind(All_dat.tp.spec$`Agaricia lamarcki`$ALRR, 
                                All_dat.tp.spec$`Agaricia lamarcki`$ASAP_rate, 
                                All_dat.tp.spec$`Agaricia lamarcki`$RLRR, 
                                All_dat.tp.spec$`Agaricia lamarcki`$SAP_recov, 
                                All_dat.tp.spec$`Agaricia lamarcki`$IRR))
rec.mds.tot.AL = metaMDS(rec_matrix.AL, k = 2, distance = 'euclidean')

data.AL<- as.data.frame(rec_matrix.AL)
head(data.AL)
ncol(data.AL)
head(data.AL)
d_matrix.AL=data.matrix(data.AL)
head(d_matrix.AL)
d_veg.AL=vegdist(d_matrix.AL, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

rec_dist.depth.AL = betadisper(d_veg.AL, group =  All_dat.tp.spec$`Agaricia lamarcki`$sitetype) #looking at variance of the data among site type
anova(rec_dist.depth.AL) #Nope
rec_dist.timepoint.AL = betadisper(d_veg.AL, group =  All_dat.tp.spec$`Agaricia lamarcki`$timepoint) #looking at variance of the data
anova(rec_dist.timepoint.AL) #nope

All_dat.total.spec$`Agaricia lamarcki` = na.omit(All_dat.total.spec$`Agaricia lamarcki`) 
rec_matrix.tot.AL = as.matrix(cbind(All_dat.total.spec$`Agaricia lamarcki`$ALRR, 
                                    All_dat.total.spec$`Agaricia lamarcki`$ASAP_rate, 
                                    All_dat.total.spec$`Agaricia lamarcki`$RLRR, 
                                    All_dat.total.spec$`Agaricia lamarcki`$SAP_recov, 
                                    All_dat.total.spec$`Agaricia lamarcki`$IRR))
rec.mds.AL = metaMDS(rec_matrix.tot.AL, k = 2, distance = 'euclidean')

data.AL<- as.data.frame(rec_matrix.tot.AL)
head(data.AL)
ncol(data.AL)
head(data.AL)
d_matrix.AL=data.matrix(data.AL)
head(d_matrix.AL)
d_veg.AL=vegdist(d_matrix.AL, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

rec_dist.depth.AL = betadisper(d_veg.AL, group =   All_dat.total.spec$`Agaricia lamarcki`$sitetype) #looking at variance of the data among site type
anova(rec_dist.depth.AL) #Nope

#Orbicella Franksi
All_dat.tp.spec$`Orbicella franksi` = na.omit(All_dat.tp.spec$`Orbicella franksi`)
rec_matrix.OFRA = as.matrix(cbind(All_dat.tp.spec$`Orbicella franksi`$ALRR, 
                                  All_dat.tp.spec$`Orbicella franksi`$ASAP_rate, 
                                  All_dat.tp.spec$`Orbicella franksi`$RLRR, 
                                  All_dat.tp.spec$`Orbicella franksi`$SAP_recov, 
                                  All_dat.tp.spec$`Orbicella franksi`$IRR))
rec.mds.OFRA = metaMDS(rec_matrix, k = 2, distance = 'euclidean')

data.OFRA<- as.data.frame(rec_matrix.OFRA)
head(data.OFRA)
ncol(data.OFRA)
head(data.OFRA)
d_matrix.OFRA=data.matrix(data.OFRA)
head(d_matrix.OFRA)
d_veg.OFRA=vegdist(d_matrix.OFRA, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

rec_dist.depth.OFRA = betadisper(d_veg.OFRA, group =  All_dat.tp.spec$`Orbicella franksi`$sitetype) #looking at variance of the data among site type
anova(rec_dist.depth.OFRA) #Interaction!
rec_dist.timepoint.OFRA = betadisper(d_veg.OFRA, group =  All_dat.tp.spec$`Orbicella franksi`$timepoint) #looking at variance of the data
anova(rec_dist.timepoint.OFRA) #no differences through time

simp.1 = simper(rec_matrix, group = list(All_dat$species,All_dat$sitetype,All_dat$timepoint))
summary(simp.1)  

########Density Distributions##############

####Total Values
ggplot(All_dat, aes(x = RLRR, fill = sitetype))+
  geom_density()+
  facet_wrap(species ~ site)

Dens.RLRR = ggplot(All_dat.total, aes(x = RLRR, fill = sitetype))+
  geom_density(aes(x = RLRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"))+
  facet_grid(rows = vars(species))+
  theme(strip.text.y = element_text(size = 20))+
  xlab(label = "% of lesion/day")+
  ylab(label = " ")+
  ggtitle(label = "Relative Recovery")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20))+
  theme(text = element_text(size = 20),
        legend.justification=c(1,1),
        legend.position=c(.82, 1),
        legend.background = element_rect(size = 0.5, colour = 'black'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))+
  theme(text = element_text(size = 20))

Dens.ALRR = ggplot(All_dat.total, aes(x = ALRR, fill = sitetype))+
  geom_density(aes(x = ALRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual("Depth Category", values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"), guide = FALSE)+
  facet_grid(rows = vars(species))+
  theme(strip.text.y = element_text(size = 20))+
  xlab(label = "cm^2/day")+
  ylab(label = " ")+
  ggtitle(label = "Absolute Recovery")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20))

Dens.IRR = ggplot(All_dat.total, aes(x = IRR, fill = sitetype))+
  geom_density(aes(x = IRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual(values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"), guide = FALSE)+
  facet_grid(rows = vars(species))+
  theme(strip.text.y = element_text(size = 20))+
  xlab(label = "color value/day")+
  ylab(label = " ")+
  ggtitle(label = "Pigmentation Recovery")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20))
  

grid.arrange(Dens.ALRR,Dens.RLRR,Dens.IRR, ncol = 3)

SCPD = data.frame(All_dat[All_dat$site == "SCPD",])
SCPS= data.frame(All_dat[All_dat$site == "SCPS",])
SCP = rbind(SCPD, SCPS)
SCP$species = revalue(SCP$species, c("AL" = "A. Lamarcki","OFRA" = "O. Franksi"))
ggplot(SCP, aes(x = RLRR, fill = site))+
  geom_density(aes(x = RLRR, y = ..scaled.., fill = site), alpha = 3/4)+
  scale_fill_discrete(name = "Location", labels = c("South Capella 35m","South Capella 23m"))+
  facet_wrap(~species)+
  xlab(label = "Relative Lesion Recovery Rate (% recovered/day)")+
  ylab(label = " ")+
  ggtitle(label = "Density Curves of Relative Recovery Rate at South Capella")+
  theme(plot.title = element_text(hjust = 0.5))




#####Timepoint x Species
ggplot(All_dat.tp, aes(x = RLRR, fill = sitetype))+
  geom_density(aes(x = RLRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual("Depth Category",values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"))+
  facet_grid(species~timepoint)+
  xlab(label = "Relative Lesion Recovery Rate (% of lesion/day)")+
  ylab(label = " ")+
  ggtitle(label = "Density Curves of Relative Recovery Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20), strip.text.y = element_text(size = 20), strip.text.x = element_text(size =20))

ggplot(All_dat.tp, aes(x = ALRR, fill = sitetype))+
  geom_density(aes(x = ALRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual("Depth Category",values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"))+
  facet_grid(species~timepoint)+
  xlab(label = "Lesion Recovery Rate (cm^2/day)")+
  ylab(label = " ")+
  ggtitle(label = "Density Curves of Absolute Recovery Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20), strip.text.y = element_text(size = 20), strip.text.x = element_text(size =20))

ggplot(All_dat.tp, aes(x = IRR, fill = sitetype))+
  geom_density(aes(x = IRR, y = ..scaled.., fill = sitetype), alpha = 3/4)+
  scale_fill_manual("Depth Category",values = c("Mesophotic Site" = "red2", "Shallow Site" = "deepskyblue1"))+
  facet_grid(species~timepoint)+
  xlab(label = "Pigmentation Recovery Rate (color value/day)")+
  ylab(label = " ")+
  ggtitle(label = "Density Curves of Pigmentation Recovery Rate")+
  theme(plot.title = element_text(hjust = 0.5, size = 20), text = element_text(size = 20), strip.text.y = element_text(size = 20), strip.text.x = element_text(size =20))



