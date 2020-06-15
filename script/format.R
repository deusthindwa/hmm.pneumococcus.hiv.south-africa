#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#rename columns for household acquisition rates
phirst.es <- rename(phirst.es, c("iid"="Table.1A","age"="Age","hiv"="HIV","hh_hiv"="Household.adult.HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

#rename columns for household acquisition probabilities
phirst.es <- rename(phirst.es, c("iid"="Table.1B","age"="Age","hiv"="HIV","hh_hiv"="Household.adult.HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

#rename columns for community acquisition rates
phirst.es <- rename(phirst.es, c("iid"="Table.1C","age"="Age","hiv"="HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

#rename columns for community acquisition probabilities
phirst.es <- rename(phirst.es, c("iid"="Table.1D","age"="Age","hiv"="HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

#rename columns for average carriage duration by ART
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,artv,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2A","age"="Age","hiv"="HIV","artv"="ART","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#rename columns for average carriage duration by ABX
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,abxcat,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2B","age"="Age","hiv"="HIV","abxcat"="ABX","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#rename columns for average carriage duration
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2C","age"="Age","hiv"="HIV","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#rename columns for probability of carriage clearance
phirst.es <- rename(phirst.es, c("iid"="Table.2D","age"="Age","hiv"="HIV","clear.est"="Estimate","Lclear.est"="Lower.bound","Uclear.est"="Upper.bound"))

