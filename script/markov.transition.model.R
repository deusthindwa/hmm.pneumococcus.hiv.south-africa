#written by Deus Thindwa
#14/05/2019

DDHP.packages <-c("foreign","tidyverse","janitor","readstata13","rstan","rethinking")
lapply(DDHP.packages, library, character.only=TRUE)

#load male questionnaire csv
phirst.cohort <-as_tibble(read.dta13(""))
