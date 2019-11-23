#Written by Deus Thindwa
#Are HIV+ adult PNC carriers more likely to transmit pneumococci than HIV- adult PNC carriers?
#Continuous-time time-inhomogeneous hidden Markov modelling study, PhD chapter 1.
#13/08/2019

#---------------load required packages into memory
phirst.packages <-c("tidyverse","plyr","msm","timetk","gridExtra","curl","dplyr","minqa","lubridate", "magrittr","data.table","parallel")
lapply(phirst.packages, library, character.only=TRUE)

#---------------load a phirst dataset
phirst <- sp_msm
phirst$state <- if_else(phirst$state==1,2L,1L)
phirst$dys <- phirst$wks*3.5

#---------------show transition frequency in state table
statetable.msm(state,iid,data=phirst)

#---------------initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.1), 
                  c(0.1,0.0))
rownames(matrix.Q) <- c("Clear","Carrying")
colnames(matrix.Q) <- c("Clear","Carrying")

#---------------fit a simpler model without misclassification
p.model1<-msm(state~dys, subject=iid, data=phirst, 
              qmatrix=matrix.Q, 
              #censor=c(3,4), censor.states=list(c(1,2),c(1,2)), 
              #covariates=~age+hiv,
              pci=7,
              control=list(trace=1,REPORT=1), cl=0.95)
printnew.msm(p.model1)

#---------------assessment of a model without misclassification
#(a) observed vs expected
dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model1, mintime=0,legend.pos=c(0,100),lwd.obs=3,lwd.exp=3,xlab="Day",ylab="Prevalence%")
#legend(0, 100, legend=c("Simulated data", "fitted to HHM"),col=c("blue","white"), lty=1:2, cex=1, lwd=3)

#(b) pearson-type good of fit
options(digits=2)
pearson.msm(p.model1, timegroups=2)

#(c) likelihood surface
surface.msm(p.model1,params=c(1,3),np=2,type="filled.contour")

#---------------fit Markov model1 with covariates seperately
qmatrix.msm(p.model1, covariates=list(hiv=0))
qmatrix.msm(p.model1, covariates=list(hiv=1))
qmatrix.msm(p.model1, covariates="mean")

qmatrix.msm(p.model1, covariates=list(age=0))
qmatrix.msm(p.model1, covariates=list(age=1))
qmatrix.msm(p.model1, covariates="mean")

#---------------fit Markov model1 with covariates in combination
qmatrix.msm(p.model1, covariates=list(hiv=1, age=0))
qmatrix.msm(p.model1, covariates=list(hiv=1, age=1))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=1))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=0))

#---------------extract transition probability from fitted Markov model
pmatrix.msm(p.model1, t=15, ci="boot", cl=0.95, B=2)

#---------------initiate emission matrix E
matrix.E <- rbind(c(1.00,0.00), 
                  c(0.05,0.95))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carrying")

#---------------fit a complex model with misclassification
p.model2 <- msm(state~dys, subject=iid, data=phirst,
                qmatrix=matrix.Q, 
                ematrix=matrix.E,
                #covariates= ~age+hiv,
                #censor=c(3,4), censor.states=list(c(1,2),c(1,2)),
                est.initprobs=T,
                control=list(trace=1,REPORT=1), cl=0.95)
                #opt.method="bobyqa")
printnew.msm(p.model2)

#---------------assessment of a model without misclassification
#(a) observed vs expected
dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model2, mintime=0,legend.pos=c(0,100),lwd.obs=3,lwd.exp=3, col.obs="blue",xlab="Times (days)",ylab="Prevalence%")
#legend(0, 100, legend=c("Simulated data", "Model fitting"),col=c("blue","red"), lty=1,2)

#(b) pearson-type good of fit
options(digits=2)
pearson.msm(p.model2, timegroups=2)

#(c) likelihood surface
surface.msm(p.model2,params=c(1,2),np=50,type="filled.contour")

#---------------fit Markov model2 with covariates seperately
qmatrix.msm(p.model2, covariates=list(hiv=0))
qmatrix.msm(p.model2, covariates=list(hiv=1))
qmatrix.msm(p.model2, covariates="mean")

qmatrix.msm(p.model2, covariates=list(age=0))
qmatrix.msm(p.model2, covariates=list(age=1))
qmatrix.msm(p.model2, covariates="mean")

#---------------fit Markov model2 with covariates in combination
qmatrix.msm(p.model2, covariates=list(hiv=0, age=0))
qmatrix.msm(p.model2, covariates=list(hiv=1, age=0))
qmatrix.msm(p.model2, covariates=list(hiv=0, age=1))
qmatrix.msm(p.model2, covariates=list(hiv=1, age=1))

#---------------other important parameters
#transition probability matrices extracted from the msm() fit function
pmatrix.msm(p.model2, t=10, ci="boot", cl=0.95, B=10)

#average period in a single stay in a state (sojourn)
sojourn.msm(p.model2)

#Probability that each state is next
pnext.msm(p.model2)

#forecasted total length of time spent in each trasient state
totlos.msm(p.model2, tot=84)

#expected time until Markov process first enters a carrying state (hitting time)
efpt.msm(p.model2, tostate=2)

#expected number of visit to a state
envisits.msm(p.model2, ci="boot", cl=0.95, B=1)

#---------------plot the fit Markov model1 with covariates in combination
hmm_plots <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/markov.chain.model.pneumococcus.hiv.rsa/master/data/hmm_plots.csv"))
dev.off()
ggplot(hmm_plots, aes(x=Age, y=Intensity, group=HIV, color=HIV)) + 
  scale_color_manual(values=c('#999999','#E69F00')) + 
  geom_point(size=3,position=position_dodge(width=0.3),stat="identity") +
  geom_errorbar(aes(ymin=Lintensity, ymax=Uintensity), width=0.2,size=1,position=position_dodge(width=0.3),stat="identity") +
  facet_grid(.~State, scales="free_y") +
  ylim(c(0,50)) +
  theme_bw() +
  ylab("Weekly transition rate (%)") +
  xlab("Age group") +
  theme(strip.text.x = element_text(size = 11, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  guides(color=guide_legend(title="")) +
  theme(legend.position = c(0.8,0.85), legend.text = element_text(size = 11), legend.title = element_text(face="bold", size=11)) 

