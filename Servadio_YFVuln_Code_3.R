
for (i in 1:20) { cat('\n') }
for (i in 1:20) { cat('#####\n\n') }


date()

setwd('~/Desktop/YF/Vuln')
source('Anonymize_Muns.R')

rm(list = ls())

# setwd('/Volumes/easystore/MayoDesk/YellowFever/Vuln') # setwd('Z:/YellowFever/Vuln')


library(Hmisc)
library(tictoc)
library(MASS)
library(sensitivity)
library(car)

create.output = 'no'

# dat = read.csv('Data_for_analyses.csv', stringsAsFactors = FALSE)
dat.yf.year = read.csv('Servadio_YFVuln_Data_2.csv', stringsAsFactors = FALSE)

dat.yf.year$EcoReg[which(dat.yf.year$EcoReg %in% c(2, 9, 13, 14))] = 10

#
#
# Full cleaned data are available
#
#

###
# Vulnerability model: 2000-2016
###

vuln.cuts = c(0, 10)
dat.yf.year$Vuln = 0
for (i in 1:length(vuln.cuts)) { dat.yf.year$Vuln = dat.yf.year$Vuln + (dat.yf.year$Inc > vuln.cuts[i]) }

# Select Dengue control municipalities
den.muns = sample(unique(dat.yf.year$Mun_ID[dat.yf.year$AnyDen == 1 & dat.yf.year$AnyYF == 0]), 
	length(unique(dat.yf.year$Mun_ID[dat.yf.year$AnyYF == 1])), replace = FALSE)
dat.yf.year$AnyDen[which(dat.yf.year$Mun_ID %in% den.muns)] = 1.5

dat = dat.yf.year[dat.yf.year$AnyYF == 1 | dat.yf.year$AnyDen == 1.5, ]
dat = dat[dat$Year < 2017, ]

# Variable selection: most closely fitting data

X = data.frame('Rainpct' = dat$Rainpct, 'Rainpct2' = dat$Rainpct^2, 
	'Temp' = dat$MeanTmp, 'Temp2' = dat$MeanTmp^2, 
	'PopDen' = dat$PopDen, 'Elevation' = dat$Elevation, 
	'EcoReg' = dat$EcoReg, 'NDVI' = dat$NDVI, 'Drainage' = dat$Drainage, 
	'n.b.Cases' = 1 * (dat$n.b.Cases > 0), 
	'n.50.Cases' = 1 * (dat$n.50.Cases > 0), 'Vuln' = dat$Vuln)
X$EcoReg = as.factor(X$EcoReg)
X$n.b.Cases = as.factor(X$n.b.Cases)
X$n.50.Cases = as.factor(X$n.50.Cases)

n.var = ncol(X)-1

covars = matrix(NA, nrow = (2^n.var)-1, ncol = n.var)
for (i in 1:nrow(covars)) {
	covars[i, ] = as.numeric(intToBits(i))[1:ncol(covars)]
}

covars = covars[-which(covars[ , 2] > 0 & covars[ , 1] == 0), ]
covars = covars[-which(covars[ , 4] > 0 & covars[ , 3] == 0), ]
covars = covars[-which(covars[ , 10] > 0 & covars[ , 11] > 0), ]

cl.auc = function(pred.vals, y) {
	cutoffs = seq(0.001, 0.999, by = 0.001)
	cutoffs = seq(min(pred), max(pred[ , 2]) + max(pred[ , 3]) + 0.01, len = 5000)
	
	coords.01 = coords.12 = matrix(NA, nrow = length(cutoffs), ncol = 2)
	
	for (cc in 1:length(cutoffs)) {
		
		is.12 = 1 * ((pred.vals[ , 2] + pred.vals[ , 3]) >= cutoffs[cc])
		is.2 = 1 * (pred.vals[ , 3] >= cutoffs[cc])
		
		is.0 = 1 * (pred.vals[ , 1])
		
		tpr.12 = sum(is.12 == 1 & X$Vuln %in% c(1, 2)) / sum(X$Vuln %in% c(1, 2))
		tpr.2 = sum(is.2 == 1 & X$Vuln == 2) / sum(X$Vuln == 2)
		
		fpr.12 = sum(is.12 == 1 & X$Vuln == 0) / sum(X$Vuln == 0)
		fpr.2 = sum(is.2 == 1 & X$Vuln < 2) / sum(X$Vuln < 2)
		
		coords.01[cc, ] = c(fpr.2, tpr.2)
		coords.12[cc, ] = c(fpr.12, tpr.12)
		
	}
	
	auc.01 = sum(
		(coords.01[-nrow(coords.01), 1] - coords.01[-1, 1]) *
		(0.5 * (coords.01[-nrow(coords.01), 2] + coords.01[-1, 2]))
	)
	auc.12 = sum(
		(coords.12[-nrow(coords.01), 1] - coords.12[-1, 1]) *
		(0.5 * (coords.12[-nrow(coords.12), 2] + coords.12[-1, 2]))
	)
	
	# plot(coords.01, xlim = c(0, 1), ylim = c(0, 1), type ='l'); abline(0, 1, col = 'gray')
	# plot(coords.12, xlim = c(0, 1), ylim = c(0, 1), pch = 20); abline(0, 1, col = 'gray')
	
	return(c('auc.01' = auc.01, 'auc.12' = auc.12))
	
}

mods = list()

tic()

for (i in 1:nrow(covars)) {
	mods[[i]] = list()
	mods[[i]][[1]] = names(X)[which(covars[i, ] > 0)]
	
	temp = X[ , c(which(covars[i, ] > 0), ncol(X))]
	mod = polr(as.factor(Vuln) ~ ., data = temp, method = 'logistic')
	mods[[i]][[2]] = mod
	mods[[i]][[3]] = suppressMessages(summary(mod))

	pred = fitted(mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))

	mods[[i]][[4]] = c(sum(score[dat$Vuln == 0]) / sum(dat$Vuln == 0), 
		sum(score[dat$Vuln == 1]) / sum(dat$Vuln == 1), 
		sum(score[dat$Vuln == 2]) / sum(dat$Vuln == 2))
	
	mods[[i]][[5]] = (2 * sum(covars[i, ])) - (2 * sum(log(score)))
	
	mods[[i]][[6]] = cl.auc(pred, X$Vuln)

	names(mods[[i]]) = c('Covariates', 'Model', 'Summary', 'Fit', 'AIC', 'AUC')
}

toc()

if (create.output == 'yes') { save(mods, file = '/mods_2000-16.Rdata') }

scores = NULL
for (i in 1:length(mods)) {
	scores = rbind(scores, mods[[i]]$Fit)
}

aics = NULL
for (i in 1:length(mods)) {
	aics = c(aics, mods[[i]]$AIC)
}

final.score = apply(scores, 1, sum)

aucs = NULL
for (i in 1:length(mods)) {
	aucs = rbind(aucs, mods[[i]]$AUC)
}
aucs = cbind(aucs, apply(aucs, 1, sum))

best = which(aics == min(aics))
best = which(final.score == max(final.score))
best = which(aucs[ , 3] == max(aucs[ , 3]))
mods[[best]]

ff = fitted(mods[[best]]$Model)
ff = cbind(dat$Mun_ID, ff)
ff = as.data.frame(ff)
names(ff) = c('MunCode', 'P0_1', 'P1_1', 'P2_1')
for (i in 2:4) { ff[ , i] = as.numeric(as.character(ff[ , i]))}
ff2 = aggregate(ff[ , -1], by = list('MunCode' = ff$MunCode), FUN = mean)
if (create.output == 'yes') { write.csv(ff2[!duplicated(ff2$MunCode), ], 
	'Vuln_fit_2000-16_ag.csv', row.names = FALSE) } 



#
# sensitivity to cut points
#

sens.scores = NULL; sens.aics = NULL; sens.aucs = NULL
for (q in c(1:35)) {

	vuln.cuts = c(0, q)
	dat$Vuln = 0
	for (i in 1:length(vuln.cuts)) { dat$Vuln = dat$Vuln + (dat$Inc > vuln.cuts[i]) }

	X$Vuln = dat$Vuln

	new.mod = polr(as.factor(Vuln) ~ ., data = X[ , c(which(covars[best, ] > 0), ncol(X))], method = 'logistic')

	pred = fitted(new.mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))
	new.aic = (2 * length(new.mod$coefficients)) - (2 * sum(log(score)))
	
	new.auc = cl.auc(pred, dat$Vuln)

	sens.scores = rbind(sens.scores, c(sum(score[X$Vuln == 0]) / sum(X$Vuln == 0), 
		sum(score[X$Vuln == 1]) / sum(X$Vuln == 1), 
		sum(score[X$Vuln == 2]) / sum(X$Vuln == 2)))
		
	sens.aics = c(sens.aics, new.aic)
	
	sens.aucs = rbind(sens.aucs, new.auc)
}

if (create.output == 'yes') {

	pdf('New_out/Sensitivity_cutpoints_1.pdf', height = 6, width = 9)

	plot(c(1:35), apply(sens.scores, 1, sum), pch = 20, xlab = 'Cut point', 
		ylab = 'Fit Score')
		
	plot(c(1:35), sens.aics, pch = 20, xlab = 'Cut point', ylab = 'AIC')
	
	plot(c(1:35), apply(sens.aucs, 1, sum), pch = 20, xlab = 'Cut point', ylab = 'AUC')

	dev.off()
	
}

# Fitting Vulnerability

dat.new = dat.yf.year #[dat.yf.year$AnyYF == 0, ]
dat.new = dat.new[dat.new$Year < 2017, ]

X.new = data.frame('Rainpct' = dat.new$Rainpct, 'Rainpct2' = dat.new$Rainpct^2, 
	'Temp' = dat.new$MeanTmp, 'Temp2' = dat.new$MeanTmp^2, 
	'PopDen' = dat.new$PopDen, 'Elevation' = dat.new$Elevation, 'Drainage' = dat.new$Drainage, 
	'EcoReg' = dat.new$EcoReg, 'NDVI' = dat.new$NDVI, 
	'n.50.Cases' = 1 * (dat.new$n.50.Cases> 0), 'n.b.Cases' = 1 * (dat.new$n.b.Cases > 0))
X.new$EcoReg = as.factor(X.new$EcoReg)
X.new$n.b.Cases = as.factor(X.new$n.b.Cases)
X.new$n.50.Cases = as.factor(X.new$n.50.Cases)
X.new$Vuln = as.factor((dat.new$Incidence > 0) + (dat.new$Incidence > 10))

prediction_2000_2016 = predict(mods[[best]]$Model, X.new, type = 'prob', interval = 'confidence')

prediction_2000_2016 = as.data.frame(cbind(dat.new$Mun_ID, 
	dat.new$Year, prediction_2000_2016))
names(prediction_2000_2016) = c('MunCode', 'Year', 'P0', 'P1', 'P2')

prediction_2000_2016$P0 = as.numeric(as.character(prediction_2000_2016$P0))
prediction_2000_2016$P1 = as.numeric(as.character(prediction_2000_2016$P1))
prediction_2000_2016$P2 = as.numeric(as.character(prediction_2000_2016$P2))

if(create.output == 'yes') { save(prediction_2000_2016, file = 'New_out/Vuln_Pred_2000_2016.Rdat') } 




# GSUA

source('my_sobol.R')

a = mysobol2(ff, X[ , which(covars[best, ] > 0)], best)

same = which(substr(a$mnames, 1, 1) == substr(a$mnames, 3, 3))

a[[2]][same]
a[[2]][same] / sum(a[[2]][same])

a[[3]][same]
a[[3]][same] / sum(a[[3]][same])

a[[4]][same]
a[[4]][same] / sum(a[[4]][same])

a

# First order for each category
for (i in 1:length(same)) {
	cat(sum(a$first.0[grep(i, a$mnames)] / sum(a$first.0)), '\n')
}
for (i in 1:length(same)) {
	cat(sum(a$first.1[grep(i, a$mnames)] / sum(a$first.1)), '\n')
}
for (i in 1:length(same)) {
	cat(sum(a$first.2[grep(i, a$mnames)] / sum(a$first.2)), '\n')
}










###
# Vulnerability model: 2017
###

rm(list = c('a', 'aics', 'covars', 'dat', 'dat.new', 
	'ff', 'ff2', 'final.score', 'i', 'mod', 'mod.v', 'mods', 
	'n.var', 'new.aic', 'new.mod', 'pred', 'q', 'same', 'score', 
	'scores', 'sens.scores', 'sens.aics', 
	'temp', 'vuln.cuts', 'X', 'X.new'))

dat = dat.yf.year[dat.yf.year$AnyYF == 1 | dat.yf.year$AnyDen == 1.5, ]
dat = dat[dat$Year == 2017, ]

# Variable selection: most closely fitting data

X = data.frame('Rainpct' = dat$Rainpct, 'Rainpct2' = dat$Rainpct^2, 
	'Temp' = dat$MeanTmp, 'Temp2' = dat$MeanTmp^2, 
	'PopDen' = dat$PopDen, 'Elevation' = dat$Elevation, 
	'EcoReg' = dat$EcoReg, 'NDVI' = dat$NDVI, 'Drainage' = dat$Drainage, 
	'n.b.Cases' = 1 * (dat$n.b.Cases > 0), 
	'n.50.Cases' = 1 * (dat$n.50.Cases > 0), 'Vuln' = dat$Vuln)
X$EcoReg = as.factor(X$EcoReg)
X$n.b.Cases = as.factor(X$n.b.Cases)
X$n.50.Cases = as.factor(X$n.50.Cases)

n.var = ncol(X)-1

covars = matrix(NA, nrow = (2^n.var)-1, ncol = n.var)
for (i in 1:nrow(covars)) {
	covars[i, ] = as.numeric(intToBits(i))[1:ncol(covars)]
}

covars = covars[-which(covars[ , 2] > 0 & covars[ , 1] == 0), ]
covars = covars[-which(covars[ , 4] > 0 & covars[ , 3] == 0), ]
covars = covars[-which(covars[ , 10] > 0 & covars[ , 11] > 0), ]

mods = list()


for (i in 1:nrow(covars)) {
	mods[[i]] = list()
	mods[[i]][[1]] = names(X)[which(covars[i, ] > 0)]
	
	temp = X[ , c(which(covars[i, ] > 0), ncol(X))]
	mod = polr(as.factor(Vuln) ~ ., data = temp, method = 'logistic')
	mods[[i]][[2]] = mod
	mods[[i]][[3]] = suppressMessages(summary(mod))

	pred = fitted(mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))

	mods[[i]][[4]] = c(sum(score[dat$Vuln == 0]) / sum(dat$Vuln == 0), 
		sum(score[dat$Vuln == 1]) / sum(dat$Vuln == 1), 
		sum(score[dat$Vuln == 2]) / sum(dat$Vuln == 2))
	
	mods[[i]][[5]] = (2 * sum(covars[i, ])) - (2 * sum(log(score)))
	
	mods[[i]][[6]] = cl.auc(pred, X$Vuln)

	names(mods[[i]]) = c('Covariates', 'Model', 'Summary', 'Fit', 'AIC', 'AUC')
}


if (create.output == 'yes') { save(mods, file = 'mods_2017.Rdata') }

scores = NULL
for (i in 1:length(mods)) {
	scores = rbind(scores, mods[[i]]$Fit)
}

aics = NULL
for (i in 1:length(mods)) {
	aics = c(aics, mods[[i]]$AIC)
}

aucs = NULL
for (i in 1:length(mods)) {
	aucs = rbind(aucs, mods[[i]]$AUC)
}
aucs = cbind(aucs, apply(aucs, 1, sum))

final.score = apply(scores, 1, sum)

best = which(final.score == max(final.score))
best = which(aics == min(aics))
best = which(aucs[ , 3] == max(aucs[ , 3]))

mods[[best]]

ff = fitted(mods[[best]]$Model)
ff = cbind(dat$Mun_ID, ff)
ff = as.data.frame(ff)
names(ff) = c('MunCode', 'P0_1', 'P1_1', 'P2_1')
ff2 = aggregate(ff[ , -1], by = list('MunCode' = ff$MunCode), FUN = mean)
if (create.output == 'yes') { write.csv(ff2[!duplicated(ff2$MunCode), ], 
	'Vuln_fit_2017_ag.csv', row.names = FALSE) } 



# VIF

temp = X[ , c(which(covars[which(aics == min(aics)), ] > 0), ncol(X))]
mod.v = lm(scale(Vuln) ~ ., data = temp)
vif(mod.v)

#
# sensitivity to cut points
#

sens.scores = NULL; sens.aics = NULL; sens.aucs = NULL
for (q in c(1:35)) {

	vuln.cuts = c(0, q)
	dat$Vuln = 0
	for (i in 1:length(vuln.cuts)) { dat$Vuln = dat$Vuln + (dat$Inc > vuln.cuts[i]) }

	X$Vuln = dat$Vuln

	new.mod = polr(as.factor(Vuln) ~ ., data = X[ , c(which(covars[best, ] > 0), ncol(X))], method = 'logistic')

	pred = fitted(new.mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))
	new.aic = (2 * length(new.mod$coefficients)) - (2 * sum(log(score)))

	sens.scores = rbind(sens.scores, c(sum(score[X$Vuln == 0]) / sum(X$Vuln == 0), 
		sum(score[X$Vuln == 1]) / sum(X$Vuln == 1), 
		sum(score[X$Vuln == 2]) / sum(X$Vuln == 2)))
		
	new.auc = cl.auc(pred, X$Vuln)
		
	sens.aics = c(sens.aics, new.aic)
	sens.aucs = c(sens.aucs, sum(new.auc))
}

if (create.output == 'yes') {

	pdf('New_out/Sensitivity_cutpoints_2.pdf', height = 6, width = 9)

	plot(c(1:35), apply(sens.scores, 1, sum), pch = 20, xlab = 'Cut point', 
		ylab = 'Fit Score')
		
	plot(c(1:35), sens.aics, pch = 20, xlab = 'Cut point', ylab = 'AIC')
	
	plot(c(1:35), sens.aucs, pch = 20, xlab = 'Cut point', ylab = 'AUC')

	dev.off()
	
}

# Fitting Vulnerability

dat.new = dat.yf.year #[dat.yf.year$AnyYF == 0, ]
dat.new = dat.new[dat.new$Year == 2017, ]

X.new = data.frame('Rainpct' = dat.new$Rainpct, 'Rainpct2' = dat.new$Rainpct^2, 
	'Temp' = dat.new$MeanTmp, 'Temp2' = dat.new$MeanTmp^2, 
	'PopDen' = dat.new$PopDen, 'Elevation' = dat.new$Elevation, 'Drainage' = dat.new$Drainage, 
	'EcoReg' = dat.new$EcoReg, 'NDVI' = dat.new$NDVI, 
	'n.50.Cases' = 1 * (dat.new$n.50.Cases > 0), 'n.b.Cases' = 1 * (dat.new$n.b.Cases > 0))
X.new$EcoReg = as.factor(X.new$EcoReg)
X.new$n.b.Cases = as.factor(X.new$n.b.Cases)
X.new$n.50.Cases = as.factor(X.new$n.50.Cases)
X.new$Vuln = as.factor((dat.new$Incidence > 0) + (dat.new$Incidence > 10))

prediction_2017 = predict(mods[[best]]$Model, X.new, type = 'prob', interval = 'confidence')

prediction_2017 = as.data.frame(cbind(dat.new$Mun_ID, 
	dat.new$Year, prediction_2017))
names(prediction_2017) = c('MunCode', 'Year', 'P0', 'P1', 'P2')

prediction_2017$P0 = as.numeric(as.character(prediction_2017$P0))
prediction_2017$P1 = as.numeric(as.character(prediction_2017$P1))
prediction_2017$P2 = as.numeric(as.character(prediction_2017$P2))

if(create.output == 'yes') { save(prediction_2017, file = 'New_out/Vuln_Pred_2017.Rdat') } 




# GSUA

source('my_sobol.R')

a = mysobol2(ff, X[ , which(covars[best, ] > 0)], best)

same = which(substr(a$mnames, 1, 1) == substr(a$mnames, 3, 3))

a[[2]][same]
a[[2]][same] / sum(a[[2]][same])

a[[3]][same]
a[[3]][same] / sum(a[[3]][same])

a[[4]][same]
a[[4]][same] / sum(a[[4]][same])

a

for (i in 1:length(same)) {
	cat(sum(a$first.0[grep(i, a$mnames)] / sum(a$first.0)), '\n')
}

for (i in 1:length(same)) {
	cat(sum(a$first.1[grep(i, a$mnames)] / sum(a$first.1)), '\n')
}

for (i in 1:length(same)) {
	cat(sum(a$first.2[grep(i, a$mnames)] / sum(a$first.2)), '\n')
}



