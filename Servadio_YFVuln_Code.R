
rm(list = ls())

# setwd('/Volumes/easystore/MayoDesk/YellowFever/Vuln') # setwd('Z:/YellowFever/Vuln')

library(Hmisc)
library(tictoc)
library(MASS)
library(sensitivity)
library(car)

create.output = 'no'

# dat = read.csv('Data_for_analyses.csv', stringsAsFactors = FALSE)
dat.yf.year = read.csv('~/Desktop/Servadio_YFVuln_Data_key.csv', stringsAsFactors = FALSE)

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

dat = dat.yf.year[dat.yf.year$AnyYF == 1, ]
dat = dat[dat$Year < 2017, ]

# Variable selection: most closely fitting data

X = data.frame('Eco56' = dat$Eco56, 'Rainpct' = dat$Rainpct, 'Rainpct2' = dat$Rainpct^2, 
	'PopDen' = dat$PopDen, 'MEAN.DrDens' = dat$MEAN.DrDens, 
	'Elevation_m' = dat$Elevation_m, 'n.50.Cases' = 1 * (dat$n.50.Cases > 0), 
	'n.b.Cases' = 1 * (dat$n.b.Cases > 0), 'Vuln' = dat$Vuln)
X$Eco56 = as.factor(X$Eco56)
X$n.b.Cases = as.factor(X$n.b.Cases)
X$n.50.Cases = as.factor(X$n.50.Cases)

n.var = ncol(X)-1

covars = matrix(NA, nrow = (2^n.var)-1, ncol = n.var)
for (i in 1:nrow(covars)) {
	covars[i, ] = as.numeric(intToBits(i))[1:ncol(covars)]
}

covars = covars[-which(covars[ , 3] > 0 & covars[ , 2] == 0), ]
covars = covars[-which(covars[ , 7] > 0 & covars[ , 8] > 0), ]

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

	names(mods[[i]]) = c('Covariates', 'Model', 'Summary', 'Fit')
}

toc()

scores = NULL
for (i in 1:length(mods)) {
	scores = rbind(scores, mods[[i]]$Fit)
}

final.score = apply(scores, 1, sum)

ff = fitted(mods[[95]]$Model)
ff = cbind(dat$Mun_ID, ff)
ff = as.data.frame(ff)
names(ff) = c('MunCode', 'P0_1', 'P1_1', 'P2_1')
ff2 = aggregate(ff[ , -1], by = list('MunCode' = ff$MunCode), FUN = mean)
if (create.output == 'yes') { write.csv(ff2[!duplicated(ff2$MunCode), ], 
	'Vuln_fit_2000-16_ag.csv', row.names = FALSE) } 



# VIF

temp = X[ , c(which(covars[95, ] > 0), ncol(X))]
mod.v = lm(scale(Vuln) ~ ., data = temp[ , -3])
vif(mod.v)

#
# sensitivity to cut points
#

sens.scores = NULL
for (q in c(1:35)) {

	vuln.cuts = c(0, q)
	dat$Vuln = 0
	for (i in 1:length(vuln.cuts)) { dat$Vuln = dat$Vuln + (dat$Inc > vuln.cuts[i]) }

	X$Vuln = dat$Vuln

	new.mod = polr(as.factor(Vuln) ~ ., data = X[ , -8], method = 'logistic')

	pred = fitted(new.mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))

	sens.scores = rbind(sens.scores, c(sum(score[X$Vuln == 0]) / sum(X$Vuln == 0), 
		sum(score[X$Vuln == 1]) / sum(X$Vuln == 1), 
		sum(score[X$Vuln == 2]) / sum(X$Vuln == 2)))
}

if (create.output == 'yes') {

	pdf('Sensitivity_cutpoints_1.pdf', height = 6, width = 9)

	plot(c(1:35), apply(sens.scores, 1, sum), pch = 20, xlab = 'Cut point', 
		ylab = 'Fit Score')

	dev.off()
	
}

# Fitting Vulnerability

dat.new = dat.yf.year #[dat.yf.year$AnyYF == 0, ]
dat.new = dat.new[dat.new$Year < 2017, ]
X.new = data.frame('Eco56' = dat.new$Eco56, 
	'Rainpct' = dat.new$Rainpct, 
	'Rainpct2' = dat.new$Rainpct^2, 'PopDen' = dat.new$PopDen, 
	'MEAN.DrDens' = dat.new$MEAN.DrDens, 'Elevation_m' = dat.new$Elevation_m, 
	'n.50.Cases' = 1 * (dat.new$n.50.Cases > 0), 
	'n.b.Cases' = 1 * (dat.new$n.b.Cases > 0))
X.new$Eco56 = as.factor(X.new$Eco56)
X.new$n.b.Cases = as.factor(X.new$n.b.Cases)
X.new$n.50.Cases = as.factor(X.new$n.50.Cases)
X.new$Vuln = as.factor((dat.new$Incidence > 0) + (dat.new$Incidence > 10))

prediction_2000_2016 = predict(mods[[95]]$Model, X.new, type = 'prob', interval = 'confidence')

prediction_2000_2016 = as.data.frame(cbind(dat.new$Mun_ID, 
	dat.new$Year, prediction_2000_2016))
names(prediction_2000_2016) = c('MunCode', 'Year', 'P0', 'P1', 'P2')

prediction_2000_2016$P0 = as.numeric(as.character(prediction_2000_2016$P0))
prediction_2000_2016$P1 = as.numeric(as.character(prediction_2000_2016$P1))
prediction_2000_2016$P2 = as.numeric(as.character(prediction_2000_2016$P2))

if(create.output == 'yes') { save(prediction_2000_2016, file = 'Vuln_Pred_2000_2016.Rdat') } 


# GSUA

source('my_sobol.R')

tic()
a = mysobol2(ff, X[ , -c(8, 9)], 95)
toc()

a[[2]][c(1, 8, 14, 19, 23, 26, 28)]
a[[2]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[2]][c(1, 8, 14, 19, 23, 26, 28)])

a[[3]][c(1, 8, 14, 19, 23, 26, 28)]
a[[3]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[3]][c(1, 8, 14, 19, 23, 26, 28)])

a[[4]][c(1, 8, 14, 19, 23, 26, 28)]
a[[4]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[4]][c(1, 8, 14, 19, 23, 26, 28)])

a

sum(a$first.2[grep('7', a$mnames)] / sum(a$first.2))

# Breaking down into two logits

X2 = X
X2$PopDen = X2$PopDen / max(X2$PopDen)
X2$Elevation_m = X2$Elevation_m / max(X2$Elevation_m)

m1 = glm((1*(Vuln > 0)) ~ ., data = X2[ , -8], family = binomial(link = 'logit'))
m2 = glm((1*(Vuln > 1)) ~ ., data = X2[ , -8], family = binomial(link = 'logit'))

r1 = sample(c(1:nrow(X)), ceiling(nrow(X)/2), replace = FALSE)
r2 = c(1:nrow(X))[-r1]

s1 = sobol(model = m1, X1 = X2[r1, -c(8, 9)], X2 = X2[r2, -c(8, 9)], order = 4, nboot = 100)
s2 = sobol(model = m2, X1 = X2[r1, -c(8, 9)], X2 = X2[r2, -c(8, 9)], order = 4, nboot = 100)






###
# Vulnerability model: 2017
###

rm(list = c('a', 'b', 'covars', 'dat', 'dat.new', 'i', 'mod', 'mods', 
	'n.var', 'new.mod', 'pred', 'q', 'score', 'scores', 'sens.scores', 
	'temp', 'vuln.cuts', 'x', 'X', 'X.new', 'y'))



dat = dat.yf.year[dat.yf.year$AnyYF == 1, ]
dat = dat[dat$Year == 2017, ]

vuln.cuts = c(0, 10)
dat$Vuln = 0
for (i in 1:length(vuln.cuts)) { dat$Vuln = dat$Vuln + (dat$Inc > vuln.cuts[i]) }


# Variable selection: most closely fitting data

X = data.frame('Eco56' = dat$Eco56, 'Rainpct' = dat$Rainpct, 'Rainpct2' = dat$Rainpct^2, 
	'PopDen' = dat$PopDen, 'MEAN.DrDens' = dat$MEAN.DrDens, 
	'Elevation_m' = dat$Elevation_m, 'n.50.Cases' = 1 * (dat$n.50.Cases > 0), 
	'n.b.Cases' = 1 * (dat$n.b.Cases > 0), 'Vuln' = dat$Vuln)
X$Eco56 = as.factor(X$Eco56)
X$n.b.Cases = as.factor(X$n.b.Cases)
X$n.50.Cases = as.factor(X$n.50.Cases)

n.var = ncol(X)-1

covars = matrix(NA, nrow = (2^n.var)-1, ncol = n.var)
for (i in 1:nrow(covars)) {
	covars[i, ] = as.numeric(intToBits(i))[1:ncol(covars)]
}

covars = covars[-which(covars[ , 3] > 0 & covars[ , 2] == 0), ]
covars = covars[-which(covars[ , 7] > 0 & covars[ , 8] > 0), ]

mods = list()

tic()

for (i in 1:nrow(covars)) {
	mods[[i]] = list()
	mods[[i]][[1]] = names(X)[which(covars[i, ] > 0)]
	
	temp = X[ , c(which(covars[i, ] > 0), 9)]
	mod = polr(as.factor(Vuln) ~ ., data = temp, method = 'logistic')
	mods[[i]][[2]] = mod
	mods[[i]][[3]] = suppressMessages(summary(mod))

	pred = fitted(mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))

	mods[[i]][[4]] = c(sum(score[dat$Vuln == 0]) / sum(dat$Vuln == 0), 
		sum(score[dat$Vuln == 1]) / sum(dat$Vuln == 1), 
		sum(score[dat$Vuln == 2]) / sum(dat$Vuln == 2))

	names(mods[[i]]) = c('Covariates', 'Model', 'Summary', 'Fit')
}

toc()

scores = NULL
for (i in 1:length(mods)) {
	scores = rbind(scores, mods[[i]]$Fit)
}

a = apply(scores, 1, sum)

which(a == max(a))

ff = fitted(mods[[143]]$Model)
ff = cbind(dat$Mun_ID, ff)
ff = as.data.frame(ff)
names(ff) = c('MunCode', 'P0_2', 'P1_2', 'P2_2')

if (create.output == 'yes') { write.csv(ff, 'Vuln_fit_2017.csv', row.names = FALSE) }


#VIF
temp = X[ , c(which(covars[143, ] > 0), 9)]
mod.v = lm(scale(Vuln) ~ ., data = temp)
vif(mod.v)

mod.v = lm(scale(Vuln) ~ ., data = temp[ , -3])
vif(mod.v)


#
# sensitivity to cut points
#

sens.scores = NULL
for (q in c(1:35)) {

	vuln.cuts = c(0, q)
	dat$Vuln = 0
	for (i in 1:length(vuln.cuts)) { dat$Vuln = dat$Vuln + (dat$Inc > vuln.cuts[i]) }

	X$Vuln = dat$Vuln

	new.mod = polr(as.factor(Vuln) ~ ., data = X[ , -7], method = 'logistic')

	pred = fitted(new.mod)

	score = (pred[ , 1] * (dat$Vuln == 0)) + (pred[ , 2] * (dat$Vuln == 1)) + 
		(pred[ , 3] * (dat$Vuln == 2))

	sens.scores = rbind(sens.scores, c(sum(score[X$Vuln == 0]) / sum(X$Vuln == 0), 
		sum(score[X$Vuln == 1]) / sum(X$Vuln == 1), 
		sum(score[X$Vuln == 2]) / sum(X$Vuln == 2)))
}

if (create.output == 'yes') { pdf('Sensitivity_cutpoints_2.pdf', height = 6, width = 9) }

plot(c(1:35), apply(sens.scores, 1, sum), pch = 20, xlab = 'Cut point', 
	ylab = 'Fit Score')

if (create.output == 'yes') { dev.off() }
	

# Fitting Vulnerability

dat.new = dat.yf.year #[dat.yf.year$AnyYF == 0, ]
dat.new = dat.new[dat.new$Year == 2017, ]
X.new = data.frame('Eco56' = dat.new$Eco56, 
	'Rainpct' = dat.new$Rainpct, 
	'Rainpct2' = dat.new$Rainpct^2, 'PopDen' = dat.new$PopDen, 
	'MEAN.DrDens' = dat.new$MEAN.DrDens, 'Elevation_m' = dat.new$Elevation_m, 
	'n.50.Cases' = 1 * (dat.new$n.50.Cases > 0), 
	'n.b.Cases' = 1 * (dat.new$n.b.Cases > 0))
X.new$Eco56 = as.factor(X.new$Eco56)
X.new$n.b.Cases = as.factor(X.new$n.b.Cases)
X.new$n.50.Cases = as.factor(X.new$n.50.Cases)

prediction_2017 = predict(mods[[143]]$Model, X.new, type = 'prob')

prediction_2017 = as.data.frame(cbind(dat.new$Mun_ID, dat.new$GID_2, 
	dat.new$Year, prediction_2017))
names(prediction_2017) = c('MunCode', 'GID_2', 'Year', 'P0', 'P1', 'P2')

prediction_2017$P0 = as.numeric(as.character(prediction_2017$P0))
prediction_2017$P1 = as.numeric(as.character(prediction_2017$P1))
prediction_2017$P2 = as.numeric(as.character(prediction_2017$P2))


if (create.output == 'yes') { save(prediction_2017, file = 'Vuln_Pred_2017.Rdat') }

# GSUA

source('my_sobol.R')
source('Z:/ProgressBar.R')

tic()
a = mysobol2(ff, X[ , -c(7, 9)], 143)
toc()

a[[2]][c(1, 8, 14, 19, 23, 26, 28)]
a[[2]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[2]][c(1, 8, 14, 19, 23, 26, 28)])

a[[3]][c(1, 8, 14, 19, 23, 26, 28)]
a[[3]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[3]][c(1, 8, 14, 19, 23, 26, 28)])

a[[4]][c(1, 8, 14, 19, 23, 26, 28)]
a[[4]][c(1, 8, 14, 19, 23, 26, 28)] / sum(a[[4]][c(1, 8, 14, 19, 23, 26, 28)])

a

sum(a$first.2[grep('7', a$mnames)] / sum(a$first.2))


# Breaking down into two logits

X2 = X
X2$PopDen = X2$PopDen / max(X2$PopDen)
X2$Elevation_m = X2$Elevation_m / max(X2$Elevation_m)

m1 = glm((1*(Vuln > 0)) ~ ., data = X2[ , -7], family = binomial(link = 'logit'))
m2 = glm((1*(Vuln > 1)) ~ ., data = X2[ , -7], family = binomial(link = 'logit'))

r1 = sample(c(1:nrow(X)), ceiling(nrow(X)/2), replace = FALSE)
r2 = c(1:nrow(X))[-r1]

s1 = sobol(model = m1, X1 = X2[r1, -c(7, 9)], X2 = X2[r2, -c(7, 9)], order = 4, nboot = 100)
s2 = sobol(model = m2, X1 = X2[r1, -c(7, 9)], X2 = X2[r2, -c(7, 9)], order = 4, nboot = 100)

