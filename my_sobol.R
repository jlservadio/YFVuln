logit = function(x) {
	return(log(x / (1 - x)))
}

library(boot)

mysobol = function(outcome, inputs) {

	for (q in 1:ncol(inputs)) { inputs[ , q] = as.numeric(as.character(inputs[ , q])) }
	inputs = as.matrix(inputs)

	first.orders.0 = first.orders.1 = first.orders.2 = NULL
	for (i in 1:ncol(inputs)) {
		progress.bar(i, 1, ncol(inputs), ncol(inputs))
		a = min(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)
		b = max(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)

		means.0 = means.1 = means.2 = NULL

		for (j in 1:100) {

			new = inputs
			new[ , i] = runif(1, a, b)

			# new.fit = predict(mods[[95]]$Model, new, type = 'prob')
			new.fit.0 = new %*% as.matrix(coef(mods[[95]]$Model))

			new.fit = cbind(rep(0, nrow(new)), rep(0, nrow(new)), 
				exp(new.fit.0 - mods[[95]]$Model$zeta[2]) / 
					(1 + exp(new.fit.0 - mods[[95]]$Model$zeta[2]))
			)
			new.fit[ , 2] = (exp(new.fit.0 - mods[[95]]$Model$zeta[2]) / 
				(1 + exp(new.fit.0 - mods[[95]]$Model$zeta[2]))) - 
				(exp(new.fit.0 - mods[[95]]$Model$zeta[1]) / 
				(1 + exp(new.fit.0 - mods[[95]]$Model$zeta[1])))
			new.fit[ , 1] = 1 - new.fit[ , 2] - new.fit[ , 3]

			means.0 = c(means.0, mean(new.fit[ , 1]))
			means.1 = c(means.1, mean(new.fit[ , 2]))
			means.2 = c(means.2, mean(new.fit[ , 3]))
		}

		v.m.0 = var(means.0); v.m.1 = var(means.1); v.m.2 = var(means.2)

		first.orders.0 = c(first.orders.0, v.m.0 / var(outcome[ , 2]))
		first.orders.1 = c(first.orders.1, v.m.1 / var(outcome[ , 3]))
		first.orders.2 = c(first.orders.2, v.m.2 / var(outcome[ , 4]))

	}

	out = list('first.0' = first.orders.0, 'first.1' = first.orders.1, 
		'first.2' = first.orders.2)

	return(out)

}



mysobol2 = function(outcome, inputs, best) {

	for (q in 1:ncol(inputs)) { inputs[ , q] = as.numeric(as.character(inputs[ , q])) }
	inputs = as.matrix(inputs)

	first.orders.0 = first.orders.1 = first.orders.2 = NULL
	names = NULL
	for (i in 1:ncol(inputs)) {
		progress.bar(i, 1, ncol(inputs), ncol(inputs))

		for (ii in i:ncol(inputs)) { 
			a1 = min(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)
			b1 = max(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)

			a2 = min(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)
			b2 = max(as.numeric(as.character(inputs[ , i])), na.rm = TRUE)

			means.0 = means.1 = means.2 = NULL
			
			for (j in 1:1000) {

				new = inputs
				new[ , i] = runif(1, a1, b1)
				new[ , ii] = runif(1, a2, b2)

				# new.fit = predict(mods[[best]]$Model, new, type = 'prob')
				new.fit.0 = new %*% as.matrix(coef(mods[[best]]$Model))

				new.fit = cbind(rep(0, nrow(new)), rep(0, nrow(new)), 
					exp(new.fit.0 - mods[[best]]$Model$zeta[2]) / 
						(1 + exp(new.fit.0 - mods[[best]]$Model$zeta[2]))
				)

				new.fit = cbind(rep(0, nrow(new)), rep(0, nrow(new)), 
					inv.logit(new.fit.0 - mods[[best]]$Model$zeta[2])
				)
				new.fit[ , 2] = 
					inv.logit(new.fit.0 - mods[[best]]$Model$zeta[2]) - 
					inv.logit(new.fit.0 - mods[[best]]$Model$zeta[1])
				new.fit[ , 1] = 1 - new.fit[ , 2] - new.fit[ , 3]

				# new.fit[ , 2] = (exp(new.fit.0 - mods[[95]]$Model$zeta[2]) / 
				# 	(1 + exp(new.fit.0 - mods[[95]]$Model$zeta[2]))) - 
				# 	(exp(new.fit.0 - mods[[95]]$Model$zeta[1]) / 
				# 	(1 + exp(new.fit.0 - mods[[95]]$Model$zeta[1])))
				# new.fit[ , 1] = 1 - new.fit[ , 2] - new.fit[ , 3]

				means.0 = c(means.0, mean(new.fit[ , 1]))
				means.1 = c(means.1, mean(new.fit[ , 2]))
				means.2 = c(means.2, mean(new.fit[ , 3]))
			}

			v.m.0 = var(means.0); v.m.1 = var(means.1); v.m.2 = var(means.2)

			first.orders.0 = c(first.orders.0, v.m.0)
			first.orders.1 = c(first.orders.1, v.m.1)
			first.orders.2 = c(first.orders.2, v.m.2)
			names = c(names, paste(i, ii, sep = '_'))
		}
	}

	out = list('mnames' = names, 'first.0' = first.orders.0, 
		'first.1' = first.orders.1, 'first.2' = first.orders.2)

	return(out)

}


