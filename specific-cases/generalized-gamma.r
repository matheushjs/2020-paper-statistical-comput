require(ggamma)

#a_list = c(0.5, 1, 10);
b_list = c(0.5, 1, 10);
k_list = c(0.5, 1, 10, 200);
n_list = c(5, 10, 20, 50);

#for(a in a_list)
for(b in b_list)
for(k in k_list)
for(n in n_list){
	data = rggamma(n=n, a=1, b=b, k=k);

	c_old = min(data) - sd(data) * sqrt(log(log(n)) / (2*n) );
	c     = min(data) - sd(data) * sqrt( -log(0.01/2) / (2*n) );
	c = max(0, c);
	#print(c);

	newdata = data - c;
	#print(newdata);

	likelihood = function(p){
		allLogs = log(dggamma(newdata, a=1, b=p[1], k=p[2]));

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning("Low amount (<5) of warnings.", call.=FALSE);
		}

		-2 * sum(allLogs);
	}

	initParams = rbind(
		c(0.5, 0.5),
		c(0.5,   1),
		c(0.5,   2),
		c(  1, 0.5),
		c(  1,   1),
		c(  1,   2),
		c(  2, 0.5),
		c(  2,   1),
		c(  2,   2)
	);
	allResults = list();
	allValues = rep(0, nrow(initParams));
	for(j in 1:nrow(initParams)){
		allResults[[j]] = optim(initParams[j,], likelihood, method="L-BFGS", lower=c(0, 0, 0));
		allValues[j] = allResults[[j]]$value;
	}

	bestResult_c = allResults[[which.min(allValues)]];
	par = bestResult_c$par;
	#print(c(a, b, k));
	#print(par);

	kl_divergence_f = function(x){
		probs1 = rep(0, length(x));
		probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k) / (1 - pggamma(c, a=1, b=b, k=k));

		probs2 = dggamma(x - c, a=1, b=par[1], k=par[2]);
		result = probs1 * log(probs1 / probs2);

		# probs2 is never zero
		# when probs1 is zero, the convention is KL = 0, so we enforce it here
		result[!is.finite(result)] = 0;
		result
	}

	kl_divergence_c = integrate(kl_divergence_f, c, Inf);
	print("------------------------------------------------------------------------------------------");
	print(paste("KL Divergence: ", kl_divergence_c$value, "  |  likelihood: ", bestResult_c$value, "  |  convergence: ", bestResult_c$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2]));


	likelihood = function(p){
		allLogs = log(dggamma(data - p[3], a=1, b=p[1], k=p[2]));

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning("Low amount (<5) of warnings.", call.=FALSE);
		}

		-2 * sum(allLogs);
	}

	allResults = list();
	allValues = rep(0, nrow(initParams));
	for(j in 1:nrow(initParams)){
		allResults[[j]] = optim(c(initParams[j,], 0), likelihood, method="L-BFGS", lower=c(0, 0, 0), control=list(maxit=1000));
		allValues[j] = allResults[[j]]$value;
	}

	bestResult_inf = allResults[[which.min(allValues)]];
	par = bestResult_inf$par;
	#print(par);

	kl_divergence_f = function(x){
		probs1 = rep(0, length(x));
		probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k) / (1 - pggamma(par[3], a=1, b=b, k=k));

		probs2 = dggamma(x - par[3], a=1, b=par[1], k=par[2]);
		result = probs1 * log(probs1 / probs2);

		# probs2 is never zero
		# when probs1 is zero, the convention is KL = 0, so we enforce it here
		result[!is.finite(result)] = 0;
		result
	}

	kl_divergence_inf = integrate(kl_divergence_f, c, Inf);
	print(paste("KL Divergence: ", kl_divergence_inf$value, "  |  likelihood: ", bestResult_inf$value, "  |  convergence: ", bestResult_inf$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3]));
	print("------------------------------------------------------------------------------------------");
}


#data = rggamma(n=N, a=1.55, b=0.61, k=10.22) + 20;
