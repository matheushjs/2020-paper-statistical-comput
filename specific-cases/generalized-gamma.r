require(ggamma)

setOption("width", 150);

#a_list = c(0.5, 1, 10);
b_list = c(0.7, 1, 10);
k_list = c(0.7, 1, 10, 20);
n_list = c(5, 10, 20, 50);

#set.seed(1);

with.infer.c = function(){
	result = NULL;

	#for(a in a_list)
	for(b in b_list)
	for(k in k_list)
	for(n in n_list){
		# Generate synthetic data
		data = rggamma(n=n, a=1, b=b, k=k);

		# Get the estimator for populational minimum (or low quantile)
		c = min(data) - sd(data) * sqrt(log(log(n)) / (2*n) );

		# We will count number of calls to the likelihood function
		count = 0;

		likelihood = function(p){
			count <<- count + 1;
			allLogs = log(dggamma(data - c, a=p[1], b=p[2], k=p[3]));

			problems = which(!is.finite(allLogs))
			allLogs[problems] = log(1e-300); # Very bad result merely to force optim to continue optimizing

			-2 * sum(allLogs);
		}

		# Try multiple initial parameters
		initParams = rbind(
			#c(1, 0.5, 0.5),
			#c(1, 0.5,   1),
			#c(1, 0.5,   2),
			#c(1,   1, 0.5),
			c(1,   1,   1)
			#c(1,   1,   2),
			#c(1,   2, 0.5),
			#c(1,   2,   1),
			#c(1,   2,   2)
		);

		allResults = list();                    # Return values from optim()
		allValues = rep(0, nrow(initParams));   # Optimized -2l (l = log likelihood) values
		estim_counts = NULL;                    # Number of calls to likelihood() during optim()
		for(j in 1:nrow(initParams)){
			allResults[[j]] = optim(initParams[j,], likelihood, method="L-BFGS", lower=c(0, 0, 0));
			estim_counts = c(estim_counts, count);
			count = 0;
			allValues[j] = allResults[[j]]$value;
		}

		bestResult_c = allResults[[which.min(allValues)]];
		par = bestResult_c$par;
		estim_par = par;
		#print(c(a, b, k));
		#print(par);

		# Measures divergence between two densities
		divergence_cb = function(x){
			# Sometimes integrate() calls this function with negative argument, so we gotta prepare for that...
			probs1 = rep(0, length(x));
			if(c > 0){
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k) / (1 - pggamma(c, a=1, b=b, k=k));
			} else {
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k);
			}

			z = x - c;
			probs2 = rep(0, length(x));
			if(c > 0){
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]);
			} else {
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]) / (1 - pggamma(-c, a=par[1], b=par[2], k=par[3]));
			}

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_c = integrate(divergence_cb, c, Inf);
		print("------------------------------------------------------------------------------------------");
		print(paste("N: ", n, "  |  KL: ", divergence_c$value, "  |  likelihood: ", bestResult_c$value, "  |  convergence: ", bestResult_c$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3]));

		# Repeat everything, but now for the case where we infer c
		count = 0;

		likelihood = function(p){
			count <<- count + 1;
			allLogs = log(dggamma(data - p[4], a=p[1], b=p[2], k=p[3]));

			problems = which(!is.finite(allLogs))
			allLogs[problems] = log(1e-300); # Very bad result merely to force optim to continue optimizing

			-2 * sum(allLogs);
		}

		allResults = list();
		allValues = rep(0, nrow(initParams));
		inf_counts = NULL;
		for(j in 1:nrow(initParams)){
			allResults[[j]] = optim(c(initParams[j,], c), likelihood, method="L-BFGS", lower=c(0, 0, 0, -Inf), control=list(maxit=1000));
			inf_counts = c(inf_counts, count);
			count = 0;
			allValues[j] = allResults[[j]]$value;
		}

		bestResult_inf = allResults[[which.min(allValues)]];
		par = bestResult_inf$par;
		inf_par = par;

		divergence_cb = function(x){
			# Sometimes integrate() calls this function with negative argument, so we gotta prepare for that...
			probs1 = rep(0, length(x));
			if(par[4] > 0){
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k) / (1 - pggamma(par[4], a=1, b=b, k=k));
			} else {
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k);
			}

			z = x - par[4];
			probs2 = rep(0, length(x));
			if(par[4] > 0){
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]);
			} else {
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]) / (1 - pggamma(-par[4], a=par[1], b=par[2], k=par[3]));
			}

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_inf = integrate(divergence_cb, par[4], Inf);
		print(paste("N: ", n, "  |  KL: ", divergence_inf$value, "  |  likelihood: ", bestResult_inf$value, "  |  convergence: ", bestResult_inf$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3], par[4]));
		print("------------------------------------------------------------------------------------------");

		result = rbind(result, c(1, b, k, n,
					   estim_par[1], estim_par[2], estim_par[3], c, divergence_c$value, bestResult_c$value, 2*bestResult_c$value + 6, sum(estim_counts),
					   inf_par[1], inf_par[2], inf_par[3], inf_par[4], divergence_inf$value, bestResult_inf$value, 2*bestResult_inf$value + 8, sum(inf_counts)));
		colnames(result) = c("orig_a", "orig_b", "orig_k", "n",
							 "estim_a", "estim_b", "estim_k", "estim_c", "estim_kl", "estim_loglik", "estim_aic", "estim_itercount",
							 "inf_a", "inf_b", "inf_k", "inf_c", "inf_kl", "inf_loglik", "inf_aic", "inf_itercount");
	}

	return(result);
}

without.infer.c = function(){
	result = NULL;

	#for(a in a_list)
	for(b in b_list)
	for(k in k_list)
	for(n in n_list){
		# Generate synthetic data
		data = rggamma(n=n, a=1, b=b, k=k);

		# Get the estimator for populational minimum (or low quantile)
		c = min(data) - sd(data) * sqrt(log(log(n)) / (2*n) );

		# We will count number of calls to the likelihood function
		count = 0;

		likelihood = function(p){
			count <<- count + 1;
			allLogs = log(dggamma(data - c, a=p[1], b=p[2], k=p[3]));

			problems = which(!is.finite(allLogs))
			allLogs[problems] = log(1e-300); # Very bad result merely to force optim to continue optimizing

			-2 * sum(allLogs);
		}

		# Try multiple initial parameters
		initParams = rbind(
			#c(1, 0.5, 0.5),
			#c(1, 0.5,   1),
			#c(1, 0.5,   2),
			#c(1,   1, 0.5),
			c(1,   1,   1)
			#c(1,   1,   2),
			#c(1,   2, 0.5),
			#c(1,   2,   1),
			#c(1,   2,   2)
		);

		allResults = list();                    # Return values from optim()
		allValues = rep(0, nrow(initParams));   # Optimized -2l (l = log likelihood) values
		estim_counts = NULL;                    # Number of calls to likelihood() during optim()
		for(j in 1:nrow(initParams)){
			allResults[[j]] = optim(initParams[j,], likelihood, method="L-BFGS", lower=c(0, 0, 0));
			estim_counts = c(estim_counts, count);
			count = 0;
			allValues[j] = allResults[[j]]$value;
		}

		bestResult_c = allResults[[which.min(allValues)]];
		par = bestResult_c$par;
		estim_par = par;
		#print(c(a, b, k));
		#print(par);

		# Measures divergence between two densities
		divergence_cb = function(x){
			# Sometimes integrate() calls this function with negative argument, so we gotta prepare for that...
			probs1 = rep(0, length(x));
			if(c > 0){
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k) / (1 - pggamma(c, a=1, b=b, k=k));
			} else {
				probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k);
			}

			z = x - c;
			probs2 = rep(0, length(x));
			if(c > 0){
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]);
			} else {
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]) / (1 - pggamma(-c, a=par[1], b=par[2], k=par[3]));
			}

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_c = integrate(divergence_cb, c, Inf);
		print("------------------------------------------------------------------------------------------");
		print(paste("N: ", n, "  |  KL: ", divergence_c$value, "  |  likelihood: ", bestResult_c$value, "  |  convergence: ", bestResult_c$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3]));

		# Repeat everything, but now for the case where we infer c
		count = 0;

		likelihood = function(p){
			count <<- count + 1;
			allLogs = log(dggamma(data, a=p[1], b=p[2], k=p[3]));

			problems = which(!is.finite(allLogs))
			allLogs[problems] = log(1e-300); # Very bad result merely to force optim to continue optimizing

			-2 * sum(allLogs);
		}

		allResults = list();
		allValues = rep(0, nrow(initParams));
		inf_counts = NULL;
		for(j in 1:nrow(initParams)){
			allResults[[j]] = optim(c(initParams[j,]), likelihood, method="L-BFGS", lower=c(0, 0, 0), control=list(maxit=1000));
			inf_counts = c(inf_counts, count);
			count = 0;
			allValues[j] = allResults[[j]]$value;
		}

		bestResult_inf = allResults[[which.min(allValues)]];
		par = bestResult_inf$par;
		inf_par = par;

		divergence_cb = function(x){
			# Sometimes integrate() calls this function with negative argument, so we gotta prepare for that...
			probs1 = rep(0, length(x));
			probs1[x >= 0] = dggamma(x[x >= 0], a=1, b=b, k=k);

			probs2 = rep(0, length(x));
			probs2[x >= 0] = dggamma(x[x >= 0], a=par[1], b=par[2], k=par[3]);

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_inf = integrate(divergence_cb, 0, Inf);
		print(paste("N: ", n, "  |  KL: ", divergence_inf$value, "  |  likelihood: ", bestResult_inf$value, "  |  convergence: ", bestResult_inf$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3]));
		print("------------------------------------------------------------------------------------------");

		result = rbind(result, c(1, b, k, n,
					   estim_par[1], estim_par[2], estim_par[3], c, divergence_c$value, bestResult_c$value, 2*bestResult_c$value + 6, sum(estim_counts),
					   inf_par[1], inf_par[2], inf_par[3], divergence_inf$value, bestResult_inf$value, 2*bestResult_inf$value + 8, sum(inf_counts)));
		colnames(result) = c("orig_a", "orig_b", "orig_k", "n",
							 "estim_a", "estim_b", "estim_k", "estim_c", "estim_kl", "estim_loglik", "estim_aic", "estim_itercount",
							 "inf_a", "inf_b", "inf_k", "inf_kl", "inf_loglik", "inf_aic", "inf_itercount");
	}

	return(result);
}

result = without.infer.c();
#result = with.infer.c();

#data = rggamma(n=N, a=1.55, b=0.61, k=10.22) + 20;
