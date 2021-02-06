require(ggamma)

setOption("width", 150);

b_list = c(0.5, 0.7, 0.9, 1, 1.5, 2, 5, 10, 20, 100);
k_list = c(0.5, 0.7, 0.9, 1, 1.5, 2, 5, 10, 20, 100);
n_list = c(5, 10, 20, 50);

set.seed(8396941); # generated randomly

with.infer.c = function(dislodge = 50, reducedInitParams=F,
						inferenceInitC = c("real", "estimator", "zero")){
	result = NULL;

	if(is.vector(inferenceInitC))
		inferenceInitC = inferenceInitC[1];

	#for(a in a_list)
	for(b in b_list)
	for(k in k_list)
	for(n in n_list){
		# Generate synthetic data
		data = rggamma(n=n, a=1, b=b, k=k) + dislodge;

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
		if(reducedInitParams){
			initParams = rbind(
				c(1,   1,   1)
			);
		} else {
			initParams = rbind(
				c(1, 0.5, 0.5),
				c(1, 0.5,   2),
				c(1,   2, 0.5),
				c(1,   2,   2)
			);
		}

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
			z = x - dislodge;
			if(c > dislodge){
				probs1[z >= 0] = dggamma(z[z >= 0], a=1, b=b, k=k) / (1 - pggamma(c - dislodge, a=1, b=b, k=k));
			} else {
				probs1[z >= 0] = dggamma(z[z >= 0], a=1, b=b, k=k);
			}

			z = x - c;
			probs2 = rep(0, length(x));
			if(c > dislodge){
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]);
			} else {
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]) / (1 - pggamma(dislodge - c, a=par[1], b=par[2], k=par[3]));
			}

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_c = integrate(divergence_cb, max(c, dislodge), Inf, stop.on.error=F);
		error = FALSE;
		error = error || (divergence_c$message != "OK");
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

		if(inferenceInitC == "real"){
			initC = dislodge;
		} else if(inferenceInitC == "estimator"){
			initC = c;
		} else if(inferenceInitC == "zero"){
			initC = 0;
		} else {
			stop("Wrong value of inferenceInitC");
		}

		allResults = list();
		allValues = rep(0, nrow(initParams));
		inf_counts = NULL;
		for(j in 1:nrow(initParams)){
			allResults[[j]] = optim(c(initParams[j,], initC), likelihood, method="L-BFGS", lower=c(0, 0, 0, -Inf), control=list(maxit=1000));
			inf_counts = c(inf_counts, count);
			count = 0;
			allValues[j] = allResults[[j]]$value;
		}

		bestResult_inf = allResults[[which.min(allValues)]];
		par = bestResult_inf$par;
		inf_par = par;

		divergence_cb = function(x){
			# Sometimes integrate() calls this function with negative argument, so we gotta prepare for that...
			z = x - dislodge;
			probs1 = rep(0, length(x));
			if(par[4] > dislodge){
				probs1[z >= 0] = dggamma(z[z >= 0], a=1, b=b, k=k) / (1 - pggamma(par[4] - dislodge, a=1, b=b, k=k));
			} else {
				probs1[z >= 0] = dggamma(z[z >= 0], a=1, b=b, k=k);
			}

			z = x - par[4];
			probs2 = rep(0, length(x));
			if(par[4] > dislodge){
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]);
			} else {
				probs2[z >= 0] = dggamma(z[z >= 0], a=par[1], b=par[2], k=par[3]) / (1 - pggamma(dislodge-par[4], a=par[1], b=par[2], k=par[3]));
			}

			result = abs(probs1 - probs2);

			result[!is.finite(result)] = 0;
			result
		}

		divergence_inf = integrate(divergence_cb, max(par[4], dislodge), Inf, stop.on.error=F);
		error = error || (divergence_inf$message != "OK");

		print(paste("N: ", n, "  |  KL: ", divergence_inf$value, "  |  likelihood: ", bestResult_inf$value, "  |  convergence: ", bestResult_inf$convergence, "  |  estimate: ", c, "  |  params: ", par[1], par[2], par[3], par[4]));
		print("------------------------------------------------------------------------------------------");

		result = rbind(result, c(1, b, k, n, error,
					   estim_par[1], estim_par[2], estim_par[3], c, divergence_c$value, bestResult_c$value, bestResult_c$value + 6, sum(estim_counts),
					   inf_par[1], inf_par[2], inf_par[3], inf_par[4], divergence_inf$value, bestResult_inf$value, bestResult_inf$value + 8, sum(inf_counts)));
		colnames(result) = c("orig_a", "orig_b", "orig_k", "n", "int_err",
							 "estim_a", "estim_b", "estim_k", "estim_c", "estim_kl", "estim_loglik", "estim_aic", "estim_itercount",
							 "inf_a", "inf_b", "inf_k", "inf_c", "inf_kl", "inf_loglik", "inf_aic", "inf_itercount");
	}

	return(result);
}

allResultsList = list();
N = 10;

for(i in 1:N){
	result = with.infer.c();
	allResultsList[[i]] = result;
}

allResults = NULL;
for(i in 1:N){
	single = allResultsList[[i]];
	allResults = rbind(allResults, single);
}

# Simple confidence interval based on t-student distribution
plus.minus = function(data){
	q = qt(0.95, df=length(data) - 1);
	return( q * sd(data) / sqrt(length(data)) );
}

generate.report = function(allResults){
	idx5  = allResults[,"n"] == 5;
	idx10 = allResults[,"n"] == 10;
	idx20 = allResults[,"n"] == 20;
	idx50 = allResults[,"n"] == 50;

	A = allResults[,"estim_kl"];
	B = allResults[,"inf_kl"];
	print(paste("Number of times our estimator is better (Absolute Difference Divergence): ",
			sum(A < B) / length(B) * 100 , "%",
			" (",
			sum(A[idx5] < B[idx5]) / length(B[idx5]) * 100,
			sum(A[idx10] < B[idx10]) / length(B[idx10]) * 100,
			sum(A[idx20] < B[idx20]) / length(B[idx20]) * 100,
			sum(A[idx50] < B[idx50]) / length(B[idx50]) * 100,
			" for N = 5, 10, 20, 50)"), sep="");

	idx = which(allResults[,"orig_b"] * allResults[,"orig_k"] > 1);

	A = allResults[,"estim_kl"];
	B = allResults[,"inf_kl"];
	print(paste("For two-tailed cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_kl"];
		B = slice[,"inf_kl"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_kl"];
	B = allResults[,"inf_kl"];
	print(paste("For one-tailed cases: ",
			sum(A[-idx] < B[-idx]) / length(B[-idx]) * 100 , "%",
			" in ", length(B[-idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[-idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_kl"];
		B = slice[,"inf_kl"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_kl"];
	B = allResults[,"inf_kl"];
	idx = which(allResults[,"orig_k"] == 1);
	print(paste("For Weibull (k=1) cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_kl"];
		B = slice[,"inf_kl"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_kl"];
	B = allResults[,"inf_kl"];
	idx = which(allResults[,"orig_b"] == 1);
	print(paste("For gamma (b=1) cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_kl"];
		B = slice[,"inf_kl"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	cat("\n");

	A = allResults[,"estim_aic"];
	B = allResults[,"inf_aic"];
	print(paste("Number of times our estimator is better (Tentative Akaike Information Criterion): ",
			sum(A < B) / length(B) * 100 , "%",
			" (",
			sum(A[idx5] < B[idx5]) / length(B[idx5]) * 100,
			sum(A[idx10] < B[idx10]) / length(B[idx10]) * 100,
			sum(A[idx20] < B[idx20]) / length(B[idx20]) * 100,
			sum(A[idx50] < B[idx50]) / length(B[idx50]) * 100,
			" for N = 5, 10, 20, 50)"), sep="");

	idx = which(allResults[,"orig_b"] * allResults[,"orig_k"] > 1);

	A = allResults[,"estim_aic"];
	B = allResults[,"inf_aic"];
	print(paste("For two-tailed cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_aic"];
		B = slice[,"inf_aic"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_aic"];
	B = allResults[,"inf_aic"];
	print(paste("For one-tailed cases: ",
			sum(A[-idx] < B[-idx]) / length(B[-idx]) * 100 , "%",
			" in ", length(B[-idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[-idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_aic"];
		B = slice[,"inf_aic"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_aic"];
	B = allResults[,"inf_aic"];
	idx = which(allResults[,"orig_k"] == 1);
	print(paste("For Weibull (k=1) cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_aic"];
		B = slice[,"inf_aic"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_aic"];
	B = allResults[,"inf_aic"];
	idx = which(allResults[,"orig_b"] == 1);
	print(paste("For gamma (b=1) cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_aic"];
		B = slice[,"inf_aic"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	cat("\n");

	A = allResults[,"estim_itercount"];
	B = allResults[,"inf_itercount"];
	print(paste("By inferring c we needed ", sum(A) / sum(B) * 100, "% more evaluations of the likelihood function"));
	cat("\n");


	A = allResults[,"estim_c"];
	B = allResults[,"inf_c"];
	print(paste("Our estimator was below the inferred c in ", sum(A < B) / length(B) * 100, "% of the trials."));

	idx = which(allResults[,"orig_b"] * allResults[,"orig_k"] > 1);

	print(paste("For two-tailed cases: ",
			sum(A[idx] < B[idx]) / length(B[idx]) * 100 , "%",
			" in ", length(B[idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_c"];
		B = slice[,"inf_c"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	A = allResults[,"estim_c"];
	B = allResults[,"inf_c"];
	print(paste("For one-tailed cases: ",
			sum(A[-idx] < B[-idx]) / length(B[-idx]) * 100 , "%",
			" in ", length(B[-idx]), " trials"), sep="");
	cat("\t => ");
	for(N in c(5, 10, 20, 50)){
		slice = allResults[-idx,];
		slice = slice[slice[,"n"] == N,];
		A = slice[,"estim_c"];
		B = slice[,"inf_c"];
		cat(sum(A < B) / length(B) * 100, " ");
	}
	cat(" for N = 5, 10, 20, 50\n");

	cat("\n");
}

allResultsList2 = list();
N = 10;

for(i in 1:N){
	result = with.infer.c(reducedInitParams=T);
	allResultsList2[[i]] = result;
}

allResults2 = NULL;
for(i in 1:N){
	single = allResultsList2[[i]];
	allResults2 = rbind(allResults2, single);
}

# Check for possible errors, and remove them
notErrors = allResults[,"int_err"] != 1;
print(paste("There were ", sum(!notErrors), " errors"));

notErrors2 = allResults2[,"int_err"] != 1;
print(paste("There were ", sum(!notErrors2), " errors"));

cat("===== FULL LIST OF PARAMETERS =====\n");
generate.report(allResults[notErrors,]);
cat("===== REDUCED LIST OF PARAMETERS =====\n");
generate.report(allResults2[notErrors2,]);

# Need this filtering before performing the following calculations
allNotErrors = notErrors && notErrors2;
allResults  = allResults[allNotErrors,];
allResults2 = allResults2[allNotErrors,];

A = allResults2[,"estim_loglik"] - allResults[,"estim_loglik"];
B = allResults2[,"inf_loglik"] - allResults[,"inf_loglik"];

print(paste("Estimating c had an average worsening of ", mean(A), "+-", plus.minus(A), " points of log likelihood"));
print(paste("Inferring c had an average worsening of ",  mean(B), "+-", plus.minus(B), " points of log likelihood"));
