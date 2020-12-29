require(ggamma)

a_list = c(0.5, 1, 10);
b_list = c(0.5, 1, 10);
k_list = c(0.5, 1, 10);
n_list = c(10, 20, 50);

for(a in a_list)
for(b in b_list)
for(k in k_list)
for(n in n_list){
	data = rggamma(n=n, a=a, b=b, k=k);

	c = min(data) - sd(data) * sqrt(log(log(n)) / (2*n) );
	#print(c);

	newdata = data - c;
	#print(newdata);

	likelihood = function(p){
		allLogs = log(dggamma(newdata, a=p[1], b=p[2], k=p[3]));

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning("Low amount (<5) of warnings.", call.=FALSE);
		}

		-2 * sum(allLogs);
	}

	result = optim(c(1, 1, 1), likelihood, method="L-BFGS", lower=c(0, 0, 0));
	par = result$par;
	#print(c(a, b, k));
	#print(par);

	kl_divergence_f = function(x){
		probs1 = dggamma(x, a=a, b=b, k=k) / (1 - pggamma(c, a=a, b=b, k=k));
		probs2 = dggamma(x - c, a=par[1], b=par[2], k=par[3]);
		result = probs1 * log(probs1 / probs2);

		# probs2 is never zero
		# when probs1 is zero, the convention is KL = 0, so we enforce it here
		result[!is.finite(result)] = 0;
		result
	}

	kl_divergence = integrate(kl_divergence_f, c, Inf);
	print(kl_divergence);
}


#data = rggamma(n=N, a=1.55, b=0.61, k=10.22) + 20;
