require(elfDistr);

N = 10;
tries = 10;

df = NULL;

R.dkwcwg = function(x, alpha, beta, gamma, a, b, log=F){
	aux1 = (gamma*x)^beta;
	aux2 = exp(-(gamma*x)^beta);
	aux3 = alpha + aux2*(1-alpha);
	cdf_g = alpha*(1 - aux2) / aux3;
	pdf_g = (alpha*beta*gamma*aux1/(gamma*x)*aux2) / (aux3*aux3);
	logresult = log(a) + log(b) + log(pdf_g) + (a-1)*log(cdf_g) + (b-1)*log(1 - cdf_g^a);

	if(log){
		return(logresult);
	} else {
		return(exp(logresult));
	}
}

R.pkwcwg = function(x, alpha, beta, gamma, a, b, log=F){
	result = 1 - ( 1 - alpha^a * (
				(1-exp(-(gamma*x)^beta)) / (alpha + (1-alpha)*exp(-(gamma*x)**beta))
			)^a
		)^b;
	if(log) {
		return(log(result));
	} else {
		return(result);
	}
}

R.qkwcwg = function(p, alpha, beta, gamma, a, b, lower.tail=T){
	aux = (1 - (1-p)^(1/b))^(1/a);
	result = log((alpha + (1-alpha)*aux)/(alpha*(1-aux)))^(1/beta) / gamma;
	if(lower.tail) {
		return(result);
	} else {
		return(1 - result);
	}
}

R.rkwcwg = function(n, alpha, beta, gamma, a, b){
	R.qkwcwg(runif(n=n, 0, 1), alpha, beta, gamma, a, b);
}

for(alpha in c(0.1, 0.5, 1))
for(beta in c(0.5, 1, 1.001, 1.1, 2, 20, 50)) # Weibull shape
for(gamma in c(1/beta, 1/(2*beta), beta, 5*beta, 1/20, 1)) # Weibull scale
for(a in c(0.05, 0.8, 1, 5))
for(b in c(0.5, 1, 10)){ # 5 and 20 good too
#for(a in c(1))
#for(b in c(15)){

	min.median = qkwcwg(1 - (1 - 0.5)**(1/N), alpha, beta, gamma, a, b);

	# 1% and 5% quantiles of the minimum distribution
	quantile05 = qkwcwg(1 - (1 - 0.05)**(1/N), alpha, beta, gamma, a, b);
	quantile01 = qkwcwg(1 - (1 - 0.01)**(1/N), alpha, beta, gamma, a, b);

	# Degenerate model, we don't want it
	if(min.median < 1e-5)
		next;

	data = matrix(rkwcwg(n=N*tries, alpha, beta,gamma, a, b), nrow=tries);

	# Filter NA/Infinite elements
	failureRate = sum(!is.finite(data)) / tries;
	if(failureRate > 0.05*N)
		cat("High failure rate: ", failureRate);

	mins   = apply(data, 1, function(row) min(row[is.finite(row)]));
	stddev = apply(data, 1, function(row) sd(row[is.finite(row)]));
	means  = apply(data, 1, function(row) mean(row[is.finite(row)]));

	data.mean = mean(data);

	estim = list();
	estim[[1]] = mins
	#estim[[2]] = mins - mins * (stddev / means) / log10(N);
	#estim[[3]] = mins - stddev / N;
	#estim[[4]] = mins - stddev * sqrt(log(log(N)) / (2*N));
	#estim[[5]] = mins - stddev * sqrt(-log(0.05/2) / (2*N));
	estim[[2]] = mins - mins * (stddev / means) / log10(N);
	estim[[3]] = mins - mins * (stddev / means) / N;
	estim[[4]] = mins - mins * (stddev / means) * sqrt(log(log(N)) / (2*N));
	estim[[5]] = mins - mins * (stddev / means) * sqrt(-log(0.05/2) / (2*N));

	# We will get:
	# 1) The mean distance from the estimator to the 01 and 05 min-quantiles
	# 2) The mean cdf on the estimator

	dist05 = list();
	dist01 = list();
	distMin = list();
	cdf = list();
	for(i in 1:5){
		dist05[[i]] = (estim[[i]] - quantile05) / mean(data);
		dist01[[i]] = (estim[[i]] - quantile01) / mean(data);
		distMin[[i]] = (mins - estim[[i]]) / mean(data);

		# Negative values should give cdf = 0 for the KW-CWG
		aux = pkwcwg(estim[[i]], alpha, beta, gamma, a, b);
		aux[which(estim[[i]] <= 0)] = 0;
		cdf[[i]] = aux;
	}

	means.dist05  = rep(0, length(estim));
	means.dist01  = rep(0, length(estim));
	means.distMin = rep(0, length(estim));
	means.cdf     = rep(0, length(estim));
	for(i in 1:length(estim)){
		means.dist05[i] = mean(dist05[[i]]);
		means.dist01[i] = mean(dist01[[i]]);
		means.distMin[i] = mean(distMin[[i]]);
		means.cdf[i]    = mean(cdf[[i]]);
	}

	names = c(
		"alpha", "beta", "gamma", "a", "b",
		"data.stddev", "data.mean",
		paste("dist05-c", 1:5, sep=""),
		paste("dist01-c", 1:5, sep=""),
		paste("distMin-c", 1:5, sep=""),
		paste("cdf", 1:5, sep="")
	);
	
	row = c(alpha, beta, gamma, a, b,
			sd(data), mean(data), # yes we do sd() and mean() over whole matrix
			means.dist05, means.dist01, means.distMin, means.cdf);
	print(names);
	print(row);

	df = rbind(df, row);
	colnames(df) = names;
}

print(df);
