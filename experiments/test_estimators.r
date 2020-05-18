require(elfDistr);

N = 100;
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

	data = matrix(rkwcwg(n=N*tries, alpha, beta,gamma, a, b), nrow=tries);
	mins = apply(data, 1, min);
	stddev = apply(data, 1, sd);
	means = apply(data, 1, mean);

	estim1 = mins - mins * (stddev / means) / log10(N);
	estim2 = mins - stddev / N;
	estim3 = mins - stddev * sqrt(log(log(N)) / (2*N));
	estim4 = mins - stddev * sqrt(-log(0.05/2) / (2*N));

	# rkwcwg often results in Inf. Here we take the maximum of those that are finite.
	estim1 = max(estim1[is.finite(estim1)]);
	estim2 = max(estim2[is.finite(estim2)]);
	estim3 = max(estim3[is.finite(estim3)]);
	estim4 = max(estim4[is.finite(estim4)]);

	aux = c(max(0, estim1), max(0, estim2), max(0, estim3), max(0, estim4));
	cdfs = pkwcwg(aux, alpha, beta, gamma, a, b);

	# 1% and 5% quantiles of the minimum distribution
	quantile05 = qkwcwg(1 - (1 - 0.05)**(1/N), alpha, beta, gamma, a, b);
	quantile01 = qkwcwg(1 - (1 - 0.01)**(1/N), alpha, beta, gamma, a, b);

	#min.max = qkwcwg(1 - (1 - 0.99)**(1/N), alpha, beta, gamma, a, b);
	#min.mean = integrate(function(x) n*(1 - pkwcwg(x, alpha, beta, gamma, a, b))**(n-1)*dkwcwg(x, alpha, beta, gamma, a, b), lower=0, upper=min.max)$value;

	names = c("alpha", "beta", "gamma", "a", "b", "min.median", "data.stddev", "data.mean", "min.05", "min.01", "c1", "c2", "c3", "c4", "cdf1", "cdf2", "cdf3", "cdf4");
	row = c(alpha, beta, gamma, a, b, min.median, sd(data), mean(data), quantile05, quantile01, estim1, estim2, estim3, estim4, cdfs);
	print(names);
	print(row);

	df = rbind(df, row);
	colnames(df) = names;
}

print(df);
