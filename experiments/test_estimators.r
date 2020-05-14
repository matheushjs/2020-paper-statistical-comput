require(elfDistr);

N = 100;
tries = 10;

df = NULL;

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

	quantile05 = qkwcwg(0.05, alpha, beta, gamma, a, b);
	quantile01 = qkwcwg(0.01, alpha, beta, gamma, a, b);

	#min.max = qkwcwg(1 - (1 - 0.99)**(1/N), alpha, beta, gamma, a, b);
	#min.mean = integrate(function(x) n*(1 - pkwcwg(x, alpha, beta, gamma, a, b))**(n-1)*dkwcwg(x, alpha, beta, gamma, a, b), lower=0, upper=min.max)$value;

	names = c("alpha", "beta", "gamma", "a", "b", "min.median", "data.stddev", "data.mean", "data.05", "data.01", "c1", "c2", "c3", "c4", "cdf1", "cdf2", "cdf3", "cdf4");
	row = c(alpha, beta, gamma, a, b, min.median, sd(data), mean(data), quantile05, quantile01, estim1, estim2, estim3, estim4, cdfs);
	print(names);
	print(row);

	df = rbind(df, row);
	colnames(df) = names;
}

print(df);
