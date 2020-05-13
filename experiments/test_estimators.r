require(elfDistr);

N = 100;
tries = 10;

df = NULL;

for(alpha in c(0.2, 0.5, 0.8))
for(beta in c(0.2, 0.5, 1, 2, 5, 10, 20))
for(gamma in c(0.2, 0.5, 1, 2, 5, 10, 20))
for(a in c(0.2, 0.5, 0.8, 1, 1.2, 2, 5))
for(b in c(0.2, 0.5, 0.8, 1, 1.2, 2, 5)){

	min.median = qkwcwg(1 - (1 - 0.5)**(1/N), alpha, beta, gamma, a, b);
	min.cdf    = pkwcwg(min.median, alpha, beta, gamma, a, b);

	data = matrix(rkwcwg(n=N*tries, alpha, beta,gamma, a, b), nrow=tries);
	mins = apply(data, 1, min);
	stddev = apply(data, 1, sd);
	means = apply(data, 1, mean);

	estim1 = mins - mins * (stddev / means) / log10(N);
	estim2 = mins - stddev / N;
	estim3 = mins - stddev * sqrt(log(log(N)) / (2*n));
	estim4 = mins - stddev * sqrt(-log(0.05/2) / (2*n));

	# rkwcwg often results in Inf. Here we take the maximum of those that are finite.
	estim1 = max(0, estim1[is.finite(estim1)])
	estim2 = max(0, estim2[is.finite(estim2)])
	estim3 = max(0, estim3[is.finite(estim3)])
	estim4 = max(0, estim4[is.finite(estim4)])

	cdfs = pkwcwg(c(estim1, estim2, estim3, estim4), alpha, beta, gamma, a, b);

	#min.max = qkwcwg(1 - (1 - 0.99)**(1/N), alpha, beta, gamma, a, b);
	#min.mean = integrate(function(x) n*(1 - pkwcwg(x, alpha, beta, gamma, a, b))**(n-1)*dkwcwg(x, alpha, beta, gamma, a, b), lower=0, upper=min.max)$value;

	names = c("alpha", "beta", "gamma", "a", "b", "min.median", "min.cdf", "data.stddev", "data.mean", "c1", "c2", "c3", "c4", "cdf1", "cdf2", "cdf3", "cdf4");
	row = c(alpha, beta, gamma, a, b, min.median, min.cdf, sd(data), mean(data), estim1, estim2, estim3, estim4, cdfs);
	print(names);
	print(row);

	df = rbind(df, row);
	colnames(df) = names;
}

print(df);
