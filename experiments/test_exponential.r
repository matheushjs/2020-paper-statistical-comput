require(colorspace);

ALL.N = c(10, 25, 50, 75, 100, 200);
#ALL.N = 100*c(100, 200);
ITERATIONS = 200;

RATE = 1

plotMeans = NULL;
plotConf  = NULL;
df = NULL;
for(N in ALL.N){
	for(idx in 1:ITERATIONS){
		min.median = qexp(1 - (1 - 0.5)**(1/N), rate=RATE);

		# 1% and 5% quantiles of the minimum distribution
		quantile05 = qexp(1 - (1 - 0.05)**(1/N), rate=RATE);
		quantile01 = qexp(1 - (1 - 0.05)**(1/N), rate=RATE);

		data = list();
		data$samples = rexp(rate=RATE, n=N);
		data$min = min(data$samples);
		data$stddev = sd(data$samples);
		data$mean   = mean(data$samples);

		estim = list();
		estim[[1]] = data$min;
		estim[[2]] = data$min - data$min * (data$stddev / data$mean) / log10(N);
		estim[[3]] = data$min - data$min * (data$stddev / data$mean) / N;
		estim[[4]] = data$min - data$min * (data$stddev / data$mean) * sqrt(log(log(N)) / (2*N));
		estim[[5]] = data$min - data$min * (data$stddev / data$mean) * sqrt(-log(0.05/2) / (2*N));

		dist01 = rep(0, 5);
		for(i in 1:5){
			dist01[i] = (estim[[i]] - quantile01);
		}

		#print(dist01);

		lik = rep(0, 5);
		for(i in 1:5){
			samples = data$samples - estim[[i]];

			lik.f = function(param){
				-sum(dexp(samples, rate=param[1], log=TRUE));
			}

			result = optim(c(1/mean(samples)), lik.f, method="L-BFGS", upper=c(Inf), lower=c(0));
			lik[[i]] = -result$value;
		}

		likBase = sum(dexp(samples, rate=RATE, log=TRUE));

		#print(lik);

		df = rbind(df, c(N, dist01, lik, likBase));
	}

	slice = df[df[,1] == N,];
	means = colMeans(slice);
	slice[,7:11] = slice[,7:11] - slice[,12];
	stds  = apply(slice, 2, function(col) sd(col));

	#delta = stds;
	#delta = qt(0.05, df=ITERATIONS) * stds / ITERATIONS;
	delta = qnorm(0.01) * stds / ITERATIONS;
	#delta = apply(slice[,-1], 2, function(col){
	#	b = boot(col, function(data, idx) mean(data[idx]), R=100, sim="balanced");
	#	ci = boot.ci(b, type="norm");
	#	print(b);
	#	print(ci);
	#	0.05;
	#});

	plotMeans = rbind(plotMeans, means[-1]);
	plotConf  = rbind(plotConf, delta[-1]);
}

colnames(df) = c("N", paste("dist", 1:5, sep=""), paste("lik", 1:5, sep=""), "likBase");
#print(df);

graphics.off();
dev.new(width=0.8*8, height=0.8*5.4);

#print(plotMeans);
base = plotMeans[,11];
par(lwd=3, mar=c(3, 2, 0.2, 0.2));
palette(qualitative_hcl(palette="Dark 3", n=5));
plot(ALL.N, plotMeans[,6]   - base, type="l", ylim=c(-0.05, 1), col=1, log="x", xlab="", ylab="");
lines(ALL.N, plotMeans[,7]  - base, col=2, lty=2);
lines(ALL.N, plotMeans[,8]  - base, col=3, lty=3);
lines(ALL.N, plotMeans[,9]  - base, col=4, lty=4);
lines(ALL.N, plotMeans[,10] - base, col=5, lty=5);

for(i in 6:10){
	arrows(x0=ALL.N, y0=plotMeans[,i] - plotConf[,i] - base, y1=plotMeans[,i] + plotConf[,i] - base, angle=90, code=3, length=0, col="#00000099");
}

title(xlab="sample size", line=2);

legend("bottom", c(expression(c1), expression(c2), expression(c3), expression(c4), expression(c5)), col=1:5, lwd=3, lty=1:5, seg.len=4, bg="#FFFFFFBB");
#savePlot("test-exponential-likelihood.png");


dev.new(width=0.8*8, height=0.8*5);

#print(plotMeans);
#base = plotMeans[,11];
par(lwd=3, mar=c(3, 2, 0.2, 0.2));
palette(qualitative_hcl(palette="Dark 3", n=5));
plot(ALL.N, plotMeans[,1], type="l", ylim=c(0, 0.1), col=1, log="x", xlab="", ylab="");
lines(ALL.N, plotMeans[,2], col=2, lty=2);
lines(ALL.N, plotMeans[,3], col=3, lty=3);
lines(ALL.N, plotMeans[,4], col=4, lty=4);
lines(ALL.N, plotMeans[,5], col=5, lty=5);

for(i in 1:5){
	arrows(x0=ALL.N, y0=plotMeans[,i] - plotConf[,i], y1=plotMeans[,i] + plotConf[,i], angle=90, code=3, length=0);
}

title(xlab="sample size", line=2);

legend("topright", c(expression(c1), expression(c2), expression(c3), expression(c4), expression(c5)), col=1:5, lwd=3, lty=1:5, seg.len=4, bg="#FFFFFFBB");
