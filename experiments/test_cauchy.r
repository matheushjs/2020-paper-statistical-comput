require(colorspace);
require(boot);

ALL.N = c(10, 25, 50, 75, 100, 200, 500, 1000, 2000, 5000);
#ALL.N = 100*c(100, 200);
ITERATIONS = 200;

LOCATION = 0;
SCALE    = 1;

plotMeans = NULL;
plotConf  = NULL;
df = NULL;
for(N in ALL.N){
	for(idx in 1:ITERATIONS){
		min.median = qcauchy(1 - (1 - 0.5)**(1/N), location=LOCATION, scale=SCALE);

		# 1% and 5% quantiles of the minimum distribution
		quantile05 = qcauchy(1 - (1 - 0.05)**(1/N), location=LOCATION, scale=SCALE);
		quantile01 = qcauchy(1 - (1 - 0.05)**(1/N), location=LOCATION, scale=SCALE);

		data = list();
		data$samples = rcauchy(location=LOCATION, scale=SCALE, n=N);
		data$min = min(data$samples);
		data$stddev = sd(data$samples);
		data$mean   = mean(data$samples);

		estim = list();
		estim[[1]] = data$min;
		estim[[2]] = data$min - abs(data$min) * abs(data$stddev / data$mean) / log10(N);
		estim[[3]] = data$min - abs(data$min) * data$stddev / N;
		estim[[4]] = data$min - abs(data$min) * data$stddev * sqrt(log(log(N)) / (2*N));
		estim[[5]] = data$min - abs(data$min) * data$stddev * sqrt(-log(0.05/2) / (2*N));

		dist01 = rep(0, 5);
		for(i in 1:5){
			dist01[i] = pcauchy(estim[[i]]);
		}

		#print(rbind(estim, dist01));

		lik = rep(0, 5);
		for(i in 1:5){
			samples = data$samples - estim[[i]];

			lik.f = function(param){
				r = -prod(dcauchy(samples, location=param[1], scale=param[2]));
				r;
			}

			result = optim(c(LOCATION, SCALE), lik.f, method="Nelder");
			lik[[i]] = -result$value;
		}

		likBase = sum(dcauchy(samples, location=LOCATION, scale=SCALE, log=TRUE));

		#print(lik);

		df = rbind(df, c(N, dist01, lik, likBase));
	}

	slice = df[df[,1] == N,];
	means = colMeans(slice);
	stds  = apply(slice, 2, function(col) sd(col));

	delta = stds / 10;
	#delta = qt(0.05, df=ITERATIONS) * stds / ITERATIONS;
	#delta = apply(slice[,-1], 2, function(col){
	#	b = boot(col, function(data, idx) mean(data[idx]), R=100, sim="balanced");
	#	ci = boot.ci(b, type="norm");
	#	print(b);
	#	print(ci);
	#	0.05;
	#});
	print(delta[-1]);

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

title(xlab="sample size", line=2);

legend("bottom", c(expression(textstyle(min)), expression(c1), expression(c2), expression(c3), expression(c4)), col=1:5, lwd=3, lty=1:5, seg.len=4, bg="#FFFFFFBB");
#savePlot("test-exponential-likelihood.png");


dev.new(width=0.8*8, height=0.8*5);

#print(plotMeans);
#base = plotMeans[,11];
par(lwd=3, mar=c(3, 3+0.4, 0.2, 0.2));
palette(qualitative_hcl(palette="Dark 3", n=5));
plot(ALL.N, plotMeans[,1], type="l", ylim=c(0, 0.1), col=1, log="x", xlab="", ylab="");
lines(ALL.N, plotMeans[,2], col=2, lty=2);
lines(ALL.N, plotMeans[,3], col=3, lty=3);
lines(ALL.N, plotMeans[,4], col=4, lty=4);
lines(ALL.N, plotMeans[,5], col=5, lty=5);

for(i in 1:5){
	arrows(x0=ALL.N, y0=plotMeans[,i] - plotConf[,i], y1=plotMeans[,i] + plotConf[,i], angle=90, code=3, length=0.02);
}

title(xlab="sample size", line=2);
title(ylab=expression(F(paste(phantom(i), widehat(c), phantom(i)))), line=2);

legend("topright", c(expression(textstyle(min)), expression(c1), expression(c2), expression(c3), expression(c4)), col=1:5, lwd=3, lty=1:5, seg.len=4, bg="#FFFFFFBB");
