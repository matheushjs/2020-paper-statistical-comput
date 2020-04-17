require(viridis);

cols = viridis(100);

graphics.off();
dev.new(width=8, height=8);

mydist = function(x){
	dgamma(x - 30, shape=15, scale=0.8);
}

x = seq(0, 100, length=2000);
plot(x, mydist(x), type="l", lwd=3, col=2, ylab="density");
abline(v=30, col=1, lwd=3, lty=2);

# Now we find the closes gamma near that function
difference = function(params){
	sum(abs(dgamma(x, shape=params[1], scale=params[2]) - mydist(x)));
}

#params = optim(c(15, 0.8), difference, method="L-BFGS", lower=c(0, 0))$par;
#print(params);
params = c(147.4261944, 0.2829535); # What the above optimization finds

lines(x, dgamma(x, shape=params[1], scale=params[2]), lwd=3, lty=3, col=4);

