require(viridis);

mycolors = c(
	"#000000FF",
	"#FF0000FF",
	3,
	"#0000FFFF",
	"#00FFFFFF",
	"#FF00FFFF",
	"#FFFF00FF",
	"#999999FF"
);

mylwd = c(2, 5);
mylty = function(idx){
	l = list(1, "32");
	l[[idx]];
}

graphics.off();
scale=0.6;
dev.new(width=scale*12, height=scale*8);

mydist = function(x){
	dgamma(x - 6*pi, shape=3, scale=10);
}

x = seq(0, 100, length=2000);
plot(x, mydist(x), type="l", lwd=mylwd[2], col=mycolors[2], xlab="x", ylab="density", xlim=c(0, 100), ylim=c(0, 0.035), lty=mylty(2));
abline(v=6*pi, col=1, lwd=2, lty="12");

# Now we find the closest gamma near that function
difference = function(params){
	sum(abs(dgamma(x, shape=params[1], scale=params[2]) - mydist(x)));
}

#params = optim(c(1.2, 5), difference, method="L-BFGS", lower=c(0, 0))$par;
#print(params);
params = c(8.243800, 5.647846); # What the above optimization finds

lines(x, dgamma(x, shape=params[1], scale=params[2]), lwd=mylwd[1], lty=mylty(1), col=mycolors[1]);

legend("topright", c("real dist.", "inferred dist."), col=mycolors[1:2], lty=1:2, lwd=3);

text(17.8, 0.10, "c = 20");
arrows(18, 0.096, 19.9, 0.090, length=0.1);
