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
mylty = c(1, "32");

graphics.off();
dev.new(width=12, height=8);

mydist = function(x){
	dgamma(x - 20, shape=15, scale=0.8);
}

x = seq(0, 70, length=2000);
plot(x[x>20], mydist(x[x>20]), type="l", lwd=2, col=mycolors[1], ylab="density", xlim=c(0, 70));
abline(v=20, col=1, lwd=2, lty="12");

# Now we find the closes gamma near that function
difference = function(params){
	sum(abs(dgamma(x, shape=params[1], scale=params[2]) - mydist(x)));
}

#params = optim(c(10, 0.8), difference, method="L-BFGS", lower=c(0, 0))$par;
#print(params);
params = c(94.4493562, 0.3364162); # What the above optimization finds

lines(x, dgamma(x, shape=params[1], scale=params[2]), lwd=5, lty="32", col=mycolors[2]);

legend("topleft", c("real dist.", "inferred dist."), col=mycolors[1:2], lty=1:2, lwd=3);
