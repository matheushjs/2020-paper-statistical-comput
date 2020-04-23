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

# Now we find the closest gamma near that function
difference = function(params){
	sum(abs(dgamma(x, shape=params[1], scale=params[2]) - mydist(x)));
}

#params = optim(c(1.2, 5), difference, method="L-BFGS", lower=c(0, 0))$par;
#print(params);
params = c(8.243800, 5.647846); # What the above optimization finds

#data = rgamma(n=20, shape=params[1], scale=params[2]);
#print(data);
data = c(45.42578,37.14454,31.91805,78.38157,51.13203,76.95387,34.06098,34.73172,26.69495,67.43486,98.01783,61.88397,41.99400,34.95689,51.49814,56.23078,55.54928,44.69236,36.25360,59.33179)

# Histogram
hist(data, freq=FALSE, border=F, col="#A0A0A070", xlab="x", ylab="density", xlim=c(0, 100), ylim=c(0, 0.035));

# Inferred dist
lines(x[x>6*pi], mydist(x[x>6*pi]), type="l", lwd=mylwd[2], col=mycolors[2], lty=mylty(2));
abline(v=6*pi, col=1, lwd=2, lty="12");

# Real dist
lines(x, dgamma(x, shape=params[1], scale=params[2]), lwd=mylwd[1], lty=mylty(1), col=mycolors[1]);

legend("topright", c("real dist.", "inferred dist.", "experimental data"), col=c(mycolors[1:2], "#A0A0A070"), lty=c(1:2, 1), lwd=c(3, 3, 15), box.lwd=0);

text(18, 0.033, "shifted\norigin", pos=4);
arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("formulation-fig1.png");
