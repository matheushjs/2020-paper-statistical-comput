require(viridis);
require(colorspace);
require(magick);

mycolors = c(
	"#000000",
	"red",
	"blue",
	"#FF00FFFF",
	"#FFFF00FF",
	"#999999FF"
);

mylwd = c(3, 5, 4);
mylty = c(1, 2, 3);

graphics.off();
scale=0.6;
dev.new(width=scale*12, height=scale*8);

realParams = c(8.243800, 5.647846);

mydist = function(x, shape, scale, c=6*pi){
	dgamma(x - c, shape=shape, scale=scale);
}

x = seq(0, 400, length=5000);

# Now we find the closest gamma near that function
difference = function(params){
	slice = x[x > 6.0*pi];
	#print(log(mydist(slice, shape=params[1], scale=params[2])) * dgamma(slice, shape=realParams[1], scale=realParams[2]));
	#print(params);
	#print(-sum(log(mydist(slice, shape=params[1], scale=params[2])) * dgamma(slice, shape=realParams[1], scale=realParams[2])));
	-sum(log(mydist(slice, shape=params[1], scale=params[2])) * 1.15*dgamma(slice, shape=realParams[1], scale=realParams[2]));
}

#params = optim(c(3, 9), difference, method="L-BFGS", lower=c(1e-3, 1e-3))$par;
#print(params);
params = c(2.590555, 10.892933); # What the above optimization finds

#data = rgamma(n=10, shape=realParams[1], scale=realParams[2]);
#print(data);
data = c(45.42578,37.14454,31.91805,78.38157,51.13203,76.95387,34.06098,34.73172,26.69495,67.43486,98.01783,61.88397,41.99400,34.95689,51.49814,56.23078,55.54928,44.69236,36.25360,59.33179)

# Histogram
hist(data, freq=FALSE, border=F, col="#A0A0A070", xlab="", ylab="", xlim=c(0, 100), ylim=c(0, 0.035), axes=FALSE);
axis(1, pos=0);
axis(2, pos=0);
title(xlab="x", ylab="density", line=1);

# Inferred dist
lines(x[x>6*pi], dgamma(x[x>6*pi], shape=realParams[1], scale=realParams[2]) * pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower=FALSE)**-1 * 1.05, type="l", lwd=mylwd[2], col=mycolors[2], lty=mylty[2]);
abline(v=6*pi, col=1, lwd=2, lty="12");

# Real dist
#lines(x[x > 6*pi], dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[3], lty=mylty[3], col=mycolors[3]);
lines(x, dgamma(x, shape=realParams[1], scale=realParams[2]), lwd=mylwd[1], lty=mylty[1], col=mycolors[1]);
#lines(x[x > 6*pi], 1.15*dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[3], lty=mylty[3], col=mycolors[1]);
#lines(x, dgamma(x, shape=realParams[1], scale=realParams[2]), lwd=mylwd[1], lty=mylty[1], col=mycolors[3]);

#legend("topright", c("real dist.", "inferred dist.", "trunc. real dist.", "experimental data"), col=c(mycolors[1:3], "#A0A0A070"), lty=c(mylty[1:3], 1), lwd=c(mylwd[1:3]/1.8, 15), box.lwd=0);
legend("topright", c("real dist.", "truncated dist.", "experimental data"), col=c(mycolors[1:2], "#A0A0A070"), lty=c(mylty[1:2], 1), lwd=c(mylwd[1:2]/1.8, 15), box.lwd=0);

text(18, 0.033, expression(group("[", list(c, infinity), ")")), pos=4);
arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("fig1.png");
system("convert fig1.png -crop 624x316+40+72 fig1.png");

img = image_read("fig1.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);

dev.new(width=scale*12, height=scale*8);

hist(data, freq=FALSE, border=F, col="#A0A0A070", xlim=c(0, 100), ylim=c(0, 0.035), axes=FALSE, xlab="", ylab="");
axis(1, pos=0);
axis(2, pos=0);
title(xlab="x", ylab="density", line=1);

lines(x[x>6*pi], mydist(x[x>6*pi], params[1], params[2]), type="l", lwd=mylwd[1], col=mycolors[3], lty=mylty[1]);
abline(v=6*pi, col=1, lwd=2, lty="12");

lines(x[x > 6*pi], dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[2], lty=mylty[2], col=mycolors[2]);

legend("topright", c("trunc. real dist.", "inferred dist.", "experimental data"), col=c(mycolors[2:3], "#A0A0A070"), lty=mylty[2:1], lwd=c(mylwd[2:1]/1.8, 15), box.lwd=0);
#text(18, 0.033, "shifted\norigin", pos=4);
#arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("fig2.png");
system("convert fig2.png -crop 624x316+40+72 fig2.png");

img = image_read("fig2.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);



dev.new(width=scale*12, height=scale*8);

hist(data, freq=FALSE, border=F, col="#A0A0A070", xlim=c(0, 100), ylim=c(0, 0.045), axes=FALSE, xlab="", ylab="");
axis(1, pos=0);
axis(2, pos=0);
title(xlab="x", ylab="density", line=1);

lines(x, dgamma(x, shape=realParams[1], scale=realParams[2]), lwd=mylwd[1], lty=mylty[1], col=mycolors[1]);

n = length(data);
CV = sd(data) / mean(data);
est1 = min(data);
est2 = min(data) - abs(min(data)) * (CV / log10(n));
est3 = min(data) - abs(min(data)) * (1 / n);
est4 = min(data) - abs(min(data)) * 1 * sqrt(log(log(n)) / (2*n));
est5 = min(data) - abs(min(data)) * 1 * sqrt(-log(0.05/2) / (2*n));

y = 0.02
lwd = 4
cols = viridis(6)[1:5];
mycex = 2;
segments(est1, 0, y1=y, col=cols[1], lwd=lwd, lty="12");
points(est1, y, cex=mycex, col=cols[1], pch=19);
#text(est1, y+0.002, "est1");

y = y + 0.003
segments(est3, 0, y1=y, col=cols[2], lwd=lwd, lty="12");
points(est3, y, cex=mycex, col=cols[2], pch=19);
#text(est3, y+0.002, "est3");

y = y + 0.003
segments(est4, 0, y1=y, col=cols[3], lwd=lwd, lty="12");
points(est4, y, cex=mycex, col=cols[3], pch=19);
#text(est4, y+0.002, "est4");

y = y + 0.003
segments(est2, 0, y1=y, col=cols[4], lwd=lwd, lty="12");
points(est2, y, cex=mycex, col=cols[4], pch=19);
#text(est2, y+0.002, "est2");

y = y + 0.003
segments(est5, 0, y1=y, col=cols[5], lwd=lwd, lty="12");
points(est5, y, cex=mycex, col=cols[5], pch=19);
#text(est5, y+0.002, "est5");

legends=list(
	"sample\nmin",
	expression(over(CV, log[10](n))),
	expression(over(1,n)),
	expression(sqrt(over(log(log(n)), 2*n))),
	expression(sqrt(over(-log(nu/2), 2*n)))
);
legends = legends[c(1,3,4,2,5)];

myx = rev(c(0.5, 25, 45, 70, 80)) + 5;
for(i in 1:5){
	legend(myx[i], 0.039, yjust=0, legend=legends[[i]], pch=19, col=cols[i], box.lwd=0, pt.cex=2);
}

savePlot("fig3.png");
system("convert fig3.png -crop 624x316+40+72 fig3.png");

img = image_read("fig3.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);




dev.new(width=scale*12, height=scale*8);

plot(c(-1,-1), xlim=c(0, 100), ylim=c(0, 0.055), axes=FALSE, xlab="", ylab="");
axis(1, pos=0);
axis(2, pos=0);
title(xlab="x", ylab="density", line=1);

plotFuncDiff = function(x, y1, y2, ...){
	forwardPart  = cbind(x, y1);
	backwardPart = cbind(rev(x), rev(y2));
	all = rbind(forwardPart, backwardPart, c(x[1], y1[1]));
	polygon(all, ...);
}

#lines(x[x>6*pi], mydist(x[x>6*pi], params[1], params[2]), type="l", lwd=mylwd[2], col=mycolors[2], lty=mylty[2]);
#abline(v=6*pi, col=1, lwd=2, lty="12");

#lines(x[x > 6*pi], dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[1], lty=mylty[1], col=mycolors[1]);

f = function(c, raiseY=0, ...){
	difference = function(params){
		slice = x[x > c];
		-sum(log(mydist(slice, shape=params[1], scale=params[2], c=c)) * dgamma(slice, shape=realParams[1], scale=realParams[2]));
	}

	mockData = rgamma(n=2000, shape=realParams[1], scale=realParams[2]) - c;
	initShape = 2;
	initScale = mean(mockData[mockData > 0]) / initShape;

	params = optim(c(initShape, initScale), difference, method="L-BFGS", lower=c(1e-3, 1e-3))$par;
	print(params);

	# abline(v=c, lty="13", lwd=2);
	plotFuncDiff(
		x[x > c],
		mydist(x[x > c], params[1], params[2], c=c) + raiseY,
		dgamma(x[x > c], shape=realParams[1], scale=realParams[2]) / pgamma(c, shape=realParams[1], scale=realParams[2], lower.tail=FALSE) + raiseY,
		...
	);

	segments(c, 0, c, raiseY, lwd=3, lty="12", col="#00000077");

	# Estimate area under curve
	y1 = mydist(x[x > c], params[1], params[2], c=c);
	y2 = dgamma(x[x > c], shape=realParams[1], scale=realParams[2]) / pgamma(c, shape=realParams[1], scale=realParams[2], lower.tail=FALSE);
	diffs = abs(y1 - y2);
	diffs = diffs[2:length(diffs)];
	
	sum(diff(x[x > c]) * diffs);
}

allC = c(25, 20, 15, 10, 7.5, 5);
allCols = viridis(length(allC)+1, alpha=0.8)[1:length(allC)];
areas = seq_along(allC);
for(i in 1:length(allC)){
	areas[i] = f(allC[i], 0.004*(length(allC) - i), col=allCols[i]);
}

legend("topright", paste("c = ", allC, "	area = ", round(areas, digits=3)), col=allCols, pch=15, pt.cex=2.5, box.lwd=0);
#text(18, 0.033, "shifted\norigin", pos=4);
#arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("fig4.png");
system("convert fig4.png -crop 624x316+40+72 fig4.png");

img = image_read("fig4.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);




dev.new(width=scale*12, height=scale*8);

#hist(data, freq=FALSE, border=F, col="#A0A0A070", xlim=c(0, 100), ylim=c(0, 0.045), axes=FALSE, xlab="", ylab="");

x = seq(0, 45, length=2000);
plot(x, dgamma(x, shape=160, scale=0.2), lwd=mylwd[1], lty=mylty[1], col=mycolors[1], axes=FALSE, xlab="", ylab="", type="l", xlim=c(20, 45));
axis(1, pos=0);
axis(2, pos=20);
title(xlab="x", ylab="density", line=1);

q = 0.05
# I will just recycle the names here... sorry I'm lazy
est1 = qgamma(1 - (1 - q)**(1/10), shape=160, scale=0.2);
est2 = qgamma(1 - (1 - q)**(1/20), shape=160, scale=0.2);
est3 = qgamma(1 - (1 - q)**(1/50), shape=160, scale=0.2);
est4 = qgamma(1 - (1 - q)**(1/100), shape=160, scale=0.2);
est5 = qgamma(1 - (1 - q)**(1/200), shape=160, scale=0.2);

y = 0.05
lwd = 4
cols = viridis(6)[1:5];
mycex = 2;
segments(est1, 0, y1=y, col=cols[1], lwd=lwd, lty="12");
points(est1, y, cex=mycex, col=cols[1], pch=19);
#text(est1, y+0.02, "10");

y = y + 0.008
segments(est2, 0, y1=y, col=cols[2], lwd=lwd, lty="12");
points(est2, y, cex=mycex, col=cols[2], pch=19);
#text(est2, y+0.02, "20");

y = y + 0.008
segments(est3, 0, y1=y, col=cols[3], lwd=lwd, lty="12");
points(est3, y, cex=mycex, col=cols[3], pch=19);
#text(est3, y+0.02, "50");

y = y + 0.008
segments(est4, 0, y1=y, col=cols[4], lwd=lwd, lty="12");
points(est4, y, cex=mycex, col=cols[4], pch=19);
#text(est4, y+0.02, "100");

y = y + 0.008
segments(est5, 0, y1=y, col=cols[5], lwd=lwd, lty="12");
points(est5, y, cex=mycex, col=cols[5], pch=19);
text(est5-1, y+0.01, "5%");

q = 0.01
# I will just recycle the names here... sorry I'm lazy
est1 = qgamma(1 - (1 - q)**(1/10), shape=160, scale=0.2);
est2 = qgamma(1 - (1 - q)**(1/20), shape=160, scale=0.2);
est3 = qgamma(1 - (1 - q)**(1/50), shape=160, scale=0.2);
est4 = qgamma(1 - (1 - q)**(1/100), shape=160, scale=0.2);
est5 = qgamma(1 - (1 - q)**(1/200), shape=160, scale=0.2);

y = 0.02
lwd = 4
cols = viridis(6)[1:5];
mycex = 2;
segments(est1, 0, y1=y, col=cols[1], lwd=lwd, lty="12");
points(est1, y, cex=mycex, col=cols[1], pch=19);
#text(est1, y+0.02, "10");

y = y + 0.008
segments(est2, 0, y1=y, col=cols[2], lwd=lwd, lty="12");
points(est2, y, cex=mycex, col=cols[2], pch=19);
#text(est2, y+0.02, "20");

y = y + 0.008
segments(est3, 0, y1=y, col=cols[3], lwd=lwd, lty="12");
points(est3, y, cex=mycex, col=cols[3], pch=19);
#text(est3, y+0.02, "50");

y = y + 0.008
segments(est4, 0, y1=y, col=cols[4], lwd=lwd, lty="12");
points(est4, y, cex=mycex, col=cols[4], pch=19);
#text(est4, y+0.02, "100");

y = y + 0.008
segments(est5, 0, y1=y, col=cols[5], lwd=lwd, lty="12");
points(est5, y, cex=mycex, col=cols[5], pch=19);
text(est5-1, y+0.01, "1%");

legend("topright", c("200", "100", "50", "20", "10"), title="sample size", yjust=0, pch=19, col=rev(cols), box.lwd=0, pt.cex=2);

savePlot("fig5.png");
system("convert fig5.png -crop 624x316+40+72 fig5.png");

img = image_read("fig5.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);
