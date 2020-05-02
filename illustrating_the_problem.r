require(viridis);
require(magick);

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

mylwd = c(2, 5, 4);
mylty = c(1, 2, 3);

graphics.off();
scale=0.6;
dev.new(width=scale*12, height=scale*8);

realParams = c(8.243800, 5.647846);

mydist = function(x, shape, scale){
	dgamma(x - 6*pi, shape=shape, scale=scale);
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
lines(x[x>6*pi], mydist(x[x>6*pi], params[1], params[2]), type="l", lwd=mylwd[2], col=mycolors[2], lty=mylty[2]);
abline(v=6*pi, col=1, lwd=2, lty="12");

# Real dist
lines(x[x > 6*pi], dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[3], lty=mylty[3], col=mycolors[3]);
lines(x, dgamma(x, shape=realParams[1], scale=realParams[2]), lwd=mylwd[1], lty=mylty[1], col=mycolors[1]);
#lines(x[x > 6*pi], 1.15*dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[3], lty=mylty[3], col=mycolors[1]);
#lines(x, dgamma(x, shape=realParams[1], scale=realParams[2]), lwd=mylwd[1], lty=mylty[1], col=mycolors[3]);

legend("topright", c("real dist.", "inferred dist.", "trunc. real dist.", "experimental data"), col=c(mycolors[1:3], "#A0A0A070"), lty=c(mylty[1:3], 1), lwd=c(mylwd[1:3]/1.8, 15), box.lwd=0);

text(18, 0.033, "shifted\norigin", pos=4);
arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("formulation-fig1.png");
system("convert formulation-fig1.png -crop 624x316+40+72 formulation-fig1.png");

img = image_read("formulation-fig1.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);



dev.new(width=scale*12, height=scale*8);

hist(data, freq=FALSE, border=F, col="#A0A0A070", xlim=c(0, 100), ylim=c(0, 0.035), axes=FALSE, xlab="", ylab="");
axis(1, pos=0);
axis(2, pos=0);
title(xlab="x", ylab="density", line=1);

lines(x[x>6*pi], mydist(x[x>6*pi], params[1], params[2]), type="l", lwd=mylwd[2], col=mycolors[2], lty=mylty[2]);
abline(v=6*pi, col=1, lwd=2, lty="12");

lines(x[x > 6*pi], dgamma(x[x > 6*pi], shape=realParams[1], scale=realParams[2]) / pgamma(6*pi, shape=realParams[1], scale=realParams[2], lower.tail=FALSE), lwd=mylwd[1], lty=mylty[1], col=mycolors[1]);

legend("topright", c("trunc. real dist.", "inferred dist.", "experimental data"), col=c(mycolors[1:2], "#A0A0A070"), lty=mylty[1:2], lwd=c(mylwd[1:2]/1.8, 15), box.lwd=0);
text(18, 0.033, "shifted\norigin", pos=4);
arrows(19, 0.031, 35, 0.031, length=0.1);

savePlot("formulation-fig2.png");
system("convert formulation-fig2.png -crop 624x316+40+72 formulation-fig2.png");

img = image_read("formulation-fig2.png");
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

savePlot("formulation-fig3.png");
system("convert formulation-fig3.png -crop 624x316+40+72 formulation-fig3.png");

img = image_read("formulation-fig3.png");
img = image_convert(img, type="grayscale");
dev.new();
plot(img);


