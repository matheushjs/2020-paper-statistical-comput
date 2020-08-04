require(elfDistr);
require(viridis);
require(magick);
require(stringr);
require(colorspace);

TRIALS = 10;

graphics.off();

# We make a ton of histograms and then later we select which of them we want

cols = qualitative_hcl("Set 2", n=TRIALS); # inferno(TRIALS+1)[1:TRIALS];
for(i in 1:TRIALS){
	scale = 0.4;
	dev.new(width=scale*12, height=scale*6);

	alpha = runif(n=1, 0.1, 1);
	beta  = sample(c(runif(n=1, 0.1, 1), runif(n=1, 1, 20)))[1];
	gamma = sample(c(runif(n=1, 0.1, 1), runif(n=1, 1, 20)))[1];
	a     = sample(c(runif(n=1, 0.1, 1), runif(n=1, 1, 20)))[1];
	b     = sample(c(runif(n=1, 0.2, 1), runif(n=1, 1, 20)))[1];

	print(c(alpha, beta, gamma, a, b));
	data = rkwcwg(n=100, alpha, beta, gamma, a, b);
	#par(bg="cyan");
	hist(data, border=FALSE, col=cols[i], axes=F);
	axis(1, pos=0, mgp=c(0,0.5,0));

	filename = paste("hist-", i, ".png", sep="");
	savePlot(filename);
	system(str_replace_all("convert ## -crop 78%x39%+75%+74% ##", "##", filename));

	img = image_read(filename);
	img = image_convert(img, type="grayscale");
	dev.new();
	plot(img);
}
