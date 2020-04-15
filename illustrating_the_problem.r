
mydist = function(x){
	dgamma(x - 10, shape=5, scale=5);
}

x = seq(0, 100, length=2000);
plot(x, mydist(x), type="l");

# Now we find the closes gamma near that function
difference = function(params){
	sum(abs(dgamma(x, shape=params[1], scale=params[2]) - mydist(x)));
}

params = optim(c(5, 5), difference, method="L-BFGS", lower=c(0, 0))$par;
print(params);

lines(x, dgamma(x, shape=params[1], scale=params[2]));
