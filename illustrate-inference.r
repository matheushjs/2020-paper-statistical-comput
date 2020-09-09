
if(FALSE){
	N = 20;
	data = rgamma(n=N, shape=2.3, scale=2);
	c.init = min(data) - sd(data)/N;
	likelihood = function(p) -sum(log(dgamma(data - p[3], shape=p[1], scale=p[2])));
	result = optim(
		par=c(1, 1, c.init), fn=likelihood);

	print(result);
	par = result$par;
	x = seq(-1, 15, length=1000);
	plot(x, dgamma(x - par[3], shape=par[1], scale=par[2]));
	hist(data, add=T, prob=T);
} else {
	N = 20;
	data = rgamma(n=N, shape=2.3, scale=2);
	likelihood = function(p){
		q = qgamma(1 - (1 - 0.5)^(1/N), shape=p[1], scale=p[2]);
		-sum(log(dgamma(data - q, shape=p[1], scale=p[2])));
	}
	result = optim(
		par=c(1, 1), fn=likelihood);

	print(result);
	par = result$par;
	x = seq(-1, 15, length=1000);
	plot(x, dgamma(x - qgamma(1 - (1 - 0.5)**(1/N), shape=par[1], scale=par[2]), shape=par[1], scale=par[2]));
	hist(data, add=T, prob=T);
}
