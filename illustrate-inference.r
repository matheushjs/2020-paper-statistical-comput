
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
