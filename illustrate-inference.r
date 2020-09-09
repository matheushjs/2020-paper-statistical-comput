

data = rgamma(n=20, shape=2.3, scale=2);
likelihood = function(p) - sum(log(dgamma(data - p[3], shape=p[1], scale=p[2])));
result = optim(par=c(1.1, 1.1, min(data) - 1/length(data)), fn=likelihood);

print(result);
par = result$par;
x = seq(-1, 15, length=1000);
plot(x, dgamma(x - par[3], shape=par[1], scale=par[2]));
hist(data, add=T, prob=T);
