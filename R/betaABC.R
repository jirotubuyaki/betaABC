library(abc);

sigma.hat <<- 0.3;
gmm_number <<- 1;
iteration_limit <<- 10000;
iteration <<- 1500;
u <<- array();
theta <<- matrix(0, nrow=iteration_limit, ncol=gmm_number);
theta_after <<- matrix(0, nrow=0, ncol=gmm_number);
diff <- 7000;
a <<- 1;
b <<- 1;

mu_1 <<- c(1,10,13,19,23);
sigma_1 <<- sqrt(c(1,9,5,8,3));
weight_1 <<- c(0.3,0.05,0.05,0.5,0.1);

mu_1 <- c(1,20);
sigma_1 <- sqrt(c(1,1));
weight_1 <- c(0.1,0.9);
data.obs <<- rmixnorm(iteration,mu_1,sigma_1,weight_1);
sortlist <<- order(data.obs);
data.obs <<- data.obs[sortlist];
print(data.obs);
mean_X <<- mean(data.obs);
max_data <<- data.obs[iteration];
min_data <<- data.obs[iteration];
max_X <<- 50;
min_X <<- -50;

x <- 1;
y <- 1;
z <- 1;
model <- function(iteration){
	while(x <= iteration_limit){
		u[x] <<- rbetanorm(a,b);
		theta_sample <- setPrior(1,u[x]);
		data <- setMixModel(theta_sample);
		sortlist <<- order(data);
		data.sample <<- data[sortlist];
		theta[x,] <<- c(theta_sample[1:gmm_number]);
		y <<- 1;
		sum_diff <<- 0;
		while(y <= iteration){
			sum_diff <<- sum_diff + abs(data.obs[y] - data.sample[y]);
			y <<- y + 1;
		}
		if(sum_diff <= diff){
			print("accept");
			print(sum_diff);
			print(x);
			theta_after <<- rbind(theta_after,theta[x,]);
			
			z <<- z + 1;
			if(u[x] > mean_X){
				a <<- a + abs((u[x]-mean_X) / abs(max_data -mean_X));
			}
			else{
				b <<- b  + abs((mean_X-u[x]) / abs(min_data -mean_X));
			}
		}
		else{
			if(u[x] > mean_X){
				b <<- b + abs((u[x]-mean_X) / abs(max_data -mean_X));
			}
			else{
				a <<- a  + abs((mean_X-u[x]) / abs(min_data -mean_X));
			}
		}
		x <<- x + 1;
	}

	#data <- sapply(theta, setNormData);
	#ss.sim <- dataCell(data,min_X,max_X);
	#return(list(theta = theta, ss = ss.sim));
}
rbetanorm <- function(shape1,shape2){
	x <- rbeta(1,shape1,shape2,ncp=0);
	x <- x * (max(data.obs) - min(data.obs)) + min(data.obs);
	return(x);
}
setPrior <- function(n,u){
	X <- rnorm(gmm_number, u, 4);
	return(X);
}
setMixModel <- function(mu){
  Y <-sample(length(1 / gmm_number), n, replace=TRUE, prob=(1 / gmm_number));
  X <-rnorm(iteration, mu[Y], sigma.hat);
  return(X);
}
dbetanorm <- function(x,shape1,shape2){
	dbeta(x,shape1,shape2,ncp=0,log =FALSE);
}
rmixnorm <- function(m,mu, sigma,weight){
  Y <-sample(length(weight), m, replace=TRUE, prob=weight);
  X <-rnorm(m, mu[Y], sigma[Y]);
  return(X);

}

# Simulate theta and  summary statistics.
sim <- model(iteration);


# Rejection algorithm (tolerance = 1%).
#abc.rej <- abc(target = ss.obs, param = matrix(sim$theta, dimnames = list(NULL, "theta")),sumstat = sim$ss, tol = 0.01, method = "rejection")
range.x <- seq(-10, 200, x <- 0.05)
oldpar <- par(mfrow = c(3,4));
hist(data.obs,breaks = 100,freq=TRUE,xlim=c(min(data.obs),max(data.obs)),main = "Observed Data",col="royalblue",border="royalblue");
hist(u, breaks = 100,freq = TRUE,xlim=c(min(data.obs),max(data.obs)),main = "Super Prior Distridution",col="royalblue",border="royalblue");
x <- 1;
while(x <= gmm_number){
	hist(theta[1:iteration_limit,x], breaks = 100,freq=TRUE,xlab = paste(expression(theta),x),main = "Prior Distridution",xlim=c(min(theta[,1]),max(theta[,1])),col="royalblue",border="royalblue");
	x <- x + 1;
}
x <- 1;
while(x <= gmm_number){
	hist(theta_after[1:z-1,x], breaks = 100,freq=TRUE,xlab = paste(expression(theta),x),main = "Posterior Distribution",xlim=c(min(theta[,1]),max(theta[,1])),col="royalblue",border="royalblue")
	x <- x + 1;
}
#oldpar <- par(mfrow = c(1,1));
#library(MASS);
#library(KernSmooth);
#x <- theta_after[1:z-1,1];
#y <- theta_after[1:z-1,2];
#f1 <- bkde2D(cbind(x, y), bandwidth=c(width.SJ(x), width.SJ(y)))
#persp(f1$fhat, phi = 30, theta = 20, d = 5)
