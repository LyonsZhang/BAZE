function [X,Y,beta]=gen_simdata_ind(n,p,gamma,b,sigmaX,sigma)

	%X is always N(0,sigmaX^2)
	%Y = bX + N(0,sigma^2)
	%n samples, p predictors.
    
    beta=gamma.*b;
	X = sigmaX*randn(n,p);
    epsilon = sigma*randn(n,1);
  	Y=X*beta'+ epsilon;
end