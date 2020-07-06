function [X,Y,beta]=gen_simdata_cor(n,p,gamma,b,sigmaX,Xcor,sigma)

	%X is always N(0,sigmaX^2)
	%Y = bX + N(0,sigma^2)
	%n samples, p predictors.
    
    beta=gamma.*b;
	X = mvnrnd(zeros(1,p),diag(sigmaX)*Xcor*diag(sigmaX),n);
    epsilon = sigma*randn(n,1);
  	Y=X*beta'+ epsilon;
end