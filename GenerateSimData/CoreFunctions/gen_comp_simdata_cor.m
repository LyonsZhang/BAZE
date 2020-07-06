function [X,beta,epsilon]=gen_comp_simdata_cor(n,gamma,b,theta,sigmaX,Xcor,sigma)

	%X is always N(0,sigmaX^2)
	%Y = bX + N(0,sigma^2)
	%n samples, p predictors.
    
    beta=gamma.*b;
%     theta=zeros(1,p);
%     theta(1:5)=log(0.5*p)*ones(1,5);
	X = mvnrnd(theta,diag(sigmaX)*Xcor*diag(sigmaX),n);
    epsilon = sigma*randn(n,1);
  	%Y=X*beta'+ epsilon;
end