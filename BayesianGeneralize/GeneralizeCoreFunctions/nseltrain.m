function res=nseltrain(j,indices,X,Y,n,p,gammatrue,nburnin,niter,T,a,Q,tau,nu,omega,seed,stand,cutoff)
        test=(indices==j);
        train=~test;
        Xtest=X(test,:);
        Ytest=Y(test);
        ntest=length(Ytest);
        Xtrain=X(train,:);
        Ytrain=Y(train);
        ntrain=n-ntest;
        nop=floor(min(ntrain,p)/2);
        gamma=gibbsgamma(nburnin,niter,p,nop,Ytrain, Xtrain, T, a, Q, ntrain,tau,nu,omega,seed,true,stand,false);

        freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
        selectindx=freq>cutoff;
        
        nvs=sum(selectindx==1);
        fp=sum((gammatrue'==0) & (selectindx==1));
        tn=sum((gammatrue'==0) & (selectindx==0));
        fn=sum((gammatrue'==1) & (selectindx==0));
        tp=sum((gammatrue'==1) & (selectindx==1));
        FPR= fp/(fp+tn);
        TPR= tp/(tp+fn);

        temp=1:p;
        xindx=temp(selectindx);
        %prediction conditional on selection
        Xpre=Xtrain(:,xindx);
        Tpre=T(:,xindx);
        Apre=Xpre'*Xpre+tau^(-2)*(Tpre'*Tpre);
        invApre=Apre\eye(length(xindx));
        tembeta=invApre*Xpre'*Ytrain;
        %betapre=stand.Sy*diag(1./stand.Sx(xindx))*tembeta;
        Ypre=stand.muy+stand.Sy*Xtest(:,xindx)*tembeta;
        mse=1/ntest*(Ytest.*stand.Sy+stand.muy-Ypre)'*(Ytest.*stand.Sy+stand.muy-Ypre);  
        res=[nvs,FPR,TPR,mse];
end