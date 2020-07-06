function [FPR,TPR,sensitivity,specificity,precision,accuracy,AUC]=rocauc(gammatrue,gamma,nburnin,niter,sn)


freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;

%%set up initial values
cutoff=0:sn:1;
sensitivity=zeros(1,1/sn+1);
specificity=zeros(1,1/sn+1);
precision=zeros(1,1/sn+1);
accuracy=zeros(1,1/sn+1);
FPR=zeros(1,1/sn+1);
TPR=zeros(1,1/sn+1);


for i=1:(1/sn+1)
    gammabin=freq>cutoff(i);

    FP = sum(gammatrue==0 & gammabin==1);
    FN = sum(gammatrue==1 & gammabin==0);
    TP = sum(gammatrue==1 & gammabin==1);
    TN = sum(gammatrue==0 & gammabin==0);

    sensitivity(i)= TP/(TP+FN);
    specificity(i)= TN/(TN+FP);
    precision(i)= TP/(TP+FP);
    accuracy(i)= (TP+TN)/(TP+TN+FP+FN);
    FPR(i)=FP/(FP+TN);
    TPR(i)=TP/(TP+FN);
end

AUC=trapz(FPR,TPR);

end

