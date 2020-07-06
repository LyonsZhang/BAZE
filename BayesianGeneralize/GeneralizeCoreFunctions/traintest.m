function [MSEpre,MSEtest,average,sem,FPR,FNR,FPN,FNN,sensitivity,specificity,precision,accuracy,etafinal,roc]=traintest(k,nburnin,niter,p,Y, X,T, beta, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand,pl,true_index,pathname,figname,doroc)

%%%%do not print iterations
display=false;
%plot=false;

if iscolumn(beta)
    beta=beta';
end

indices=crossvalind('Kfold',n,k);
MSEtest=zeros(1,k);
MSEpre=zeros(1,k);
%MSEpre=zeros(1,niter);
FPR= zeros(1,k);
FNR=zeros(1,k);
FPN=zeros(1,k);
FNN=zeros(1,k);
sensitivity=zeros(1,k);
specificity=zeros(1,k);
precision=zeros(1,k);
accuracy=zeros(1,k);
etafinal=zeros(k,p);
if doroc
    sn=0.05;
    rFPR=zeros(k,1/sn+1);
    rTPR=zeros(k,1/sn+1);
    rsensitivity=zeros(k,1/sn+1);
    rspecificity=zeros(k,1/sn+1);
    rprecision=zeros(k,1/sn+1);
    raccuracy=zeros(k,1/sn+1);
    rAUC=zeros(1,k);
end

parfor i=1:k
    test=(indices==i);
    train=~test;
    Xtest=X(test,:);
    Ytest=Y(test);
    ntest=length(Ytest);
    Xtrain=X(train,:);
    Ytrain=Y(train);
    ntrain=n-ntest;
    nop=floor(min(ntrain,p)/2);
    signal=round(var(Xtrain*beta'),2);
    noise=round(var(Ytrain-Xtrain*beta'),2);
    [gamma,betahat,MSE]=gibbsgamma(nburnin,niter,p,nop,Ytrain, Xtrain, T, a, Q, ntrain,tau,nu,omega,seed,predict,stand,display);
    seq=(nburnin+1):(nburnin+niter);
    MSEtemp=arrayfun(@(x) 1/(ntest)*(Ytest.*stand.Sy-(Xtest*diag(stand.Sx))*betahat(:,x))'*(Ytest.*stand.Sy-(Xtest*diag(stand.Sx))*betahat(:,x)),seq);
    MSEtest(i)=mean(MSEtemp); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%plot selected gamma%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
    top=freq(freq>cutoff);
    temp=1:p;
    xindx=temp(freq>cutoff);
    
    if pl
    plot(true_index,0.1*ones(size(true_index)),'s','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.800000 0.500000 0.200000])
    hold on
    plot(temp,freq,'ko')
    hold on
    stem(xindx,top,'-r')
    hold off
    xlabel('Variable')
    ylabel('Frequency')
    % Set the figure properties
    fig = figure(1);
    fig.Resize = 'off';
    fig.PaperUnits = 'inches';
    fig.Units = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0, 0, 10, 5];
    fig.PaperSize = [10, 5];
    fig.Position = [0.1, 0.1, 9.9, 4.9];
    fig.InvertHardcopy = 'off';

    ax = gca;
    ax.TickLabelInterpreter = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.FontName = 'LaTeX';
    ax.Title.Interpreter = 'LaTeX';
    ax.XLabel.Interpreter = 'LaTeX';
    ax.YLabel.Interpreter = 'LaTeX';
    ax.XTick = 0:100:1000;
    ax.Box = 'off';
    ax.LineWidth = 1.5;
    ax.FontSize = 16;

    %use that when you save
    figfile1 = fullfile(pathname, ['selected' num2str(i) figname '_signal' num2str(signal) 'noise' num2str(noise) '.eps']);
    figfile2 = fullfile(pathname, ['selected' num2str(i) figname '_signal' num2str(signal) 'noise' num2str(noise) '.pdf']);
    saveas(gcf,figfile1,'epsc')
    saveas(gcf,figfile2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%plot MSE convergence%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(0:(nburnin+niter),MSE,'-o','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5)
    xlabel('Iteration')
    ylabel('MSE')

    % Set the figure properties
    fig = figure(1);
    fig.Resize = 'off';
    fig.PaperUnits = 'inches';
    fig.Units = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0, 0, 10, 5];
    fig.PaperSize = [10, 5];
    fig.Position = [0.1, 0.1, 9.9, 4.9];
    fig.InvertHardcopy = 'off';

    ax = gca;
    ax.TickLabelInterpreter = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.FontName = 'LaTeX';
    ax.Title.Interpreter = 'LaTeX';
    ax.XLabel.Interpreter = 'LaTeX';
    ax.YLabel.Interpreter = 'LaTeX';
    ax.XTick = 0:1000:(nburnin+niter);
    ax.Box = 'off';
    ax.LineWidth = 1.5;
    ax.FontSize = 16;
    %use that when you save
    figfile1 = fullfile(pathname, ['MSE' num2str(i) figname '_signal' num2str(signal) 'noise' num2str(noise) '.eps']);
    figfile2 = fullfile(pathname, ['MSE' num2str(i) figname '_signal' num2str(signal) 'noise' num2str(noise) '.pdf']);
    saveas(gcf,figfile1,'epsc')
    saveas(gcf,figfile2)
    end
    
    %prediction conditional on selection
    Xpre=Xtrain(:,xindx);
    Tpre=T(:,xindx);
    Apre=Xpre'*Xpre+tau^(-2)*(Tpre'*Tpre);
    invApre=Apre\eye(length(xindx));
    tembeta=invApre*Xpre'*Ytrain;
    %betapre=stand.Sy*diag(1./stand.Sx(xindx))*tembeta;
    Ypre=stand.muy+stand.Sy*Xtest(:,xindx)*tembeta;
    MSEpre(i)=1/ntest*(Ytest.*stand.Sy+stand.muy-Ypre)'*(Ytest.*stand.Sy+stand.muy-Ypre);  
    
    %selection accuracy
    trueindx=(beta~=0);
    selectindx=(freq>cutoff);
    FP=sum((trueindx==0) & (selectindx==1));
    TN=sum((trueindx==0) & (selectindx==0));
    FN=sum((trueindx==1) & (selectindx==0));
    TP=sum((trueindx==1) & (selectindx==1));
    FPR(i)= FP/(FP+TN);
    FNR(i)=FN/(TP+FN);
    FPN(i)=FP;
    FNN(i)=FN;
    sensitivity(i)=TP/(TP+FN);
    specificity(i)=TN/(TN+FP);
    precision(i)=TP/(TP+FP);
    accuracy(i)=(TP+TN)/(TP+TN+FP+FN);
    
    if doroc
        [rFPR(i,:),rTPR(i,:),rsensitivity(i,:),rspecificity(i,:),rprecision(i,:),raccuracy(i,:),rAUC(i)]=rocauc(trueindx,gamma,nburnin,niter,sn);       
    end
    
    %coefficient distance
    betafinal=zeros(1,p);
    Xri=X(:,xindx);
    Tri=T(:,xindx);
    Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
    invAri=Ari\eye(size(Xri,2));
    betafinal(xindx)=stand.Sy*diag(1./stand.Sx(xindx))*invAri*X(:,xindx)'*Y;
    etafinal(i,:)=betafinal;
    %coef_l2(i)=sqrt(sum(abs(etafinal'-beta).^2));
end

average=mean(MSEtest);
sem=std(MSEtest)/sqrt(k);

if doroc
    roc = struct('rFPR',rFPR,'rTPR',rTPR,'rsensitivity',rsensitivity,'rspecificity',rspecificity,'rprecision',rprecision,'raccuracy',raccuracy,'rAUC',rAUC);
end

end