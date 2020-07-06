function [MSEpre,MSEtest,etafinal]=traintestreal0(nburnin,niter,p,Y, X,T, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand,pl,pathname,figname)

%%%%do not print iterations
display=false;
%plot=false;

%indices=crossvalind('Kfold',n,k);
[train,test]=crossvalind('HoldOut',n,0.25);

Xtest=X(test,:);
Ytest=Y(test);
ntest=length(Ytest);
Xtrain=X(train,:);
Ytrain=Y(train);
ntrain=n-ntest;

%standadize data
[Xtrain,stand1.mux,stand1.Sx] = zscore(log(Xtrain));
[Ytrain,stand1.muy,stand1.Sy] = zscore(Ytrain);
[Xtest,stand2.mux,stand2.Sx] = zscore(log(Xtest));
[Ytest,stand2.muy,stand2.Sy] = zscore(Ytest);

%specify penalty matrix
c=100;
T1=[eye(p);c*ones(1,p)]*diag(1./stand1.Sx);

% MSEtest=0;
% MSEpre=0;
%MSEpre=zeros(1,niter);

% etafinal=zeros(1,p);


nop=floor(min(ntrain,p)/2);
trainsignal=Ytrain'*Ytrain/ntrain;
testsignal=Ytest'*Ytest/ntest;

[gamma,betahat,MSE]=gibbsgamma(nburnin,niter,p,nop,Ytrain, Xtrain, T, a, Q, ntrain,tau,nu,omega,seed,predict,stand,display);
seq=(nburnin+1):(nburnin+niter);
MSEtemp=arrayfun(@(x) 1/(ntest)*(Ytest.*stand.Sy-(Xtest*diag(stand.Sx))*betahat(:,x))'*(Ytest.*stand.Sy-(Xtest*diag(stand.Sx))*betahat(:,x)),seq);
MSEtest=mean(MSEtemp); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%plot selected gamma%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
top=freq(freq>cutoff);
temp=1:p;
xindx=temp(freq>cutoff);

if pl
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
figfile1 = fullfile(pathname, ['selected' figname '_trainsignal' num2str(trainsignal) 'testsignal' num2str(testsignal) '.eps']);
figfile2 = fullfile(pathname, ['selected' figname '_trainsignal' num2str(trainsignal) 'testsignal' num2str(testsignal) '.pdf']);
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
figfile1 = fullfile(pathname, ['MSE' figname '_trainsignal' num2str(trainsignal) 'testsignal' num2str(testsignal) '.eps']);
figfile2 = fullfile(pathname, ['MSE' figname '_trainsignal' num2str(trainsignal) 'testsignal' num2str(testsignal) '.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)
end

%prediction conditional on selection
Xpre=Xtrain(:,xindx);
Tpre=T1(:,xindx);
Apre=Xpre'*Xpre+tau^(-2)*(Tpre'*Tpre);
invApre=Apre\eye(length(xindx));
tembeta=invApre*Xpre'*Ytrain;
%betapre=stand.Sy*diag(1./stand.Sx(xindx))*tembeta;

Ypre=stand2.muy+stand2.Sy*Xtest(:,xindx)*tembeta;
MSEpre=1/ntest*(Ytest.*stand2.Sy+stand2.muy-Ypre)'*(Ytest.*stand2.Sy+stand2.muy-Ypre);  

%coefficient distance
betafinal=zeros(1,p);
Xri=Xtrain(:,xindx);
Tri=T1(:,xindx);
Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
invAri=Ari\eye(size(Xri,2));
betafinal(xindx)=stand1.Sy*diag(1./stand1.Sx(xindx))*invAri*Xtrain(:,xindx)'*Ytrain;
etafinal=betafinal;
%coef_l2(i)=sqrt(sum(abs(etafinal'-beta).^2));


% average=mean(MSEtest);
% sem=std(MSEtest)/sqrt(k);

end