%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%load simulated data
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated\')
%%%set path
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
mkdir('results');
cd('results')

Z=csvread('train_otu.csv');
Y=csvread('train_bmi.csv');
Ztest=csvread('test_otu.csv');
Ytest=csvread('test_bmi.csv');

%standadize data
[X,stand.mux,stand.Sx] = zscore(log(Z));
[Y,stand.muy,stand.Sy] = zscore(Y);
[Xtest,stand1.mux,stand1.Sx] = zscore(log(Ztest));
[Ytest,stand1.muy,stand1.Sy] = zscore(Ytest);

[n,p]=size(X); %n is number of samples
 % p is number of predictors

%specify penalty matrix
c=100;
T=[eye(p);c*ones(1,p)]*diag(1./stand.Sx);

%%%read in phylogenetic tree similarity matrix
Q_full = csvread('sortcorr.csv');
save('Qorg.mat','Q_full')
%Q_full(Q_full<0.6)=0;
Q_full(Q_full>0.9)=0.9;
Q_full=Q_full+0.1*eye(p);
Q_full=Q_full\eye(p);
save('Qinv.mat','Q_full')
Q_full=Q_full-diag(diag(Q_full));
Q_full(Q_full>1)=1;
Q_full(Q_full<-1)=-1;
%Q_full=Q_full-Q_full;
sqrt(sum(sum(Q_full>0)))
Q=sparse(Q_full);
%Q=zeros(p);
histogram(Q)
sum(sum(Q,2))/p
HeatMap(Q_full)

% imagesc(Q_full)
% colormap(gray)
imagesc(abs(Q_full)>0.1)
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar



%set the number of burn-in steps and iterations
nburnin=20000;
niter=5000;

%initialize gamma and setsort
nop=floor(n/2);

a0=-7.8;
a=a0*ones(1,p);

tau=1;
nu=0;
omega=0;
seed=2020;

tic
[gamma,betahat,MSE,nselect,Yhat]=gibbsgamma(nburnin,niter,p,nop,Y, X,T, a, Q, n,tau,nu,omega,seed,true,stand,true);
toc

cutoff=0.5;

freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
top=freq(freq>cutoff);
temp=1:p;
xindx=temp(freq>cutoff);

betafinal=zeros(1,p);
Xri=X(:,xindx);
Tri=T(:,xindx);
Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
invAri=Ari\eye(size(Xri,2));
tembeta=invAri*X(:,xindx)'*Y;
betafinal(xindx)=stand.Sy.*(1./stand.Sx(xindx)).*tembeta';
Yest=stand.muy+stand.Sy*Xri*tembeta;
Yobs=Y.*stand.Sy+stand.muy;
intcept=stand.muy-stand.Sy*stand.mux(xindx)./stand.Sx(xindx)*tembeta;
%%%%%%%%%%%%%%
% betafinal=zeros(1,p);
% Xri=Xtest(:,xindx);
% Tri=T(:,xindx);
% Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
% invAri=Ari\eye(size(Xri,2));
% tembeta=invAri*Xtest(:,xindx)'*Ytest;
% %betafinal(xindx)=stand1.Sy.*(1./stand1.Sx(xindx)).*tembeta';
Ytestest=stand1.muy+stand1.Sy*Xtest(:,xindx)*tembeta;
Ytestobs=Ytest.*stand1.Sy+stand1.muy;
% Ytestest=stand1.muy+log(Ztest)*betafinal'*stand1.Sy;
% Ytestest=intcept+log(Ztest)*betafinal';
% Ytestobs=Ytest.*stand1.Sy+stand1.muy;

1/size(Ytestest,1)*(Ytestobs-Ytestest)'*(Ytestobs-Ytestest)

plot(temp,freq,'ko')
hold on
stem(xindx,top,'-r')
hold off
xlabel('Variable')
ylabel('Frequency')

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
%ax.XTick = 0:100:1000;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);


plot(Ytestobs,Ytestest,'o','MarkerSize',10)
refline(1,0)
xlabel('Observed BMI')
ylabel('Fitted BMI')

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
%ax.XTick = 0:100:1000;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);

%save('realestimation.mat','xindx','gamma','Yobs','Yest','betahat','MSE','Ytestobs','Ytestest')

%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
