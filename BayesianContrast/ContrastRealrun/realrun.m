%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianContrast\ContrastCoreFunctions\')
%%%load simulated data
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated\')
%%%set path
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
mkdir('results');
cd('results')

Z=csvread('otu.csv');
Y=csvread('bmi.csv');

%standadize data
pid=size(Z,2);
[X,Ts,stand,Xs]=trans_comp(Z,'log','additive',pid,true);
%choose from log, logit, additive, centered, isometric
[n,p]=size(X); %n is number of samples
 % p is number of predictors

%dimension degenerate
[Y,stand.muy,stand.Sy] = zscore(Y);
T=Ts*diag(1./stand.Sx(1:p));

%%%read in phylogenetic tree similarity matrix
Q_full = csvread('sortcorr.csv');
Q_full(:,pid)=[];
Q_full(pid,:)=[];
Q_full(Q_full>0.9)=0.9;
Q_full=Q_full+0.1*eye(p);
Q_full=Q_full\eye(p);
Q_full=Q_full-diag(diag(Q_full));
Q_full(Q_full>1)=1;
Q_full(Q_full<-1)=-1;

sqrt(sum(sum(Q_full>0)))
Q=sparse(Q_full);
%Q=zeros(p);
histogram(Q)
sum(sum(Q,2))/p


%set the number of burn-in steps and iterations
nburnin=15000;
niter=5000;

%initialize gamma and setsort
nop=floor(n/2);

a0=-8;
a=a0*ones(1,p);

tau=1;
nu=0;
omega=0;
seed=100;

tic
[gamma,betahat,MSE,nselect,Yhat,etahat]=gibbsgamma(nburnin,niter,p,nop,Y, X, T,Xs,Ts, a, Q, n,tau,nu,omega,seed,true,stand,true);
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


plot(Yobs,Yest,'o','MarkerSize',10)
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

%save('realestimation.mat',xindx,gamma,Yobs,Yest,betahat,MSE)

%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianContrast\ContrastCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')