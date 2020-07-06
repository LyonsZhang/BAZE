%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
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
[X,stand.mux,stand.Sx] = zscore(log(Z));
[Y,stand.muy,stand.Sy] = zscore(Y);

[n,p]=size(X); %n is number of samples
 % p is number of predictors

%specify penalty matrix
c=100;
T=[eye(p);c*ones(1,p)]*diag(1./stand.Sx);

%%%read in phylogenetic tree similarity matrix
Q_full = csvread('sortcorr.csv');
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

ss=50;
al=-15;
au=0;
a0=al:(au-al)/ss:au;

tau=1;
nu=0;
omega=0;
seed=100;


cutoff=0.5;

mse=zeros(1,ss+1);
nvs=zeros(1,ss+1);

parfor i=1:(ss+1)
    a=a0(i)*ones(1,p);   
    [gamma,betahat,MSE]=gibbsgamma(nburnin,niter,p,nop,Y, X, T, a, Q, n,tau,nu,omega,seed,true,stand,false);
    freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
    selectindx=freq>cutoff;
    nvs(i)=sum(selectindx==1);
    mse(i)=mean(MSE((nburnin+1):(nburnin+niter)));
end


figname='realrun_';
%get the directory of your input files:
pathname = pwd;

save('stable.mat','a0','mse')


plot(a0,mse,'o-r','MarkerSize',12,'LineWidth',2)
xlabel('Shrinkage parameter a')
ylabel('MSE')

% Set the figure properties
fig = figure(1);
fig.Resize = 'off';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 10, 7];
fig.PaperSize = [10, 7];
fig.Position = [0.1, 0.1, 9.9, 6.9];
fig.InvertHardcopy = 'off';

ax = gca;
ax.TickLabelInterpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
%set(gca,'XColor','none')
%box off
set(gcf, 'Color', [1,1,1]);

%use that when you save
figfile1 = fullfile(pathname, [figname 'stablemse.eps']);
figfile2 = fullfile(pathname, [figname 'stablemse.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)


plot(a0,nvs,'o-r','MarkerSize',12,'LineWidth',2)
xlabel('Shrinkage parameter a')
ylabel('Number of variables selected')

% Set the figure properties
fig = figure(1);
fig.Resize = 'off';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 10, 7];
fig.PaperSize = [10, 7];
fig.Position = [0.1, 0.1, 9.9, 6.9];
fig.InvertHardcopy = 'off';

ax = gca;
ax.TickLabelInterpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);

%use that when you save
figfile1 = fullfile(pathname, [figname 'stablenvs.eps']);
figfile2 = fullfile(pathname, [figname 'stablenvs.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)


%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')