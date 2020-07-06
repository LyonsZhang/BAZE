%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%load simulated data
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
%%%set path
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
mkdir('results\simulation');
cd('results\simulation')

rng(2018)

load str2_snr10.mat

%standadize data
[X,stand.mux,stand.Sx] = zscore(log(Z));
[Y,stand.muy,stand.Sy] = zscore(Y);

%specify penalty matrix
c=100;
T=[eye(p);c*ones(1,p)]*diag(1./stand.Sx);

%set the number of burn-in steps and iterations
nburnin=10000;
niter=5000;

%initialize gamma and setsort
nop=floor(n/2);

sc=100;
sl=0;
su=3;
ss=sl:(su-sl)/sc:su;

a0=-10;
a=a0*ones(1,p);

tau=1;
%nu=0;
%omega=0;
seed=100;

Q=0.002*(ones(p,p)-eye(p));
for i=180:20:380
    for j=(i+20):20:400
        Q(i,j)=4;
    end
end

for i=580:20:780
    for j=(i+20):20:800
        Q(i,j)=4;
    end
end

for i=445:459
    for j=(i+1):460
        Q(i,j)=4;
    end
end

for i=945:959
     for j=(i+1):960
         Q(i,j)=4;
     end
end

for i=45:59
     for j=(i+1):60
         Q(i,j)=4;
     end
end

Q=(Q+Q')/2;

cutoff=0.5;

FP=zeros(1,sc+1);
TN=zeros(1,sc+1);
FN=zeros(1,sc+1);
TP=zeros(1,sc+1);
FPR=zeros(1,sc+1);
TPR=zeros(1,sc+1);
mse=zeros(1,sc+1);
nvs=zeros(1,sc+1);
coef_l2=zeros(1,sc+1);
alpha0=zeros(1,sc+1);
beta0=zeros(1,sc+1);

parfor i=1:(sc+1)
    nu=ss(i);
    omega=ss(i);
    [gamma,betahat,MSE]=gibbsgamma(nburnin,niter,p,nop,Y, X, T, a, Q, n,tau,nu,omega,seed,true,stand,false);
    
    freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
    
    temp=1:p;
    xindx=temp(freq>cutoff);
    betafinal=zeros(p,1);
    Xri=X(:,xindx);
    Tri=T(:,xindx);
    Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
    invAri=Ari\eye(size(Xri,2));
    betafinal(xindx)=stand.Sy*diag(1./stand.Sx(xindx))*invAri*X(:,xindx)'*Y;
    coef_l2(i)=sqrt(sum(abs(betafinal-beta).^2));
    SSE=(Y.*stand.Sy-(X*diag(stand.Sx))*betafinal)'*(Y.*stand.Sy-(X*diag(stand.Sx))*betafinal);
    alpha0(i)=(n+nu)/2;
    beta0(i)=(SSE+nu*omega)/2;
    
    
    selectindx=freq>cutoff;
    nvs(i)=sum(selectindx==1);
    fp=sum((gammatrue'==0) & (selectindx==1));
    tn=sum((gammatrue'==0) & (selectindx==0));
    fn=sum((gammatrue'==1) & (selectindx==0));
    tp=sum((gammatrue'==1) & (selectindx==1));
    FPR(i)= fp/(fp+tn);
    TPR(i)= tp/(tp+fn);
    FP(i)=fp;
    TN(i)=tn;
    FN(i)=fn;
    TP(i)=tp;
    mse(i)=mean(MSE((nburnin+1):(nburnin+niter)));
end

figname='str2_snr10_';
%get the directory of your input files:
pathname = pwd;

save('stable.mat','a0','FPR','TPR','FP','TN','FN','TP','mse')

%%%ROC curve
patch = fill([0,0,1,1], [0,1,1,0],[128 193 219]./255);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
hold on
plot(linspace(0,1),linspace(0,1),'LineWidth',5,'color',[52 148 186]./255);
hold on;
plot(FPR,TPR,'*-r','MarkerSize',12,'LineWidth',5)
xlabel('FPR')
ylabel('TPR')

% Set the figure properties
fig = figure(1);
fig.Resize = 'off';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 8, 8];
fig.PaperSize = [8, 8];
fig.Position = [0.1, 0.1, 7.9, 7.9];
fig.InvertHardcopy = 'off';

ax = gca;
ax.TickLabelInterpreter = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.FontName = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.XTick = 0:0.2:1;
ax.YTick = 0:0.2:1;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
box off

%use that when you save
figfile1 = fullfile(pathname, [figname 'stableROC.eps']);
figfile2 = fullfile(pathname, [figname 'stableROC.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

hold off
%%%%%stable curve
plot(a0,TPR,'o-r','MarkerSize',12,'LineWidth',2)
xlabel('Shrinkage parameter a')
ylabel('TPR')

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
box off

%use that when you save
figfile1 = fullfile(pathname, [figname 'stableTPR.eps']);
figfile2 = fullfile(pathname, [figname 'stableTPR.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

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
set(gca,'XColor','none')
box off


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

plot(ss,nvs,'o-r','MarkerSize',12,'LineWidth',2)
xlabel('Prior parameters \nu and \omega')
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
% ax.TickLabelInterpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% ax.FontName = 'LaTeX';
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);

%use that when you save
figfile1 = fullfile(pathname, [figname 'nuomega_nvs.eps']);
figfile2 = fullfile(pathname, [figname 'nuomega_nvs.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

plot(ss,round(beta0./alpha0,4),'o-r','MarkerSize',12,'LineWidth',2)
xlabel('Prior parameters \nu and \omega')
ylabel('Mean value of \sigma^2 posterior')

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
% ax.TickLabelInterpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% ax.FontName = 'LaTeX';
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);

%use that when you save
figfile1 = fullfile(pathname, [figname 'nuomega_meanposterior.eps']);
figfile2 = fullfile(pathname, [figname 'nuomega_meanposterior.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

tem=round((beta0*2-ss.^2)/n,4);
plot(ss,tem,'o-r','MarkerSize',12,'LineWidth',2)
ylim([-0.1,0.1])
xlabel('Prior parameters \nu and \omega')
ylabel('The mean squared errors')

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
% ax.TickLabelInterpreter = 'LaTeX';
% ax.TickLabelInterpreter = 'LaTeX';
% ax.FontName = 'LaTeX';
% ax.Title.Interpreter = 'LaTeX';
% ax.XLabel.Interpreter = 'LaTeX';
% ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
set(gcf, 'Color', [1,1,1]);

%use that when you save
figfile1 = fullfile(pathname, [figname 'nuomega_MSE.eps']);
figfile2 = fullfile(pathname, [figname 'nuomega_MSE.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)


%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')