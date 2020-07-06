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

load str2_snr1.mat

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

ss=100;
al=-13;
au=-10;
a0=al:(au-al)/ss:au;

tau=1;
nu=0;
omega=0;
seed=100;

Q=0.002*(ones(p,p)-eye(p));
for i=180:20:380
    for j=(i+20):20:400
        Q(i,j)=4;
    end
end

for i=350:10:440
    for j=(i+10):10:450
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

% FP=zeros(1,ss+1);
% TN=zeros(1,ss+1);
% FN=zeros(1,ss+1);
% TP=zeros(1,ss+1);
FPR=zeros(1,ss+1);
TPR=zeros(1,ss+1);
mse=zeros(1,ss+1);
nvs=zeros(1,ss+1);

k=5;
tn=10;
indtemp=zeros(n,tn);
for i=1:tn
    indtemp(:,i)=crossvalind('Kfold',n,k);
end

parfor i=1:(ss+1)
    a=a0(i)*ones(1,p);
    res=arrayfun(@(t) arrayfun(@(j) nseltrain(j,indtemp(:,t),X,Y,n,p,gammatrue,nburnin,niter,T,a,Q,tau,nu,omega,seed,stand,cutoff),1:k,'UniformOutput',false),1:tn,'UniformOutput',false);
    temp=median(cell2mat(cat(2,res{:})'),1);
    nvs(i)=temp(1);
    TPR(i)=temp(2);
    FPR(i)=temp(3);
    mse(i)=temp(4);
end


figname='str2_snr1_';
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
set(gcf, 'Color', [1,1,1]);

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
set(gcf, 'Color', [1,1,1]);


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
%set(gca,'XColor','none')
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