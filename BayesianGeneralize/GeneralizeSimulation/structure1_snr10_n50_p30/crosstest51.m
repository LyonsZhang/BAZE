%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%load simulated data
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
%%%set path
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
mkdir('results\crosstest');
cd('results\crosstest')

rng(2018)

load str1_snr10_n50_p30.mat

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
%nop=floor(n/2);

a0=-10;
a=a0*ones(1,p);

tau=1;
nu=0;
omega=0;
seed=100;

Q=0.000*ones(p,p);

% for i=1:7
%     for j=(i+1):8
%         Q(i,j)=4;
%     end
% end

%% cross validation
Q=(Q+Q')/2;
cutoff=0.5;
predict=true;
k=10;

pl=true;
if pl
%get the directory of your input files:
pathname = pwd;
figname='sim51'; 
else
    pathname=[];
    figname= [];  
end

rn=10;
MSEall=zeros(rn,k);
MSEpreall=zeros(rn,k);

FPR_all=zeros(rn,k);
FNR_all=zeros(rn,k);
FPN_all=zeros(rn,k);
FNN_all=zeros(rn,k);
sensitivity_all=zeros(rn,k);
specificity_all=zeros(rn,k);
precision_all=zeros(rn,k);
accuracy_all=zeros(rn,k);

coef_l1_all=zeros(rn,k);
coef_l2_all=zeros(rn,k);
coef_linf_all=zeros(rn,k);

%%%%%roc curve
doroc=true;
sn=0.05;
if doroc
    rFPR=zeros(rn*k,1/sn+1);
    rTPR=zeros(rn*k,1/sn+1);
    rsensitivity=zeros(rn*k,1/sn+1);
    rspecificity=zeros(rn*k,1/sn+1);
    rprecision=zeros(rn*k,1/sn+1);
    raccuracy=zeros(rn*k,1/sn+1);
    rAUC=zeros(rn,k);
end

tic
for repeat=1:rn
    if doroc
    [MSEpre,MSEtest,average,sem,FPR,FNR,FPN,FNN,sensitivity,specificity,precision,accuracy,etafinal,roc]=traintest(k,nburnin,niter,p,Y, X, T,beta, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand,pl,true_index,pathname,figname,doroc);
    rFPR(((repeat-1)*k+1):(repeat*k),:)=roc.rFPR;
    rTPR(((repeat-1)*k+1):(repeat*k),:)=roc.rTPR;
%     rsensitivity(((repeat-1)*k+1):(repeat*k),:)=roc.rsensitivity;
%     rspecificity(((repeat-1)*k+1):(repeat*k),:)=roc.rspecificity;
%     rprecision(((repeat-1)*k+1):(repeat*k),:)=roc.rprecision;
%     raccuracy(((repeat-1)*k+1):(repeat*k),:)=roc.raccuracy;
    rAUC(repeat,:)=roc.rAUC; 
    else
        [MSEpre,MSEtest,average,sem,FPR,FNR,FPN,FNN,sensitivity,specificity,precision,accuracy,etafinal]=traintest(k,nburnin,niter,p,Y, X, T,beta, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand, pl,true_index,pathname,figname,doroc);
    end
    MSEall(repeat,:)=MSEtest;
    MSEpreall(repeat,:)=MSEpre;
    FPR_all(repeat,:)=FPR;
    FNR_all(repeat,:)=FNR;
    FPN_all(repeat,:)=FPN;
    FNN_all(repeat,:)=FNN;
    sensitivity_all(repeat,:)=sensitivity;
    specificity_all(repeat,:)=specificity;
    precision_all(repeat,:)=precision;
    accuracy_all(repeat,:)=accuracy;
    
    coef_l1_all(repeat,:)=sum(abs(etafinal-beta'),2);
    coef_l2_all(repeat,:)=sqrt(sum(abs(etafinal-beta').^2,2))';
    coef_linf_all(repeat,:)=max(abs(etafinal-beta'),[],2);
end
toc

varCell{1} =MSEall;
varCell{2} = FPR_all;
varCell{3}= FNR_all;
varCell{4} = FPN_all;
varCell{5}= FNN_all;
varCell{6}=sensitivity_all;
varCell{7}=specificity_all;
varCell{8}=precision_all;
varCell{9}=accuracy_all;

varCell{10}=coef_l1_all;
varCell{11}=coef_l2_all;
varCell{12}=coef_linf_all;
varCell{13}=MSEpreall;

meanall=cellfun(@mean2,varCell)
meanstd=cellfun(@std2,varCell)/sqrt(k*rn)

varCell{14}=[meanall;meanstd];

if doroc
    pathname= pwd;
    FPRm=mean(rFPR,1);
    FPRs=std(rFPR,1)/sqrt(k*rn);
    TPRm=mean(rTPR,1);
    TPRs=std(rTPR,1)/sqrt(k*rn);
    AUCm=mean2(rAUC);
    AUCs=std2(rAUC)/sqrt(k*rn);
    xlswrite(fullfile(pathname,'roc.xlsx'),FPRm,'Sheet1');
    xlswrite(fullfile(pathname,'roc.xlsx'),TPRm,'Sheet2');
    x_vector = [FPRm+FPRs, fliplr(FPRm-FPRs)];
    patch = fill(x_vector, [TPRm+TPRs,fliplr(TPRm-TPRs)],[128 193 219]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
    hold on
    plot(linspace(0,1),linspace(0,1));
    hold on;
    plot(FPRm, TPRm, 'color', [52 148 186]./255, ...
        'LineWidth', 2);
    xlabel('False positive rate') 
    ylabel('True positive rate')
    title('ROC curve')
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
    %ax.XTick = 0:1000:(nburnin+niter);
    ax.Box = 'off';
    ax.LineWidth = 1.5;
    ax.FontSize = 16;
    %use that when you save
    figfile1 = fullfile(pathname, ['ROC_signal' num2str(signal) 'noise' num2str(noise) '.eps']);
    figfile2 = fullfile(pathname, ['ROC_signal' num2str(signal) 'noise' num2str(noise) '.pdf']);
    saveas(gcf,figfile1,'epsc')
    saveas(gcf,figfile2)
end

filename = 'crosstest.xlsx';
for i=1:14
    sheetname=strcat('Sheet',num2str(i));
    xlswrite(filename,varCell{i},sheetname);
end

%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')