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
[X,Ts,stand,Xs]=trans_comp(Z,'log','centered',pid,true);
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
Q_full(Q_full<0.6)=0;
Q_full= Q_full - diag(diag(Q_full));

sqrt(sum(sum(Q_full>0)))
Q=sparse(Q_full);
%Q=zeros(p);
histogram(Q)
sum(sum(Q,2))/p

%set the number of burn-in steps and iterations
nburnin=10000;
niter=5000;

%initialize gamma and setsort
nop=floor(n/2);

a0=-11;
a=a0*ones(1,p);

tau=1;
nu=0;
omega=0;
seed=100;


cutoff=0.5;
predict=true;
k=10;

pl=false;
if pl
%get the directory of your input files:
pathname = pwd;
figname='realcontrast'; 
else
    pathname=[];
    figname= [];  
end

rn=10;
MSEall=zeros(rn,k);
MSEpreall=zeros(rn,k);

beta_all=zeros(rn*k,p+1);

tic
for repeat=1:rn
    if pl
    [MSEpre,MSEtest,average,sem,etafinal]=traintestreal(k,nburnin,niter,p,Y, X, T,Xs,Ts, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand,pl,pathname,figname);
    else
        [MSEpre,MSEtest,average,sem,etafinal]=traintestreal(k,nburnin,niter,p,Y, X, T,Xs,Ts, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand, pl,pathname,figname);
    end
    MSEall(repeat,:)=MSEtest;
    MSEpreall(repeat,:)=MSEpre;
    
    beta_all(((repeat-1)*k+1):(repeat*k),:)=etafinal;
end
toc

varCell{1} =MSEall;
varCell{2}=MSEpreall;

varCell{3}=beta_all;


meanall=cellfun(@mean2,varCell)
meanstd=cellfun(@std2,varCell)/sqrt(k*rn)

varCell{4}=[meanall;meanstd];


filename = 'crosstest_centered.xlsx';
for i=1:4
    sheetname=strcat('Sheet',num2str(i));
    xlswrite(filename,varCell{i},sheetname);
end

%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianContrast\ContrastCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')