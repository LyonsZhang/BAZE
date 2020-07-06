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
[Xall,stand.mux,stand.Sx] = zscore(log(Z));
[Yall,stand.muy,stand.Sy] = zscore(Y);
X=Z;

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
nburnin=22000;
niter=3000;

%initialize gamma and setsort
nop=floor(n/2);

a0=-8;
a=a0*ones(1,p);

tau=1;
nu=0;
omega=0;
seed=2020;


cutoff=0.5;
predict=true;

pl=false;
if pl
%get the directory of your input files:
pathname = pwd;
figname='realgeneralize'; 
else
    pathname=[];
    figname= [];  
end

rn=100;
MSEall=zeros(rn,1);
MSEpreall=zeros(rn,1);

beta_all=zeros(rn,p);

tic
parfor repeat=1:rn
    if pl
    [MSEpre,MSEtest,etafinal]=crossvalidatereal0(nburnin,niter,p,Y, X, T, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand,pl,pathname,figname);
    else
        [MSEpre,MSEtest,etafinal]=crossvalidatereal0(nburnin,niter,p,Y, X, T, a, Q, n,tau,nu,omega,seed,predict,cutoff,stand, pl,pathname,figname);
    end
    MSEall(repeat)=MSEtest;
    MSEpreall(repeat)=MSEpre;
    
    beta_all(repeat,:)=etafinal;
end
toc

varCell{1} =MSEall;
varCell{2}=MSEpreall;

varCell{3}=beta_all;


meanall=cellfun(@mean2,varCell)
meanstd=cellfun(@std2,varCell)/sqrt(k*rn)

varCell{4}=[meanall;meanstd];

tt=(1:p);
newindx=tt(mean(beta_all>0)>0.7)


Ztrain=csvread('train_otu.csv');
Ytrain=csvread('train_bmi.csv');
Ztest=csvread('test_otu.csv');
Ytest=csvread('test_bmi.csv');

%standadize data
[Xtrain,stand1.mux,stand1.Sx] = zscore(log(Ztrain));
[Ytrain,stand1.muy,stand1.Sy] = zscore(Ytrain);
[Xtest,stand2.mux,stand2.Sx] = zscore(log(Ztest));
[Ytest,stand2.muy,stand2.Sy] = zscore(Ytest);
%specify penalty matrix
c=100;
Ttrain=[eye(p);c*ones(1,p)]*diag(1./stand1.Sx);

Xri=Xtrain(:,newindx);
Tri=Ttrain(:,newindx);
Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
invAri=Ari\eye(size(Xri,2));
tembeta=invAri*Xtrain(:,newindx)'*Ytrain;
Ytestest1=stand2.muy+stand2.Sy*Xtest(:,newindx)*tembeta;
Ytestobs1=Ytest.*stand2.Sy+stand2.muy;

1/size(Ytestest1,1)*(Ytestobs1-Ytestest1)'*(Ytestobs1-Ytestest1)


filename = 'test.xlsx';
for i=1:4
    sheetname=strcat('Sheet',num2str(i));
    xlswrite(filename,varCell{i},sheetname);
end


%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianGeneralize\GeneralizeCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')