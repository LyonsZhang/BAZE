addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\CoreFunctions')
cd('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
%simulate data
rng(2018)

n=100; %number of samples
p=1000;  % number of predictors
snr=5;
gammatrue=zeros(1,p);

% index of covariate
%true_index=5:8;
true_index=[1:3 6:8];
gammatrue(true_index)=1;

b=zeros(1,p);
%b(gammatrue==1)=b0;
b(1:3)=[1 -0.8 0.6];
b(6:8)=[-1.5 -0.5 1.2];

%rng(2018)

sigmaX=1; % std of X
sigma=1/snr*mean(abs(b(b~=0))); % noise variance.

theta=zeros(1,p);
theta(1:5)=log(0.5*p)*ones(1,5);
%theta([45:60 445:460 945:960])=log(0.25*p)*ones(1,48);

[Xorg, beta,epsilon]=gen_comp_simdata_ind(n,p,gammatrue,b,theta,sigmaX,sigma);

temp=exp(2*Xorg);
Z=bsxfun(@rdivide,temp,sum(temp,1));
X=log(Z);
Y=X*beta'+epsilon;
beta=beta';
gammatrue=gammatrue';

coeffs=table(beta,gammatrue);

%signal=round(std(log(Z)*beta'),2)
%noise=round(std(Y-log(Z)*beta'),2)

signal=mean(abs(b(b~=0)))
noise=sigma

type='structure1';

datafile1 = fullfile([type '_snr' num2str(round(snr)) num2str(n) '_p' num2str(p) 'simdata.xlsx']);
datafile2 = fullfile([type '_snr' num2str(round(snr)) num2str(n) '_p' num2str(p) 'simdata.mat']);

xlswrite(datafile1,Y,1)
xlswrite(datafile1,X,2)
writetable(coeffs,datafile1,'Sheet',3)
save str1_snr5_n100_p1000.mat

