addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\CoreFunctions')
cd('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
%simulate data
%rng(2018)

n=100; %number of samples
p=1000;  % number of predictors
gammatrue=zeros(1,p);

% index of covariate
%true_index=5:8;
true_index=[180:20:400,580:20:800];
gammatrue(true_index)=1;

%b0=2;
% k=11;
% samtemp=round(0.5+1.5*rand(1,k),2);
% select=randsample(k,floor(k/2));
% samtemp(select)=-samtemp(select);
% b1=[samtemp -sum(samtemp)];
% 
% samtemp=round(0.5+1.5*rand(1,k),2);
% select=randsample(k,floor(k/2));
% samtemp(select)=-samtemp(select);
% b2=[samtemp -sum(samtemp)];
b1=[0.8800   -1.3900    1.0400    1.2100   -1.8600   -1.3400    1.7600   -0.9900    0.6900   -0.5400    1.3500   -0.8100];
b2=[-1.4100   -1.1500    0.5100   -1.9500    1.9300   -0.8500   -1.6600    1.4800    1.8700    0.7200    0.6700   -0.1600];

b=zeros(1,p);
%b(gammatrue==1)=b0;
b(true_index(1:2:length(true_index)))=b1;
b(true_index(2:2:length(true_index)))=b2;

%rng(2018)

sigmaX=ones(1,p); % std of X
sigma=0.11758; % noise variance.

Xcor=zeros(p);

for i=180:20:380
    for j=(i+20):20:400
        Xcor(i,j)=0.75-0.0015*abs(i-j);
    end
end

for i=580:20:780
    for j=(i+20):20:800
        Xcor(i,j)=0.75-0.0015*abs(i-j);
    end
end

for i=445:459
    for j=(i+1):460
        Xcor(i,j)=0.4-0.02*abs(i-j);
    end
end

for i=945:959
    for j=(i+1):960
        Xcor(i,j)=0.4-0.02*abs(i-j);
    end
end

Xcor=Xcor+Xcor'+eye(p);

theta=zeros(1,p);
theta([180:20:400 580:20:800])=log(0.5*p)*ones(1,24);
theta([45:60 445:460 945:960])=log(0.25*p)*ones(1,48);

[Xorg, beta,epsilon]=gen_comp_simdata_cor(n,gammatrue,b,theta,sigmaX,Xcor,sigma);

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
snr=signal/noise

type='structure2';

datafile1 = fullfile([type '_snr' num2str(round(snr)) 'simdata.xlsx']);
datafile2 = fullfile([type '_snr' num2str(round(snr)) 'simdata.mat']);

xlswrite(datafile1,Y,1)
xlswrite(datafile1,X,2)
writetable(coeffs,datafile1,'Sheet',3)
save str2_snr10.mat

