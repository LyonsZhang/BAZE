%%%load functions
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianContrast\ContrastCoreFunctions\')
%%%load simulated data
addpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')
%%%set path
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
mkdir('results\simulation');
cd('results\simulation')

rng(2018)

load str2_snr0d1.mat

%standadize data
pid=p;
[X,Ts,stand,Xs]=trans_comp(Z,'log','additive',pid,true);
%choose from log, logit, additive, centered, isometric
p=p-1;
%dimension degenerate
[Y,stand.muy,stand.Sy] = zscore(Y);
T=Ts*diag(1./stand.Sx(1:p));

%set the number of burn-in steps and iterations
nburnin=15000;
niter=5000;

%initialize gamma and setsort
nop=floor(n/2);

a0=-10;
a=a0*ones(1,p);

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

tic
[gamma,betahat,MSE,nselect,Yhat,etahat]=gibbsgamma(nburnin,niter,p,nop,Y, X, T,Xs,Ts, a, Q, n,tau,nu,omega,seed,true,stand,true);
toc



%%%%%%%%%%%Describe blocks%%%%%%
x=[180:20:400];                  %#initialize x array
y1=0.6+zeros(1,length(x));                      %#create first curve
y2=0.8+zeros(1,length(x));                   %#create second curve
X1=[x,fliplr(x)];                %#create continuous x value array for plotting
Y1=[y1,fliplr(y2)];              %#create y values for out and then back

x=[580:20:800];                  %#initialize x array
y1=0.6+zeros(1,length(x));                      %#create first curve
y2=0.8+zeros(1,length(x));                   %#create second curve
X2=[x,fliplr(x)];                %#create continuous x value array for plotting
Y2=[y1,fliplr(y2)];  %#create y values for out and then back

x=[445:460];                  %#initialize x array
y1=0.6+zeros(1,length(x));                      %#create first curve
y2=0.8+zeros(1,length(x));                   %#create second curve
X3=[x,fliplr(x)];                %#create continuous x value array for plotting
Y3=[y1,fliplr(y2)];  %#create y values for out and then back


x=[945:960];                  %#initialize x array
y1=0.6+zeros(1,length(x));                      %#create first curve
y2=0.8+zeros(1,length(x));                   %#create second curve
X4=[x,fliplr(x)];                %#create continuous x value array for plotting
Y4=[y1,fliplr(y2)];              %#create y values for out and then back

x=[45:60];                    %#initialize x array
y1=0.6+zeros(1,length(x));                      %#create first curve
y2=0.8+zeros(1,length(x));                   %#create second curve
X5=[x,fliplr(x)];                %#create continuous x value array for plotting
Y5=[y1,fliplr(y2)];              %#create y values for out and then back



%%%%%%%%%%%plot figure%%%%%%%%%%%

figname='sim51';

cutoff=0.5;

freq=sum(gamma((nburnin+1):(nburnin+niter),:))/niter;
top=freq(freq>cutoff);
temp=1:p;
xindx=temp(freq>cutoff);

etafinal=zeros(1,p);
%betafinal=zeros(1,p+1);
Xri=X(:,xindx);
Tri=T(:,xindx);
Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
invAri=Ari\eye(size(Xri,2));
etafinal(xindx)=stand.Sy*diag(1./stand.Sx(xindx))*invAri*X(:,xindx)'*Y;
betafinal=Ts*etafinal';
coef_l2=sqrt(sum(abs(betafinal-beta).^2));
nsel=sum(xindx~=0);
betasigma=1/(n+nu-2)*(MSE(nburnin+niter)*(n-nsel)+nu*omega).*diag(Ts([xindx true],xindx)*invAri*Ts([xindx true],xindx)');
%diag(1./(stand.Sx(xindx).^2))

plot(true_index,0.1*ones(size(true_index)),'s','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.800000 0.500000 0.200000])
hold on
%%%%add shaded area of priors%%%%%

fill(X1,Y1,[.466 .674 .188],'FaceAlpha', 0.1);
fill(X2,Y2,[.466 .674 .188],'FaceAlpha', 0.1);
fill(X3,Y3,'m','FaceAlpha', 0.1);
fill(X4,Y4,'m','FaceAlpha', 0.1);
fill(X5,Y5,'m','FaceAlpha', 0.1);
hold on
plot(temp,freq,'ko')
hold on
plot(xindx,top,'rp','MarkerSize',12,'MarkerFaceColor','red')
hold on
stem(xindx,top,'-r')
hold off
xlabel('Variable')
ylabel('Frequency')

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
ax.XTick = 0:100:1000;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;

%get the directory of your input files:
pathname = pwd;
%use that when you save
figfile1 = fullfile(pathname, [figname '_signal' num2str(signal) 'noise' num2str(noise) 'selected.eps']);
figfile2 = fullfile(pathname, [figname '_signal' num2str(signal) 'noise' num2str(noise) 'selected.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

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
ax.XTick = 0:1000:(nburnin+niter);
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 16;
%use that when you save
figfile1 = fullfile(pathname, [figname '_signal' num2str(signal) 'noise' num2str(noise) 'MSE.eps']);
figfile2 = fullfile(pathname, [figname '_signal' num2str(signal) 'noise' num2str(noise) 'MSE.pdf']);
saveas(gcf,figfile1,'epsc')
saveas(gcf,figfile2)

close

xindx
sum(betafinal(xindx))
coef_l2

%%%unload functions
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\BayesianContrast\ContrastCoreFunctions\')
%%%unload simulated data
rmpath('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\GenerateSimData\')