function [gamma,betahat,MSE,nselect,Yhat]=gibbsgamma(nburnin,niter,p,nop,Y, X,T, a, Q, n,tau,nu,omega,seed,predict,stand,display)  
    
    % initialize the random number generator
    rng(seed);
    
    %initialize index, nop is the number of initial 1's
    %sample with replacement, p more than nop
    index=randsample(1:p,nop);
    index=sort(index);
    nop=length(index);
    
    %initialize gamma
    gamma=zeros(nburnin+niter+1,p);
    gamma(1,index)=1;
    if isempty(stand)
        invSx=eye(p);
        Sy=1;
        Yobs=Y;
        %Xobs=X;
    else
        invSx=diag(1./stand.Sx);
        Sy=stand.Sy;
        Yobs=Y.*stand.Sy+stand.muy;
        %Xobs=X*diag(stand.Sx)+stand.mux;
    end
    
    %setsort=[sort(index),zeros(1,nsort-nop)];
    
    %generate a sequence of index for update
    %rng('shuffle')
    proposeindx=randsample(1:p,nburnin+niter+1,true);
    
    % initialize Ari and keep
    Xri=X(:,index);
    Tri=T(:,index);
    Ari=Xri'*Xri+tau^(-2)*(Tri'*Tri);
    invAri=Ari\eye(nop);
    Lri=chol(invAri,'lower');
    
    keep.resAi=Y'*Y-Y'*Xri*invAri*Xri'*Y;
    keep.sqrtdetinvAi=sum(log(diag(Lri)));
    
    if predict
        %initialize beta
        betahat=zeros(p,nburnin+niter+1);
        tembeta=invAri*Xri'*Y;
        betahat(index,1)=Sy*invSx(index,index)*tembeta;
        Yhat=zeros(n,nburnin+niter+1);
        Yhat(:,1)=stand.muy+Sy*X(:,index)*tembeta;
        MSE=zeros(1,nburnin+niter+1);
        %MSE(1)=1/(n-nop)*(Yobs-Yhat(:,1))'*(Yobs-Yhat(:,1));
        MSE(1)=1/n*(Yobs-Yhat(:,1))'*(Yobs-Yhat(:,1));
        nselect=zeros(1,nburnin+niter+1);
        nselect(1)=nop;
    end
    
    %rng('shuffle')
    if display
        k=1;
        disp('Gibbs Sampling is starting and will print every 500 iterations.')
    end
    for i=1:(nburnin+niter)      
        gamma(i+1,:)=gamma(i,:);
        betahat(:,i+1)=zeros(p,1);
        flag=any(index==proposeindx(i));
        if flag
            %% propsed variable is inside existing index
            indxtemp=index(index~=proposeindx(i));
            Xri=X(:,indxtemp);
            Xi=X(:,proposeindx(i));
            XIi=[Xri Xi];
            Tri=T(:,indxtemp);
            Ti=T(:,proposeindx(i));
            TIi=[Tri Ti];
            invAritemp=invAri;
            idx = repmat({':'}, ndims(invAri), 1);
            tn=length(index);
            seq=1:tn;

            %test example
%             X=[1,2,3,6,8;4,3,2,7,5;6,7,7,0,4];
%             A=X'*X+diag(ones(1,5))
%             T=inv(A);
%             T(1:4,1:4)- T(1:4,5)*T(5,5)^(-1)*T(5,1:4)
%             inv(A(1:4,1:4))
            
            if proposeindx(i)==max(index)
                idx{1}=seq;
                idx{2}=seq;
            elseif proposeindx(i)==min(index)
                    idx{1}=[2:tn 1];
                    idx{2}=[2:tn 1];
            else
                ti=seq(index==proposeindx(i));
                idx{1}=[1:(ti-1) (ti+1):tn ti];
                idx{2}=[1:(ti-1) (ti+1):tn ti];
            end
            invAitemp=invAri(idx{:});
            %swap the ti variable to the last one of the list
            invAri=invAitemp(1:(tn-1),1:(tn-1))-invAitemp(1:(tn-1),tn)*invAitemp(tn,tn)^(-1)*invAitemp(tn,1:(tn-1));
%            Ari=Xri'*Xri+tau^(-2)*eye(length(indxtemp));
%            invAri=Ari\eye(length(indxtemp));
            [F,keep] = BayesFactor(Y,Xri,Xi,XIi,Tri,Ti,TIi,invAri,n,tau,nu,omega,flag,keep);
            pgammai1=exp(a(proposeindx(i))+Q(proposeindx(i),indxtemp)*gamma(i,indxtemp)');
            
            pcond=1/(1+F^(-1)/pgammai1);
            newgamma=binornd(1,pcond);
            gamma(i+1,proposeindx(i))=newgamma;
            if (newgamma==0) && (tn>1)
                %remove the proposed variable
                index=indxtemp;
                keep.resAi=keep.resAri;
                keep.sqrtdetinvAi=keep.sqrtdetinvAri;
                if predict
                    tembeta=invAri*Xri'*Y;
                    betahat(index,i+1)=Sy*invSx(index,index)*tembeta;
                    Yhat(:,i+1)=stand.muy+Sy*X(:,index)*tembeta;
                    %MSE(i+1)=1/(n-tn+1)*(Yobs-Yhat(:,i+1))'*(Yobs-Yhat(:,i+1));
                    MSE(i+1)=1/n*(Yobs-Yhat(:,i+1))'*(Yobs-Yhat(:,i+1));
                    nselect(i+1)=tn-1;
                end
            elseif (newgamma~=0) || (tn==1)
                invAri=invAritemp;
                if predict
                    betahat(:,i+1)=betahat(:,i);
                    Yhat(:,i+1)=Yhat(:,i);
                    MSE(i+1)=MSE(i);
                    nselect(i+1)=nselect(i);
                end
            end
            
        else
            %% proposed variable is outside existing index
            indxtemp=index;
            Xri=X(:,indxtemp);
            Xi=X(:,proposeindx(i));
            XIi=[Xri Xi];
            Tri=T(:,indxtemp);
            Ti=T(:,proposeindx(i));
            TIi=[Tri Ti];
            [F,keep,invAi] = BayesFactor(Y,Xri,Xi,XIi,Tri,Ti,TIi,invAri,n,tau,nu,omega,flag,keep);
            pgammai1=exp(a(proposeindx(i))+Q(proposeindx(i),indxtemp)*gamma(i,indxtemp)');
            
            pcond=1/(1+F^(-1)/pgammai1);
            newgamma=binornd(1,pcond);
            gamma(i+1,proposeindx(i))=newgamma;
            if newgamma==1
                %include the proposed variable
                index=sort([indxtemp proposeindx(i)]);
                invAri=invAi;
                idx = repmat({':'}, ndims(invAri), 1);
                tn=length(index);
                seq=1:tn;
                if proposeindx(i)>max(indxtemp)
                    idx{1}=seq;
                    idx{2}=seq;
                elseif proposeindx(i)<min(indxtemp)
                    idx{1}=[tn 1:(tn-1)];
                    idx{2}=[tn 1:(tn-1)];
                else
                    ti=seq(index==proposeindx(i));
                    idx{1}=[1:(ti-1) tn ti:(tn-1)];
                    idx{2}=[1:(ti-1) tn ti:(tn-1)];
                end
                invAri=invAri(idx{:});
                %put the last variable to the new index list with order
                if predict
                    tembeta=invAri*X(:,index)'*Y;
                    betahat(index,i+1)=Sy*invSx(index,index)* tembeta;
                    Yhat(:,i+1)=stand.muy+Sy*X(:,index)*tembeta;
                    %MSE(i+1)=1/(n-tn)*(Yobs-Yhat(:,i+1))'*(Yobs-Yhat(:,i+1));
                    MSE(i+1)=1/n*(Yobs-Yhat(:,i+1))'*(Yobs-Yhat(:,i+1));
                    nselect(i+1)=tn;
                end
            else
                keep.resAi=keep.resAri;
                keep.sqrtdetinvAi=keep.sqrtdetinvAri;
                if predict
                    betahat(:,i+1)=betahat(:,i);
                    Yhat(:,i+1)=Yhat(:,i);
                    MSE(i+1)=MSE(i);
                    nselect(i+1)=nselect(i);
                end
            end
        end
        if display
            if mod(k,500)==0
            disp(['iteration is ',num2str(k)])
            end
        k=k+1;
        end
    end
    if display
        disp('Gibbs Sampling is ending and will print the frequency of select variables.')
        sum(gamma(nburnin:(nburnin+niter),:),1)
    end
end
