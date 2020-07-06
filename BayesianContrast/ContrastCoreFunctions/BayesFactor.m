function [F,keep,invAi] = BayesFactor(Y,Xri,Xi,XIi,Tri,Ti,TIi,invAri,n,tau,nu,omega,flag,keep)
%    digits(50)
    %Ai=Xri'*Xri+tau^(-2)*eye(pri);
    keep=keep;
    Lri=chol(invAri,'lower');
    
    sqrtdetinvAri=sum(log(diag(Lri)));
    resAri=Y'*Y-Y'*Xri*invAri*Xri'*Y;
         

     if flag

         %copy invAi from previous step
         sqrtdetinvAi=keep.sqrtdetinvAi;
         resAi=keep.resAi;
         
         keep.resAri=resAri;
         keep.sqrtdetinvAri=sqrtdetinvAri;
         logratiodetT=sum(log(diag(chol(Tri'*Tri))))-sum(log(diag(chol(TIi'*TIi))));

     else
         %update invAi from invAri
        Srii=Xri'*Xi+tau^(-2)*(Tri'*Ti);
        sii=Xi'*Xi+tau^(-2)*(Ti'*Ti);
        
        %calculate inverse Ai^(-1) given Ari^(-1)
        v=sqrt(1/(sii*(1-Srii'*invAri*Srii/sii)))*invAri*Srii;
        A11=invAri+v*v';
        %A112=invAri-Srii/sii*Srii';
        A12=-A11*Srii/sii;
        A21=-1/sii*Srii'*A11;
        A22=1/sii+1/sii*Srii'*A11*Srii/sii;
        invAi=[A11,A12;A21,A22];
        %eig(invAi);
        resAi=Y'*Y-Y'*XIi*invAi*XIi'*Y;

        
        
        %calculate determinant of Ai^(-1) and Ari^(-1)
        %Rank 1 update to Cholesky factorization
        tilLri=cholupdate(Lri',v)';
        %tilLri=chol(A11,'lower');
        opts.LT = true;
        Lrii=linsolve(tilLri,A12,opts);
        lii=sqrt(A22-Lrii'*Lrii);
        Li=[tilLri,zeros(size(tilLri,1),1);Lrii',lii];
        sqrtdetinvAi=sum(log(diag(Li)));
        
        
        %keep quantities for future steps
        keep.resAri=resAri;
        keep.sqrtdetinvAri=sqrtdetinvAri;

        keep.resAi=resAi;
        keep.sqrtdetinvAi=sqrtdetinvAi;
        logratiodetT=-sum(log(diag(chol(Tri'*Tri))))+sum(log(diag(chol(TIi'*TIi))));
       
    end
    
    %calculate the Bayes factor bf=P(Y|gammai=0,gamma_-i)/P(Y|gammai=1,gamma_-i)
         
    F=-log(tau)+(sqrtdetinvAri-sqrtdetinvAi)+logratiodetT+(n+nu)/2*log((nu*omega+resAri)/(nu*omega+resAi));

    F=exp(F);
end
