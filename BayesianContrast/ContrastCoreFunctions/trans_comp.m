function [out,T,stand,Xs]=trans_comp(X,self_type,contrast_type,position,standard)

if strcmp(self_type,'log')
    if ~standard
        Xs=log(X);
        stand=[];
    else
        [Xs,stand.mux,stand.Sx] = zscore(log(X));
    end
elseif strcmp(self_type,'logit')
    if ~standard
        Xs=log(X./(1-X));
        stand=[];
    else
        [Xs,stand.mux,stand.Sx] = zscore(log(X./(1-X)));
    end       
else
    error(' "self_type" should be one of the two: "log" or "logit" ')
end

switch contrast_type
    case 'additive'
        T=add_contrast(size(Xs,2),position);
        out=Xs*T;
    case 'centered'
        T=center_contrast(size(Xs,2),position);
        out=Xs*T;
    case 'isometric'
        T=helm(size(Xs,2),position);
        out=Xs*T;
    otherwise
        error(' "contrast_type" should be one of the three: "additive" or "centered" or "isometric" ')
end


end

