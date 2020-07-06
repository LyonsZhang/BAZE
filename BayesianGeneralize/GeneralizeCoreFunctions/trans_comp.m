function out=trans_comp(X,self_type,contrast_type,position)

if strcmp(self_type,'log')
    X=log(X);
elseif strcmp(self_type,'logit')
    X=log(X./(1-X));
else
    error(' "self_type" should be one of the two: "log" or "logit" ')
end

switch contrast_type
    case 'additive'
        out=X*add_contrast(size(X,2),position);
    case 'centered'
        out=X*center_contrast(size(X,2),position);
    case 'isometric'
        out=X*helm(size(X,2),position);
    otherwise
        error(' "contrast_type" should be one of the three: "additive" or "centered" or "isometric" ')
end


end

