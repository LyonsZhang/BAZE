function mat=helm(n,position)

mat =zeros(n - 1, n);
i = 2:n;
r = 1./sqrt(i.*(i - 1));

for j= 1:(n - 1) 
    mat(j, 1:(j + 1)) = [r(j)*ones(1,j), -j * r(j)];
end
mat=mat';

if position==n
    mat=[mat(2:n,:);mat(1,:)];
elseif position>1 || position <n
    mat=[mat(2:position,:);mat(1,:);mat((position+1):n,:)];
end

end