function mat=add_contrast(n,position)

if position ==1
    mat=[-1*ones(1,n-1);eye(n-1)];
elseif position==n
    mat=[eye(n-1);-1*ones(1,n-1)];
else
    mat=[eye(position-1),zeros(position-1,n-position);-1*ones(1,n-1);zeros(n-position,position-1),eye(n-position)];
end

end