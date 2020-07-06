function mat=center_contrast(n,position)

if position ==1
    mat=[-1/n*ones(1,n-1);eye(n-1)-1/n];
elseif position==n
    mat=[eye(n-1)-1/n;-1/n*ones(1,n-1)];
else
    mat=[eye(position-1)-1/n,-1/n*ones(position-1,n-position);-1/n*ones(1,n-1);-1/n*ones(n-position,position-1),eye(n-position)-1/n];
end

end