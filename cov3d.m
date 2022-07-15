function x = cov3d(x)
% Based on Matlab's cov, version 5.16.4.10

[m,n,p] = size(x);
if m == 1
    x = zeros(n,n,p,class(x));
else
    x = bsxfun(@minus,x,sum(x,1)/m);
    for i = 1:p
        xi = x(:,:,i);
        x(:,:,i) = xi'*xi;
    end
    x = x/(m-1);
end