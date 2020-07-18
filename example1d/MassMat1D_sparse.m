%-------Assemble sparse mass matrix--------%

function [B] = MassMat1D_sparse(x)
%----------x : node coordinate vector-----------------%

if (size(x,1)==1) 
    x = x';
end;

n = length(x) -1;
h = x(2:end)-x(1:end-1);

I_sparse = reshape([ repmat([n 1:n-1],2,1) ; repmat(1:n,2,1) ],[],1);
J_sparse = reshape(repmat([n 1:n-1;(1:n)],2,1),[],1);
A_sparse = reshape(repmat([1/3 1/6 1/6 1/3]',1,n) .* repmat(h',4,1),[],1);

B = sparse(I_sparse,J_sparse,A_sparse,n,n);

return;