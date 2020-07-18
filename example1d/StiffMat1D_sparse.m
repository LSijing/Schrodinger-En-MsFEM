%-------Assemble sparse stiffness matrix--------%
%-------mid-point integral---------------%

function [B] = StiffMat1D_sparse(x,a)
%----------x : node coordinate vector-----------------%
%----------a : a function-----------------------------%

if (size(x,1)==1) 
    x = x';
end;
n = length(x) -1;
xc = (x(2:end) + x(1:end-1)) /2;
h = x(2:end)-x(1:end-1);

I_sparse = reshape([ repmat([n 1:n-1],2,1) ; repmat(1:n,2,1) ],[],1);
J_sparse = reshape(repmat([n 1:n-1;(1:n)],2,1),[],1);
amid = arrayfun(a,xc);

A_sparse = reshape(repmat([1 -1 -1 1]',1,n) .* repmat((amid./h)',4,1),[],1);

B = sparse(I_sparse,J_sparse,A_sparse,n,n);

return;
