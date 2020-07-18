%-------Assemble sparse potential matrix--------%
%-------second order quatrature rule --------%

function [C] = PotMat1D_sparse(x,f)
%----------x : node coordinate vector-----------------%

if (size(x,1)==1) 
    x = x';
end;

n = length(x) -1;
xc = (x(2:end) + x(1:end-1)) /2;
h = x(2:end)-x(1:end-1);
fc = arrayfun(f,xc);

I_sparse = reshape([ repmat([n 1:n-1],2,1) ; repmat(1:n,2,1) ],[],1);
J_sparse = reshape(repmat([n 1:n-1;(1:n)],2,1),[],1);
C_sparse = reshape(repmat([1/3 1/6 1/6 1/3]',1,n) .* repmat(fc' .* h',4,1),[],1);

C = sparse(I_sparse,J_sparse,C_sparse,n,n);

return;