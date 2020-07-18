%------assemble sparse load vector-----%
%------trapezoidal rule------------------%

function [v] = LoadVec1D_sparse(x,f);
%-------f : right hand side function---------------------%

if (size(x,1)==1) 
    x = x';
end;
n = length(x) -1;
h = x(2:end)-x(1:end-1);
farray = arrayfun(f,x);

I_sparse = reshape([n 1:n-1;1:n],[],1);
J_sparse = ones(2*n,1);
V_sparse = reshape([(farray(1:end-1).*h/2)'; (farray(2:end).*h/2)'],[],1);

v = sparse(I_sparse,J_sparse,V_sparse,n,1);


