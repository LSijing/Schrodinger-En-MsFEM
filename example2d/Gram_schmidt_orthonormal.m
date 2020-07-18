

function [V_orthobasis] = Gram_schmidt_orthonormal(V_basis,M);
%% Gram schmidt orthogonalization
% V_basis : each column is a FEM basis coefficients
% M       : mass matrix of fine FEM space Vh

for k = 1:size(V_basis,2)
    V_basis(:,k) = V_basis(:,k) - V_basis(:,1:k-1)*(V_basis(:,1:k-1)'*(M*V_basis(:,k)));
    V_basis(:,k) = V_basis(:,k) / sqrt(V_basis(:,k)'*M*V_basis(:,k));
end;

V_orthobasis = V_basis;