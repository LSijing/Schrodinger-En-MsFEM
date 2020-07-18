

function [I_peroid] = Index_transform_to_periodicBC(I,n);
%% transform (n+1)*(n+1) indices to n*n indices

% 1d index ---> 2d index
j = mod(I-1,n+1)+1;
i = ceil((I/(n+1)));

% [1:n+1] ----> [1:n,1]
j = mod(j-1,n)+1;
i = mod(i-1,n)+1;

% 2d index ---> 1d index
I_peroid = (i-1)*n + j;