
%% ------assemble mass matrix-------%
%% 

function [A] = MassMat2D_sparse(p,t)

np = size(p,1);
nt = size(t,1);
n = round(sqrt(np))-1;

Isparse_stiff = reshape(repmat(t,1,3)',[],1);
Jsparse_stiff = reshape(repmat(reshape(t',[],1),1,3)',[],1);
Isparse_stiff = Index_transform_to_periodicBC(Isparse_stiff,n);
Jsparse_stiff = Index_transform_to_periodicBC(Jsparse_stiff,n);

x1_elem = p(t(:,1),1);
x2_elem = p(t(:,2),1);
x3_elem = p(t(:,3),1);
y1_elem = p(t(:,1),2);
y2_elem = p(t(:,2),2);
y3_elem = p(t(:,3),2);
area = polyarea([x1_elem';x2_elem';x3_elem'],[y1_elem';y2_elem';y3_elem'])';

AK = repmat([2 1 1 1 2 1 1 1 2],nt,1);

Asparse_stiff = AK.*repmat(area,1,9).*repmat(1/12,nt,9);
Asparse_stiff = reshape(Asparse_stiff',[],1);
A = sparse(Isparse_stiff,Jsparse_stiff,Asparse_stiff,n^2,n^2);

