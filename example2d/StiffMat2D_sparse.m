
%% ------assemble stiffness matrix-------%
%% equation : - Div( a * Gradient (u) ) = f
% variational form : integral{a * Grad(phi_i) * sum{Grad(phi_j)*u_j} } = integral{f * phi_i}

function [A] = StiffMat2D_sparse(a,p,t)

np = size(p,1);
n = round(sqrt(np))-1;
nt = size(t,1);

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
xc_elem = (x1_elem + x2_elem + x3_elem) /3;
yc_elem = (y1_elem + y2_elem + y3_elem) /3;
abar = arrayfun(a,xc_elem,yc_elem);
b_elem = [y2_elem-y3_elem  y3_elem-y1_elem  y1_elem-y2_elem] ./ repmat(area,1,3) /2;
c_elem = [x2_elem-x3_elem  x3_elem-x1_elem  x1_elem-x2_elem] ./ repmat(area,1,3) /2;
AK = [b_elem(:,1).*b_elem(:,1)+c_elem(:,1).*c_elem(:,1) ...,
    b_elem(:,2).*b_elem(:,1)+c_elem(:,2).*c_elem(:,1) ...,
    b_elem(:,3).*b_elem(:,1)+c_elem(:,3).*c_elem(:,1) ...,
      b_elem(:,1).*b_elem(:,2)+c_elem(:,1).*c_elem(:,2) ...,
    b_elem(:,2).*b_elem(:,2)+c_elem(:,2).*c_elem(:,2) ...,
    b_elem(:,3).*b_elem(:,2)+c_elem(:,3).*c_elem(:,2) ...,
      b_elem(:,1).*b_elem(:,3)+c_elem(:,1).*c_elem(:,3) ...,
    b_elem(:,2).*b_elem(:,3)+c_elem(:,2).*c_elem(:,3) ...,
    b_elem(:,3).*b_elem(:,3)+c_elem(:,3).*c_elem(:,3)];

Asparse_stiff = AK.*repmat(area,1,9).*repmat(abar,1,9);
Asparse_stiff = reshape(Asparse_stiff',[],1);
A = sparse(Isparse_stiff,Jsparse_stiff,Asparse_stiff,n^2,n^2);

