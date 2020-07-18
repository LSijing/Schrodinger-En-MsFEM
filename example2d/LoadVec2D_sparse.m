
%% ------assemble right hand side-------%
%  equation : - Div( a * Gradient (u) ) = f
%  variational form : integral{a * Grad(phi_i) * sum{Grad(phi_j)*u_j} } = integral{f * phi_i}
%  vector 'LoadVec2D_sparse' here is right side 'int{f * phi_i}' in above

%%
function [b] = LoadVec2D_sparse(p,t,f)

np = size(p,1);
nt = size(t,1);
n = round(sqrt(np))-1;


x1_elem = p(t(:,1),1);
x2_elem = p(t(:,2),1);
x3_elem = p(t(:,3),1);
y1_elem = p(t(:,1),2);
y2_elem = p(t(:,2),2);
y3_elem = p(t(:,3),2);
area = polyarea([x1_elem';x2_elem';x3_elem'],[y1_elem';y2_elem';y3_elem'])';
xc_elem = (x1_elem + x2_elem + x3_elem) /3;
yc_elem = (y1_elem + y2_elem + y3_elem) /3;


Isparse_loadv = reshape(t',[],1);
Jsparse_loadv = ones(nt*3,1);
Isparse_loadv = Index_transform_to_periodicBC(Isparse_loadv,n);
% (center on three edges)
x_ec = [(x1_elem+x2_elem)/2  (x2_elem+x3_elem)/2  (x3_elem+x1_elem)/2];
y_ec = [(y1_elem+y2_elem)/2  (y2_elem+y3_elem)/2  (y3_elem+y1_elem)/2];
fbar = ( arrayfun(f,x_ec(:,1),y_ec(:,1)) + arrayfun(f,x_ec(:,2),y_ec(:,2)) + arrayfun(f,x_ec(:,3),y_ec(:,3)) )/3;

bK = fbar .* area /3;
bsparse_loadv = reshape(repmat(bK,1,3)',[],1);
b = sparse(Isparse_loadv,Jsparse_loadv,bsparse_loadv,n^2,1);

