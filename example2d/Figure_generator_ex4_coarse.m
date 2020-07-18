clear all;


%% load problem settings and fine solution 
line_width = 0.5;
marker_size = 6;

epsilon = 1/8;
E0 = 20;
T = 1;
N_fine = 3*2^7;
x_fine = ((0:N_fine)/N_fine) * 1 -0;
deltat = 1/2^18;

gap = 1;
N = 96;
filename1 = sprintf('Ex4_coarseEnMs_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename1);
u_fine = u_EnMs;

[X_fine,Y_fine] = ndgrid((0:N_fine)/N_fine,(0:N_fine)/N_fine);
p = [X_fine(:) Y_fine(:)];
tri = delaunay(p);
e = freeBoundary(triangulation(tri,p));


% Robin matrix, Load vector, Robin vector (because mesh, Bd Cond, RHS are fixed)

A = StiffMat2D_sparse(DiffCoef,p,tri);
M = MassMat2D_sparse(p,tri);
H1 = A * 2 / epsilon^2 + M;
P1 = PotMat2D_sparse(p,tri,PotentialV1);
P2 = PotMat2D_sparse(p,tri,PotentialV2);

%  坐标向量含头(x=0)不含尾(x=1)
X_fine = X_fine(1:end-1,1:end-1);
Y_fine = Y_fine(1:end-1,1:end-1);
p = [X_fine(:) Y_fine(:)];

%% plot the potential

figure(7);
v1 = reshape(arrayfun(PotentialV1,p(:,1),p(:,2)),size(X_fine));
%trimesh(tri,p(:,1),p(:,2),arrayfun(PotentialV1,p(:,1),p(:,2)));
%trisurf(tri,p(:,1),p(:,2),arrayfun(PotentialV1,p(:,1),p(:,2)),'edgecolor','none');
surf(X_fine,Y_fine,v1,'edgecolor','none');
xlabel('x');
ylabel('y');
colormap parula;
colorbar;
view(0,90);


%% load coarse mesh solution and compute final time errors
wave_fine_end = u_fine(:,end);

Num = 4;
H_coarse = zeros(Num,1);
N_coarse = [48, 32, 24, 16];
gap = 1;

L2err_FEM     =  zeros(Num,1);
L2err_Ms =  zeros(Num,1);
L2err_EnMs   =  zeros(Num,1);
H1err_FEM     =  zeros(Num,1);
H1err_Ms =  zeros(Num,1);
H1err_EnMs   =  zeros(Num,1);


for i = 1:Num
filename2 = sprintf('Ex4_coarseFEMMs_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N_coarse(i),round(1/deltat),gap);
load(filename2);

wave_FEM_end = u_FEM(:,end);
wave_Ms_end = u_Ms(:,end);


H_coarse(i) = 1 / N ;

L2err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * M * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
L2err_Ms(i)  =  sqrt(diag( (wave_Ms_end - wave_fine_end)' * M * (wave_Ms_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
H1err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * H1 * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
H1err_Ms(i)  =  sqrt(diag( (wave_Ms_end - wave_fine_end)' * H1 * (wave_Ms_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));


filename3 = sprintf('Ex4_coarseEnMs_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N_coarse(i),round(1/deltat),gap);
load(filename3);
wave_EnMs_end = u_EnMs(:,end);

L2err_EnMs(i)    =  sqrt(diag( (wave_EnMs_end - wave_fine_end)' * M * (wave_EnMs_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
H1err_EnMs(i)    =  sqrt(diag( (wave_EnMs_end - wave_fine_end)' * H1 * (wave_EnMs_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
end;




%% errors

figure(1);
hold on;
plot(H_coarse,abs(L2err_FEM(:,end)),'k--o');
plot(H_coarse,abs(L2err_Ms(:,end)),'r--*');
plot(H_coarse,abs(L2err_EnMs(:,end)),'b--^');
xlim([1e-2 ,  1e-1]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('H');
ylabel('relative L^2 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
%grid on;
box on;

figure(2);
hold on;
plot(H_coarse,abs(H1err_FEM(:,end)),'k--o');
plot(H_coarse,abs(H1err_Ms(:,end)),'r--*');
plot(H_coarse,abs(H1err_EnMs(:,end)),'b--^');
xlim([1e-2 ,  1e-1]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('H');
ylabel('relative H^1 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
%grid on;
box on;



%%  load a coarse mesh solution and compute time-series errors
Mass_fine_series = diag(u_fine' * M * u_fine);
Energy_fine_series = zeros(size(t_series));
for i = 1:length(t_series)
    Energy_fine_series(i) = u_fine(:,i)' * ( A+P1+ P2*PotentialWt(t_series(i)) ) * u_fine(:,i);
end;

N = round(1/epsilon) * 4;
filename5 = sprintf('Ex4_coarseFEMMs_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename5);

L2err_Ms_series = sqrt(diag( (u_Ms - u_fine)' * M * (u_Ms - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
H1err_Ms_series = sqrt(diag( (u_Ms - u_fine)' * H1 * (u_Ms - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
Mass_Ms_series = diag(u_Ms' * M * u_Ms);
Energy_Ms_series = zeros(size(t_series));
for i = 1:length(t_series)
    Energy_Ms_series(i) = u_Ms(:,i)' * ( A+P1+ P2*PotentialWt(t_series(i)) ) * u_Ms(:,i);
end;

filename6 = sprintf('Ex4_coarseEnMs_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename6);

L2err_EnMs_series = sqrt(diag( (u_EnMs - u_fine)' * M * (u_EnMs - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
H1err_EnMs_series = sqrt(diag( (u_EnMs - u_fine)' * H1 * (u_EnMs - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
Mass_EnMs_series = diag( u_EnMs' * M * u_EnMs );
Energy_EnMs_series = zeros(size(t_series));
for i = 1:length(t_series)
    Energy_EnMs_series(i) = u_EnMs(:,i)' * ( A+P1+ P2*PotentialWt(t_series(i)) ) * u_EnMs(:,i);
end;


figure(101);
hold on;
plot(t_series,abs(L2err_Ms_series),'r--');
plot(t_series,abs(L2err_EnMs_series),'b-');
xlabel('t');
ylabel('relative L^2 error');
leg1 = legend('MsFEM','En-MsFEM','location','northwest');
set(leg1,'fontsize',11);
box on;

figure(102);
hold on;
plot(t_series,abs(H1err_Ms_series),'r--');
plot(t_series,abs(H1err_EnMs_series),'b-');
xlabel('t');
ylabel('relative H^1 error');
leg1 = legend('MsFEM','En-MsFEM','location','northwest');
set(leg1,'fontsize',11);
box on;


%% plot some density/engdensity profiles
idx = 65;
density_fine = abs(u_fine(:,idx)).^2;
engdensity_fine = epsilon^2/2 * abs(reshape(gradient(reshape(u_fine(:,idx),size(X_fine)),1/N_fine),[],1)).^2 +  ...,
                        (arrayfun(PotentialV1,p(:,1),p(:,2)) + arrayfun(PotentialV2,p(:,1),p(:,2))*PotentialWt(t_series(idx)) ).*density_fine;
density_Ms = abs(u_Ms(:,idx)).^2;
engdensity_Ms = epsilon^2/2 * abs(reshape(gradient(reshape(u_Ms(:,idx),size(X_fine)),1/N_fine),[],1)).^2 +  ...,
                        (arrayfun(PotentialV1,p(:,1),p(:,2)) + arrayfun(PotentialV2,p(:,1),p(:,2))*PotentialWt(t_series(idx)) ).*density_Ms;
density_EnMs = abs(u_EnMs(:,idx)).^2;
engdensity_EnMs = epsilon^2/2 * abs(reshape(gradient(reshape(u_EnMs(:,idx),size(X_fine)),1/N_fine),[],1)).^2 +  ...,
                        (arrayfun(PotentialV1,p(:,1),p(:,2)) + arrayfun(PotentialV2,p(:,1),p(:,2))*PotentialWt(t_series(idx)) ).*density_EnMs;




figure(300+idx);
climit = [0 ,  max( [max(density_fine),max(density_Ms),max(density_EnMs)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_Ms,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_EnMs,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
subplot('position',[0.65 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.94 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);


figure(400+idx);
climit = [min( [min(engdensity_fine),min(engdensity_Ms),min(engdensity_EnMs)] ) ,  max( [max(engdensity_fine),max(engdensity_Ms),max(engdensity_EnMs)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_Ms,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_EnMs,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
subplot('position',[0.65 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.94 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);


figure(500+idx);
climit = [min([min(density_Ms-density_fine),min(density_EnMs-density_fine)]) ,  max( [max(density_Ms-density_fine),max(density_EnMs-density_fine)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_Ms-density_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_EnMs-density_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.64 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);


figure(600+idx);
climit = [min( [min(engdensity_Ms-engdensity_fine),min(engdensity_EnMs-engdensity_fine)] ) ,  max( [max(engdensity_Ms-engdensity_fine),max(engdensity_EnMs-engdensity_fine)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_Ms-engdensity_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_EnMs-engdensity_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.64 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);
