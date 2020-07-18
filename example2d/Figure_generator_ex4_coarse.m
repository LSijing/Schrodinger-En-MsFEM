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
N = 64;
filename1 = sprintf('Ex4_coarseFEMcombo_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename1);
u_fine = u_combo;

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
%{
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
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'ex4_potentialV1','-dpdf');
close;
%}

%}

%% load coarse mesh solution and compute final time errors
wave_fine_end = u_fine(:,end);

Num = 4;
H_coarse = zeros(Num,1);
N_coarse = [48, 32, 24, 16];
gap = 1;

L2err_FEM     =  zeros(Num,1);
L2err_gamblet =  zeros(Num,1);
L2err_combo   =  zeros(Num,1);
H1err_FEM     =  zeros(Num,1);
H1err_gamblet =  zeros(Num,1);
H1err_combo   =  zeros(Num,1);


for i = 1:Num
filename2 = sprintf('Ex4_coarseFEMgamblet_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N_coarse(i),round(1/deltat),gap);
load(filename2);

wave_FEM_end = u_FEM(:,end);
wave_gamblet_end = u_gamblet(:,end);


H_coarse(i) = 1 / N ;

L2err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * M * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
L2err_gamblet(i)  =  sqrt(diag( (wave_gamblet_end - wave_fine_end)' * M * (wave_gamblet_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
H1err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * H1 * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
H1err_gamblet(i)  =  sqrt(diag( (wave_gamblet_end - wave_fine_end)' * H1 * (wave_gamblet_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));


filename3 = sprintf('Ex4_coarseFEMcombo_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N_coarse(i),round(1/deltat),gap);
load(filename3);
wave_combo_end = u_combo(:,end);

L2err_combo(i)    =  sqrt(diag( (wave_combo_end - wave_fine_end)' * M * (wave_combo_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
H1err_combo(i)    =  sqrt(diag( (wave_combo_end - wave_fine_end)' * H1 * (wave_combo_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
end;




%% errors

figure(1);
hold on;
plot(H_coarse,abs(L2err_FEM(:,end)),'k--o');
plot(H_coarse,abs(L2err_gamblet(:,end)),'r--*');
plot(H_coarse,abs(L2err_combo(:,end)),'b--^');
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
plot(H_coarse,abs(H1err_gamblet(:,end)),'r--*');
plot(H_coarse,abs(H1err_combo(:,end)),'b--^');
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
filename5 = sprintf('Ex4_coarseFEMgamblet_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename5);

L2err_gamblet_series = sqrt(diag( (u_gamblet - u_fine)' * M * (u_gamblet - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
H1err_gamblet_series = sqrt(diag( (u_gamblet - u_fine)' * H1 * (u_gamblet - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
Mass_gamblet_series = diag(u_gamblet' * M * u_gamblet);
Energy_gamblet_series = zeros(size(t_series));
for i = 1:length(t_series)
    Energy_gamblet_series(i) = u_gamblet(:,i)' * ( A+P1+ P2*PotentialWt(t_series(i)) ) * u_gamblet(:,i);
end;

filename6 = sprintf('Ex4_coarseFEMcombo_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
load(filename6);

L2err_combo_series = sqrt(diag( (u_combo - u_fine)' * M * (u_combo - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
H1err_combo_series = sqrt(diag( (u_combo - u_fine)' * H1 * (u_combo - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
Mass_combo_series = diag( u_combo' * M * u_combo );
Energy_combo_series = zeros(size(t_series));
for i = 1:length(t_series)
    Energy_combo_series(i) = u_combo(:,i)' * ( A+P1+ P2*PotentialWt(t_series(i)) ) * u_combo(:,i);
end;


figure(101);
hold on;
plot(t_series,abs(L2err_gamblet_series),'r--');
plot(t_series,abs(L2err_combo_series),'b-');
xlabel('t');
ylabel('relative L^2 error');
leg1 = legend('MsFEM','En-MsFEM','location','northwest');
set(leg1,'fontsize',11);
box on;

figure(102);
hold on;
plot(t_series,abs(H1err_gamblet_series),'r--');
plot(t_series,abs(H1err_combo_series),'b-');
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
density_gamblet = abs(u_gamblet(:,idx)).^2;
engdensity_gamblet = epsilon^2/2 * abs(reshape(gradient(reshape(u_gamblet(:,idx),size(X_fine)),1/N_fine),[],1)).^2 +  ...,
                        (arrayfun(PotentialV1,p(:,1),p(:,2)) + arrayfun(PotentialV2,p(:,1),p(:,2))*PotentialWt(t_series(idx)) ).*density_gamblet;
density_combo = abs(u_combo(:,idx)).^2;
engdensity_combo = epsilon^2/2 * abs(reshape(gradient(reshape(u_combo(:,idx),size(X_fine)),1/N_fine),[],1)).^2 +  ...,
                        (arrayfun(PotentialV1,p(:,1),p(:,2)) + arrayfun(PotentialV2,p(:,1),p(:,2))*PotentialWt(t_series(idx)) ).*density_combo;




figure(300+idx);
climit = [0 ,  max( [max(density_fine),max(density_gamblet),max(density_combo)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_gamblet,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_combo,size(X_fine)),'edgecolor','none');
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
climit = [min( [min(engdensity_fine),min(engdensity_gamblet),min(engdensity_combo)] ) ,  max( [max(engdensity_fine),max(engdensity_gamblet),max(engdensity_combo)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_gamblet,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_combo,size(X_fine)),'edgecolor','none');
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
climit = [min([min(density_gamblet-density_fine),min(density_combo-density_fine)]) ,  max( [max(density_gamblet-density_fine),max(density_combo-density_fine)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_gamblet-density_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(density_combo-density_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.64 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);


figure(600+idx);
climit = [min( [min(engdensity_gamblet-engdensity_fine),min(engdensity_combo-engdensity_fine)] ) ,  max( [max(engdensity_gamblet-engdensity_fine),max(engdensity_combo-engdensity_fine)] ) ];
subplot('position',[0.05 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_gamblet-engdensity_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
caxis(climit);
%view(0,90);
subplot('position',[0.35 0.15 0.25 0.8]);
surf(X_fine,Y_fine,reshape(engdensity_combo-engdensity_fine,size(X_fine)),'edgecolor','none');
xlabel('x');
ylabel('y');
%view(0,90);
caxis(climit);
colorbar('position',[0.64 0.15 0.025 0.8]);
set(gcf,'position',[100 100 1000 300]);
