clear all;

%% plotting errors

%% load problem settings and fine solution 
line_width = 0.5;
marker_size = 6;

epsilon = 1/32;
E0 = 20;
T = 1;
N_fine = 3*2^15;
x_fine = ((0:N_fine)/N_fine) * 1 -0;
deltat = 1/2^20;

gap = 1;
filename1 = sprintf('Ex1_finesolution_eps1over%d_E0%d_T%.3f_h1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,round(1/deltat),gap);
load(filename1);

% FEM fine mesh stiffness, mass, Potential Matrix; 
A = StiffMat1D_sparse(x_fine,DiffCoef);
M = MassMat1D_sparse(x_fine);
H1 = A * 2 / epsilon^2 + M;
P1 = PotMat1D_sparse(x_fine,PotentialV1);
P2 = PotMat1D_sparse(x_fine,PotentialV2);

x_fine = x_fine(2:end);

density_fine = abs(u_fine).^2;
engdensity_fine = epsilon^2/2 * abs( ( (u_fine([2:end,1],:) - u_fine([end,1:end-1],:)) /(2/N_fine) ).^2) ...,
    + (repmat(arrayfun(PotentialV1,x_fine)',1,n_snap+1) + repmat(arrayfun(PotentialV2,x_fine)',1,n_snap+1) .* repmat(PotentialWt(t_series),N_fine,1) ) .* density_fine;



%% load coarse mesh solution and compute final time errors
wave_fine_end = u_fine(:,end);
density_fine_end = density_fine(:,end);
engdensity_fine_end = engdensity_fine(:,end);

Num = 6;
H_coarse = zeros(Num,1);
N_coarse = [3*2^8, 2^9, 3*2^7, 2^8, 3*2^6, 2^7, 3*2^5, 2^6, 3*2^4, 2^5];
gap = 4;
L2err_FEM     =  zeros(Num,1);
L2err_gamblet =  zeros(Num,1);
L2err_combo   =  zeros(Num,1);
H1err_FEM     =  zeros(Num,1);
H1err_gamblet =  zeros(Num,1);
H1err_combo   =  zeros(Num,1);
density_L2err_FEM     =  zeros(Num,1);
density_L2err_gamblet =  zeros(Num,1);
density_L2err_combo   =  zeros(Num,1);
engdensity_L2err_FEM     =  zeros(Num,1);
engdensity_L2err_gamblet =  zeros(Num,1);
engdensity_L2err_combo   =  zeros(Num,1);


for i = 1:Num
filename2 = sprintf('Ex1_coarsesol_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N_coarse(i),round(1/deltat),gap);
load(filename2);

wave_FEM_end = u_FEM(:,end);
density_FEM_end = abs(wave_FEM_end).^2;
engdensity_FEM_end = epsilon^2/2 * abs( ( (wave_FEM_end([2:end,1]) - wave_FEM_end([end,1:end-1])) /(2/N_fine) ).^2) ...,
    + (arrayfun(PotentialV1,x_fine)' + arrayfun(PotentialV2,x_fine)' .* repmat(PotentialWt(t_series(end)),N_fine,1) ) .* density_FEM_end;

wave_gamblet_end = u_gamblet(:,end);
density_gamblet_end = abs(wave_gamblet_end).^2;
engdensity_gamblet_end = epsilon^2/2 * abs( ( (wave_gamblet_end([2:end,1]) - wave_gamblet_end([end,1:end-1])) /(2/N_fine) ).^2) ...,
    + (arrayfun(PotentialV1,x_fine)' + arrayfun(PotentialV2,x_fine)' .* repmat(PotentialWt(t_series(end)),N_fine,1) ) .* density_gamblet_end;

wave_combo_end = u_combo(:,end);
density_combo_end = abs(wave_combo_end).^2;
engdensity_combo_end = epsilon^2/2 * abs( ( (wave_combo_end([2:end,1]) - wave_combo_end([end,1:end-1])) /(2/N_fine) ).^2) ...,
    + (arrayfun(PotentialV1,x_fine)' + arrayfun(PotentialV2,x_fine)' .* repmat(PotentialWt(t_series(end)),N_fine,1) ) .* density_combo_end;

H_coarse(i) = 1 / N1 ;
L2err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * M * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
L2err_gamblet(i)  =  sqrt(diag( (wave_gamblet_end - wave_fine_end)' * M * (wave_gamblet_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
L2err_combo(i)    =  sqrt(diag( (wave_combo_end - wave_fine_end)' * M * (wave_combo_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * M * wave_fine_end  ));
H1err_FEM(i)      =  sqrt(diag( (wave_FEM_end - wave_fine_end)' * H1 * (wave_FEM_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
H1err_gamblet(i)  =  sqrt(diag( (wave_gamblet_end - wave_fine_end)' * H1 * (wave_gamblet_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
H1err_combo(i)    =  sqrt(diag( (wave_combo_end - wave_fine_end)' * H1 * (wave_combo_end - wave_fine_end)  )) ./ sqrt(diag( wave_fine_end' * H1 * wave_fine_end  ));
density_L2err_FEM(i)      =  sqrt(diag( (density_FEM_end - density_fine_end)' * M * (density_FEM_end - density_fine_end)  )) ./ sqrt(diag( density_fine_end' * M * density_fine_end  ));
density_L2err_gamblet(i)  =  sqrt(diag( (density_gamblet_end - density_fine_end)' * M * (density_gamblet_end - density_fine_end)  )) ./ sqrt(diag( density_fine_end' * M * density_fine_end  ));
density_L2err_combo(i)    =  sqrt(diag( (density_combo_end - density_fine_end)' * M * (density_combo_end - density_fine_end)  )) ./ sqrt(diag( density_fine_end' * M * density_fine_end  ));
engdensity_L2err_FEM(i)      =  sqrt(diag( (engdensity_FEM_end - engdensity_fine_end)' * M * (engdensity_FEM_end - engdensity_fine_end)  )) ./ sqrt(diag( engdensity_fine_end' * M * engdensity_fine_end  ));
engdensity_L2err_gamblet(i)  =  sqrt(diag( (engdensity_gamblet_end - engdensity_fine_end)' * M * (engdensity_gamblet_end - engdensity_fine_end)  )) ./ sqrt(diag( engdensity_fine_end' * M * engdensity_fine_end  ));
engdensity_L2err_combo(i)    =  sqrt(diag( (engdensity_combo_end - engdensity_fine_end)' * M * (engdensity_combo_end - engdensity_fine_end)  )) ./ sqrt(diag( engdensity_fine_end' * M * engdensity_fine_end  ));
end;





%% final time errors


figure(1);
hold on;
plot(H_coarse,abs(L2err_FEM),'k--o');
plot(H_coarse,abs(L2err_gamblet),'r--*');
plot(H_coarse,abs(L2err_combo),'b--^');
xlim([H_coarse(1)/2 ,  H_coarse(end)*2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
%set(gca,'Ygrid','on');
xlabel('H');
ylabel('relative L^2 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
box on;

figure(2);
hold on;
plot(H_coarse,abs(H1err_FEM),'k--o');
plot(H_coarse,abs(H1err_gamblet),'r--*');
plot(H_coarse,abs(H1err_combo),'b--^');
xlim([H_coarse(1)/2 ,  H_coarse(end)*2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('H');
ylabel('relative H^1 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
%grid on;
box on;

figure(3);
hold on;
plot(H_coarse,density_L2err_FEM,'k--o');
plot(H_coarse,density_L2err_gamblet,'r--*');
plot(H_coarse,density_L2err_combo,'b--^');
xlim([H_coarse(1)/2 ,  H_coarse(end)*2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('H');
ylabel('relative L^2 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
%grid on;
box on;

figure(4);
hold on;
plot(H_coarse,engdensity_L2err_FEM,'k--o');
plot(H_coarse,engdensity_L2err_gamblet,'r--*');
plot(H_coarse,engdensity_L2err_combo,'b--^');
xlim([H_coarse(1)/2 ,  H_coarse(end)*2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('H');
ylabel('relative L^2 error');
leg1 = legend('FEM','MsFEM','En-MsFEM','location','southeast');
set(leg1,'fontsize',11);
%grid on;
box on;



%% load a coarse mesh solution and compute time-series errors

N1 = round(1/epsilon)*8;
filename3 = sprintf('Ex1_coarsesol_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N1,round(1/deltat),gap);
load(filename3);

density_FEM = abs(u_FEM).^2;
engdensity_FEM = epsilon^2/2 * abs( ( (u_FEM([2:end,1],:) - u_FEM([end,1:end-1],:)) /(2/N_fine) ).^2) ...,
    + (repmat(arrayfun(PotentialV1,x_fine)',1,n_snap+1) + repmat(arrayfun(PotentialV2,x_fine)',1,n_snap+1) .* repmat(PotentialWt(t_series),N_fine,1) ) .* density_FEM;

density_gamblet = abs(u_gamblet).^2;
engdensity_gamblet = epsilon^2/2 * abs( ( (u_gamblet([2:end,1],:) - u_gamblet([end,1:end-1],:)) /(2/N_fine) ).^2) ...,
    + (repmat(arrayfun(PotentialV1,x_fine)',1,n_snap+1) + repmat(arrayfun(PotentialV2,x_fine)',1,n_snap+1) .* repmat(PotentialWt(t_series),N_fine,1) ) .* density_gamblet;

density_combo = abs(u_combo).^2;
engdensity_combo = epsilon^2/2 * abs( ( (u_combo([2:end,1],:) - u_combo([end,1:end-1],:)) /(2/N_fine) ).^2) ...,
    + (repmat(arrayfun(PotentialV1,x_fine)',1,n_snap+1) + repmat(arrayfun(PotentialV2,x_fine)',1,n_snap+1) .* repmat(PotentialWt(t_series),N_fine,1) ) .* density_combo;


L2err_gamblet_series = sqrt(diag( (u_gamblet - u_fine)' * M * (u_gamblet - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
L2err_combo_series = sqrt(diag( (u_combo - u_fine)' * M * (u_combo - u_fine)  )) ./ sqrt(diag( u_fine' * M * u_fine  ));
H1err_gamblet_series = sqrt(diag( (u_gamblet - u_fine)' * H1 * (u_gamblet - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
H1err_combo_series = sqrt(diag( (u_combo - u_fine)' * H1 * (u_combo - u_fine)  )) ./ sqrt(diag( u_fine' * H1 * u_fine  ));
density_L2err_gamblet_series = sqrt(diag( (density_gamblet - density_fine)' * M * (density_gamblet - density_fine)  )) ./ sqrt(diag( density_fine' * M * density_fine  ));
density_L2err_combo_series = sqrt(diag( (density_combo - density_fine)' * M * (density_combo - density_fine)  )) ./ sqrt(diag( density_fine' * M * density_fine  ));
engdensity_L2err_gamblet_series = sqrt(diag( (engdensity_gamblet - engdensity_fine)' * M * (engdensity_gamblet - engdensity_fine)  )) ./ sqrt(diag( engdensity_fine' * M * engdensity_fine  ));
engdensity_L2err_combo_series = sqrt(diag( (engdensity_combo - engdensity_fine)' * M * (engdensity_combo - engdensity_fine)  )) ./ sqrt(diag( engdensity_fine' * M * engdensity_fine  ));


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

figure(103);
hold on;
plot(t_series,abs(density_L2err_gamblet_series),'r--');
plot(t_series,abs(density_L2err_combo_series),'b-');
xlabel('t');
ylabel('relative L^2 error');
leg1 = legend('MsFEM','En-MsFEM','location','northwest');
set(leg1,'fontsize',11);
box on;

figure(104);
hold on;
plot(t_series,abs(engdensity_L2err_gamblet_series),'r--');
plot(t_series,abs(engdensity_L2err_combo_series),'b-');
xlabel('t');
ylabel('relative L^2 error');
leg1 = legend('MsFEM','En-MsFEM','location','northwest');
set(leg1,'fontsize',11);
box on;
