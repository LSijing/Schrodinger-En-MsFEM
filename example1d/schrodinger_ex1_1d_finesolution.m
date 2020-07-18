clear all;

%% Schrodinger equation (time dependent potential)  i*epsilon*u_t = -epsilon^2 /2 u_xx + V1 * u + (V2(x)*U(t)) * u
%%   Fine solution generator

%% problem setting

epsilon = 1/32;
DiffCoef = @(x) epsilon^2/2;
PotentialV1 = @(x) cos(2*pi*x/epsilon);
E0 = 20;
PotentialV2 = @(x) E0 * x;
PotentialWt = @(t) sin(2*pi*t); %exp(2*sin(2*pi*t))-1;
sigma1 = 0.2;
IniFunc = @(x)  exp( -(x-0.5)^2 / (4*sigma1^2) ) / sqrt((sqrt(2*pi)*sigma1));


%% spatial setting
% FEM fine mesh 
N_fine = 3*2^15;
x_fine = ((0:N_fine)/N_fine) * 1 -0;

% FEM fine mesh stiffness, mass, Potential Matrix; 
A = StiffMat1D_sparse(x_fine,DiffCoef);
M = MassMat1D_sparse(x_fine);
P1 = PotMat1D_sparse(x_fine,PotentialV1);
P2 = PotMat1D_sparse(x_fine,PotentialV2);

%  initial Vector;
u_ini = M \ (LoadVec1D_sparse(x_fine,IniFunc));
%  坐标向量含尾(x=1)不含头(x=0)
x_fine = x_fine(2:end);

%% temporal setting, Crank-Nicolson
T = 1;
deltat = 1/(2^20);
t = (0:round(T/deltat))*deltat;
% 粗时间网格（指标数组）
gap = 1;
coarse_t_idx = 1:gap:length(t);
% 存snapshots的数量（和指标数组,时间数组）
n_snap = 2^6;
snap_t_idx = 1:floor((length(coarse_t_idx)-1)/n_snap):length(coarse_t_idx);
n_snap = length(snap_t_idx) - 1;
t_series = t(coarse_t_idx(snap_t_idx));

%% saving data
filename1 = sprintf('Ex1_finesolution_eps1over%d_E0%d_T%.3f_h1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,round(1/deltat),gap);
save(filename1,'epsilon','DiffCoef','PotentialV1','E0','PotentialV2','PotentialWt','sigma1','IniFunc','N_fine','-v7.3');
save(filename1,'t_series','T','deltat','gap','snap_t_idx','coarse_t_idx','t','n_snap','-append');

%% reference exact solution
time_record_fine = 0;

u_fine = zeros(length(x_fine),n_snap+1);
u_fine(:,1) = u_ini;
u_evolve = u_ini;

epMplusdtAP1 = 1i * epsilon * M + deltat/2 * (A + P1);
epMminusdtAP1 = 1i * epsilon * M - deltat/2 * (A + P1);
dtP2 = deltat/2 * P2;

for i = 2:length(coarse_t_idx)
    t_f = t(coarse_t_idx(i-1):coarse_t_idx(i));
    for i_f = 2:length(t_f)
        tic;
        u_evolve = (epMminusdtAP1 - dtP2 * PotentialWt(t_f(i_f)) ) \ ( (epMplusdtAP1 + dtP2 * PotentialWt(t_f(i_f-1)) )* u_evolve );
        time_record_fine = time_record_fine + toc;
    end;
    if ismember(i,snap_t_idx)
        q_i = find(snap_t_idx==i);
        u_fine(:,q_i) = u_evolve;
        fprintf('%d\n',q_i);
    end;
end;

save(filename1,'time_record_fine','u_fine','-append');
