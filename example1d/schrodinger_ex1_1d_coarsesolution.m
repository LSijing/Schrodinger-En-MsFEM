clear all;

%% Schrodinger equation (time dependent potential)  i*epsilon*u_t = -epsilon^2 /2 u_xx + V1 * u + (V2(x)*U(t)) * u
%%   Coarse solution generator

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
gap = 4;
coarse_t_idx = 1:gap:length(t);
% 存snapshots的数量（和指标数组,时间数组）
n_snap = 2^6;
snap_t_idx = 1:floor((length(coarse_t_idx)-1)/n_snap):length(coarse_t_idx);
n_snap = length(snap_t_idx) - 1;
t_series = t(coarse_t_idx(snap_t_idx));

%% computing coarse solutions
N1_array = log([3*2^7, 2^8, 3*2^6, 2^7, 3*2^5, 2^6, 3*2^4, 2^5]) /log(2);
N2_array = max(N1_array-3,1);
time_record = zeros(3,length(N1_array));
time_description = {'coarse FEM','coarse MsFEM','coarse En MsFEM'};



for j = 1:length(N1_array)
    u_FEM = zeros(length(x_fine),n_snap+1);
    u_FEM(:,1) = u_ini;
    u_Ms = zeros(length(x_fine),n_snap+1);
    u_Ms(:,1) = u_ini;
    u_EnMs = zeros(length(x_fine),n_snap+1);
    u_EnMs(:,1) = u_ini;

    N1 = round(2^(N1_array(j)));
    N2 = round(2^(N2_array(j)));
    [Phi,Psi] = New_basis_optimization(A+P1,M,N_fine,N1,'periodic',12);
    
    % coarse FEM
    U_FEM = (Phi' * M * Phi) \ (Phi' * M * u_ini);
    M_FEM = Phi' * M * Phi;
    A1_FEM = Phi' * (A + P1) * Phi;
    A2_FEM = Phi' * P2 * Phi;
    
    % coarse MsFEM
    U_Ms = (Psi' * M * Psi) \ (Psi' * M * u_ini);
    M_Ms = Psi' * M * Psi;
    A1_Ms = Psi' * (A + P1) * Psi;
    A2_Ms = Psi' * P2 * Psi;

    % coarse En MsFEM
    [Phi2,Psi2] = New_basis_optimization(A+P1+P2,M,N_fine,N2,'periodic',12);
    Psi_EnMs = [Psi,Psi2];
    Psi_EnMs = Gram_schmidt_orthonormal(Psi_EnMs,M);

    U_EnMs = (Psi_EnMs' * M * Psi_EnMs) \ (Psi_EnMs' * M * u_ini);
    M_EnMs = Psi_EnMs' * M * Psi_EnMs;
    A1_EnMs = Psi_EnMs' * (A + P1) * Psi_EnMs;
    A2_EnMs = Psi_EnMs' * P2 * Psi_EnMs;
    
    %% time evolution 
    t_c = t(coarse_t_idx);
    q_i = 0; % 计数器
    for i = 1:length(coarse_t_idx)
        if (i>1)
            % coarse solution (evolution)
            tau = t_c(i) - t_c(i-1);
            tic;
            U_FEM = (1i * epsilon * M_FEM - tau/2 * A1_FEM - tau/2 * A2_FEM * PotentialWt(t_c(i)) ) \ ( (1i * epsilon * M_FEM + tau/2 * A1_FEM + tau/2 * A2_FEM * PotentialWt(t_c(i-1)) ) * U_FEM );
            time_record(1,j) = time_record(1,j) + toc;
            tic;
            U_Ms = (1i * epsilon * M_Ms - tau/2 * A1_Ms - tau/2 * A2_Ms * PotentialWt(t_c(i)) ) \ ( (1i * epsilon * M_Ms + tau/2 * A1_Ms + tau/2 * A2_Ms * PotentialWt(t_c(i-1)) ) * U_Ms );
            time_record(2,j) = time_record(2,j) + toc;
            tic;
            U_EnMs = (1i * epsilon * M_EnMs - tau/2 * A1_EnMs - tau/2 * A2_EnMs * PotentialWt(t_c(i)) ) \ ( (1i * epsilon * M_EnMs + tau/2 * A1_EnMs + tau/2 * A2_EnMs * PotentialWt(t_c(i-1)) ) * U_EnMs );                
            time_record(3,j) = time_record(3,j) + toc;
        end;
        %% error snapshots
        if ismember(i,snap_t_idx)
            q_i = find(snap_t_idx==i);
            
            u_FEM(:,q_i) = Phi * U_FEM;
            u_Ms(:,q_i) = Psi * U_Ms;
            u_EnMs(:,q_i) = Psi_EnMs * U_EnMs;
            
            if (mod(q_i-1,round(n_snap/8))==0)
                fprintf('%d,',q_i);
            end;
            
        end;
    end;
    filename1 = sprintf('Ex1_coarsesol_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N1,round(1/deltat),gap);
    save(filename1,'time_record','time_description','N1','N2','-v7.3');
    save(filename1,'u_FEM','u_Ms','u_EnMs','-append');
    fprintf('### N1=%d, N2=%d \n',N1,N2);
end;

