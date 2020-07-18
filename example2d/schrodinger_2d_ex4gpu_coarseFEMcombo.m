clear all;

%% Schrodinger equation (time dependent potential)  i*epsilon*u_t = -epsilon^2 /2 u_xx + V1 * u + (V2(x)*U(t)) * u
%%   COARSE solution generator

%% problem setting
% i * epsilon * u_t = -epsilon^2 /2 u_xx + V1 * u + (V2(x)*U(t)) * u
% periodic boundary condition ！！！！(0,dx,2dx,3dx,...,(n-1)dx=1-dx)^2 (除掉1点算上0点)
% Gaussian initial function
epsilon = 1/8;
DiffCoef = @(x,y) epsilon^2/2;
PotentialV1 = @(x,y)  (   (sin(2*pi*x/epsilon)+1) * (cos(2*pi*y/epsilon)) ) * xor(x<0.5,y<0.5) + ...,
                      (   (sin(2*pi*x*(1/epsilon-2))) * (cos(2*pi*y*(1/epsilon-2))+1) ) * (1-xor(x<0.5,y<0.5)) ;
E0 = 20;
PotentialV2 = @(x,y) E0 * (x + y);
PotentialWt = @(t)  sin(2*pi*t); %exp(2*sin(2*pi*t)) -1;
sigma1 = 0.2;
IniFunc = @(x,y)  exp( - ((x-0.5)^2+(y-0.5)^2) / (4*sigma1^2) ) / (sqrt(2*pi)*sigma1);

%% spatial setting
% FEM fine mesh, stiffness, mass, Potential Matrix; 
N_fine = 3*2^7;
h = 1/N_fine;
[X_fine,Y_fine] = ndgrid((0:N_fine)/N_fine,(0:N_fine)/N_fine);
p = [X_fine(:) Y_fine(:)];
tri = delaunay(p);
e = freeBoundary(triangulation(tri,p));

% Robin matrix, Load vector, Robin vector (because mesh, Bd Cond, RHS are fixed)

A = StiffMat2D_sparse(DiffCoef,p,tri);
M = MassMat2D_sparse(p,tri);
P1 = PotMat2D_sparse(p,tri,PotentialV1);
P2 = PotMat2D_sparse(p,tri,PotentialV2);

%  initial Vector; (这里是FEM空间L2投影)
u_ini = M \ LoadVec2D_sparse(p,tri,IniFunc);

%  坐标向量含头(x=0)不含尾(x=1)
X_fine = X_fine(1:end-1,1:end-1);
Y_fine = Y_fine(1:end-1,1:end-1);
p = [X_fine(:) Y_fine(:)];


%% temporal setting, Crank-Nicolson
T = 1;
deltat = 1/(2^18);
t = (0:round(T/deltat))*deltat;
% 粗时间网格（指标数组）
gap = 1;
coarse_t_idx = 1:gap:length(t);
% 存snapshots的数量（和指标数组,时间数组）
n_snap = 2^6;
snap_t_idx = 1:floor((length(coarse_t_idx)-1)/n_snap):length(coarse_t_idx);
n_snap = length(snap_t_idx) - 1;
t_series = t(coarse_t_idx(snap_t_idx));


%% Various reduced basis solution

N_array = [64, 48, 32, 24, 16, 12];

for j = N_array
    j1 = find(N_array==j);
    N = j;
    N2 = round(j/4);
    time_record = 0;
    [Psi,Phi] = New_basis_gamblet_oldversion1(A+P1,M,N_fine,N,12);
    
    
    % combo
    [Psi2,Phi2] = New_basis_gamblet_oldversion1(A+P1+P2,M,N_fine,N2,12);
    Psi_combo = full([Psi,Psi2]);
    clear Psi2 Phi2 Phi Psi;
    Psi_combo = Gram_schmidt_orthonormal(Psi_combo,M);
    
    Ab_combo = gpuArray( Psi_combo' * (1i * epsilon * M - deltat/2 * (A + P1)) * Psi_combo );
    Af_combo = gpuArray( Psi_combo' * (1i * epsilon * M + deltat/2 * (A + P1)) * Psi_combo );
    P2_combo = gpuArray( Psi_combo' * P2 * Psi_combo );

    % initial data
    U_combo = gpuArray( (Psi_combo' * M * Psi_combo) \ (Psi_combo' * M * u_ini) );
    u_combo = zeros(N_fine^2,n_snap+1);
    u_combo(:,1) = Psi_combo * gather(U_combo);
    
    t_c = t(coarse_t_idx);
    for i = 1:length(coarse_t_idx)
        if (i>1)
            % coarse solution (evolution)
            tau = t_c(i) - t_c(i-1);
            t_m = ( t_c(i) + t_c(i-1) ) / 2;
            tic;
            U_combo = ( Ab_combo - tau/2 * P2_combo * PotentialWt(t_m) ) \ ( ( Af_combo + tau/2 * P2_combo * PotentialWt(t_m) )*U_combo );
            time_record = time_record + toc;
        end;
        if (ismember(i,snap_t_idx))
            q_i = find(snap_t_idx==i);
            u_combo(:,q_i) = Psi_combo * gather(U_combo);
            fprintf('%d  , N = %d ,  %.2f(hr.) \n',q_i,N,time_record/3600);
        end;
    end;

    fprintf('\n');
    clear Psi_combo;
    
    filename1 = sprintf('Ex4_coarseFEMcombo_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
    save(filename1,'epsilon','DiffCoef','PotentialV1','E0','PotentialV2','PotentialWt','sigma1','IniFunc','N_fine','-v7.3');
    save(filename1,'t_series','T','deltat','gap','snap_t_idx','coarse_t_idx','t','n_snap','-append');
    save(filename1,'N','N2','u_combo','time_record','-append');

end;
