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
time_record = zeros(2,length(N_array));

for j = N_array
    j1 = find(N_array==j);
    N = j;
    N2 = round(j/4);
    [Psi,Phi] = New_basis_gamblet_oldversion1(A+P1,M,N_fine,N,12);
    Psi = full(Psi);
    
    % coarse FEM
    Ab_FEM = gpuArray( Phi' * (1i * epsilon * M - deltat/2 * (A + P1)) * Phi );
    Af_FEM = gpuArray( Phi' * (1i * epsilon * M + deltat/2 * (A + P1)) * Phi );
    P2_FEM = gpuArray( Phi' * P2 * Phi );
    
    % gamblet 
    Ab_gamblet = gpuArray( (Psi' * (1i * epsilon * M - deltat/2 * (A + P1))) * Psi );
    Af_gamblet = gpuArray( (Psi' * (1i * epsilon * M + deltat/2 * (A + P1))) * Psi );
    P2_gamblet = gpuArray( (Psi' * P2) * Psi );
    
    
    % initial data
    U_FEM = gpuArray( (Phi' * M * Phi) \ (Phi' * M * u_ini) );
    U_gamblet = gpuArray( (Psi' * M * Psi) \ (Psi' * M * u_ini) );
    u_FEM = zeros(N_fine^2,n_snap+1);
    u_FEM(:,1) = Phi * gather(U_FEM);
    u_gamblet = zeros(N_fine^2,n_snap+1);
    u_gamblet(:,1) = Psi * gather(U_gamblet);
    
    t_c = t(coarse_t_idx);
    for i = 1:length(coarse_t_idx)
        if (i>1)
            % coarse solution (evolution)
            tau = t_c(i) - t_c(i-1);
            t_m = ( t_c(i) + t_c(i-1) ) / 2;
            tic;
            [U_FEM,ccc1,ccc2,ccc3,ccc4] = gmres(Ab_FEM - tau/2 * P2_FEM * PotentialWt(t_m) , ( Af_FEM + tau/2 * P2_FEM * PotentialWt(t_m) )*U_FEM , 10 , gap/2*1e-12);
            time_record(1,j1) = time_record(1,j1) + toc;
            tic;
            U_gamblet = ( Ab_gamblet - tau/2 * P2_gamblet * PotentialWt(t_m) ) \ ( ( Af_gamblet + tau/2 * P2_gamblet * PotentialWt(t_m) )*U_gamblet );
            time_record(2,j1) = time_record(2,j1) + toc;
        end;
        if (ismember(i,snap_t_idx))
            q_i = find(snap_t_idx==i);
            u_FEM(:,q_i) = Phi * gather(U_FEM);
            u_gamblet(:,q_i) = Psi * gather(U_gamblet);
            fprintf('%d   ,%.2f(hr.), %.2f(hr.),\n',q_i,time_record(1,j1)/3600,time_record(2,j1)/3600);
        end;
    end;

    fprintf('\n');
    clear Psi Phi;
    
    filename1 = sprintf('Ex4_coarseFEMgamblet_eps1over%d_E0%d_T%.3f_h1over%d_H1over%d_dt1over%d_gap%d.mat',round(1/epsilon),E0,T,N_fine,N,round(1/deltat),gap);
    save(filename1,'epsilon','DiffCoef','PotentialV1','E0','PotentialV2','PotentialWt','sigma1','IniFunc','N_fine','-v7.3');
    save(filename1,'t_series','T','deltat','gap','snap_t_idx','coarse_t_idx','t','n_snap','-append');
    save(filename1,'N','N2','u_FEM','u_gamblet','time_record','-append');

end;
