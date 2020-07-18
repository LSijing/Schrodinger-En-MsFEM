
function [Psi,Phi] = New_basis_optim(A,M,N_fine,N_coarse,opt2);

% standard FEM basis (coarse grid)
Phi = sparse(N_fine^2,N_coarse^2);
[X_fine,Y_fine] = ndgrid((0:N_fine)/N_fine,(0:N_fine)/N_fine);
X_fine = X_fine(1:end-1,1:end-1);
n = round(N_fine/N_coarse);
for j = 1:N_coarse
    for i = 1:N_coarse
        if ((i==1) && (j==1))
            temp = zeros(size(X_fine));
            for k = 0:n-1
                for m = 0: n-k-1
                    temp(i*n+1+k,j*n+1+m) = 1 - (k+m)/n;
                    temp(i*n+1-k,j*n+1-m) = 1 - (k+m)/n;
                end;
                for m = 1:n-1
                    temp(i*n+1+k,j*n+1-m) = 1 - max(k,m)/n;
                    temp(i*n+1-k,j*n+1+m) = 1 - max(k,m)/n;
                end;
            end;
        else
            temp = reshape(Phi(:,1),N_fine,N_fine);
            temp = [temp(n*(N_coarse-i+1)+1:end,:) ; temp(1:n*(N_coarse-i+1),:)];
            temp = [temp(:,n*(N_coarse-j+1)+1:end) , temp(:,1:n*(N_coarse-j+1))];
        end;
        Phi(:,(j-1)*(N_coarse)+i) = reshape(temp,[],1);
    end;
end;

if (opt2==0)
    Psi = 0;
    return;
end;

% gamblet basis
Aeq = Phi' * M;
% basis
Psi = sparse(size(Phi,1),size(Phi,2));

options = optimoptions('quadprog','Display','off');

tic;
parfor i = 1:size(Psi,2)
    beq = zeros(size(Phi,2),1);
    beq(i) = 1;
    [Psi(:,i),fval,exitflag,output,Lag_multiplier] = quadprog(A,zeros(N_fine^2,1),[],[],Aeq,beq,[],[],[],options);
    %fprintf('finish constructing basis %d\n',i);
end;
toc;


