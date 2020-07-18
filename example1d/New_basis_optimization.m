
function [Phi,Psi] = New_basis_optimization(A,M,N_fine,N_coarse,opt,opt2);


% standard FEM basis (coarse grid)
switch opt
    case 'periodic'
    Phi = sparse(N_fine,N_coarse);
    for i = 1:N_coarse
        if (i==1)
            n = round(N_fine/N_coarse);
            Phi(1:n,i) = (1:n)'/n;
            Phi(n:2*n,i) = (n:-1:0)'/n;
        else
            Phi(:,i) = [Phi(end-n+1:end,i-1);Phi(1:end-n,i-1)];
        end;
    end;

    case 'zeroboundary'
        Phi = sparse(N_fine+1,N_coarse-1);
        for i = 1:N_coarse-1
            if (i==1)
                n = round(N_fine/N_coarse);
                Phi(1:n+1,i) = (0:n)'/n;
                Phi(n+1:2*n+1,i) = (n:-1:0)'/n;
            else
                Phi(:,i) = [Phi(end-n+1:end,i-1);Phi(1:end-n,i-1)];
            end;
        end;
end;


if (opt2==0)
    Psi = 0;
    return;
end;

% multiscale basis
Aeq = Phi' * M;
% basis
Psi = sparse(size(Phi,1),size(Phi,2));

options = optimoptions('quadprog','Display','off');
parfor i = 1:size(Psi,2)
    beq = zeros(size(Phi,2),1);
    beq(i) = 1;
    %fprintf('finish constructing basis %d\n',i);
    [Psi(:,i),fval,exitflag,output,Lag_multiplier] = quadprog(A,zeros(size(Phi,1),1),[],[],Aeq,beq,[],[],[],options);
end;

