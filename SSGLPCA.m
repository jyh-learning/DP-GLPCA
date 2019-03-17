function [P,obj,err,iter,Q] = SSGLPCA(L,D,lambda,k,opts)% lambda is alpha


% 

tol = 1e-4; 
max_iter = 500;
rho = 1.05;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


n = size(L,1);
P = zeros(n);
Q = P;
Y = P;

iter = 0;
for iter = 1 : max_iter
    Pk = P;
    Qk = Q;
  
    [P]=olver_H_sub(D,Q,L,Y,mu,lambda);
    % update Q
    temp = P+Y/mu;
    temp = (temp+temp')/2;
    Q = project_fantope(temp,k);
    
    dY = P-Q;
    chgP = max(max(abs(Pk-P)));
    chgQ = max(max(abs(Qk-Q)));
    chg = max([chgP chgQ max(abs(dY(:)))]);
    if DEBUG        
        if iter == 1 || mod(iter, 10) == 0
            obj = trace(P'*L)+lambda*norm(Q(:),1);
            err = norm(dY,'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
    DQ=D.*Q;
    obj(iter) = trace(P'*L)+lambda*norm(DQ(:),1);
    disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
end
% obj = trace(P'*L)+lambda*norm(Q(:),1);
err = norm(dY,'fro');