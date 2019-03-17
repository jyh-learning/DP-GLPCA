function [Dold,obj]=DP(W,S,P,T,DPpara)
%%%%%%% function readme %%%%%%%%
% the DP function solves the pairwise constraints propagation probelm,
% i.e, min 0.5**trace(D'*L*D)+alpha*trace(S'*D)+0.5*beta*norm(P.*(D-T))
% the first term is the smoothness tern to guarantee the propagation
% the second term is coding term to guarantee the correctness
% and the last term guarantees the consistency  of the Pairwise information

% W the local similarity matrix for 
% S the dense similarity matrix for codeing
% A the diagonal matrix for W
% L the Laplacian matrix for W 
% P the supervisory information position matric
% T the pairwise supervisory information matrix
%%%%%%% parameter setting
A=diag(sum(W,2));
L=A-W;
alpha=DPpara.alpha;
beta=DPpara.beta;
maxiter=DPpara.maxiter;

Dold=rand(size(A));
Dold=(Dold+Dold')/2;
% Dnew=zeros(size(Dold));
obj(1)=0.5*trace(Dold'*L*Dold)+alpha*trace(S'*Dold)+0.5*beta*norm(P.*(Dold-T));
betaPT=beta*(P.*T);
alphaS=alpha*S;
for iter=1:maxiter
    
%     Dnew=Dold.*((W*Dold+beta*(P.*T))./(A*Dold+alpha*S+beta*(P.*Dold)));
    Dnew=Dold.*((W*Dold+betaPT)./(A*Dold+alphaS+beta*(P.*Dold)+eps));
    obj(iter+1)=0.5*trace(Dnew'*L*Dold)+alpha*trace(S'*Dnew)+0.5*beta*norm(P.*(Dnew-T));
    disp(['the ',num2str(iter),' iteration. obj value: ',num2str(obj(iter+1))]);
%     if max(max(abs(Dnew-Dold)))<10^-3 || abs(obj(iter)-obj(iter+1))<10^-3
%         break;
%     end
    if max(max(abs(Dnew-Dold)))<10^-3 
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dold=(Dnew+Dnew')/2;
    % or
%     Dold=Dnew;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



