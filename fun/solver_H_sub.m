function [P]=solver_H_sub(D,Q,L,Y,mu,lambda)
P=zeros(size(D));

temp=Q-(Y+L)/mu;
P(D==0)=temp(D==0);
temp=(mu./(2*lambda*D+mu)).*Q-(Y+L)/mu;
P(D~=0)=temp(D~=0);