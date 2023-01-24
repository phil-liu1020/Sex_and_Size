function [mres,alphares]=Evol_Dynamics_Mass_Alpha_Many_Realisations(Nrealz,m0,alpha0,f0,A,M,T,C,beta,mu,delta,NEVOL)

% Obtains multiple realisations of the simulation for the coevolutionary dynamics between mass and \alpha.
% Parameters: Nrealz - number of realisations

m=zeros(Nrealz,0);
alpha=zeros(Nrealz,0);


parfor i=1:Nrealz
    
[m(i,:),alpha(i,:)]=Evol_Dynamics_Mass_Alpha(m0,alpha0,f0,A,M,T,C,beta,mu,delta,NEVOL);


end

assignin('base','m',m)
assignin('base','alpha',alpha)

