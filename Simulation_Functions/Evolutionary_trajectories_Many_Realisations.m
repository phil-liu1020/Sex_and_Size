function [m,alpha,m2,alpha2]=Evolutionary_trajectories_Many_Realisations(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits )
                                                         
% Obtains multiple realisations of the simulation for the coevolutionary dynamics between mass and \alpha in a system undergoing bet-hedging in a stochastically switching environment.
% Parameters: Nrealz - number of realisations
%             f0 - initial frequency of rare mutant

m=zeros(number_of_realisations,0);
alpha=zeros(number_of_realisations,0);
m2=zeros(number_of_realisations,0);
alpha2=zeros(number_of_realisations,0);

if plasticity==1

    parfor i=1:number_of_realisations
    
    [~,~,m(i,:),alpha(i,:),m2(i,:),alpha2(i,:)]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity);
    end

assignin('base','m',m)
assignin('base','alpha',alpha)
assignin('base','m2',m2)
assignin('base','alpha2',alpha2)

else

    parfor i=1:number_of_realisations
    
    [~,~,m(i,:),alpha(i,:),~,~]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity)
  
    end

assignin('base','m',m)
assignin('base','alpha',alpha)


end




