function [StrainData_m,StrainData_alpha,m,alpha,m2,alpha2]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits )

if number_of_realisations>1

[m,alpha,m2,alpha2]=Evolutionary_trajectories_Many_Realisations(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );    

StrainData_m=0;
StrainData_alpha=0;

else

[StrainData_m,StrainData_alpha,m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity);

end

% This is the daddy function.