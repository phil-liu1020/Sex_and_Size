
% fig 4b
% initialising parameters

A=100; M=10; T=0.1; C=0.7; betaB=4; betaG=0.01; lambdaGB=1/6; lambdaBG=13/222; mu=0.000035; delta=0.007; NEVOL=3500; f0=0.002; m0=2; alpha0=0.6;  switching_environments=1; plasticity=0; return_traits=0; Nrealz=25; alphamax=1000;


[~,~,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );

