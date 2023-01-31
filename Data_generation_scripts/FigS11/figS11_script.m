% fig S11
% initialising parameters

A=100; M=10; T=0.1; C=0.7; beta1=4; beta2=0.01; lambda21=1/6; lambda12=13/222; mu=0.00035; delta=0.007; NEVOL=4000; f0=0.002; m0=2; alpha0=0.6; 
switching_environments=1; plasticity=0; return_traits=1; number_of_realisations=1; alphamax=1000;

addpath(genpath('Sex_and_Size-main'))

[StrainsData_m,StrainsData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );



cd ..
cd ..


save('Data_files\FigS11\m.mat','m');
save('Data_files\FigS11\alpha.mat','alpha');
save('Data_files\FigS11\StrainsData_m.mat','StrainsData_m');
save('Data_files\FigS11\StrainsData_alpha.mat','StrainsData_alpha');
