% fig S13
% initialising parameters

A=100; M=1; T=1; C=0.99; beta1=1; beta2=1; lambda21=0; lambda12=0; mu=0.001; delta=0.01; NEVOL=2500; f0=0.05; m0=0.25; alpha0=1.5; 
switching_environments=0; plasticity=0; return_traits=1; number_of_realisations=1; alphamax=1000;

cd ..
cd ..
addpath(genpath('Simulation_Functions'))
cd Data_generation_scripts/FigS13

[StrainsData_m,StrainsData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );



cd ..
cd ..


save('Data_files\FigS13\m.mat','m');
save('Data_files\FigS13\alpha.mat','alpha');
save('Data_files\FigS13\StrainsData_m.mat','StrainsData_m');
save('Data_files\FigS13\StrainsData_alpha.mat','StrainsData_alpha');
