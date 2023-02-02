% fig S9
% initialising parameters

A=100; M=10; T=0.1; C=0.35; beta1=3; beta2=0.01; lambda12=1/4; lambda21=67/532; mu=0.00035; delta=0.007; NEVOL=4000; f0=0.002; m0=2; alpha0=0.4; 
switching_environments=1; plasticity=0; return_genotypes=1; number_of_realisations=1; alphamax=1000;

cd ..
cd ..
addpath(genpath('Simulation_Functions'))
cd Data_generation_scripts/FigS9

[genotypesData_m,genotypesData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_genotypes );



cd ..
cd ..


save('Data_files\FigS9\m.mat','m');
save('Data_files\FigS9\alpha.mat','alpha');
save('Data_files\FigS9\genotypesData_m.mat','genotypesData_m');
save('Data_files\FigS9\genotypesData_alpha.mat','genotypesData_alpha');
