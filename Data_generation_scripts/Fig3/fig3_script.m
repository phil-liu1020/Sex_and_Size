% fig 3
% initialising parameters

A=100; M=1; T=1; C=0.6; beta1=1; beta2=1; lambda12=0; lambda21=0; mu=0.001; delta=0.01; NEVOL=3500; f0=0.002; m0=1.5; alpha0=0.6;  switching_environments=0; plasticity=0; return_genotypes=1; number_of_realisations=1; alphamax=1000;

cd ..
cd ..

addpath(genpath('Simulation_Functions'))

cd Data_generation_scripts/Fig3

[genotypeData_m,genotypeData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_genotypes );


cd ..
cd ..


save('Data_files\Fig3\m.mat','m');
save('Data_files\Fig3\alpha.mat','alpha');
save('Data_files\Fig3\genotypeData_m.mat','genotypeData_m');
save('Data_files\Fig3\genotypeData_alpha.mat','genotypeData_alpha');
