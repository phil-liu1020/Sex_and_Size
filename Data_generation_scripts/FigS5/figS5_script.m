% fig S5
% initialising parameters

A=100; M=1; T=1; C=0.6; beta1=1; beta2=1; lambda21=0; lambda12=0; mu=0.0005; delta=0.005; NEVOL=2000; f0=0.002; m0=0.35; alpha0=0; 
switching_environments=0; plasticity=0; return_genotypes=1; number_of_realisations=1; alphamax=0.05;

cd ..
cd ..
addpath(genpath('Simulation_Functions'))
cd Data_generation_scripts/FigS5

[genotypesData_m,genotypesData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_genotypes );

% I need to define what return_genotypes means exactly!!!!!!



cd ..

cd ..

save('Data_files\FigS5\m.mat','m');
save('Data_files\FigS5\alpha.mat','alpha');
save('Data_files\FigS5\genotypesData_m.mat','genotypesData_m');
save('Data_files\FigS5\genotypesData_alpha.mat','genotypesData_alpha');
