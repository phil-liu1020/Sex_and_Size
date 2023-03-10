% fig S9
% initialising parameters

A=100; M=1; T=1; C=0.99; beta1=1; beta2=1; lambda21=0; lambda12=0; mu=0.001; delta=0.01; NEVOL=2000; f0=0.05; m0=0.25; alpha0=1.5; 
switching_environments=0; plasticity=0; return_traits=1; number_of_realisations=1; alphamax=1000;


[StrainsData_m,StrainsData_alpha,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );


addpath(genpath('Sex_and_Size-main'))

cd ..
cd ..


save('Data_files\FigS9\m.mat','m');
save('Data_files\FigS9\alpha.mat','alpha');
save('Data_files\FigS9\StrainsData_m.mat','StrainsData_m');
save('Data_files\FigS9\StrainsData_alpha.mat','StrainsData_alpha');