% fig 4b
% initialising parameters

A=100; M=10; T=0.1; C=0.7; beta1=4; beta2=0.01; lambda21=0.0005/6; lambda12=13/444000; mu=0.0005; delta=0.005; NEVOL=3500; f0=0.002; m0=2; alpha0=0.6; 
switching_environments=1; plasticity=0; return_traits=0; number_of_realisations=25; alphamax=1000;

cd ..
cd ..
cd ..

addpath('Simulation_Functions')

cd Data_generation_scripts/Fig4/panel_b


[~,~,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );


cd ..

cd ..

cd ..

save('Data_files\Fig4\panel_b\m_FRTE.mat','m');
save('Data_files\Fig4\panel_b\alpha_FRTE.mat','alpha');
