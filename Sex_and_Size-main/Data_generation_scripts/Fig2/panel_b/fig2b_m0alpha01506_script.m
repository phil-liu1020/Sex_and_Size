 %fig 2b (m(0),alpha(0))=(1.5,0.6)
% initialising parameters

A=100; M=1; T=1; C=0.6; beta1=1; beta2=1; lambda21=0; lambda12=0; mu=0.0005; delta=0.005; NEVOL=5500; f0=0.002; m0=1.5; alpha0=0.6; 
switching_environments=0; plasticity=0; return_traits=0; number_of_realisations=25; alphamax=1000;                   


addpath(genpath('Sex_and_Size-main'))

[~,~,m,alpha,~,~]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,C,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, plasticity, return_traits );


cd ..

cd ..

cd ..

save('Data_files\Fig2\panel_b\m_1506.mat','m');
save('Data_files\Fig2\panel_b\alpha_1506.mat','alpha');
