% Instruction for FigS15 panel a)

cd ..
cd ..
addpath(genpath('Sex_and_Size-main'))
cd Sex_and_Size-main/Figure_generation_scripts

A=100; M=1; T=1; C=0.5;  beta2=2; lambda12=0.05; lambda21=3/140; mu=0.002; NEVOL=500; f0=0.002; delta=0.02; alphamax=1000; switching_environments=1;return_genotypes=0; plasticity=1; 
alpha0=0; 

figS15output= Plotting_FigS15_a(alpha0,A,M,T,C,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_genotypes,plasticity);

plot(figS15output,'o')
