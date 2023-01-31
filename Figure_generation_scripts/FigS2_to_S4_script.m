% Instructions for producing FigS2 to FigS4

% FigS2

addpath(genpath('Sex_and_Size-main'))

cd C:\Users\liuph\Downloads\Sex_and_Size-main\Simulation_Functions

mstrains=[0,0.4,0.395]; alphastrains=[0,0.1,0.1]; S=2; beta=1; C=0.6; G=4000; T=1; A=100; M=1; f=[0.998,0.002];

[g,~] = Plotting_Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f);

plot(g(2,:))

% FigS3

addpath(genpath('Sex_and_Size-main'))

cd C:\Users\liuph\Downloads\Sex_and_Size-main\Simulation_Functions

mstrains=[0,0.4,0.4]; alphastrains=[0,0.1,0.105]; S=2; beta=1; C=0.6; G=15000; T=1; A=100; M=1; f=[0.998,0.002];

[g,~] = Plotting_Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f);

plot(g(2,:))

% FigS4

addpath(genpath('Sex_and_Size-main'))

cd C:\Users\liuph\Downloads\Sex_and_Size-main\Simulation_Functions

mstrains=[0,0.3,0.305]; alphastrains=[0,0.1,0.1]; S=2; beta1=4; beta2=0.01; lambda12=13/222; lambda21=1/6; C=0.7; G=600; T=0.1; A=100; M=10; f=[0.998,0.002];

Invtraj=Plotting_Invasion_dynamics_SwitchingEnv(mstrains,alphastrains,S,beta1,beta2,lambda21,lambda12,C,G,T,A,M,f,stochastic);

plot(Invtraj)