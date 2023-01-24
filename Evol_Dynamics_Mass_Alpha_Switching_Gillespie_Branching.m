function [StrainData_m,StrainData_alpha,m,alpha]=Evol_Dynamics_Mass_Alpha_Switching_Gillespie_Branching(m0,alpha0,A,M,T,C,betaB,betaG,lambdaGB,lambdaBG,mu,NEVOL,f0,delta)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over NEVOL mutations for a system undergoing bet-hedging in a
% stochastically switching environment. Here, the outputs StrainData_m and StrainData_alpha gives the set of all strains in one simulation of the coevoultionary dynamics. This enables us to
% visualise evolutionary branching in both mass and \alpha.

% As is the case for the code "Evol_Dynamics_Mass_Alpha_EvolBranching", branching is visualised by plotting up the set of all strains using the
% functions "Evol_Branching_plots" or "Evol_Branching_plots_Coevolution".


% Input Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             betaB - resistance to survival in bad environment
%             betaG - resistance to survival in good environment
%             lambdaGB - rate of switching from good to bad environment
%             lambdaBG - rate of switching from bad to good environment
% Output Parameters:
%             StrainData_m and StrainData_alpha are cell arrays storing information on the mass and fusion rate of each strain at a given time during the simulation. They both store the frequency of each strain at a given time in the simulation.
%             m and alpha are vectors showing the mean mass and fusion rate of the population respectively following each mutation event.

Good=1; Bad=0;
StrainData_m=cell(NEVOL,2);
StrainData_alpha=cell(NEVOL,2);

% This block of code is for initialisation.
%-------------------------------------------------------------
Nmut=1;
Nev=1;
if randn<0
mstrains=[0,m0,m0+delta*sign(randn)];
alphastrains=[0,alpha0,alpha0];
else
mstrains=[0,m0,m0];
alphastrains=[0,alpha0,alpha0+delta*sign(randn)];
end

f=[1-f0,f0];
S=2; %initialising S
mcur=0;
alphacur=0;


randomG=geornd(   (mu+lambdaGB)*ones(1,2500000)      )+ones(1,2500000) ;
randomB=geornd(   (mu+lambdaBG)*ones(1,2500000)      )+ones(1,2500000) ;
%-------------------------------------------------------------

while Nmut<NEVOL

 if Good==1
   tnext=randomG(Nev);


   if rand<mu/(mu+lambdaGB)
       f=f(:); 
       [mstrains,alphastrains,f]=Evol_Dynamics_Mass_Alpha_forBetHedging(mstrains,alphastrains,f,f0,A,M,T,C,betaG,tnext,delta,1);
       S=length(mstrains)-1;

       % This block of code calculates the mean mass and alpha of the
       % population before the next mutation occurs.
       %-------------------------------------------------------------
       for i=1:length(f)        
       mcur=mcur+(f(i)/sum(f))*mstrains(i+1);
       alphacur=alphacur+(f(i)/sum(f))*alphastrains(i+1);
       end       
       m(Nmut)=mcur;
       alpha(Nmut)=alphacur;
       mcur=0;
       alphacur=0;
       %-------------------------------------------------------------
       
       StrainData_m{Nmut,1}=f;
       StrainData_m{Nmut,2}=mstrains(2:end);

       StrainData_alpha{Nmut,1}=f;
       StrainData_alpha{Nmut,2}=alphastrains(2:end);
       
       
       Nmut=Nmut+1;       
 
       fprintf('Processing %d...',Nmut);
   else
       f=Invasion_dynamics(mstrains,alphastrains,S,betaG,C,tnext,T,A,M,f);
       
       Good=0; Bad=1;
       
   end
  
Nev=Nev+1;
 else
   tnext=randomB(Nev);

   if rand<mu/(mu+lambdaBG)
       f=f(:); 
       [mstrains,alphastrains,f]=Evol_Dynamics_Mass_Alpha_forBetHedging(mstrains,alphastrains,f,f0,A,M,T,C,betaB,tnext,delta,1);
       S=length(mstrains)-1;

       % This block of code calculates the mean mass and alpha of the
       % population before the next mutation occurs.
       %-------------------------------------------------------------
       for i=1:length(f)        
       mcur=mcur+(f(i)/sum(f))*mstrains(i+1);
       alphacur=alphacur+(f(i)/sum(f))*alphastrains(i+1);
       end       
       m(Nmut)=mcur;
       alpha(Nmut)=alphacur;
       mcur=0;
       alphacur=0;
       %-------------------------------------------------------------

       StrainData_m{Nmut,1}=f;
       StrainData_m{Nmut,2}=mstrains(2:end);

       StrainData_alpha{Nmut,1}=f;
       StrainData_alpha{Nmut,2}=alphastrains(2:end);
       
       Nmut=Nmut+1;

       fprintf('Processing %d...',Nmut);
   else
       f=Invasion_dynamics(mstrains,alphastrains,S,betaB,C,tnext,T,A,M,f);

       Good=1; Bad=0;
   end
Nev=Nev+1;
 end


end
