function [m,alpha,muttime,meanswBG,meanswGB,randomG]=Evol_Dynamics_Mass_Alpha_Switching_Gillespie(m0,alpha0,A,M,T,C,betaB,betaG,lambdaGB,lambdaBG,mu,NEVOL,f0,delta)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over NEVOL mutations for a system undergoing bet-hedging in a stochastically switching environment.
% Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             betaB - resistance to survival in bad environment
%             betaG - resistance to survival in good environment
%             lambdaGB - rate of switching from good to bad environment
%             lambdaBG - rate of switching from bad to good environment
%             NEVOL is the number of mutations to occur. NEVOL=(total number of generations)*mu. E.g if you want the system to run for 1.1x10^(7) generations and the mutation rate is mu=5x10^(-4), then NEVOL=5.5x10^(3).

% Please note that this code can also be used to simulate the coevolutionary dynamics between m and \alpha in a fixed
% environment by either setting betaG=betaB or setting lambdaGB and lambdaBG to 0.


Good=1; Bad=0;
tnextcum=0;
tswG=0;
tswB=0;

% This block of code is for initialisation.
%-------------------------------------------------------------
NswGB=1;
NswBG=1;
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
muttime=0;
switchtimeGB=0;
switchtimeBG=0;
meanswGB=0;
meanswBG=0;


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
      
       tnextcum=tnextcum+tnext;
       
       tswG=tnext;

       muttime(Nmut)=tnextcum;
       
       meanswGB(Nmut)=mean(switchtimeGB);
       meanswBG(Nmut)=mean(switchtimeBG);
       switchtimeGB=0;
       switchtimeBG=0;

       tnextcum=0;
       
       Nmut=Nmut+1;
       
       NswGB=1;
       NswBG=1;
       
       fprintf('Processing %d...',Nmut);
   else
       f=Invasion_dynamics(mstrains,alphastrains,S,betaG,C,tnext,T,A,M,f);
       
       if tswG~=0
       tswG=tswG+tnext;
       switchtimeGB(NswGB)=tswG;      
       else
       switchtimeGB(NswGB)=tnext;   
       end
      
 
       tnextcum=tnextcum+tnext;
       tswG=0;
       NswGB=NswGB+1;
       Good=0; Bad=1;
       
   end
  
Nev=Nev+1;
 else
   tnext=randomB(Nev);

   if rand<mu/(mu+lambdaBG)
       f=f(:); % forces f to be a column vector
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
       tnextcum=tnextcum+tnext;
       
       tswB=tnext;
       
       
       muttime(Nmut)=tnextcum;
       
       meanswGB(Nmut)=mean(switchtimeGB);
       meanswBG(Nmut)=mean(switchtimeBG);
       switchtimeGB=0;
       switchtimeBG=0;
       
       tnextcum=0;
       
       Nmut=Nmut+1;
       
       NswBG=1;
       NswGB=1;

       fprintf('Processing %d...',Nmut);
   else
       f=Invasion_dynamics(mstrains,alphastrains,S,betaB,C,tnext,T,A,M,f);

       
       if tswB~=0
       tswB=tswB+tnext;
       switchtimeBG(NswBG)=tswB;      
       else
       switchtimeBG(NswBG)=tnext;    
       end
       
       tnextcum=tnextcum+tnext;
       tswB=0;       
       NswBG=NswBG+1;
       Good=1; Bad=0;
   end
Nev=Nev+1;
 end


end
