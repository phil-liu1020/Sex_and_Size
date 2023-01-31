function [StrainData_m,StrainData_alpha,m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over NEVOL mutations for a system undergoing bet-hedging in a
% stochastically switching environment. Here, the outputs StrainData_m and StrainData_alpha gives the set of all strains in one simulation of the coevoultionary dynamics. This enables us to
% visualise evolutionary branching in both mass and \alpha.

% As is the case for the code "Evol_Dynamics_Mass_Alpha_EvolBranching", branching is visualised by plotting up the set of all strains using the
% functions "Evol_Branching_plots" or "Evol_Branching_plots_Coevolution".


% Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             beta1 - resistance to survival in Env1 environment
%             beta2 - resistance to survival in Env2 environment
%             lambda21 - rate of switching from Env2 to Env1 environment
%             lambda12 - rate of switching from Env1 to Env2 environment

if beta2~=beta1 && switching_environments==0
    error('No switching environments but environmental switching specified')
end

if beta2==beta1 && switching_environments==1
    error('switching environments but environmental switching not specified')
end

if lambda21==0 && lambda12==0 && switching_environments==1
    error('switching environments but environmental switching not specified')
end

if switching_environments==0 && plasticity==1
   error('No switching environments but plasticit specified')
end


% This block of code is for initialisation.
%-------------------------------------------------------------
Nmut=1;
Nev=1;
if randn<0
mstrains=[0,m0,m0+delta*sign(randn)];
alphastrains=[0,alpha0,alpha0];
else
mstrains=[0,m0,m0];
if alpha0<delta
alphastrains=[0,alpha0,alpha0+delta];    
else
alphastrains=[0,alpha0,alpha0+delta*sign(randn)];
end
end

if plasticity==1
   if randn<0
   mstrains2=[0,m0,m0+delta*sign(randn)];
   alphastrains2=[0,alpha0,alpha0];
   else
   mstrains2=[0,m0,m0];
   if alpha0<delta
   alphastrains2=[0,alpha0,alpha0+delta];
   else
   alphastrains2=[0,alpha0,alpha0+delta*sign(randn)];
   end
   end 
else
mstrains2=0;
alphastrains2=0;
end

f=[1-f0,f0];
S=2; %initialising S
mcur=0;
alphacur=0;
if plasticity==1
f2=[1-f0,f0];
S2=2; %initialising S
mcur2=0;
alphacur2=0;    
else
f2=[];   % emptty set!!! to ensure length of 0
end


random2=geornd(   (mu+lambda21)*ones(1,2500000)      )+ones(1,2500000) ;
random1=geornd(   (mu+lambda12)*ones(1,2500000)      )+ones(1,2500000) ;

StrainData_m=cell(NEVOL,2);
StrainData_alpha=cell(NEVOL,2);

Env2=1; Env1=0;

    if plasticity==1
    mstrains2=mstrains;
    alphastrains2=alphastrains;
    else
    mstrains2=0;
    alphastrains2=0;
    end

%-------------------------------------------------------------

while Nmut<NEVOL

 if Env2==1
   tnext=random2(Nev);


   if rand<mu/(mu+lambda21)
       f=f(:); f2=f2(:);
       [mstrains,alphastrains,mstrains2,alphastrains2,f,f2]=Introduce_Mutant_Plasticity(mstrains,mstrains2,alphastrains,alphastrains2,f,f2,f0,A,M,T,C,beta1,beta2,tnext,delta,plasticity,alphamax);
       S=length(mstrains)-1;

       if plasticity==1
       S2=length(mstrains2)-1;
       end

       % This block of code calculates the mean mass and alpha of the
       % population in environment 1 before the next mutation occurs.
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

       if plasticity==1
       % This block of code calculates the mean mass and alpha of the
       % population in environment 2 before the next mutation occurs.
       %-------------------------------------------------------------
       for i=1:length(f2)        
       mcur2=mcur2+(f2(i)/sum(f2))*mstrains2(i+1);
       alphacur2=alphacur2+(f2(i)/sum(f2))*alphastrains2(i+1);
       end       
       m2(Nmut)=mcur2;
       alpha2(Nmut)=alphacur2;
       mcur2=0;
       alphacur2=0;
       %-------------------------------------------------------------
       else
       m2=0; alpha2=0;
       end

       if return_traits==1
       StrainData_m{Nmut,1}=f;
       StrainData_m{Nmut,2}=mstrains(2:end);

       StrainData_alpha{Nmut,1}=f;
       StrainData_alpha{Nmut,2}=alphastrains(2:end);
       end



       
       
       Nmut=Nmut+1;       
 
       fprintf('Processing %d...',Nmut);
   else
       
       if plasticity==1
       f2=Invasion_dynamics(mstrains2,alphastrains2,S2,beta2,C,tnext,T,A,M,f2);  
      % f=Invasion_dynamics(mstrains,alphastrains,S,beta2,C,tnext,T,A,M,f);
       else
       
            f=Invasion_dynamics(mstrains,alphastrains,S,beta2,C,tnext,T,A,M,f);
       end
       
       Env2=0; Env1=1;
       
   end
  
Nev=Nev+1;
 else
   tnext=random1(Nev);

   if rand<mu/(mu+lambda12)
       f=f(:); f2=f2(:);
       [mstrains,alphastrains,mstrains2,alphastrains2,f,f2]=Introduce_Mutant_Plasticity(mstrains,mstrains2,alphastrains,alphastrains2,f,f2,f0,A,M,T,C,beta1,beta2,tnext,delta,plasticity,alphamax);
       S=length(mstrains)-1; 

       if plasticity==1
       S2=length(mstrains2)-1;
       end

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


       if plasticity==1
       % This block of code calculates the mean mass and alpha of the
       % population in environment 2 before the next mutation occurs.
       %-------------------------------------------------------------
       for i=1:length(f2)        
       mcur2=mcur2+(f2(i)/sum(f2))*mstrains2(i+1);
       alphacur2=alphacur2+(f2(i)/sum(f2))*alphastrains2(i+1);
       end       
       m2(Nmut)=mcur2;
       alpha2(Nmut)=alphacur2;
       mcur2=0;
       alphacur2=0;
       %-------------------------------------------------------------
       else
       m2=0; alpha2=0;
       end



       if return_traits==1
       StrainData_m{Nmut,1}=f;
       StrainData_m{Nmut,2}=mstrains(2:end);

       StrainData_alpha{Nmut,1}=f;
       StrainData_alpha{Nmut,2}=alphastrains(2:end);
       end




       
       Nmut=Nmut+1;

       fprintf('Processing %d...',Nmut);
   else
       if plasticity==1
       %f2=Invasion_dynamics(mstrains2,alphastrains2,S2,beta1,C,tnext,T,A,M,f2);  
       f=Invasion_dynamics(mstrains,alphastrains,S,beta1,C,tnext,T,A,M,f);
       else
       f=Invasion_dynamics(mstrains,alphastrains,S,beta1,C,tnext,T,A,M,f);
       end

       Env2=1; Env1=0;
   end
Nev=Nev+1;
 end


end