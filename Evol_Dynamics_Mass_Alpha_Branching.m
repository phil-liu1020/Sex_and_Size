function [StrainData_m,StrainData_alpha,m,alpha,mstrains,alphastrains,f]=Evol_Dynamics_Mass_Alpha_Branching(m0,alpha0,f0,A,M,T,C,beta,mu,delta,NEVOL,alphamax)

% This code simulates the coevolutionary dynamics between mass m and fusion rate \alpha over NEVOL mutations. 
% Here, the outputs StrainData_m and StrainData_alpha gives the set of all strains in one simulation of the coevoultionary dynamics. This enable us to visualise evolutionary branching in both mass and \alpha.
%
% Branching is visualised by plotting up the set of all strains using the
% functions "Evol_Branching_plots" or "Evol_Branching_plots_Coevolution".
% To visualise branching in m or alpha over time, run
% "Evol_Branching_plots(StrainData_m)" or
% "Evol_Branching_plots(StrainData_alpha)".
% To visualise branching in m and alpha simultaneously, run
% "Evol_Branching_plots_Coevolution(StrainData_m,StrainData_alpha)".

% Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             alphamax - maximum possible fusion rate imposed by energetic/environmental constraints. If we wish not to consider alphamax, set alphamax to some large value e.g. 1000. 

% This block of code is for initialisation.
%-------------------------------------------------------------
Nev=1;
mstrains=[0,m0];
alphastrains=[0,alpha0];
f=1;
S=2; %initialising S
mcur=0;
alphacur=0;
StrainData_m=cell(NEVOL,2);
StrainData_alpha=cell(NEVOL,2);
%-------------------------------------------------------------



while Nev<=NEVOL
    
       keepstrains=0;    
       r=rand;
       prob = f;
       q= sum(r >= cumsum([0;prob/sum(f)]));                                % q is the index of 
       
       if randn<0

       mstrains(length(mstrains)+1)=mstrains(q+1)+delta*sign(randn);
               if mstrains(end)<delta
               mstrains(end)=delta;    
               else
               end       
       alphastrains(length(alphastrains)+1)=alphastrains(q+1);              % IMPORTANT NOTE! the alpha value of the mutant should be the alpha value of the strain that underwent mutation in mass! Hence the "alphastrains(q+1)"
       alphamutated=0; mmutated=1;
       S=length(mstrains)-1;
       
       else
       alphastrains(length(alphastrains)+1)=alphastrains(q+1)+delta*sign(randn);
               if alphastrains(end)<0
               alphastrains(end)=0;    
               end 
               
               if alphastrains(end)>=alphamax
               alphastrains(end)=alphamax;    
               end
       mstrains(length(mstrains)+1)=mstrains(q+1);                         % IMPORTANT NOTE! the m value of the mutant should be the m value of the strain that underwent mutation in alpha! Hence the "mstrains(q+1)"
       S=length(alphastrains)-1;   
       alphamutated=1; mmutated=0;
       
       end
 
             % This block of code checks if there's any duplicate strains
             %-------------------------------------------------------------
             aa=[mstrains' alphastrains'];   
             [u,I,J] = unique(aa, 'rows', 'first');
             hasDuplicates = size(u,1) < size(aa,1);
             ixDupRows = setdiff(1:size(aa,1), I);
             dupRowValues = aa(ixDupRows,:);
             duprow=0;

             for i=1:length(aa(:,1))
             duprow(i)=isequal(aa(i,:),dupRowValues);
             end

             if hasDuplicates>0
             f=(1-f0)*f;
             f(find(duprow,1, 'first')-1)=f(find(duprow,1, 'first')-1) + f0;
             mstrains(max(ixDupRows))=[];
             alphastrains(max(ixDupRows))=[];
             S=length(mstrains)-1;
             else
             f=(1-f0)*f;
             f(end+1)=f0;
             end
             %-------------------------------------------------------------
             
             G=geornd(mu)+1;
     
                  fout=Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f);
             f=fout;
     
     
             % Block of code that Discards strains with frequency below 10^(-3)
             %-------------------------------------------------------------
             for i=1:length(f)
                if f(i)<10^(-3)
                keepstrains(i,:)=0;
                else
                keepstrains(i,:)=1;    
                end
             end

             f(keepstrains==0)=[];
             f=f/sum(f);
             mstrains([1; keepstrains]==0)=[];
             alphastrains([1; keepstrains]==0)=[];
             S=length(mstrains)-1;
             %-------------------------------------------------------------

             StrainData_m{Nev,1}=f;
             StrainData_m{Nev,2}=mstrains(2:end);

             StrainData_alpha{Nev,1}=f;
             StrainData_alpha{Nev,2}=alphastrains(2:end);
             

          
             % This block of code calculates the mean mass and alpha of the
             % population before the next mutation occurs.
             %-------------------------------------------------------------
             for i=1:length(f)        
             mcur=mcur+(f(i)/sum(f))*mstrains(i+1);
             alphacur=alphacur+(f(i)/sum(f))*alphastrains(i+1);
             end       
             m(Nev)=mcur;
             alpha(Nev)=alphacur;
             mcur=0;
             alphacur=0;
             %-------------------------------------------------------------
           
Nev=Nev+1;

end
