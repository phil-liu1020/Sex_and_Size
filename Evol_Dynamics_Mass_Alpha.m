function [mres,alphares]=Evol_Dynamics_Mass_Alpha(m0,alpha0,f0,A,M,T,C,beta,mu,delta,NEVOL)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over NEVOL mutations.
% Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             NEVOL - the number of mutations to occur. NEVOL=(total number of generations)*mu. E.g if you want the system to run for 10^(7) generations and the mutation rate is mu=10^(-3), then NEVOL=10^(4).

% This block of code is for initialisation.
%-------------------------------------------------------------
Nev=1;
mstrains=[0,m0];
alphastrains=[0,alpha0];
f=1;
S=2; %initialising S
mcur=0;
alphacur=0;
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
               else
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
     
             f=Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f);
     
     
     
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

             
             
             % This block of code calculates the mean mass and alpha of the
             % population before the next mutation occurs.
             %-------------------------------------------------------------
             for i=1:length(f)        
             mcur=mcur+(f(i)/sum(f))*mstrains(i+1);
             alphacur=alphacur+(f(i)/sum(f))*alphastrains(i+1);
             end       
             mres(Nev)=mcur;
             alphares(Nev)=alphacur;
             mcur=0;
             alphacur=0;
             %-------------------------------------------------------------
             
Nev=Nev+1;

end
