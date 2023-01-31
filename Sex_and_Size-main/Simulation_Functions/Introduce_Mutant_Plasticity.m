function [mstrains,alphastrains,mstrains2,alphastrains2,f,f2]=Introduce_Mutant_Plasticity(mstrains,mstrains2,alphastrains,alphastrains2,f,f2,f0,A,M,T,C,beta,beta2,G,delta,plasticity,alphamax)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over 1 mutations. The difference between this and "Evol_Dynamics_Mass_Alpha" is that the evolutionary dynamics can be started from any number of strains. 
% This code is used for simulating
% "Evol_Dynamics_Mass_Alpha_Switching_Gillespie".

% Parameters: mstrains and alphastrains - the mass and fusion rate of all the strains.
%             f - initial frequency of all the strains.
%             f0 - initial frequency of rare mutant.

% This block of code is for initialisation.
%-------------------------------------------------------------

if length(mstrains)~=length(alphastrains) || length(mstrains2)~=length(alphastrains2)
error('length of mstrains not equal to length of alphastrains')
end 

if length(f)~=length(alphastrains)-1
error('length of f not equal to length of alphastrains-1')
end 

if length(f2)~=length(alphastrains2)-1
error('length of f2 not equal to length of alphastrains-1')
end 

Nev=1;
S=length(mstrains)-1; %initialising S
mcur=0;
alphacur=0;
S2=length(mstrains2)-1; %initialising S
mcur2=0;
alphacur2=0;

%-------------------------------------------------------------



while Nev<=1
    
       keepstrains=0;    
       r=rand;
       prob = f;
       q= sum(r >= cumsum([0;prob/sum(f)]));                                % q is the index of 

       keepstrains2=0;    
       r2=rand;
       prob2 = f2;
       q2= sum(r2 >= cumsum([0;prob2/sum(f2)]));
       
       if randn<0                                                           % if randn<0 mutate in mass, else mutate in \alpha.

           if plasticity==1
              if rand<0.5                                                   % if rand<0.5, mutate in environment 2, else mutate in environment 1.
              mstrains2(length(mstrains2)+1)=mstrains2(q2+1)+delta*sign(randn);
                  if mstrains2(end)<delta
                  mstrains2(end)=delta;    
                  end
              alphastrains2(length(alphastrains2)+1)=alphastrains2(q2+1);
              S2=length(alphastrains2)-1;
              m1=0;                                                         % m1, alpha1, m2 and alpha2 are set to 1 if a mutation in m or alpha has occured in environment 1 or 2.
              alpha1=0;
              m2=1;
              alpha2=0;
              else
              mstrains(length(mstrains)+1)=mstrains(q+1)+delta*sign(randn);
                  if mstrains(end)<delta
                  mstrains(end)=delta;    
                  end
              alphastrains(length(alphastrains)+1)=alphastrains(q+1);
              S=length(alphastrains)-1;
              m1=1;
              alpha1=0;
              m2=0;
              alpha2=0;
              end

           else

                mstrains(length(mstrains)+1)=mstrains(q+1)+delta*sign(randn);
         
                if mstrains(end)<delta
                mstrains(end)=delta;    
                end 
                alphastrains(length(alphastrains)+1)=alphastrains(q+1);              % IMPORTANT NOTE! the alpha value of the mutant should be the alpha value of the strain that underwent mutation in mass! Hence the "alphastrains(q+1)"
                S=length(mstrains)-1;
                m1=1;
                alpha1=0;
                m2=0;
                alpha2=0;

            end
       
       else
                
           
            if plasticity==1
               if rand<0.5
               alphastrains2(length(alphastrains2)+1)=alphastrains2(q2+1)+delta*sign(randn);
                   if alphastrains2(end)<0
                   alphastrains2(end)=0;    
                   end
               mstrains2(length(mstrains2)+1)=mstrains2(q2+1);
               S2=length(mstrains2)-1;
               m1=0;
               alpha1=0;
               m2=0;
               alpha2=1;
               else
               alphastrains(length(alphastrains)+1)=alphastrains(q+1)+delta*sign(randn);
                   if alphastrains(end)<0
                   alphastrains(end)=0;    
                   end
               mstrains(length(mstrains)+1)=mstrains(q+1); 
               S=length(mstrains)-1;
               m1=0;
               alpha1=1;
               m2=0;
               alpha2=0;
               end

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
                m1=0;
                alpha1=1;
                m2=0;
                alpha2=0;
            end
       
       end
 
             if m1==1 || alpha1==1
             % This block of code checks if there's any duplicate strains
             % in environment 1. Only to be done if a mutation has occured in environment 1.
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
             end



             if m2==1 || alpha2==1
             % This block of code checks if there's any duplicate strains
             % in environment 2 if plasiticty==1. Only to be done if a mutation has occured in environment 2.
             %-------------------------------------------------------------
             aa2=[mstrains2' alphastrains2'];   
             [u2,I2,J2] = unique(aa2, 'rows', 'first');
             hasDuplicates2 = size(u2,1) < size(aa2,1);
             ixDupRows2 = setdiff(1:size(aa2,1), I2);
             dupRowValues2 = aa2(ixDupRows2,:);
             duprow2=0;

             for i=1:length(aa2(:,1))
             duprow2(i)=isequal(aa2(i,:),dupRowValues2);
             end

             if hasDuplicates2>0
             f2=(1-f0)*f2;
             f2(find(duprow2,1, 'first')-1)=f2(find(duprow2,1, 'first')-1) + f0;
             mstrains2(max(ixDupRows2))=[];
             alphastrains2(max(ixDupRows2))=[];
             S2=length(mstrains2)-1;
             else
             f2=(1-f0)*f2;
             f2(end+1)=f0;
             end
             %-------------------------------------------------------------    
             end
             
             
     f=Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f);
     if plasticity==1
     f2=Invasion_dynamics(mstrains2,alphastrains2,S2,beta2,C,G,T,A,M,f2);
     end
     
             % Block of code that Discards strains with frequency below
             % 10^(-3) in environment 1
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

             if plasticity==1
             % Block of code that Discards strains with frequency below
             % 10^(-3) in environment 2
             %-------------------------------------------------------------
             for i=1:length(f2)
                if f2(i)<10^(-3)
                keepstrains2(i,:)=0;
                else
                keepstrains2(i,:)=1;    
                end
             end

             f2(keepstrains2==0)=[];
             f2=f2/sum(f2);
             mstrains2([1; keepstrains2]==0)=[];
             alphastrains2([1; keepstrains2]==0)=[];
             S2=length(mstrains2)-1;
             %-------------------------------------------------------------
             end

             
             
             % This block of code calculates the mean mass and alpha of the
             % population before the next mutation occurs in environment 1.
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

             if plasticity==1
             % This block of code calculates the mean mass and alpha of the
             % population before the next mutation occurs in environment 2.
             %-------------------------------------------------------------
             for i=1:length(f2)        
             mcur2=mcur2+(f2(i)/sum(f2))*mstrains2(i+1);
             alphacur2=alphacur2+(f2(i)/sum(f2))*alphastrains2(i+1);
             end       
             mres2(Nev)=mcur2;
             alphares2(Nev)=alphacur2;
             mcur2=0;
             alphacur2=0;
             %-------------------------------------------------------------
             end
             
Nev=Nev+1;

end