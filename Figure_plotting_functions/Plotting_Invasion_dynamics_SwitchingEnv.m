function Invtraj=Plotting_Invasion_dynamics_SwitchingEnv(mgenotypes,alphagenotypes,S,betaB,betaG,lambdaGB,lambdaBG,C,G,T,A,M,f,stochastic)

% This function simulates the invasion dynamics of a system undergoing bet-hedging in an environment that switches FRTI. It does so over G generations.
% Parameters: G - number of generations
%             f - initial frequency of each genotype.
%             mgenotypes and alphagenotypes are the mass and fusion rates of each genotype.
%             S - number of genotypes.
%             A, M T, C and \beta are just the parameters given in the analytical part of the text.
%             betaB and betaG are the resistance to survival in the bad and good environments respectively.
%             lambdaGB and lambdaBG are the rates of switching from the good to bad and bad to good environments repsectively.
%             If you want to simulate stochastic switching, set
%             stochastic=1. If you want to simulate fixed switching, set stochastic=0.


tinc=0;                      % tinc is the time elapsed during a switching event. tinc is incremented right after a switching event.
Invtraj=zeros(1,G);          % Invtraj will represent the invasion trajectory of an invading mutant.
Good=1; Bad=0;               % The default environmental state is the good environmental state.

while tinc<=G

if Good==1

    if stochastic==1
    tnext=geornd(  lambdaGB )+1;        % tnext is the time until the next switching event. 
    else
    tnext=round(1/lambdaGB);    
    end

[g,f]=Plotting_Invasion_dynamics(mgenotypes,alphagenotypes,S,betaG,C,tnext,T,A,M,f);
Invtraj(tinc+1:tnext+tinc)=g(2,:);
tinc=tinc+tnext;   % Here is wheree tinc is incremented.
Good=0; Bad=1;
else

    if stochastic==1
    tnext=geornd(  lambdaBG )+1;
    else
    tnext=round(1/lambdaBG); 
    end


[g,f]=Plotting_Invasion_dynamics(mgenotypes,alphagenotypes,S,betaB,C,tnext,T,A,M,f);
Invtraj(tinc+1:tnext+tinc)=g(2,:);
tinc=tinc+tnext;
Good=1; Bad=0;
end

end
