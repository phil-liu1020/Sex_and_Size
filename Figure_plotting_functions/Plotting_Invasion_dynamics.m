function [g,gend] = Plotting_Invasion_dynamics(mgenotypes,alphagenotypes,S,beta,C,G,T,A,M,f)

% This function simulates then invasion dynamics over G generations.
% Parameters: G - number of generations
%             f - initial frequency of each genotype.
%             mgenotypes and alphagenotypes are the mass and fusion rates of each genotype.
%             S - number of genotypes.
%             A, M T, C and \beta are just the parameters given in the analytical part of the text.

%             The difference between this function and "Invasion_dynamics" is that this function outputs the entire invasion trajectory 
%             (i.e. the frequency at each generation from 0 to G) whereas "Invasion_dynamics" only outputs the frequency at the G-th generation.

g=zeros(S,1);

% The for loop iterates the function "Single_Generation" over G generations. This function outputs the frequency of each genotype at the G-th generation.

for i=1:G

x = Fertilisation_Kinetics(mgenotypes,alphagenotypes,f,S,A,M,T);

f=Survival_function(mgenotypes,alphagenotypes,x,S,beta,C);

g(:,i)=f;

end

gend=g(:,end);
