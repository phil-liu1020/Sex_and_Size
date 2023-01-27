function gend = Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f)

% This function simulates then invasion dynamics over G generations.
% Parameters: G - number of generations
%             f - initial frequency of each strain.
%             mstrains and alphastrains are the mass and fusion rates of each strain.
%             S - number of strains.
%             A, M T, C and \beta are just the parameters given in the analytical part of the text.

g=zeros(S,1);

% The for loop iterates the function "Single_Generation" over G generations. This function outputs the frequency of each strain at the G-th generation.

for i=1:G

x = Fertilisation_Kinetics(mstrains,alphastrains,f,S,A,M,T);

f=Survival_function(mstrains,alphastrains,x,S,beta,C);

g(:,i)=f;

end

gend=g(:,end);