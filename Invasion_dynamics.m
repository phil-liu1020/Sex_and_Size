function gend = Invasion_dynamics(mstrains,alphastrains,S,beta,C,G,T,A,M,f0)

% This function simulates then invasion dynamics over G generations.
% Parameters: G - number of generations
%             f0 - initial frequency of each strain.

g=zeros(S,1);

f=f0;

for i=1:G

x = Fertilisation_Kinetics(mstrains,alphastrains,f,S,A,M,T);

f=Single_Generation(mstrains,alphastrains,x,S,beta,C);

g(:,i)=f;

end

gend=g(:,end);