function x = Fertilisation_Kinetics(mstrains,alphastrains,f,S,A,M,T)

% Simulates the fertilisation kinetics for a fertilisation period of T.
% Input parameters: mstrains and alphastrains - the mass and fusion rate of all the strains.
%                   f - frequency of each strain.
%                   S - number of distinct strains.

if length(mstrains)~=S+1 || length(alphastrains)~=S+1
    error('length of mstrains or alphastrains incorrect')
end

if length(f)~=S
    error('length of f incorrect')
end


for i=1:S

N0(i)=( A*M*f(i) )/mstrains(i+1);          % N0(i) is the number of unfused cells of strain i.

end

[t,x] = ode45(@(t,x)unisexual_ODE_Multistrain(t,x,alphastrains,S),[0,T],[N0,zeros(1,nchoosek(S+1 ,2))]);
