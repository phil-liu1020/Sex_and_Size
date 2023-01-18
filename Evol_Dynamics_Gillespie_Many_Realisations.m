function [mres,alphares,muttime,meanswBG,meanswGB,randomG]=Evol_Dynamics_Gillespie_Many_Realisations(Nrealz,m0,alpha0,A,M,T,C,betaB,betaG,lambdaGB,lambdaBG,mu,NEVOL,f0,delta)
                                                         
% Obtains multiple realisations of the simulation for the coevolutionary dynamics between mass and \alpha in a system undergoing bet-hedging in a stochastically switching environment.
% Parameters: Nrealz - number of realisations
%             f0 - initial frequency of rare mutant

mres=zeros(Nrealz,0);
alphares=zeros(Nrealz,0);


parfor i=1:Nrealz
    
[mres(i,:),alphares(i,:),muttime(i,:),meanswBG(i,:),meanswGB(i,:),randomG(i,:)]=Evol_Dynamics_Mass_Alpha_Switching_Gillespie(m0,alpha0,A,M,T,C,betaB,betaG,lambdaGB,lambdaBG,mu,NEVOL,f0,delta);

end

assignin('base','mres',mres)
assignin('base','alphares',alphares)
assignin('base','muttime',muttime)
assignin('base','meanswBG',meanswBG)
assignin('base','meanswGB',meanswGB)

