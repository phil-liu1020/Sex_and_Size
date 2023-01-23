function m0 = Invasion_Dynamics_Mass_FRTI_FixedSw(betaB,betaG,tsG,tsB,C,alpha,A,M,T,f0,delta,m0,plots,G)

% Parameters: betaB and betaG are the resistance to survival in good and bad environments respectively.
%             tsG and tsB are the times spent in good and bad environments respectively.
%             f0 is the initial frequency of mutant strain.
%             theta1 is the change in frequency 
%             m0 is the initial mass
%             G is the number of generations we run the invasion dynamics for.

Nsw=1;       % Nsw is the number of switches, defined by the number of jumpos between the two for loops.
tg=1;        % tg represents which generation the system is in.
Ngen=1;      % Ngen is the number of generations.
Good=1; Bad=0;  % We have chosen the simulation to start in the good environment by random choice.

g=0; 

                                                          
while Ngen<G

    if Bad==1

    for i=1:tsB
    Nres0=(A*M*(1-f0))/m0; 
    Nmut0=(A*M*f0)/(m0+delta);

    [t,x] = ode45(@(t,x)unisexual_ODE(t,x,alpha),[0,T],[Nres0,Nmut0,0,0,0]);
    wmut=x(length(x(:,1)),5)*exp(-betaB/(2*(m0+delta)))*(1-C)+x(length(x(:,1)),4)*exp(-betaB/(2*m0+delta))*(1-C)+x(length(x(:,1)),2)*exp(-betaB/(m0+delta)); 
    wres=x(length(x(:,1)),3)*exp(-betaB/(2*m0))*(1-C)+x(length(x(:,1)),4)*exp(-betaB/(2*m0+delta))*(1-C)+x(length(x(:,1)),1)*exp(-betaB/m0);
    fmut=wmut./(wmut+wres);


    f0=fmut;
    tg = tg + 1;
    if tg<=length(g)
    g(tg)=f0;
    else
    g(end+1) = f0; 
    end


    Ngen=Ngen+1;
    end
    Good=1;Bad=0;
    Nsw=Nsw+1;



    elseif Good==1
   
    
    for i=1:tsG
    Nres0=(A*M*(1-f0))/m0; 
    Nmut0=(A*M*f0)/(m0+delta);

    [t,x] = ode45(@(t,x)unisexual_ODE(t,x,alpha),[0,T],[Nres0,Nmut0,0,0,0]);
    wmut=x(length(x(:,1)),5)*exp(-betaG/(2*(m0+delta)))*(1-C)+x(length(x(:,1)),4)*exp(-betaG/(2*m0+delta))*(1-C)+x(length(x(:,1)),2)*exp(-betaG/(m0+delta)); 
    wres=x(length(x(:,1)),3)*exp(-betaG/(2*m0))*(1-C)+x(length(x(:,1)),4)*exp(-betaG/(2*m0+delta))*(1-C)+x(length(x(:,1)),1)*exp(-betaG/m0);
    fmut=wmut./(wmut+wres);
 



    f0=fmut;
    tg = tg + 1;
    if tg<=length(g)
    g(tg)=f0;
    else
    g(end+1) = f0; 
    end

    Ngen=Ngen+1;
    end
    Nsw=Nsw+1;

    Bad=1; Good=0;
    end

 
   
end



if plots==1
    plot(g)
    
    xlabel('t_g')
    ylabel('f_{m}')
else

end