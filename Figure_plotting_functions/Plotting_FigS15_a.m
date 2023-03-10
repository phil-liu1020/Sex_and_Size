function figS15output= Plotting_FigS15_a(alpha0,A,M,T,C,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity)

% For a fixed value of \beta_2, this function numerically simulates evolutionary trajectories starting at
% (m(0),\alpha(0))=(\beta_1,0) and
% (m(0),\alpha(0))=(P_1\beta_1+P_2\beta_2,0) to determine whether obligate
% fusion evolves from both initial conditions given various values of \beta_1.

% Output of 2 means obligate fusion evolves given (m(0),\alpha(0))=(\beta_1,0) and (m(0),\alpha(0))=(P_1\beta_1+P_2\beta_2,0)
% Output of 1 means obligate fusion evolves given (m(0),\alpha(0))=(\beta_1,0) only.
% Output of 0 means obligate fusion evolves given neither initial conditions.

beta1=(0.08:0.08:4);

fusion_env2beta1=zeros(1,length(beta1));
fusion_env1beta1=zeros(1,length(beta1));
fusion_env2betaBH=zeros(1,length(beta1));
fusion_env1betaBH=zeros(1,length(beta1));
figS15output=zeros(1,length(beta1));

addpath(genpath('Sex_and_Size-main'))


for i=1:length(beta1)

m0=beta1(i);

if beta1(i)==beta2
    continue
end

[~,~,m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1(i),beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity);

if (abs(m2(end)-beta2/4)<2*delta && alpha2(end)>4*delta) || alpha2(end)>4*delta
fusion_env2beta1(i)=1;
end

if (abs(m(end)-beta1(i)/4)<2*delta && alpha(end)>4*delta) || alpha(end)>4*delta
fusion_env1beta1(i)=1;
end

        if beta1(i)>beta2
        P1=lambda12/(lambda12+lambda21);
        P2=lambda21/(lambda12+lambda21);    
        else
        P1=lambda21/(lambda12+lambda21);
        P2=lambda12/(lambda12+lambda21);
        end

m0=P1*beta1(i)+P2*beta2;

[~,~,m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,C,beta1(i),beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax, switching_environments, return_traits,plasticity);

if (abs(m2(end)-beta2/4)<2*delta && alpha2(end)>4*delta) || alpha2(end)>4*delta
fusion_env2betaBH(i)=1;
end

if (abs(m(end)-beta1(i)/4)<2*delta && alpha(end)>4*delta) || alpha(end)>4*delta
fusion_env1betaBH(i)=1;
end




figS15output(i)=fusion_env2beta1(i)+fusion_env1beta1(i)+fusion_env2betaBH(i)+fusion_env1betaBH(i);


end

