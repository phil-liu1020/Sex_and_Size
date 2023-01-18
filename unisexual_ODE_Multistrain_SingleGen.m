function M = unisexual_ODE_Multistrain_SingleGen(alphastrains,S)

% x is the ODE variable. x(1) to x(S) represent dN1/dt to dNS/dt, where N
% is the number of unfertilised gametes.

% The dimension of vector x is (S+nchoosek(S+1 ,2),1)

% S is the number of strains in the system.



x = zeros(S+nchoosek(S+1 ,2),1);
k=1;
M=zeros(1,3);

for i=1:S
  dx(i) = 0;

  for j=1:S
 
  dx(i)=dx(i)+(-0.5*( alphastrains(i+1)+alphastrains(j+1) )*x(i)*x(j));                   % This evaluates the expression dNi/dt=-alpha(Ni^2+\sum^{n}_{j\neq i} N_iN_j)
  
  end
  
  M(i,1)=i;
  M(i,2)=i;
  
end


for i=1:S
   for j=1:S

       if j<i
           continue
       else
       end

   dx(k+S) = 0.5*(  alphastrains(i+1)+alphastrains(j+1)  )*x(i)*x(j);                                 % This evaluates the expression dFi/dt=alpha\sum^{n}_{j\neq i} N_iN_j
   M(k+S,1)=k+S;
   M(k+S,2)=i;
   M(k+S,3)=j;
   k=k+1;
   end
end
