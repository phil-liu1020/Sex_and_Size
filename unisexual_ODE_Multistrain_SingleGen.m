function M = unisexual_ODE_Multistrain_SingleGen(alphastrains,S)

% x is the ODE variable. x(1) to x(S) represent dN1/dt to dNS/dt, where N
% is the number of unfertilised gametes.

% The dimension of vector x is (S+nchoosek(S+1 ,2),1)

% S is the number of strains in the system.
% The output M is an array with first row giving the row number of the array. The second and third rows give information on which strain each row of x (each ODE) contains. If a row of x represents a cell with strain 1 fused with strain 2, this corresponding row of M will contain a 1 and 2 in its 2nd and 3rd columns.   


k=1;
M=zeros(1,3);

for i=1:S
  
  M(i,1)=i;
  M(i,2)=i;
  
end


for i=1:S
   for j=1:S

       if j<i
           continue
       else
       end

   M(k+S,1)=k+S;
   M(k+S,2)=i;
   M(k+S,3)=j;
   k=k+1;
   end
end
