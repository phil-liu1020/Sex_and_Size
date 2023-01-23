function f = Single_Generation(mstrains,alphastrains,x,S,beta,C)

% Simulates the frequency after each successive generation.
% x is the output of the function "Fertilisation_Kinetics" which represents the number of cells of each type e.g. number of unfused cells of strain 1, number of strain 1 cells that have fused with strain 3 etc.
% \beta is the resistance to survival and C is the fusion cost as defined in the main text.

Mv=unisexual_ODE_Multistrain_SingleGen(alphastrains,S);  % writes the array Mv which is used to determine what each column of x (each ODE) represents. e.g. when S=2, x(4) represents F^{1}_{2} the number of strain 1 cells that've fused with strain 2 cells.
f=zeros(S,1);

for k=1:S          % This for loop with index k calculates the fitness of each strain i.e. w_k for k in [1,S] as defined in the text.
   wk=0;   
   Mk = Mv(any(Mv==k,2),:);        % Mk outputs an array showing which column of x (which ODE) contains strain k cells.     
      for ii=1:length(Mk(:,1))

         if Mk(ii,3)==0 || Mk(ii,2)==0       % checks if the progeny is a zygote, if so we multiply by (1-C).
         wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) );    % any rows with k in their 2nd or 3rd columns. the +1 is here because the first entry of mstrains is set to 0 for a purpose. If an entry in Mk is 0, then it adds mstrains(1).  
         else
         wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) )*(1-C);
         end
      
      end 
      w(k)=wk;                   % fitness of kth strain
      
end

for i=1:length(w)                % frequency of each strain
f(i)=w(i)/sum(w);   
end
