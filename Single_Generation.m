function f = Single_Generation(mstrains,alphastrains,x,S,beta,C)

% Simulates the frequency after each successive generation.
% x is the number of gametes of each type.
% \beta is the resistance to survival and C is the fusion cost as defined in the main text.

Mv=unisexual_ODE_Multistrain_SingleGen(alphastrains,S);  % writes the array Mv which is used to determine which ODE represents which progeny. e.g. when S=2 x(4) represents zygote containing strain s1 and s2.
f=zeros(S,1);

for k=1:S
   wk=0;   
   Mk = Mv(any(Mv==k,2),:);               % any rows with k in their 2nd or 3rd columns. the +1 is here because the first entry of mstrains is set to 0 for a purpose. If an entry in Mk is 0, then it adds mstrains(1). 
      for ii=1:length(Mk(:,1))


      if Mk(ii,3)==0 || Mk(ii,2)==0       % checks if the progeny is a zygote, if so we multiply by (1-C).
      wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) );    
      else
      wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) )*(1-C);
      end
      
      end 
      w(k)=wk;                   % fitness of kth strain
      
end
w(S+1:end)=[];

for i=1:length(w)                % frequency of each strain
f(i)=w(i)/sum(w);   
end