function f = Single_Generation(mstrains,alphastrains,x,S,beta,C)

% Simulates the frequency after each successive generation.
% x, the output of the function "Fertilisation_Kinetics", is a vector of length S+Binomial(S+1,2) that represents the number of cells of each type potentially present at the end of the fusion period. This includes S types of unfused cells (x(1),...,x(S)) (one for each genotype) and Binomial(S,2) fused cells (x(S+1),...,x(S(S+3)/2)) (one for each pair of genotypes in a fused cell).
% \beta is the resistance to survival and C is the fusion cost as defined in the text.
% mstrains and alphastrains are the mass and fusion rates of each strain.

Mv=unisexual_ODE_Multistrain_SingleGen(alphastrains,S);  % writes the array Mv, a (S+Binomial(S,2)) x 3 matrix. The first column of Mv is the enumeration of the cell type produced at the end of the fusion period (corresponding to the enumeration of x). The second and third columns are the strains contributing to fused cells recorded in the ith row. If the third column is zero, this indicates the cell is unfused
f=zeros(S,1);

for k=1:S          % This for loop with index k calculates the fitness of each strain in the system i.e. w_k for k in [1,S] as defined in the text.
   wk=0; 
   
        -------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %%%This block of code computes the absolute fitness of the k-th strain.
   
   
        Mk = Mv(any(Mv==k,2),:);        % Mk is just a subset of Mv formed of only rows that contains k-th strain. E.g. if Mv = [1 1 0;2 2 0;3 1 1;4 1 2;5 2 2], then for k=1, Mk = [1 1 0;3 1 1;4 1 2].     
           for ii=1:length(Mk(:,1))     % This for loop with index ii calculates the absolute fitness of the k-th strain.

              if Mk(ii,3)==0 || Mk(ii,2)==0       % checks if the ii-th column of Mk represents a zygote, if so we multiply by (1-C).
              wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) );    % any rows with k in their 2nd or 3rd columns. the +1 is here because the first entry of mstrains is set to 0 for a purpose. If an entry in Mk is 0, then it adds mstrains(1).  
              else
              wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mstrains(Mk(ii,2)+1) + mstrains(Mk(ii,3)+1) ) )*(1-C);
              end
       
           end 
         ----------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      w(k)=wk;                   % writes the absolute fitness of k-th strain into a vector w. 
end

for i=1:length(w)                % frequency of each strain in the subsequent generation
f(i)=w(i)/sum(w);   
end
