function sb = Evol_Branching_plots_Coevolution(CellArray,CellArray2)

% Plots the mass of all the genotypes against their respective values for
% \alpha to visualise the coevolution between mass and \alpha.

% This enables us to see whether macrogametes evolve small fusion rates and vice versa for microgametes i.e. Oogamy.

sb=zeros(1,length(CellArray));
for j=1:round(length(CellArray)/10)-1
    
   a=CellArray{10*j,2};  
   b=CellArray{10*j,1};
   cccc=CellArray2{10*j,2};
   
   
   
    % This block of code checks if there's any duplicate genotypes
             %-------------------------------------------------------------
            % aa=[a' cccc'];   
            % [u,I,J] = unique(aa, 'rows', 'first');
            % hasDuplicates = size(u,1) < size(aa,1);
            % ixDupRows = setdiff(1:size(aa,1), I);
            % dupRowValues = aa(ixDupRows,:);
            % duprow=0;

            % for i=1:length(aa(:,1))
            % duprow(i)=isequal(aa(i,:),dupRowValues);
            % end

           % if hasDuplicates>0
            % firstDupRow=find(duprow,1, 'first');
            % c=b.*duprow';
            % b(firstDupRow)=sum(c);
            % duprow(firstDupRow)=0;
            % b(duprow==1)=[];
             
           %  a(max(ixDupRows))=[];
            % S=length(a)-1;
            % end
             %-------------------------------------------------------------
   

   
   
   
   
   for i=1:length(a)   
   patch(a(i),cccc(i),b(i),'EdgeColor','none','Marker','o','MarkerFaceColor','flat');
   hold on
   end
   sb(j)=sum(b);
end

colorbar
pbaspect([1 1 1])
set(gca,'fontsize',14)
box on
xlabel('m')
ylabel('\alpha')
