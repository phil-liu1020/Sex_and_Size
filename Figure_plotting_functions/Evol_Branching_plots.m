function sb = Evol_Branching_plots(CellArray)

% Plots the mass of all the genotypes against time to visualise branching in mass.

sb=zeros(1,length(CellArray));
for j=1:round(length(CellArray)/10)-1
    
   a=CellArray{10*j,2};  
   b=CellArray{10*j,1};
   
   
   
    % This block of code checks if there's any duplicate genotypes
             %-------------------------------------------------------------
             aa=[a'];   
             [u,I,J] = unique(aa, 'rows', 'first');
             hasDuplicates = size(u,1) < size(aa,1);
             ixDupRows = setdiff(1:size(aa,1), I);
             dupRowValues = aa(ixDupRows,:);
             duprow=0;

             for i=1:length(aa(:,1))
             duprow(i)=isequal(aa(i,:),dupRowValues);
             end

             if hasDuplicates>0
             firstDupRow=find(duprow,1, 'first');
             c=b.*duprow';
             b(firstDupRow)=sum(c);
             duprow(firstDupRow)=0;
             b(duprow==1)=[];
             
             a(max(ixDupRows))=[];
             S=length(a)-1;
             end
             %-------------------------------------------------------------
   

   
   
   
   
   for i=1:length(a)   
   patch(10*j,a(i),b(i),'EdgeColor','none','Marker','o','MarkerFaceColor','flat');
   end
   sb(j)=sum(b);
end

colorbar
pbaspect([1 1/3 1])
set(gca,'fontsize',14)
box on
xlabel('\tau')
if strcmp(inputname(1),'genotypesData_m')==1
ylabel('m')
elseif strcmp(inputname(1),'genotypesData_alpha')==1
ylabel('\alpha')
else
end
