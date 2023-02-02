% Instructions for producing FigS5 to FigS14.

% First, open "Data_Files", then open the folder for the figure you'd like to plot.
% Next, load all the mat files in the folder you've opened.

% To reproduce all the panels of a given figure, execute the command

cd ..
cd ..
addpath(genpath('Sex_and_Size-main'))
cd Sex_and_Size-main/Figure_generation_scripts
  
Evol_Branching_plots(genotypesData_m);

figure

Evol_Branching_plots(genotypesData_alpha);

figure

Evol_Branching_plots_Coevolution(genotypesData_m,genotypesData_alpha);
