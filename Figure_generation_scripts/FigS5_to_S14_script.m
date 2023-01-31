% Instructions for producing FigS5 to FigS14.

% First, open "Data_Files", then open the folder for the figure you'd like to plot.
% Next, load all the mat files in the folder you've opened.

% To reproduce all the panels of a given figure, execute the command
  
Evol_Branching_plots(StrainsData_m);

figure

Evol_Branching_plots(StrainsData_alpha);

figure

Evol_Branching_plots_Coevolution(StrainsData_m,StrainsData_alpha);
