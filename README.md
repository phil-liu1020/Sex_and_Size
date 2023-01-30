Evolutionary_trajectories:
This is the function for simulating evolutionary trajectories for all case scenarios. It has the option to produce a single realisation or multiple independent realisations. This function can call either "Evolutionary_trajectories_Many_realisations" or "Evolutionary_trajectories_1realisation".
              
              Evolutionary_trajectories_Many_realisations:
              This function produces multiple realisations of the evolutionary trajectory by calling the function "Evolutionary_trajectories_1realisation". This                       function is applicable for simulating evolutionary trajectories under scenarios with plasticity, bet-hedging, fixed environment and alphamax. This does not               produce data on individual traits within each realisation of the simulation (i.e. cannot be used to visualise branching).
              
              Evolutionary_trajectories_1realisation:
              This function produces one single realisation of the evolutionary trajectory. This function is applicable for simulating trajectories for scenarios with                 plasticity, bet-hedging, fixed environment and alphamax. This fucntion has the option to output data on individual traits within each realisation of the                 simulation, allowing us to visualise branching. This function calls "Invasion_dynamics" and runs it until the next mutation/environmental switching event.
              
                                                   Invasion_dynamics:
                                                   The invasion dynamics function contains a loop that iterates over the number of generations G specified in the input                                                      parameter. G can represent the number of generations until the next mutation/environmental switching event. The                                                          loop takes as an input the mass, alpha and frequencies of each genotype and runs the "fertilisation_kinetics" function                                                    (which depends on the alpha value of each genotype) followed by the survival function (which depends on the mass of                                                      each genotype). The loop ends each iteration with the new frequency of each genotype in the population.
                                                   
                                                                    Survival_function: This function takes as input the parameters A, M, T, C, \beta as well as mass,                                                                         alpha and freqnency of each genotype in the population and calculates the frequency in the subsequent                                                                     generation.
                                                                    
                                                                    Fertilisation_kinetics:
                                                   
