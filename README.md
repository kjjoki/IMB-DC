# IMB-DC
The interactive multibundle method for constrained multi-objective nonsmooth DC optimization

MDBDC is a multiobjective double bundle solver (Fortran 95) for nonsmooth multiobjective programming by K. Joki and O. Montonen. MDBDC is able to handle problems having objective and constraint functions which can be presented as a difference of two convex (DC) functions. The method is a descent type and it generalizes the ideas of DBDC to multiobjective and constraint case. Solutions obtained are guaranteed to be weakly Pareto stationary.

The software utilizes OpenMP at each round of 'main iteration' to calculate subproblems in parallel. To turn down OpenMP, see instructions in tmpbdc.f95. In addition, there is a possibility to use simple stepsize determination after each 'main iteration'. The software uses code PLQDF1 by Prof. Ladislav Luksan to solve quadratic direction finding problem.

The software is free for academic teaching and research purposes but I ask you to refer the reference given below if you use it. To use the software modify tmdbdc.f95 and functions.f95 as needed. If you have any questions concerning the software, please contact directly the author Kaisa Joki (email: kjjoki@utu.fi).

# Codes include:        
                                                                                              
tmdbdc.f95         - Main program for MDBDC (this file)    

constants.f95      - Double precision (also some parameters)                       

bundle1.f95        - Bundle of DC component f_1                                    

bundle2.f95        - Bundle of DC component f_2                                    

functions.f95      - User-specified DC components f_1 and f_2 together with subgradients of DC components. Contains also user-specified initial values for parameters                                 

mdbdc.f95          - MDBDC method                                                  
                                                                                              
plqdf1.f           - Quadratic solver by Ladislav Luksan                           
                                                                                              
Makefile           - Makefile                                                      
                                                                                              
                                                                                
# References:                                                                        
                                                                                              
[1] Outi Montonen and Kaisa Joki: "Bundle-based descent method for nonsmooth multiobjective DC optimization with inequality constraint". Journal of Global Optimization, 72 (2018), 403–429. https://doi.org/10.1007/s10898-018-0651-0                                      
