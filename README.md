# IMB-DC
The interactive multibundle method for constrained multi-objective nonsmooth DC optimization

IMB-DC is an interactive multibundle solver (Fortran 95) for constrained nonsmooth multiobjective programming by Kaisa Joki and Outi Montonen. IMB-DC is able to handle problems having objective and constraint functions which can be presented as a difference of two convex (DC) functions. The method is a descent type and it provides an option for the user to steer the solution process. Solutions obtained are guaranteed to be weakly Pareto stationary.

The software uses code PLQDF1 by Prof. Ladislav Luksan to solve quadratic direction finding problem.

The software is free for academic teaching and research purposes but I ask you to refer the reference given below if you use it. To use the software modify timbdc.f95, imbdc.f95 and functions.f95 as needed. If you have any questions concerning the software, please contact directly the author Kaisa Joki (email: kjjoki@utu.fi).

# Codes include:        
                                                                                              
timbdc.f95         - Main program for IMB-DC    

constants.f95      - Double precision (also some parameters)                       

bundle1.f95        - Bundle of DC component f_1                                    

bundle2.f95        - Bundle of DC component f_2                                    

functions.f95      - User-specified DC components f_1 and f_2 together with subgradients of DC components. Contains also user-specified initial values for parameters                                 

imbdc.f95          - IMB-DC method and actual values for paramters                                                  
                                                                                              
plqdf1.f           - Quadratic solver by Ladislav Luksan                           
                                                                                              
Makefile           - Makefile                                                      
                                                                                              
                                                                                
# References:                                                                        
                                                                                              
[1] Outi Montonen and Kaisa Joki: "Interactive multibundle method for constrained multiobjective DC optimization". (manuscript)                                     
