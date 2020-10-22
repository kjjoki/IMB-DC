      MODULE imbdc
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                IMB-DC - THE INTERACTIVE MULTIBUNDLE METHOD FOR                   | |
        !| |                     CONSTRAINED NONSMOOTH DC OPTIMIZATION                        | | 
        !| |                                 (version 1)                                      | |
        !| |                                                                                  | |            
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    Utilizes the new version of PLQDF1 by Ladislav Luksan as a quadratic solver.      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
 
       
        USE constants, ONLY : dp    ! Double precision (i.e. accuracy)
        USE bundle1                 ! The BUNDLE of the DC component f_1
        USE bundle2                 ! The BUNDLE of the DC component f_2
        USE functions               ! Contains INFORMATION from the USER
        
        IMPLICIT NONE    
        
        EXTERNAL PLQDF1             ! The QUADRATIC SOLVER by Ladislav Luksan
        
        CONTAINS
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      new_algorithm       : The new interactive doublebundle method for           | |
        !| |                            nonsmooth multiobjective DC optimization.             | | 
        !| |                                                                                  | |       
        !| |      separate_direction  : The double bundle method to determine                 | | 
        !| |                            the descent direction for one objective.              | |
        !| |                                                                                  | |       
        !| |      common_direction    : The common search direction determination             | |       
        !| |                                                                                  | |       
        !| |      escape_procedure_d  : The algorithm guaranteeing Clarke stationarity        | |       
        !| |                            for one objective                                     | |
        !| |                                                                                  | |
        !| |      escape_procedure_H  : The algorithm guaranteeing Clarke stationarity        | |       
        !| |                            for the improvement function.                         | |
        !| |                                                                                  | |       
        !| |      quadratic_solver    : The solver for the quadratic norm minimization        | |       
        !| |                            problem. Needed in the 'Escape procedure'             | |
        !| |                                                                                  | |       
        !| |      subproblem_solver   : The solver for subproblems in the search direction    | |
        !| |                            problem. Needed in the separate_direction algorithm.  | |   
        !| |                                                                                  | |   
        !| |      subproblem_solver_const : The solver for subproblems in the search          | |
        !| |                                direction problem when constraints are used.      | |
        !| |                                Needed in the main iteration algorithm.           | |   
        !| |                                                                                  | |   
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        
    
       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                THE NEW ALGORITHM FOR NONSMOOTH DC OPTIMIZATION                 |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE new_algorithm( x_solution, f_solution, mit,&
                            & mrounds, mrounds_escape, termination, counter,  &
                            & iprint, agg_used, scale_func, t_adjust, use_reset_b1, use_reset_b2, &
                            & max_null_step0, &
                            & user_n, startpoint, number_of_obj, number_of_const, number_of_func, & 
                            & time, i_user0, counter_g, use_const_in_sep_d, scaling_const)
            ! 
            ! Solves the constrained nonsmooth multiobjective DC minimization problem
            !
            ! INPUT:  * 'mit'            : The maximum number of 'main iterations'      
            !         * 'mrounds'        : The maximum number of rounds during the descent direction determination for the objective f_i
            !         * 'mrounds_escape' : The maximum number of rounds during one 'Escape procedure' 
            !         * 'iprint'         : Specifies the print          
            !         * 'agg_used  '     : Specifies whether aggregation is used during separate descent direction determination         
            !         * 'scale_func'     : Specifies whether scaling of objective functions is used          
            !         * 't_adjust'       : Specifies whether proxmity paramter adjusting is used          
            !         * 'max_null_step0' : The maximum number of consecutive null steps          
            !         * 'use_reset_b1'   : Specifies whether bundles B1 are reset after serious step and escape procedure     
            !         * 'use_reset_b2'   : Specifies whether bundles B2 are reset after serious step and escape procedure     
            !         * 'use_const_in_sep_d'   : Specifies whether constraints are used in the individual search direction problem   
            !         * 'scaling_const'  : Specifies the scaling constraint used in the individual search direction problem   
            !         * 'user_n'         : The dimension of the problem          
            !         * 'startpoint'     : Specifies the starting point used        
            !         * 'number_of_obj'  : The number of objective functions        
            !         * 'number_of_const': The number of constants functions        
            !         * 'number_of_func' : The sum of objective and constants functions  
            !         * 'i_user0'        : The maximum number of times the user can steer the algorithm
            !
            ! OUTPUT: * 'x_solution' : The solution obtained to the minimizatin problem
            !         * 'f_solution' : The objective function values at the solution 'x_solution'
            !         * 'termination': The cause of the termination in the algorithm
            !         * 'counter'    : Gives the values of different counters for objectives
            !         * 'counter_g'  : Gives the values of different counters for constraints
            !         * 'time'       : The CPU time used
            !
            ! NOTICE: * The dimension of vectors 'x_solution' has to be 'user_n' defined by USER in MODULE functions.
            !         * The dimension of the vector 'counter' has to be DIMENSION(12,number_of_obj).
            !         * 'mit', 'mrounds' and 'mrounds_escape' have to be integers.
            !         * IF ('mit' <= 0) THEN DEFAULT value 1000 is used
            !         * IF ('mrounds' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * IF ('mrounds_escape' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * 'iprint' has to be -4, -3, -2, -1, 0, 1, 2, 3 or 4 (5). If it is NOT then DEFAULT value 1 is used. 
            !         * If 'iprint = 5' then everything is printed step by step (this is NOT recommended!!!)
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER ************************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_solution  ! the solution obtained to the problem
               
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: f_solution  ! the objective function values at the solution 'x_solution'
               REAL(KIND=dp),  INTENT(OUT) :: time                     ! the CPU time
               
               REAL(KIND=dp),  INTENT(IN) :: scaling_const             ! the scaling constraint used in the individual search direction problem
                              
               INTEGER, INTENT(INOUT) :: mit                  ! the maximum number of 'main iterations'
               INTEGER, INTENT(INOUT) :: mrounds              ! the maximum number of rounds during separate descent direction determination for the objective f_i
               INTEGER, INTENT(INOUT) :: mrounds_escape       ! the maximum number of rounds during one 'Escape procedure' for the improvement function and for the separate descent direction determination
               
               INTEGER, INTENT(OUT) :: termination        ! 1 - the stopping condition is satisfied (i.e. approximate Clarke stationarity for the improvement function)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < beta_0 in Escape procedure)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed during separate descent direction determination for the objective f_i
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_escape' of rounds executed in one 'Escape procedure'
                                                          
               INTEGER, DIMENSION(:,:), INTENT(OUT) :: counter  ! contains the values of different counteres: 
                                                                !   counter(1,1) = iter_counter:         the number of 'main iterations' executed during the new method 
                                                                !   counter(2,j) = subprob_counter:      the number of subproblems solved during all separate descent direction determinations for f_j 
                                                                !   counter(3,j) = f_counter:            the number of function values evaluated for the objective f_j (the same value holds for the DC components also)
                                                                !   counter(4,j) = subgrad1_counter:     the number of subgradients calculated for the first DC component of f_j (i.e. f_1)
                                                                !   counter(5,j) = subgrad2_counter:     the number of subgradients calculated for the second DC component of f_j (i.e f_2)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(6,1) = stop_cond_counter_d:    the number of times 'Escape procedure' is executed for the separate decsent direction determination
                                                                !   counter(7,1) = escape_f_counter_d:     the number of function values evaluated for each objective f_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                !   counter(8,1) = escape_sub_counter_d:   the number of subgradients caluculated for each objective f_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(9,1) = stop_cond_counter_H:    the number of times 'Escape procedure' is executed for the improvement function
                                                                !   counter(10,1) = escape_f_counter_H:    the number of function values evaluated for each objective f_i (in 'Escape procedure' for the improvement function)
                                                                !   counter(11,1) = escape_sub_counter_H:  the number of subgradients caluculated for each objective f_i (in 'Escape procedure' for the improvement function)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(12,1) = f_counter_line:      the number of function values evaluated for the objective f_j during line search (the same value holds for the DC components also)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                
               INTEGER, DIMENSION(:,:), INTENT(OUT) :: counter_g  ! contains the values of different counteres: 
                                                                  !   counter_g(1,j) = g_counter:              the number of function values evaluated for the constraint g_j (the same value holds for the DC components also)
                                                                  !   counter_g(2,j) = subgrad1_g_counter:     the number of subgradients calculated for the first DC component of g_j (i.e. g_1)
                                                                  !   counter_g(3,j) = subgrad2_g_counter:     the number of subgradients calculated for the second DC component of g_j (i.e g_2)
                                                                  !   counter_g(4,j) = g_counter_line:         the number of function values evaluated for the constraint g_j during line search (the same value holds for the DC components also)
                                                                  !--------------------------------------------------------------------------------------------------------------------------     
                                                                  !   counter_g(5,i) = escape_g_counter_d:       the number of function values evaluated for each constraint g_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                  !   counter_g(6,i) = escape_sub_g_counter_d:   the number of subgradients caluculated for each constraint g_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                  !--------------------------------------------------------------------------------------------------------------------------     
                                                                  !   counter_g(7,1) = escape_f_counter_H:    the number of function values evaluated for each constraint g_i (in 'Escape procedure' for the improvement function)
                                                                  !   counter_g(8,1) = escape_sub_counter_H:  the number of subgradients caluculated for each constraint g_i (in 'Escape procedure' for the improvement function)
                                                                  !-----------------------------------------------------------------------------------------------------------------------

               LOGICAL, INTENT(IN) :: agg_used          ! .TRUE. if aggragation is used during the separate descent direction determination. Otherwise .FALSE.                                                             
               LOGICAL, INTENT(IN) :: scale_func        ! .TRUE. if scaling procedure is used in the algorithm. Otherwise .FALSE.       
               LOGICAL, INTENT(IN) :: t_adjust          ! .TRUE. if the proximity parameter is adjusted in algorithm. Otherwise .FALSE.                        
               LOGICAL, INTENT(IN) :: use_reset_b1      ! .TRUE. if bundles B1 are reset after serious step and escape procedure. Otherwise .FALSE.                        
               LOGICAL, INTENT(IN) :: use_reset_b2      ! .TRUE. if bundles B2 are reset after serious step and escape procedure. Otherwise .FALSE.                        
               LOGICAL, INTENT(IN) :: use_const_in_sep_d ! .TRUE. if constraints are used in the individual search direction problem. Otherwise .FALSE.                        
                        
               INTEGER, INTENT(INOUT) :: iprint ! variable that specifies print option:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result 
                                                !   iprint = -1: basic print of final result (without the solution vector)
                                                !   iprint = 2 : extended print of final result 
                                                !   iprint = -2: extended print of final result (without the solution vector)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results
                                                !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                                !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                                !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                                !   iprint = 5 : prints each step of the bundle algorithm (i.e. everything is printed step by step) 
               
               INTEGER, INTENT(IN) :: user_n         ! The dimension of the problem
               INTEGER, INTENT(IN) :: startpoint     ! The starting point
               INTEGER, INTENT(IN) :: max_null_step0 ! The maximum number of consecutive null steps 

               INTEGER, INTENT(IN) :: number_of_obj      ! The number of objectives
               INTEGER, INTENT(IN) :: number_of_const    ! The number of constraints
               INTEGER, INTENT(IN) :: number_of_func     ! The sum of objectives and constraints
               
               INTEGER, INTENT(IN) :: i_user0           ! The maximum number of times the user can steer the algorithm
           
           !***************************** LOCAL VARIABLES ************************************  
            
             ! ** Bundles **
               TYPE(kimppu1), DIMENSION(number_of_obj) :: B1    ! The bundles B_1 for the DC components f_1 of the objectives
               TYPE(kimppu2), DIMENSION(number_of_obj) :: B2    ! The bundles B_2 for the DC components f_2 of the objectives
             !-------------------------------------------               

             ! ** Different points and function values ** 
               REAL(KIND=dp), DIMENSION(user_n) :: x_0                  ! the starting point
               REAL(KIND=dp), DIMENSION(user_n) :: x_current            ! the current iteration point (the dimension 'user_n' is the number of variables)
               REAL(KIND=dp), DIMENSION(user_n) :: x_new                ! the new iteration point (obtained from the previous 'main iteration')

               REAL(KIND=dp), DIMENSION(number_of_obj) :: f_0_all       ! the values of f_i at x_0   
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_k_all    ! the values of the DC components f1 and g1 at x_current      
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_k_all    ! the values of the DC components f2 and g2 at x_current  
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg_k_all     ! the values of the DC components f and g at x_current 
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_new_all  ! the values of the DC components f1 and g1 at x_new        
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_new_all  ! the values of the DC components f2 and g2 at x_new  
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg_new_all   ! the values of the DC components f and g at x_new 

               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg1_k_all_sep    ! the values of the DC components f1 and g1 at x_current in individual search direction determination      
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg2_k_all_sep    ! the values of the DC components f2 and g2 at x_current in individual search direction determination  
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg_k_all_sep     ! the values of the DC components f and g at x_current in individual search direction determination  

               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg1_new_all_sep  ! the values of the DC components f1 and g1 at x_new in individual search direction determination        
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg2_new_all_sep  ! the values of the DC components f2 and g2 at x_new in individual search direction determination  
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: fg_new_all_sep   ! the values of the DC components f and g at x_new in individual search direction determination 
               
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad1_k_all         ! the subgradients of the DC components f1 and g1 at x_k 
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad2_k_all         ! the subgradients of the DC components f2 and g2 at x_k 
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad1_new_all       ! the subgradients of the DC components f1 and g1 at x_new   
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad2_new_all       ! the subgradients of the DC components f2 and g2 at x_new                  
               
               REAL(KIND=dp), DIMENSION(1+number_of_const, user_n) :: grad1_k_all_sep         ! the subgradients of the DC components f1 and g1 at x_k in individual search direction determination 
               REAL(KIND=dp), DIMENSION(1+number_of_const, user_n) :: grad2_k_all_sep         ! the subgradients of the DC components f2 and g2 at x_k in individual search direction determination 
               REAL(KIND=dp), DIMENSION(1+number_of_const, user_n) :: grad1_new_all_sep       ! the subgradients of the DC components f1 and g1 at x_new in individual search direction determination   
               REAL(KIND=dp), DIMENSION(1+number_of_const, user_n) :: grad2_new_all_sep       ! the subgradients of the DC components f2 and g2 at x_new in individual search direction determination   

               
               REAL(KIND=dp), DIMENSION(number_of_func) :: f1_new       ! the values of the DC components for objectives after the separate descent direction determination
               REAL(KIND=dp), DIMENSION(number_of_func) :: f2_new       ! qthe values of the DC components for objectives after the separate descent direction determination
 
               REAL(KIND=dp), DIMENSION(user_n) :: y                    ! a help point used during stepsize determination
               REAL(KIND=dp) :: f_y                                     ! function value of f at y and at x_current
               REAL(KIND=dp) :: f1_y, f2_y                              ! function value of f_1 and f_2 at y
               REAL(KIND=dp) :: g1_y, g2_y                              ! function value of g_1 and g_2 at y
               REAL(KIND=dp) :: alpha                                   ! the linearization error
               REAL(KIND=dp) :: lin_err1, lin_err2                      ! the linearization errors
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: vf_y          ! vector of function value of objectives at y
             !-------------------------------------------               
    
             ! ** Subgradients **
               REAL(KIND=dp), DIMENSION(user_n) :: grad1, grad2   ! subgradients (the dimenison 'user_n' is the length of subgradient)       
             !-------------------------------------------               
 
             ! ** Values related to different directions ** 
               REAL(KIND=dp), DIMENSION(user_n) :: separate_d     ! the search direction obtained for a separate objective
               REAL(KIND=dp), DIMENSION(user_n) :: common_d       ! the common search direction obtained 
               REAL(KIND=dp), DIMENSION(user_n) :: dd             ! the direction used to update the bundles
               REAL(KIND=dp), DIMENSION(user_n) :: d_cause        ! the direction send to escape procedure

               REAL(KIND=dp), DIMENSION(user_n*number_of_obj) :: mD       ! the vector containing separate search directions

               REAL(KIND=dp) :: common_obj                        ! the objective function value when the common direction is determined
             !-------------------------------------------               
              
             ! ** Stepsizes **
               REAL(KIND=dp), DIMENSION(number_of_func) :: f_best    ! function value of f_1 and f_2 at y 
               REAL(KIND=dp) :: tau                                 ! A stepsize 'tau' for each objective into the common direction
               REAL(KIND=dp) :: tau_smallest                        ! The smallest stepsize 'tau' for objectives into the common direction
               REAL(KIND=dp) :: tol_tau                             ! The tolerance value for the step-size 'tau'
              
               LOGICAL :: improvement                               ! If .TRUE. all objectives decrase compared to the previous value
               LOGICAL :: inc_in_tau                                ! If .TRUE. then stepsize is increased during previous round
             !-------------------------------------------               
               
             ! ** Control of termination **
               LOGICAL :: stop_alg       ! .TRUE. if the proximal bundle algorithm can be terminated 
               LOGICAL :: step_stop      ! .TRUE. if stepsize determination can be stopped for objective
             !-------------------------------------------               

             ! ** Counters for the objective functions **
               INTEGER, DIMENSION(number_of_obj)  :: f_counter           ! the number of function values evaluated for a DC components f_1              
               INTEGER, DIMENSION(number_of_const)  :: g_counter           ! the number of function values evaluated for a DC components f_1              
               INTEGER, DIMENSION(number_of_obj)  :: f_counter_line      ! the number of function values evaluated for objectives in line search          
               INTEGER, DIMENSION(number_of_const):: g_counter_line    ! the number of function values evaluated for constraints in line search          
               INTEGER, DIMENSION(number_of_obj)  :: subgrad1_counter    ! the number of subgradients calculated for f_1 
               INTEGER, DIMENSION(number_of_const)  :: subgrad1_g_counter    ! the number of subgradients calculated for f_1 
               INTEGER, DIMENSION(number_of_obj)  :: subgrad2_counter    ! the number of subgradients calculated for f_2  
               INTEGER, DIMENSION(number_of_const)  :: subgrad2_g_counter    ! the number of subgradients calculated for f_2  
               INTEGER, DIMENSION(number_of_obj)  :: subprob_counter     ! the number of subproblems solved during separate descent direction determinations for objectives  
               
               INTEGER, DIMENSION(number_of_obj)  :: escape_counter_d   ! the number of escape procedure used for separate objectives
               INTEGER, DIMENSION(number_of_obj)  :: escape_f_counter_d    ! the number of escape procedure used for separate objectives
               INTEGER, DIMENSION(number_of_obj)  :: escape_sub_counter_d  ! the number of escape procedure used for separate objectives
               INTEGER, DIMENSION(number_of_const)  :: escape_g_counter_d    ! the number of escape procedure used for separate objectives
               INTEGER, DIMENSION(number_of_const)  :: escape_sub_g_counter_d  ! the number of escape procedure used for separate objectives
             !-------------------------------------------               
               
             ! ** Other counters for the new algorithm **
               INTEGER :: iter_counter           ! the number of rounds executed 
               
               INTEGER :: escape_counter_H       ! the number of times 'Escape procedure' is used during the algorithm for the improvement function
               INTEGER :: escape_f_counter_H     ! the number of function values evaluated for the improvement function in 'Escape procedure' during the execution of the new method
               INTEGER :: escape_sub_counter_H   ! the number of subgradients evaluated for improvement function in 'Escape procedure'  during the execution of the new method 
               INTEGER :: escape_g_counter_H     ! the number of function values evaluated for the improvement function in 'Escape procedure' during the execution of the new method
               INTEGER :: escape_sub_g_counter_H ! the number of subgradients evaluated for improvement function in 'Escape procedure'  during the execution of the new method           
             !-------------------------------------------               
             
             ! ** Improvement related stuff **
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_scale_all       ! the values of the DC components f1 and g1 at x_current  (SCALED!)    
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_scale_all       ! the values of the DC components f2 and g2 at x_current  (SCALED!)           
               REAL(KIND=dp), DIMENSION(number_of_func) :: AandB_k_all         ! the values of the components A_i and B_l at x_current                  
               
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: AandB_k_all_sep   ! the values of the components A_i and B_l at x_current in individual search direction determination 
               REAL(KIND=dp), DIMENSION(1+number_of_const) :: AandB_new_all_sep ! the values of the components A_i and B_l at new point in individual search direction determination 
               
               REAL(KIND=dp), DIMENSION(number_of_obj,number_of_const) :: AB_const_k 
               REAL(KIND=dp), DIMENSION(number_of_obj,number_of_const) :: AB_const_new 
               REAL(KIND=dp), DIMENSION(number_of_obj) :: AB_obj_k 
               REAL(KIND=dp), DIMENSION(number_of_obj) :: AB_obj_new 
               
               REAL(KIND=dp) :: H1_current, H2_current                         ! the values of the DC components of the improvement function at x_current                  
               REAL(KIND=dp), DIMENSION(number_of_obj) :: H1_current_sep, H2_current_sep ! the values of the DC components of the improvement function at x_current in individual search direction determination                 
               REAL(KIND=dp), DIMENSION(number_of_obj) :: H1_new_sep, H2_new_sep ! the values of the DC components of the improvement function at x_current in individual search direction determination                 
               REAL(KIND=dp) :: H1_new, H2_new                                 ! the values of the DC components of the improvement function at x_current                  
             !-------------------------------------------               
        
             ! ** Scaling ot data **  
               REAL(KIND=dp), DIMENSION(number_of_func) :: factor             ! the scaling factors for f_i and g_l in the improvement function            
               REAL(KIND=dp), DIMENSION(number_of_func) :: factor_new         ! the new scaling factors for f_i and g_l in the improvement function            
               REAL(KIND=dp), DIMENSION(number_of_func) :: factors_one        ! the scaling factors for f_i and g_l inside the descent direction determination  
               
               INTEGER, DIMENSION(number_of_func) :: fac_ind                  ! the scaling factor indices for f_i and g_l              
               INTEGER, DIMENSION(number_of_func) :: ero_ind                  ! the difference between scaling factor indices for f_i and g_l              
             !-------------------------------------------               
             
             ! ** Proximity parameter related **             
               REAL(KIND=dp), DIMENSION(number_of_obj) :: t                   ! the value of the proximity parameter
               REAL(KIND=dp), DIMENSION(number_of_obj) :: t_min               ! the lower bound for the proximity parameter
               REAL(KIND=dp), DIMENSION(number_of_obj) :: t_max               ! the upper bound for the proximity parameter
               
               REAL(KIND=dp) :: norm1                                         ! the norm of the subgradient for the first DC component of f_i at the current iteration point 
               REAL(KIND=dp) :: max_norm                                      ! the maximum norm of subgradient in B_2
               INTEGER, DIMENSION(number_of_obj) :: i_t                       ! the index used in the update process of the proximity parameter
             !-------------------------------------------               
  
             ! ** Parametrs needed to ryn the new method **            
               REAL(KIND=dp) :: c              ! Decrease parameter used during the descent direction determination
               REAL(KIND=dp) :: eps            ! 'eps'=the proximity measure 
               REAL(KIND=dp) :: m              ! the descent parameter used in 'main iteration'
               REAL(KIND=dp) :: r_dec, r_inc   ! 'r_dec'=the decrease parameter and 'R_inc'=the increase parameter                 

               REAL(KIND=dp) :: crit_tol_H     ! 'crit_tol_H'=the stopping tolerance for 'Escape procedure' with the improvement function
               REAL(KIND=dp) :: crit_tol_d     ! 'crit_tol_d'=the stopping tolerance for 'Escape procedure' inside the descent direction determination             
               REAL(KIND=dp) :: step_tol_H     ! Step-length tolerance used in 'Escape procedure' with the improvement function 
               REAL(KIND=dp) :: step_tol_d     ! Step-length tolerance used in 'Escape procedure' inside the descent direction determination     
               REAL(KIND=dp) :: m_escape_H     ! the descent parameter used in 'Escape procedure' with the improvement function 
               REAL(KIND=dp) :: m_escape_d     ! the descent parameter used in 'Escape procedure' inside the descent direction determination
               
               REAL(KIND=dp) :: escape_tol_H   ! the escape tolerance launching the 'Escape procedure' for the improvement function
                          
               REAL(KIND=dp) :: decrease       ! approximation of decrease  
               
               REAL(KIND=dp) :: change         ! 'change' = f_i(x_new) - f_i(x_current)  (i.e. the change in the objective function value)
               REAL(KIND=dp) :: change1        ! 'change1' = f1_i(x_new) - f1_i(x_current)  (i.e. the change in the DC component f1_i)
               REAL(KIND=dp) :: change2        ! 'change2' = f2_i(x_new) - f2_i(x_current)  (i.e. the change in the DC component f2_i)

               INTEGER :: max_null_step        ! The maximum number of consecutive null steps            
               INTEGER :: size_b1 , size_b2    ! The biggest possible size of the bundles B_1 and B_2 
             !-------------------------------------------               
           
             ! ** Different times during the execution of the algorithm **
               REAL(KIND=dp) :: start_time, finish_time                   ! start and finish CPU time of whole algorithm
               REAL(KIND=dp) :: start_direction, finish_direction         ! start and finish CPU time to calculate a search direction
               REAL(KIND=dp) :: start_round, finish_round                 ! start and finish CPU time of one 'main iteration'                
               REAL(KIND=dp) :: elapsed_time                              ! elapsed 'clock' time in seconds
               
               INTEGER :: clock_start, clock_end, clock_rate              ! start and finish 'clock' time   
             !-------------------------------------------               
               
             ! ** Help counters ** 
               INTEGER :: help_round_counter       ! the number of iteration rounds executed during separate descent direction determination
               INTEGER :: help_subprob_counter     ! the number of subproblems solved during separate descent direction determination
               INTEGER :: help_f_counter           ! the number of function values evaluated for a DC component during separate descent direction determination (same for f_1 and f_2)
               INTEGER :: help_stop_counter        ! the number of function values evaluated for a DC component during separate descent direction determination (same for f_1 and f_2)
               INTEGER :: help_subgrad1_counter    ! the number of subgradients calculated for f_1 during separate descent direction determination
               INTEGER :: help_subgrad2_counter    ! the number of subgradients calculated for f_2 during separate descent direction determination
               
               INTEGER :: help_escape_f_counter    ! the number of function values evaluated for f_i in 'Escape procedure' algorithm
               INTEGER :: help_escape_sub_counter  ! the number of subgradients evaluated for f_i in 'Escape procedure' algorithm
                
               INTEGER, DIMENSION(number_of_obj) :: direction_stop        ! the reason for stop during separate descent direction determination
                                                                          !  0  -  a new iteration point found
                                                                          !  1  -  stopping condition satisfied (Clarke stationarity)
                                                                          !  2  -  approximate stopping condition satisfied (i.e. the step-length beta* < step_tol_d)  
                                                                          !  3  -  the biggest possible number of rounds executed during descent direction determination
                                                                          !  5  -  the biggest possible number of rounds executed in the 'Escape procedure' when descent direction determination is done
              
               INTEGER :: escape_stop           ! the reason for stop during the escape procedure for the improvement function
             !-------------------------------------------               

             ! ** Used in parallellization ** 
               INTEGER :: max_threads           ! the maximum number of threads that can be used in parallellization
               INTEGER :: threads               ! the number of threads used in parallellization
               INTEGER :: max_sub_prob          ! the maximum number of subproblems solved at the same time
             !-------------------------------------------               
               
             ! ** Help values **
               REAL(KIND=dp), DIMENSION(user_n) :: vect           ! a 'help' vector
               REAL(KIND=dp) :: linerr1, linerr2                  ! help linearization errors

               INTEGER :: i, j, k         
               INTEGER :: ind
               INTEGER :: small_ind                   
               INTEGER :: line_round       ! The round in line search                          
             !-------------------------------------------        
                               
             ! ** Needed when we deterine if a separate direction is descent for each objective function ** 
               REAL(KIND=dp), DIMENSION(number_of_obj) :: v_dec_ratio  ! The vector of decrease ratios for the separate directions
               REAL(KIND=dp) :: dec_ratio                              ! decrase ratio 
               REAL(KIND=dp) :: best_dec_ratio                         ! The best decrease ratio obtained
               REAL(KIND=dp) :: tot_ratio                              ! The total decrease ratio obtained
               
               LOGICAL :: dec_direction                                ! .TRUE. then separate direction is descent for each objective
               INTEGER :: best_ind                                     ! The index of the objective for which the separate direction gives the best decrase in each objective
             !-------------------------------------------
                 
             ! ** the number of consecutive null steps **            
               INTEGER :: counter_null                                 ! The number of consecutive null steps ('Escape procedure' for the improvement function is also a null step)
  
             ! ** NOT really used **
               REAL(KIND=dp) :: delta          ! variable which can be used to relax the descent condition when separete search directions are calculated      

               REAL(KIND=dp), DIMENSION(number_of_obj,number_of_obj) :: f_for_sep_d
               LOGICAL, DIMENSION(number_of_obj) :: dec_d_list
               LOGICAL, DIMENSION(number_of_obj) :: feasible_d_list
               
               INTEGER :: command     ! Command defined by the user when escape procedure for the improvement function can be executed
               INTEGER :: i_user      ! The maximum number of times the user can steer the algorithm
               
               INTEGER :: n_const      ! The maximum number of times the user can steer the algorithm
               
               LOGICAL :: execute_escape
               LOGICAL :: correct_value
               
               
             ! ** The results file **
               CHARACTER*30 outfi/'results.txt'/
               OPEN(40,file=outfi)             
              
             ! ** END of variables definition **
             
               
               CALL cpu_time(start_time)                ! Start CPU timing     

               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time       
               
               
              ! *******************************************************************************
              !
              !                        INITIALIZATIONS
              ! 
              !********************************************************************************
                
                
              !------------------------------------------------       
              !             **  STARTING POINT **
              !------------------------------------------------
              CALL starting_point(x_0, startpoint, user_n)
       
           
               !*********************************************************************************
               !
               !     - INITIAL INTEGER VALUES OF PARAMETERs FROM USER ARE LOOKED THROUGH 
               !       AND FINAL VALUES ARE SET BASED ON THEM.
               !  
               !________________________________________________________________________________
               !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
               !***** PARAMETRES NEEDED ONLY IN SEPARATE DESCENT DIRECTION DETERMINATIONS ******
               !
               ! ***** descent parameter 'm'='m_1' *****
               SELECT CASE(user_m)
                   
                   CASE(0)
                       m = 0.2_dp
                   CASE(1)
                       m = 0.5_dp
                   CASE(2)
                       m = 0.1_dp
                   CASE(3)
                       m = 0.05_dp
                   CASE(4)  
                       m = 0.01_dp
                   CASE(5)             
                       
               END SELECT

   
               ! ***** decrease parameter 'r_2' *****               
               SELECT CASE(user_c)
                   
                   CASE(0)
                       c = 0.1_dp
                   CASE(1)
                       c = 0.05_dp
                   CASE(2)
                      
                   CASE(3)
               
               
               END SELECT
              
              
               ! ***** decrease parameter 'r_dec'='r_1' *****
               SELECT CASE(user_r_dec)
                   
                   CASE(0)
                       IF( user_n < 10) THEN 
                          r_dec = 0.75_dp
                       ELSE IF ( user_n >= 10 .AND. user_n <300 ) THEN 
                          r_dec = user_n/(user_n+5.0_dp)
                       ELSE IF ( user_n >= 300 ) THEN 
                          r_dec = 0.99_dp                    
                       END IF                      
                   CASE(1)
                          r_dec = 0.99_dp                    
                   CASE(2)
                          r_dec = 0.5_dp                    

                   CASE(3)
               
               
               END SELECT 
               
               ! ***** increase parameter 'r_inc'='R' *****
               SELECT CASE(user_r_inc)
                   
                   CASE(0)
                       IF (user_n <= 10) THEN 
                          r_inc = (10.0_dp)**6
                       ELSE IF ((user_n > 10) .AND. (user_n <= 50)) THEN 
                          r_inc = (10.0_dp)**7
                       ELSE 
                          r_inc = (10.0_dp)**7
                       END IF                  
                   CASE(1)
                       r_inc = (10.0_dp)**6
                   CASE(2)
                       r_inc = (10.0_dp)**7                
                   CASE(3)
                       r_inc = (10.0_dp)**8
                   CASE(4)
                      IF (user_n < 50) THEN 
                       r_inc = (10.0_dp)**7               
                      ELSE IF (user_n >= 50) THEN 
                       r_inc = (10.0_dp)**8                
                      END IF 
                   CASE(5)    
                      IF (user_n <= 100) THEN 
                       r_inc = (10.0_dp)**8               
                      ELSE IF (user_n > 100) THEN 
                       r_inc = (10.0_dp)**7                
                      END IF 
                   CASE(6)    
                      IF (user_n < 10) THEN 
                          r_inc = (10.0_dp)**7    
                      ELSE IF (user_n >= 10 .AND. user_n < 101) THEN 
                          r_inc = (10.0_dp)**8
                      ELSE IF (user_n > 100) THEN 
                          r_inc = (10.0_dp)**7                
                      END IF 
                   CASE(7)    
                      IF (user_n < 10) THEN 
                          r_inc = (10.0_dp)**7    
                      ELSE IF (user_n >= 10 .AND. user_n < 100) THEN 
                          r_inc = (10.0_dp)**8                    
                      ELSE IF (user_n == 100) THEN 
                          r_inc = (10.0_dp)**8                    
                      ELSE IF (user_n > 100) THEN 
                          r_inc = (10.0_dp)**6                
                      END IF  
                   CASE(8)
                      r_inc = (10.0_dp)**4   
                   CASE(9) 
                      IF (user_n < 10) THEN 
                          r_inc = (10.0_dp)**7    
                      ELSE IF (user_n >= 10 .AND. user_n < 100) THEN 
                          r_inc = (10.0_dp)**5                    
                      ELSE IF (user_n == 100) THEN 
                          r_inc = (10.0_dp)**8                    
                      ELSE IF (user_n > 100) THEN 
                          r_inc = (10.0_dp)**6                
                      END IF   
               END SELECT

               ! ***** bundle B1 size *****
               SELECT CASE(user_size_b1)
                   
                   CASE(0)
                       size_b1 = MIN((user_n+5),50)
                   CASE(1)
                       size_b1 = MIN((user_n+5),25)                   
                   CASE(2)
                       size_b1 = MIN((user_n+5),100)                   
                   CASE(3)
                       size_b1 = MIN((user_n+5),150)      
                   CASE(4)
                       size_b1 = MIN((user_n+5)*2,50)                   

               
               END SELECT


               ! ***** bundle B2 size *****
               SELECT CASE(user_size_b2)
                   
                   CASE(0)
                       size_b2 = 3
                   CASE(1)
                       size_b2 = 5
                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT


               ! ***** maximum number of consecutive null steps 'max_null_step' *****
               SELECT CASE(max_null_step0)
                   
                   CASE(0)
                      max_null_step = 5
                      
                   CASE(1)
                      max_null_step = 10
                      
                   CASE(2)
                      max_null_step = 2
               
               END SELECT 
               
               !________________________________________________________________________________
               !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
               !****************** PARAMETRES NEEDED IN STOPPING CONDITIONS ********************  
               !  
               ! ***** stopping tolerance 'crit_tol_H'='delta_2' for the escape procedure with the improvement function *****
               SELECT CASE(user_crit_tol_H)
                   
                   CASE(0)
                       crit_tol_H = (10.0_dp)**(-5)
                       
                   CASE(1)
                      IF (user_n < 10) THEN
                       crit_tol_H = (10.0_dp)**(-5)
                      ELSE 
                       crit_tol_H = (10.0_dp)**(-3)
                      END IF                      

                   CASE(2)
                      IF (user_n < 10) THEN
                       crit_tol_H = (10.0_dp)**(-5)
                      ELSE 
                       crit_tol_H = (10.0_dp)**(-2)
                      END IF    
                      
                   CASE(3)
               
               
               END SELECT

               ! ***** stopping tolerance 'crit_tol_d'='delta_1' for the escape procedure inside the descent direction determination *****
               SELECT CASE(user_crit_tol_d)
                   
                   CASE(0)
                       crit_tol_d = (10.0_dp)**(-5)
                   CASE(1)
                      IF (user_n < 10) THEN
                       crit_tol_d = (10.0_dp)**(-5)
                      ELSE 
                       crit_tol_d = (10.0_dp)**(-3)
                      END IF                      

                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT

              
               ! ***** enlargement parameter 'eps'='theta' *****
               SELECT CASE(user_eps)
                   
                   CASE(0)
                       eps = 5*(10.0_dp)**(-5)
                   CASE(1)
                   
                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT

               
               ! ***** step-size tolerance 'tol_tau'='tau' *****
               SELECT CASE(user_tol_tau)
                   
                   CASE(0)
                       IF (user_n < 101) THEN             
                          tol_tau = (10.0_dp)**(-2)
                       ELSE 
                          tol_tau = (10.0_dp)**(-3)
                       END IF                      
                   CASE(1)
                   
                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT
               
               !________________________________________________________________________________
               !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
               !*************** PARAMETRES NEEDED ONLY IN ESCAPE PROCEDURE *********************
              
               ! ***** Descent parameter 'm_escape_H'='m_2' for the escape procedure with the improvement function *****
               SELECT CASE(user_m_escape_H)
                   
                   CASE(0)
                        m_escape_H = 0.01_dp
                   CASE(1)
                   
                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT
           
               ! ***** Descent parameter 'm_escape_d'='m_2' for the escape procedure inside the descent direction determination *****
               SELECT CASE(user_m_escape_d)
                   
                   CASE(0)
                        m_escape_d = 0.01_dp
                   CASE(1)
                   
                   CASE(2)
                   
                   CASE(3)
               
               
               END SELECT 

 
               ! ***** Tolerance launching the 'Escape procedure'='gamma' with the improvement function *****
               SELECT CASE(user_escape_tol_H)
                   
                   CASE(0)
                       IF (user_n <= 50) THEN
                          escape_tol_H = (10.0_dp)**(-4)
                       ELSE IF ((user_n >50) .AND. (user_n  <= 100) ) THEN 
                          escape_tol_H = (10.0_dp)**(-4)               
                       ELSE   
                          escape_tol_H = (10.0_dp)**(-3)                             
                       END IF                      
                   CASE(1)
                       IF (user_n <= 10) THEN
                          escape_tol_H = (10.0_dp)**(-5)          
                       
                       ELSE IF ((user_n >10) .AND. (user_n  <= 100) ) THEN
                          escape_tol_H = (10.0_dp)**(-4)          
                       ELSE   
                          escape_tol_H = (10.0_dp)**(-3)                             
                       END IF                         
                   CASE(2)
                       IF (user_n <= 10) THEN
                          escape_tol_H = (10.0_dp)**(-5)          
                       
                       ELSE IF ((user_n >10) .AND. (user_n  <= 100) ) THEN
                          escape_tol_H = (10.0_dp)**(-4)          
                       ELSE   
                          escape_tol_H = (10.0_dp)**(-4)                             
                       END IF                        
                   CASE(3)
                          escape_tol_H = (10.0_dp)**(-3)                             
             
                   CASE(4)
                       IF (user_n <= 10) THEN
                          escape_tol_H = (10.0_dp)**(-6)                   
                       ELSE   
                          escape_tol_H = (10.0_dp)**(-5)                             
                       END IF                      
                   
                   CASE(5)
                       IF (user_n <= 10) THEN
                          escape_tol_H = (10.0_dp)**(-7)                   
                       ELSE   
                          escape_tol_H = (10.0_dp)**(-6)                             
                       END IF   
                   
                   
               END SELECT

            
               ! ***** Step-length tolerance 'step_tol_H'='eps2' for the escape procedure with the improvement function *****
               SELECT CASE(user_step_tol_H)
                   
                   CASE(0)
                       IF (user_n <= 50) THEN
                          step_tol_H = (10.0_dp)**(-6)
                       ELSE
                          step_tol_H = (10.0_dp)**(-5)               
                       END IF                      
                   
                   CASE(1)
                       IF (user_n <= 100) THEN
                          step_tol_H = (10.0_dp)**(-4)
                       ELSE
                          step_tol_H = (10.0_dp)**(-3)               
                       END IF                       
                   
                   CASE(2)
                       step_tol_H = (10.0_dp)**(-4)
                       
                   CASE(3)
                        step_tol_H = (10.0_dp)**(-3)
              
                   CASE(4)
                       IF (user_n <= 100) THEN
                          step_tol_H = (10.0_dp)**(-5)
                       ELSE
                          step_tol_H = (10.0_dp)**(-4)               
                       END IF     
                       
                   CASE(5)
                      IF (user_n < 10) THEN 
                          step_tol_H = (10.0_dp)**(-5) 
                      ELSE IF (user_n >= 10 .AND. user_n < 50) THEN 
                          step_tol_H = (10.0_dp)**(-4)           
                      ELSE 
                          step_tol_H = (10.0_dp)**(-3)
                      END IF     

                   CASE(6)
                      IF (user_n < 10) THEN 
                          step_tol_H = (10.0_dp)**(-4) 
                      ELSE IF (user_n >= 10 .AND. user_n < 50) THEN 
                          step_tol_H = (10.0_dp)**(-3)           
                      ELSE 
                          step_tol_H = (10.0_dp)**(-2)
                      END IF                      
                       
               END SELECT

 
               ! ***** Step-length tolerance 'step_tol_d'='eps1' for the escape procedure inside the descent direction determination *****
               SELECT CASE(user_step_tol_d)
                   
                   CASE(0)
                       step_tol_d = (10.0_dp)**(-6)    

               END SELECT


               !________________________________________________________________________________
               !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
               !*************** CONTROL COUNTERS TO STOP EXECUTION OF METHOD *******************               
               !_______________________________________________________________________________

               ! ***** maximum number of iterations 'mit' *****
               IF ( mit <= 0 ) THEN
                   mit = 1000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds'  during the separate descent direction determination *****
               IF ( mrounds <= 0 ) THEN
                   mrounds = 5000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds' during 'Escape procedure' *****
               IF ( mrounds_escape <= 0 ) THEN
                   mrounds_escape = 5000
               END IF                  
               
               ! ***** print option 'iprint' *****
               IF ( (iprint < -4 ) .OR. (iprint > 5) ) THEN   !executed if print value is wrong
                   iprint = 1
               END IF              
               !_______________________________________________________________________________
              
               ! ***** Steer option 'i_user' *****
               IF ( (i_user0 < 0) ) THEN   !executed if i_user value is wrong
                   i_user = 0
               ELSE
                   i_user = i_user0            
               END IF                 
              
              
           !_______________________________________________________________________________
           !************************ STEP 0: PARAMETER INITIALIZATION *********************    
               
               ! Are contraints used in the individual search direction problem?
               IF (use_const_in_sep_d) THEN 
                  n_const = number_of_const   ! YES
               ELSE
                  n_const = 0                 ! NO                     
               END IF 
               
               
              ! the initial scaling factors
               factor = 1.0_dp                    ! For improvement function        
               factors_one = 1.0_dp               ! For separate direction determination            
 
             ! The values of DC components f1_i and f2_i at x_0 (NO scaling)
               DO i = 1, number_of_obj                                      
                  fg1_k_all(i) = f1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = f2(x_0,problem2(i),factor(i),user_n) 
               END DO  
             
             ! The values of DC components g1_l and g2_l at x_0  (NO scaling)
               DO i = number_of_obj+1, number_of_func                       
                  fg1_k_all(i) = g1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = g2(x_0,problem2(i),factor(i),user_n) 
               END DO              
               fg_k_all = fg1_k_all - fg2_k_all
               
               IF (iprint > 1) THEN    
                   IF (number_of_const==0) THEN
                              WRITE(*,*) 'Round:', 0, 'f(x):', fg_k_all 
                   ELSE
                              WRITE(*,*) 'Round:', 0, 'f(x):', fg_k_all(1:number_of_obj), &
                                       & 'g(x):', fg_k_all(number_of_obj+1:number_of_func)
                   END IF
               ELSE IF (iprint==-1) THEN       
                              WRITE(*,*) 'Round:', 0, 'x:', x_0                
               END IF
               
               !----------------------------------------------------
               !            ** SCALING OF FUNCTIONS **
               !----------------------------------------------------
               ! Only used when dimension is > 100
               IF (user_n > 100) THEN 
                   small_ind = 1
                   DO i = 2, number_of_func                      
                      IF (ABS(fg_k_all(small_ind))> ABS(fg_k_all(i))) THEN 
                        small_ind = i
                      END IF
                   END DO       
                   
                   IF (scale_func) THEN
                    DO i = 1, number_of_func
                      IF ( ABS(fg_k_all(i)) < (10.0_dp)**(0)) THEN
                        fac_ind(i) = 0                      
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(1)) THEN
                        fac_ind(i) = 1  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(2)) THEN
                        fac_ind(i) = 2  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(3)) THEN
                        fac_ind(i) = 3  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(4)) THEN
                        fac_ind(i) = 4  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(5)) THEN
                        fac_ind(i) = 5  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(6)) THEN
                        fac_ind(i) = 6  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(7)) THEN
                        fac_ind(i) = 7  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(8)) THEN
                        fac_ind(i) = 8  
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(9)) THEN
                        fac_ind(i) = 9
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(10)) THEN
                        fac_ind(i) = 10 
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(11)) THEN
                        fac_ind(i) = 11 
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(12)) THEN
                        fac_ind(i) = 12 
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(13)) THEN
                        fac_ind(i) = 13
                      ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(14)) THEN
                        fac_ind(i) = 14                     
                      END IF    
                    END DO

                    DO i = 1, number_of_func
                      ero_ind(i) = fac_ind(i)-fac_ind(small_ind)
                       IF(ero_ind(i) > 1) THEN
                          ero_ind(i) = ero_ind(i)-1
                       END IF
                       factor(i) = (10.0_dp)**(-ero_ind(i))                   
                    END DO
                   
                   END IF
               END IF 
               !----------------------------------------------------
               !                ** END OF SCALING **
               !----------------------------------------------------
              
             ! Scaling is also used in the separate direction determination when dimension > 100
               factors_one = factor
               
               IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN
                   IF (user_n < 100) THEN
                       factors_one(number_of_func) = scaling_const * (10.0_dp)**(0)
                   ELSE
                       IF (user_n == 100) THEN 
                           factors_one(number_of_func) = scaling_const * (10.0_dp)**(1)
                       END IF                      
                   END IF                    
               END IF     
                              
             ! The current iteration point is the starting point x_0 
               x_current = x_0                         

             ! Initialization of bundles B_1 and B_2 for each objective
               DO i = 1, number_of_obj                       
                  CALL init_bundle_b1(B1(i), size_b1, user_n)    
                  CALL init_bundle_b2(B2(i), size_b2, user_n)            
               END DO 

             ! The values of DC components f1_i and f2_i at x_0
               DO i = 1, number_of_obj                       
                  fg1_k_all(i) = f1(x_0,problem1(i),factors_one(i),user_n) 
                  fg2_k_all(i) = f2(x_0,problem2(i),factors_one(i),user_n)
                
                ! The values of the objectives f_i at x_0    
                  f_0_all(i) = fg1_k_all(i) - fg2_k_all(i) 
                ! The value of the objectives f_i at x_current  
                  fg_k_all(i) = f_0_all(i)                
               END DO         
           
               IF (number_of_const>0) THEN
                 ! The values of DC components g1_l and g2_l at x_0                 
                   DO i = number_of_obj+1, number_of_func       
                      fg1_k_all(i) = g1(x_0,problem1(i),factors_one(i),user_n) 
                      fg2_k_all(i) = g2(x_0,problem2(i),factors_one(i),user_n)
                      fg_k_all(i) = fg1_k_all(i)-fg2_k_all(i)
                   END DO    
               END IF
               
             ! One function value evaluated for each DC component f1_i               
               f_counter = 1               
               
               IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN  
                   g_counter = 1                ! One more function value evaluated for the constraints                       
               ELSE
                   g_counter = 0
               END IF 
               
             ! The values of subgradients for DC components f1_i and f2_i at x_0   
               DO i = 1, number_of_obj                       
                   grad1 = subgradient_f1(x_0, problem1(i),factors_one(i),user_n)
                   grad2 = subgradient_f2(x_0, problem2(i),factors_one(i),user_n)
                   
                   IF (n_const==0) THEN 
                      CALL add_first_element_b1(B1(i), grad1, i)
                      CALL add_first_element_b2(B2(i), grad2)   
                   ELSE               
                       DO j = 1, user_n
                          grad1_k_all(i,j) = grad1(j)
                          grad2_k_all(i,j) = grad2(j)
                       END DO    
                   END IF                  
               END DO  
               
               IF (n_const>0) THEN
                ! The values of subgradients for DC components g1_l and g2_l at x_0
                   DO i = number_of_obj+1, number_of_func
                       grad1 = subgradient_g1(x_0, problem1(i),factors_one(i),user_n)
                       grad2 = subgradient_g2(x_0, problem2(i),factors_one(i),user_n)
                       DO j = 1, user_n
                          grad1_k_all(i,j) = grad1(j)
                          grad2_k_all(i,j) = grad2(j)
                       END DO
                   END DO
                END IF   

               subgrad1_counter = 1        ! one subgradient calculated for each DC component f1_i and g1_l
               subgrad2_counter = 1        ! one subgradient calculated for each DC component f2_i and g2_l 
               
               IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN 
                    subgrad1_g_counter = 1        ! one subgradient calculated for each DC component g1_l
                    subgrad2_g_counter = 1        ! one subgradient calculated for each DC component g2_l
               ELSE
                    subgrad1_g_counter = 0        ! NO subgradient calculated for each DC component g1_l
                    subgrad2_g_counter = 0        ! NO subgradient calculated for each DC component g2_l               
               END IF 
               
              IF (n_const>0) THEN
                  ! The subdradients are calculated in the individual search direction problems
                   DO i = 1, number_of_obj
                       
                     ! The value of objective in the individual search direction problem
                       fg1_k_all_sep(1) = fg1_k_all(i) 
                       fg2_k_all_sep(1) = fg2_k_all(i) 
                     ! The values of constraints in the individual search direction problem
                       IF (n_const>0) THEN  
                           DO j = 1, n_const
                               fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj) 
                               fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)             
                           END DO
                       END IF
                                          
                    ! The subgradient of the individual objctive in the search direction problem  
                       DO j = 1, user_n
                          grad1_k_all_sep(1,j) = grad1_k_all(i,j) 
                          grad2_k_all_sep(1,j) = grad2_k_all(i,j) 
                       END DO
                           
                       IF (n_const>0) THEN  
                         ! The subradients of constraints
                           DO j = 1, n_const  
                             DO k = 1, user_n
                                grad1_k_all_sep(j+1,k) = grad1_k_all(j+number_of_obj,k) 
                                grad2_k_all_sep(j+1,k) = grad2_k_all(j+number_of_obj,k)                    
                             END DO                
                           END DO                  
                       END IF

                     ! the values of components A_i and B_l for H_1 at x_0 in the individual search direction problem
                       DO j = 1, 1+n_const                 
                           AandB_k_all_sep(j) = AorB_i_sep(fg1_k_all_sep, fg2_k_all_sep, fg1_k_all_sep, fg2_k_all_sep ,j)
                       END DO
                       
                       ! A and B at the current point for the individual search direction problem are stored
                       AB_obj_k(i) = AandB_k_all_sep(1)
                       DO j = 1, n_const
                            AB_const_k(i,j) = AandB_k_all_sep(j+1)              
                       END DO   
                       
                       H1_current_sep(i) = H1_sep(AandB_k_all_sep)              ! the value of the DC component H_1 at x_0 in the individual search direction problem
                       H2_current_sep(i) = H2_sep(fg2_k_all_sep)                ! the value of the DC component H_2 at x_0 in the individual search direction problem
          
                       ind = index_H1_sep(AandB_k_all_sep)                                           ! the index of the component A_i or B_l yielding the subgradient of H_1 in the individual search direction problem
                       grad1 = subgradient_AorB_i_sep(grad1_k_all_sep, grad2_k_all_sep, ind, user_n) ! the subgradient of H_1 at x_0 in the individual search direction problem
                       grad2 = subgradient_H2_sep(grad2_k_all_sep, user_n)                           ! the subgradient of H_2 at x_0 in the individual search direction problem
                       
                       ! The first bundle element is added into the bundle B_1 and B_2 
                       ! (i.e. the one corresponding to the starting point x_0)
                       CALL add_first_element_b1(B1(i), grad1, ind)
                       CALL add_first_element_b2(B2(i), grad2)      
                       
                       DO j = 1, n_const+1
                          IF ( j/= ind) THEN 
                             grad1 = subgradient_AorB_i_sep(grad1_k_all_sep, grad2_k_all_sep, j, user_n)        ! the subgradients of A_i and B_l at x_0 in the individual search direction problem
                             alpha = 0
                             CALL add_element_b1(B1(i), grad1, alpha, j)
                          END IF
                       END DO                  
                       
                   END DO
                   
               END IF
               
               subprob_counter = 0         ! the number of 'subproblems' solved so far is zero for each objective
               
               escape_counter_d = 0        ! the number times escape procedure is executed for separate objectives
               escape_f_counter_d = 0      ! the number function values during escape procedures for separate objectives
               escape_sub_counter_d = 0    ! the number subgradients during escape procedures for separate objectives
               escape_g_counter_d = 0      ! the number function values during escape procedures for separate objectives
               escape_sub_g_counter_d = 0  ! the number subgradients during escape procedures for separate objectives
               
               f_counter_line = 0          ! The number of objective function calculations during line searches
               g_counter_line = 0          ! The number of objective function calculations during line searches
               
               iter_counter = 0            ! the number of main iterations executed so far is zero
               counter_null = 0            ! the number of consecutive null steps (Escape procedure for the improvement function is a null step also)
                  
               stop_alg = .FALSE.          ! we cannot stop the bundle algorithm

               escape_counter_H = 0        ! the number of times 'Escape procedure' is used for the improvement function                             
               escape_f_counter_H = 0      ! the initialization of objective function counter for 'Escape procedure' with the improvement function
               escape_sub_counter_H = 0    ! the initialization of subdradient counter for 'Escape procedure' with the improvement function 
               escape_g_counter_H = 0      ! the initialization of objective function counter for 'Escape procedure' with the improvement function
               escape_sub_g_counter_H = 0  ! the initialization of subdradient counter for 'Escape procedure' with the improvement function 
               
               delta = 0.0_dp
               
               !----------------------------------------------------
               !      ** Proximity parameter initialization **
               !----------------------------------------------------
               ! NOT needed at the moment
                IF (t_adjust) THEN
                  DO i = 1, number_of_obj
                    vect = give_subgrad_b1(B1(i),0)                          ! the vector \bxi_1(x_k)
                    norm1 = SQRT(DOT_PRODUCT(vect,vect))                     ! the value of the norm ||\bxi_1(x_k)||
                    max_norm = max_norm_value(B2(i))                         ! the value of the maximum subgradient norm in the bundle B_2
                    t_min(i) = (0.5_dp * r_dec * eps) / ( norm1 + max_norm ) ! the lower bound for t
                    t_max(i) = r_inc * t_min(i)                              ! the upper bound for t
                    i_t(i) = 0                                               ! index initialization
                    
                    t(i) = 0.8_dp * (t_min(i) + t_max(i))
                
                  END DO
                END IF
                
               !----------------------------------------------------
               !    ** Proximity parameter initialization END **
               !----------------------------------------------------                
                
                
! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
!               max_threads = omp_get_max_threads()
!               max_sub_prob = give_max_size_b2(B2)+1
!               threads = MIN(max_threads, max_sub_prob)
!               CALL omp_set_num_threads(threads)  
! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
               
               
               IF (iprint > 1) THEN    
                   IF (number_of_const==0) THEN
                              WRITE(*,*) 'Round:', 0, 'f(x):', fg_k_all 
                   ELSE
                              WRITE(*,*) 'Round:', 0, 'f(x):', fg_k_all(1:number_of_obj), &
                                       & 'g(x):', fg_k_all(number_of_obj+1:number_of_func)
                   END IF
               END IF              
           !_______________________________________________________________________________
           !************************ STEP 0: END ******************************************                   

           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------  
           
           
           
            DO WHILE ( ( .NOT. stop_alg) .AND. (iter_counter < mit )  )  ! is repeated until the bundle algorithm DBDC can be terminated 
                                                                         
               !_______________________________________________________________________________
               !************************ STEP 1: SEPARATE DIRECTIONS **************************
               
              ! a new iteration is started
                iter_counter = iter_counter + 1                                  
             
              ! start time for the current iteration 
                CALL cpu_time(start_round)                         
               
              ! decrease ratios for separate directions are initialized               
                v_dec_ratio = 0.0_dp            
              
              !--------------------------------------------------------------- 
              !     ** SEPARATE DESCENT DIRECTIONS FOR EACH OBJECTIVE  **
              !---------------------------------------------------------------
                DO i = 1, number_of_obj
                
                    CALL cpu_time(start_direction)
                    
                    IF (n_const==0) THEN 
                        H1_current_sep(i) = fg1_k_all(i)
                        H2_current_sep(i) = fg2_k_all(i)
                    ELSE
                        AandB_k_all_sep(1) = AB_obj_k(i)  
                        fg1_k_all_sep(1) = fg1_k_all(i)                     
                        fg2_k_all_sep(1) = fg2_k_all(i)                     
                        DO j = 1, n_const
                            AandB_k_all_sep(j+1) = AB_const_k(i,j)
                            fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj)                       
                            fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)                                        
                        END DO                      
                        
                    END IF                   
                 
                  ! Descent direction 
                    CALL separate_direction(i, factors_one(i), x_current, H1_current_sep(i), H2_current_sep(i), &   
                            &  delta, user_n, x_new, H1_new_sep(i), H2_new_sep(i), f_0_all(i), B1(i), B2(i),  &
                            & crit_tol_d, eps, m, c, r_dec, r_inc, m_escape_d, step_tol_d, mrounds, mrounds_escape, & 
                            & direction_stop(i), help_round_counter, help_stop_counter, &
                            & help_subprob_counter, help_f_counter, help_subgrad1_counter, &
                            & help_subgrad2_counter, agg_used, t(i), t_min(i), t_max(i),  &
                            & t_adjust, i_t(i), decrease, help_escape_f_counter, help_escape_sub_counter, &
                            & n_const, AandB_k_all_sep, fg1_k_all_sep, fg2_k_all_sep, factors_one(number_of_func))  


                    IF (n_const==0) THEN 
                        f1_new(i) = H1_new_sep(i)
                        f2_new(i) = H2_new_sep(i)
                    ELSE
                        f1_new(i) = f1(x_new,problem1(i),factors_one(i),user_n)
                        f2_new(i) = f2(x_new,problem1(i),factors_one(i),user_n)                
                    END IF

                    CALL cpu_time(finish_direction)
                
                    ! Different counters are updated
                    subprob_counter(i) = subprob_counter(i) + help_subprob_counter
                    f_counter(i) = f_counter(i) + help_f_counter
                    subgrad1_counter(i) = subgrad1_counter(i) + help_subgrad1_counter
                    subgrad2_counter(i) = subgrad2_counter(i) + help_subgrad2_counter 
                    
                    escape_counter_d(i) = escape_counter_d(i) + help_stop_counter 
                    escape_f_counter_d(i) = escape_f_counter_d(i) + help_escape_f_counter 
                    escape_sub_counter_d(i) = escape_sub_counter_d(i) + help_escape_sub_counter 
 
                    ! Only done if constraints are used in the individual search direction problem
                    IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN
                        g_counter = g_counter + help_f_counter
                        subgrad1_g_counter = subgrad1_g_counter + help_subgrad1_counter
                        subgrad2_g_counter = subgrad2_g_counter + help_subgrad1_counter 
                        escape_sub_g_counter_d = escape_sub_g_counter_d + help_escape_sub_counter 
                        escape_g_counter_d = escape_g_counter_d + help_escape_f_counter                         
                    END IF                 
                   
                  ! We have a descent direction
                    IF (direction_stop(i) == 0) THEN
                        separate_d = (x_new-x_current)   ! The separate direction
                        
                        dec_direction = .TRUE. 
                        dec_d_list(i) = .TRUE.
                        
                      ! The decrease ratio for the separate direction is determined
                        DO j = 1, number_of_obj
                            
                              y = x_current + separate_d

                             IF (n_const==0) THEN 
                                 IF (j==i) THEN 
                                    f_y = f1_new(i)-f2_new(i)
                                    f_for_sep_d(i,j) = f_y/factors_one(j)
                                 ELSE                              
                                    f1_y = f1(y,problem1(j),factors_one(j),user_n)
                                    f2_y = f2(y,problem2(j),factors_one(j),user_n)
                                    f_y = f1_y - f2_y
                                    f_counter(i) = f_counter(i) + 1
                                   
                                    f_for_sep_d(i,j) = f_y/factors_one(j)
                                 END IF 

                              ELSE
                                    f1_y = f1(y,problem1(j),factors_one(j),user_n)
                                    f2_y = f2(y,problem2(j),factors_one(j),user_n)
                                    f_y = f1_y - f2_y
                                    f_counter(i) = f_counter(i) + 1
                                   
                                    f_for_sep_d(i,j) = f_y/factors_one(j)
                              END IF                              

                             dec_ratio = fg1_k_all(j)-fg2_k_all(j)-f_y
                             
                             dec_ratio = dec_ratio / (ABS(fg1_k_all(j)-fg2_k_all(j))+1.0_dp)
                         
                             v_dec_ratio(i) = v_dec_ratio(i) + dec_ratio  
                             
                             IF (dec_ratio < 0.0_dp) THEN   ! The separate direction is not descent for the objective f_j
                                dec_direction = .FALSE.
                             END IF
                             
                        END DO
                        
                        IF (number_of_const>0) THEN 
                        
                           feasible_d_list(i) = .TRUE.
                           DO j = number_of_obj+1, number_of_func
                                g1_y = g1(y,problem1(j),factors_one(j),user_n)
                                g2_y = g2(y,problem2(j),factors_one(j),user_n)
                                IF (g1_y - g2_y>0.0_dp) THEN 
                                  feasible_d_list(i) = .FALSE.
                                END IF 
                               
                           END DO
                        END IF
                    ELSE
                        separate_d = 0.0_dp
                        dec_direction = .FALSE.
                        dec_d_list(i) = .FALSE.
                        feasible_d_list(i) = .FALSE.
                    END IF  
                    
                    
                  ! If the separate direction is not descent for each objective then decrease ratio is negative
                    IF (.NOT. dec_direction) THEN 
                         v_dec_ratio(i) = -1.0_dp
                    END IF 
                    
                  ! The obtained direction is stored to the correct place in the vector 'mD'                   
                    ind = (user_n)*(i-1)
                    DO j = 1, user_n
                        mD(ind+j) = separate_d(j)
                    END DO
                END DO
               
              ! Calculation of the best decrease ratio of separate directions and the index of the corresponding direction
                best_dec_ratio = -10.0_dp
                best_ind = 0
                DO i = 1, number_of_obj
                    IF (best_dec_ratio < v_dec_ratio(i)) THEN
                        best_dec_ratio = v_dec_ratio(i)
                        best_ind = i
                    END IF
                END DO

              !--------------------------------------------------------------- 
              !   ** END: SEPARATE DESCENT DIRECTIONS FOR EACH OBJECTIVE  **
              !---------------------------------------------------------------               

                
              !------------------------------------------------------- 
              !      ** A COMMON DIRECTION FOR ALL OBJECTIVES  **
              !-------------------------------------------------------       
                IF (best_dec_ratio > 0.0_dp) THEN 
                
                 ! The common descent direction is the separate direction with the largest decrease ratio                
                    ind = (user_n)*(best_ind-1)
                    DO j = 1, user_n
                        common_d(j) = mD(ind+j)
                    END DO                  
                    
                    common_obj = (10.0_dp)**6
                    
                    IF (iprint >= 3) THEN    
                        WRITE(*,*) 'Search direction is the separate direction for the objective', best_ind  
                    END IF 
                    
                ELSE 

                  ! The common direction is solved from the minimum norm problem 
                    CALL common_direction(x_current, user_n, number_of_obj, mD, 1.0_dp, common_d, common_obj)  
              
                    common_d = - common_d

                END IF 



              !------------------------------------------------------- 
              !     ** DO WE HAVE COMMON DESCENT DIRECTION ?  **
              !-------------------------------------------------------           
                  
             IF (common_obj > escape_tol_H) THEN  
                   
                 tau = 1.0_dp          ! The initial step-size
                 step_stop = .FALSE.   ! .TRUE. if the step-size determination can be stopped
                 inc_in_tau = .FALSE.  ! .TRUE. if there is incerase in the value of tau
                 f_best = fg_k_all     ! The best values of the objectives so far
                 line_round = 1        ! The number of round in step-size determination
                 
                 IF (best_dec_ratio > 0.0_dp) THEN 
                     f_counter_line = f_counter_line - 1
                     
                    ! Only done if constraints are used in the individual search direction problem
                    IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN
                          g_counter_line = g_counter_line - 1
                    END IF
                 END IF 
                 
                
                 DO WHILE (.NOT. step_stop)
                  
                     y = x_current + tau * common_d
                     improvement = .TRUE.
                     
                     DO i = 1, number_of_obj
                     
                        IF (improvement) THEN 
                            f1_y = f1(y,problem1(i),factors_one(i),user_n)
                            f2_y = f2(y,problem2(i),factors_one(i),user_n)
                            vf_y(i) = f1_y - f2_y
                            f_counter_line(i) = f_counter_line(i) + 1 
                            
                            IF (iprint >= 3 ) THEN 
                                WRITE(*,*) 'tau', tau, 'obj', problem1(i), 'f_y', vf_y(i), 'f_k', f_best(i)
                            END IF 
                            
                            IF (vf_y(i) > f_best(i) - (10.0_dp)**(-10) ) THEN ! There is at least one objective for which the step-size does not decrease the value of the objective 
                                improvement = .FALSE.   
                            END IF 
                        END IF      
                        
                        IF ((.NOT. improvement) .AND. (best_dec_ratio > 0.0_dp) .AND. (line_round==1)) THEN 
                             WRITE(*,*) 'HELP, something is wrong'
                        END IF  
                     END DO  
                     
                     IF (number_of_const>0) THEN
                         DO i = number_of_obj+1, number_of_func

                            IF (improvement) THEN                          
                                g1_y = g1(y,problem1(i),factors_one(i),user_n)
                                g2_y = g2(y,problem2(i),factors_one(i),user_n)
                                g_counter_line(i-number_of_obj) = g_counter_line(i-number_of_obj) + 1                         
                                
                                IF (iprint >= 3 ) THEN 
                                    WRITE(*,*) 'tau', tau, 'const', problem1(i), 'g_y', g1_y-g2_y, '<=', 0 ,'?'
                                END IF 
                                
                                IF (g1_y-g2_y > 0.0_dp ) THEN !We break one of the constraints
                                    improvement = .FALSE.   
                                END IF 
                            END IF 
                         END DO  
                     END IF
                     
                     IF (improvement) THEN  
                         tau_smallest = tau         
                         IF (tol_tau < tau .AND. tau < 100.0_dp) THEN 
                            IF (tau< 8.0_dp) THEN 
                              tau = 2*tau
                              inc_in_tau = .TRUE.
                            ELSE      
                              tau = tau + 5.0_dp
                              inc_in_tau = .TRUE.
                            END IF    
                         ELSE IF (tau > 100.0_dp) THEN 
                            step_stop = .TRUE.
                         ELSE IF (tau< 2*tol_tau*0.1_dp) THEN
                            tau = 0.5_dp                     
                         END IF
                        ! We know now that f_k > f_y1 - f_y2
                          f_best = vf_y
                          line_round = line_round + 1
                     ELSE
                   
                        IF (tau < tol_tau .OR. inc_in_tau ) THEN
                           step_stop = .TRUE.
                           IF (.NOT. inc_in_tau) THEN                         
                              tau_smallest = 0.0_dp
                           END IF     
                        ELSE    
                        
                          IF (line_round == 1) THEN 
                              tau = tol_tau*0.1_dp
                          ELSE
                              tau = 0.5_dp * tau    
                          END IF          
                          
                          line_round = line_round + 1 
                        
                        END IF
                     END IF
                 END DO  
         
              
            ELSE
            
              tau_smallest = 0.0_dp 
              
            END IF            

              
            ! We do not have a descent direction
              IF (tau_smallest < tol_tau) THEN   
              
               ! Escape procedure for the improveemnt function is launced
               IF (( common_obj < escape_tol_H) .OR. ( counter_null == max_null_step )) THEN 
                                
              !------------------------------------------------------- 
              !                ** ESCAPE PROCEDURE  **
              !-------------------------------------------------------          
              
                IF (iprint >= 3) THEN             
                    WRITE(*,*) 'escape', escape_counter_H, 'tau', tau_smallest, 'common_norm', common_obj             
                END IF 
                 
                 fg1_scale_all = fg1_k_all
                 fg2_scale_all = fg2_k_all  
 
               ! Only done if constraints are NOT used in the individual search direction problem
                 IF (.NOT. use_const_in_sep_d) THEN
                          g_counter = g_counter + 1    !The function value for each constraint is used from fg1_scale_all and fg2_scale_all
                 END IF
   
               ! The values of components A_i and B_l for H_1 at x_0
                 DO i = 1, number_of_func                   
                     AandB_k_all(i) = AorB_i(fg1_scale_all, fg2_scale_all, fg1_scale_all, fg2_scale_all ,i)
                 END DO
                 
               ! The improvement function H at x_0
                 H1_current = H1(AandB_k_all)              ! the value of the DC component H_1 at x_0
                 H2_current = H2(fg2_scale_all)            ! the value of the DC component H_2 at x_0                   
            
                escape_counter_H = escape_counter_H + 1
                d_cause = 1.0_dp
           
               IF (escape_counter_H<i_user+1) THEN 
                  ! The user is asked what is done 
                  WRITE(*,*)                  
                  WRITE(*,*) 'Tell what needs to be done by selecting one of the following intetegers:'
                  WRITE(*,*) '  ', 0, ' -   Execute Escape procedure (either decreases all objectives or STOPs the whole algorithm)'
                  DO i = 1,  number_of_obj
                     IF (dec_d_list(i) ) THEN 

                       ! Only done if constraints are NOT used in the individual search direction problem
                         IF (.NOT. use_const_in_sep_d) THEN
                              g_counter = g_counter + 1    !The function value for each constraint is used from fg1_scale_all and fg2_scale_all
                         END IF
                     
                         IF (feasible_d_list(i)) THEN 
                             WRITE(*,*) '  ', i,  ' -   Use the separate direction obtained for the objective ', i, & 
                             & 'with the following changes:'
                             DO j = 1, number_of_obj
                                 WRITE(*,*) '                    Objective ', problem1(j), 'changes from', & 
                                 &  fg_k_all(j)/factors_one(j), 'to', f_for_sep_d(i,j)
                             END DO   
                         END IF                          
                     END IF      
                  END DO   

                  correct_value = .FALSE.
                  DO WHILE( .NOT. correct_value)                  
                      WRITE(*,*)
                      READ *, command
                      
                      IF (command < 0) THEN 
                         WRITE(*,*) 'The integer cannot be negative! Give a new integer.'
                      ELSE IF (command > number_of_obj) THEN 
                         WRITE(*,*) 'The integer is too big! Give a new integer.'    
                      ELSE IF (command == 0) THEN
                          correct_value = .TRUE.
                      ELSE
                          IF (dec_d_list(command)) THEN 
                              correct_value = .TRUE.
                          ELSE 
                              WRITE(*,*) 'The objective ', command, ' cannot be selected. Give a new integer.'   
                          END IF 
                      END IF                      
                      
                  END DO
                  
                  IF (command == 0) THEN 
                      execute_escape = .TRUE.
                      i_user = 0
                  ELSE
                      execute_escape = .FALSE.

                   ! The obtained direction is stored to the correct place in the vector 'mD'                   
                     ind = (user_n)*(command-1)
                     DO j = 1, user_n
                         common_d(j)= mD(ind+j)
                     END DO               
                    
                     x_new = x_current + common_d
                    
                     escape_stop = 0
                  END IF     
               ELSE
                 execute_escape = .TRUE.                  
               END IF     


               IF (execute_escape) THEN             

              ! Escape procedure for the improvement function
                CALL escape_procedure_H( x_current, H1_current, H2_current, fg1_scale_all, fg2_scale_all, &
                            & escape_stop, x_new , H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            & d_cause, crit_tol_H, m_escape_H, step_tol_H, iprint, mrounds_escape, &
                            & help_round_counter, help_escape_f_counter, help_escape_sub_counter, factor, & 
                            & user_n, number_of_obj, number_of_const, number_of_func)                
                                 
                common_d = x_new - x_current
                
              ! Counter updates
                escape_f_counter_H = escape_f_counter_H + help_escape_f_counter
                escape_sub_counter_H = escape_sub_counter_H + help_escape_sub_counter
                escape_g_counter_H = escape_g_counter_H + help_escape_f_counter
                escape_sub_g_counter_H = escape_sub_g_counter_H + help_escape_sub_counter
                
                END IF
                
                IF (escape_stop == 0) THEN  ! New better point found
              !------------------------------------------------------- 
              !         ** ESCAPE PROCEDURE: SERIOUS STEP  **
              !-------------------------------------------------------                  
                  
                   IF (iprint >= 3) THEN                 
                       WRITE(*,*) 'escape: serious step'             
                   END IF 
 
                  ! The values of DC components f1_i and f2_i at x_new
                    DO i = 1, number_of_obj                       
                       fg1_new_all(i) = f1(x_new,problem1(i),factors_one(i),user_n) 
                       fg2_new_all(i) = f2(x_new,problem2(i),factors_one(i),user_n)
                       fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)                             
                    END DO         
                   
                  ! One function value evaluated for each DC component f1_i      
                    f_counter = f_counter+1         
                    
                    IF (number_of_const>0) THEN
                      ! The values of DC components g1_l and g2_l at x_new                 
                        DO i = number_of_obj+1, number_of_func       
                          fg1_new_all(i) = g1(x_new,problem1(i),factors_one(i),user_n) 
                          fg2_new_all(i) = g2(x_new,problem2(i),factors_one(i),user_n)
                          fg_new_all(i) = fg1_new_all(i)-fg2_new_all(i)                           
                        END DO 
                        IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN                        
                           g_counter = g_counter+1                                                     
                        END IF   
                    END IF

                  ! The values of subgradients for DC components f1_i and f2_i at x_0   
                    DO i = 1, number_of_obj             
                     
                       ! Reset of bundles                  
                        IF (use_reset_b1) THEN
                            CALL reset_b1(B1(i))
                        END IF
                       
                        IF (use_reset_b2) THEN
                            CALL reset_b2(B2(i))
                        END IF
                   
                       grad1 = subgradient_f1(x_new, problem1(i),factors_one(i),user_n)
                       grad2 = subgradient_f2(x_new, problem2(i),factors_one(i),user_n)
                       
                       IF (n_const==0) THEN 
                           subgrad1_counter(i) = subgrad1_counter(i) + 1
                           change1 = fg1_new_all(i) - fg1_k_all(i)
                        
                           subgrad2_counter(i) = subgrad2_counter(i) + 1
                           change2 = fg2_new_all(i) - fg2_k_all(i)             
                        
                           dd = x_new - x_current                           ! the search direction
                           change = change1 - change2                       ! the change in the objective function value                           
                    
                           CALL update_b1(B1(i), grad1, dd, change1, i)     ! bundle update for B_1
                           CALL update_b2(B2(i), grad2, dd, change2)        ! bundle update for B_2
                       ELSE               
                           subgrad1_counter(i) = subgrad1_counter(i) + 1
                           subgrad2_counter(i) = subgrad2_counter(i) + 1

                           DO j = 1, user_n
                              grad1_new_all(i,j) = grad1(j)
                              grad2_new_all(i,j) = grad2(j)
                           END DO    
                       END IF                  
                    END DO  
                    
                    IF (n_const>0) THEN
                    ! The values of subgradients for DC components g1_l and g2_l at x_0
                       DO i = number_of_obj+1, number_of_func
                           grad1 = subgradient_g1(x_new, problem1(i),factors_one(i),user_n)
                           grad2 = subgradient_g2(x_new, problem2(i),factors_one(i),user_n)                      
                           DO j = 1, user_n
                              grad1_new_all(i,j) = grad1(j)
                              grad2_new_all(i,j) = grad2(j)
                           END DO
                       END DO
                       subgrad1_g_counter = subgrad1_g_counter + 1
                       subgrad2_g_counter = subgrad2_g_counter + 1                         
                    END IF   
                      
                    IF (n_const>0) THEN
                      ! The subdradients are calculated in the individual search direction problems
                       DO i = 1, number_of_obj
                                               
                         ! The value of objective in the individual search direction problem (current iteration point)
                           fg1_k_all_sep(1) = fg1_k_all(i) 
                           fg2_k_all_sep(1) = fg2_k_all(i) 
                         ! The values of constraints in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj) 
                                   fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)             
                               END DO
                           END IF                            
                         ! The value of A in the individual search direction problem (current iteration point)
                           AandB_k_all_sep(1) = AB_obj_k(i) 
                         ! The values of B in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   AandB_k_all_sep(j+1) = AB_const_k(i,j)
                               END DO
                           END IF              
                           
                         ! The value of objective in the individual search direction problem (new iteration point)
                           fg1_new_all_sep(1) = fg1_new_all(i) 
                           fg2_new_all_sep(1) = fg2_new_all(i) 
                         ! The values of constraints in the individual search direction problem (new iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_new_all_sep(j+1) = fg1_new_all(j+number_of_obj) 
                                   fg2_new_all_sep(j+1) = fg2_new_all(j+number_of_obj)             
                               END DO
                           END IF
                       
                        ! The subgradient of the individual objctive in the individual search direction problem (new iteration point) 
                           DO j = 1, user_n
                              grad1_new_all_sep(1,j) = grad1_new_all(i,j) 
                              grad2_new_all_sep(1,j) = grad2_new_all(i,j) 
                           END DO
                               
                           IF (n_const>0) THEN  
                             ! The subradients of constraints in the individual search direction problem (new iteration point) 
                               DO j = 1, n_const  
                                 DO k = 1, user_n
                                    grad1_new_all_sep(j+1,k) = grad1_new_all(j+number_of_obj,k) 
                                    grad2_new_all_sep(j+1,k) = grad2_new_all(j+number_of_obj,k)                    
                                 END DO                
                               END DO                  
                           END IF

                         ! the values of components A and B for H_1 at x_new in the individual search direction problem
                           DO j = 1, 1+n_const                 
                               AandB_new_all_sep(j) = AorB_i_sep(fg1_new_all_sep, fg2_new_all_sep, &
                                           & fg1_new_all_sep, fg2_new_all_sep ,j)
                           END DO
                           
                           AB_obj_new(i) = AandB_new_all_sep(1)
                           DO j = 1, n_const
                               AB_const_new(i,j) = AandB_new_all_sep(j+1)              
                           END DO                              
                           
                           H1_new_sep(i) = H1_sep(AandB_new_all_sep)              ! the value of the DC component H_1 at x_0 in the individual search direction problem
                           H2_new_sep(i) = H2_sep(fg2_new_all_sep)                ! the value of the DC component H_2 at x_0 in the individual search direction problem
              
                           ind = index_H1_sep(AandB_new_all_sep)                                           ! the index of the component A_i or B_l yielding the subgradient of H_1 in the individual search direction problem
                           grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, ind, user_n) ! the subgradient of H_1 at x_0 in the individual search direction problem
                           change1 = H1_new_sep(i) - H1_current_sep(i)

                           grad2 = subgradient_H2_sep(grad2_new_all_sep, user_n)                           ! the subgradient of H_2 at x_0 in the individual search direction problem
                           change2 = H2_new_sep(i) - H2_current_sep(i)  
                        
                           dd = x_new - x_current                              ! the search direction
                           change = H1_new_sep(i) - H2_new_sep(i) 
                           change = change - (H1_current_sep(i) -H2_current_sep(i)) ! the change in the objective function value    
                    
                           CALL update_b1_imp(B1(i), grad1, ind, dd, AandB_new_all_sep, AandB_k_all_sep, &  ! bundle update for B_1
                                        & fg1_k_all_sep, fg2_k_all_sep, fg1_new_all_sep, fg2_new_all_sep, &
                                        & 1, n_const+1 )     

                           CALL update_b2(B2(i), grad2, dd, change2)                            ! bundle update for B_2
                           
                           DO j = 1, n_const+1
                              IF ( j/= ind) THEN 
                                 grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, j, user_n)        ! the subgradients of A_i and B_l at x_0 in the individual search direction problem
                                 alpha = 0.0_dp
                                 CALL add_element_b1(B1(i), grad1, alpha, j)
                              END IF
                           END DO                  
                           
                       END DO
                       
                    END IF
        
                
                    x_current = x_new                                   ! update of the current iteration point
                    fg1_k_all = fg1_new_all                             ! update of the function value f_1
                    fg2_k_all = fg2_new_all                             ! update of the function value f_2  
                    fg_k_all = fg1_k_all - fg2_k_all                
                    
                    H1_current_sep = H1_new_sep
                    H2_current_sep = H2_new_sep
                    AB_obj_k = AB_obj_new
                    AB_const_k = AB_const_new
                
                    CALL cpu_time(finish_round)                         ! The end time for the current round of the algorithm                 
                
                  ! The objective values at the new point without scaling 
                    DO i = 1, number_of_func
                       fg1_new_all(i) = fg1_new_all(i)/factors_one(i)
                       fg2_new_all(i) = fg2_new_all(i)/factors_one(i)
                       fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)
                    END DO      


                !----------------------------------------------------------------------------------------
                !----------------------------------------------------------------------------------------

                
                
                  IF (iprint > 1) THEN    
                    IF (number_of_const==0) THEN                  
                         WRITE(*,*) 'Round:', iter_counter, 'f(x):', fg_new_all, 'time', finish_round-start_round       
                    ELSE
                        WRITE(*,*) 'Round:', iter_counter, 'f(x):', fg_k_all(1:number_of_obj), &
                                       & 'g(x):', fg_k_all(number_of_obj+1:number_of_func), &
                                       & 'time', finish_round-start_round
                    END IF                       
                   ELSE IF (iprint==-1) THEN 
                              WRITE(*,*) 'Round:', iter_counter, 'x=', x_current                   
                   END IF 



              !------------------------------------------------------- 
              !      ** END: ESCAPE PROCEDURE: SERIOUS STEP  **
              !------------------------------------------------------- 
                
              ! Stopping condition satisfied (approximate Clarke stationarity)      
                ELSE IF (escape_stop == 1) THEN 
                   stop_alg = .TRUE.
                   termination = 1
                   CALL cpu_time(finish_round)                                    
              
              ! Approximate stopping condition satisfied (i.e. the step-length beta* < beta_0
                ELSE IF (escape_stop == 2) THEN 
                   stop_alg = .TRUE.
                   termination = 2
                   CALL cpu_time(finish_round)                                        
            
              ! The maximum number of rounds executed in 'Escape procedure'         
                ELSE IF (escape_stop == 5) THEN 
                   stop_alg = .TRUE.
                   termination = 5
                   CALL cpu_time(finish_round)                                        

                END IF

                
              !------------------------------------------------------- 
              !                ** END: ESCAPE PROCEDURE  **
              !------------------------------------------------------- 
                
                ELSE
 
              !------------------------------------------------------- 
              !                ** BEGIN: NULL STEP  **
              !------------------------------------------------------- 
 
                   IF (iprint >= 3) THEN
                      WRITE(*,*) 'null_Step', counter_null + 1
                   END IF 
                   
                 ! One new null step  
                   counter_null = counter_null + 1
                 
                 !------------------------------------------
                 !                POINT 1
                 !------------------------------------------
                   
                 ! The point x_new is calculated using the common direction                    
                    x_new = x_current + common_d

                  ! The values of DC components f1_i and f2_i at x_new
                    DO i = 1, number_of_obj                       
                       fg1_new_all(i) = f1(x_new,problem1(i),factors_one(i),user_n) 
                       fg2_new_all(i) = f2(x_new,problem2(i),factors_one(i),user_n)
                       fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)                             
                    END DO         
                   
                  ! One function value evaluated for each DC component f1_i      
                    f_counter = f_counter+1         
                    
                    IF (n_const>0) THEN
                      ! The values of DC components g1_l and g2_l at x_new                 
                        DO i = number_of_obj+1, number_of_func       
                          fg1_new_all(i) = g1(x_new,problem1(i),factors_one(i),user_n) 
                          fg2_new_all(i) = g2(x_new,problem2(i),factors_one(i),user_n)
                          fg_new_all(i) = fg1_new_all(i)-fg2_new_all(i)                                         
                        END DO  
                        IF (use_const_in_sep_d) THEN
                           g_counter = g_counter+1     
                        END IF                         
                    END IF

                   DO i = 1, number_of_obj
                       
                     ! Subgradients of f_1 and f_2 at the new point x_new
                       grad1 = subgradient_f1(x_new, problem1(i),factors_one(i),user_n)
                       grad2 = subgradient_f2(x_new, problem2(i),factors_one(i),user_n)   
                       
                       subgrad1_counter(i) = subgrad1_counter(i) + 1
                       subgrad2_counter(i) = subgrad2_counter(i) + 1
                       
                       IF (n_const==0) THEN 
                           linerr1 = fg1_k_all(i) - fg1_new_all(i) + DOT_PRODUCT(common_d,grad1)       ! a new linearization error for f_1
                           linerr2 = fg2_k_all(i) - fg2_new_all(i) + DOT_PRODUCT(common_d,grad2)       ! a new linearization error for f_2
                        
                           CALL add_element_b1(B1(i), grad1, linerr1, i)
                           CALL add_element_b2(B2(i), grad2, linerr2)
                       ELSE               
                           DO j = 1, user_n
                              grad1_new_all(i,j) = grad1(j)
                              grad2_new_all(i,j) = grad2(j)
                           END DO    
                       END IF                  
                   END DO
                   
                    IF (n_const>0) THEN
                    ! The values of subgradients for DC components g1_l and g2_l at x_0
                       DO i = number_of_obj+1, number_of_func
                           grad1 = subgradient_g1(x_new, problem1(i),factors_one(i),user_n)
                           grad2 = subgradient_g2(x_new, problem2(i),factors_one(i),user_n)
                           DO j = 1, user_n
                              grad1_new_all(i,j) = grad1(j)
                              grad2_new_all(i,j) = grad2(j)
                           END DO
                       END DO
                       subgrad1_g_counter = subgrad1_g_counter + 1
                       subgrad2_g_counter = subgrad2_g_counter + 1  
                    END IF                 
           
                    IF (n_const>0) THEN
                      ! The subdradients are calculated in the individual search direction problems
                       DO i = 1, number_of_obj
                                               
                         ! The value of objective in the individual search direction problem (current iteration point)
                           fg1_k_all_sep(1) = fg1_k_all(i) 
                           fg2_k_all_sep(1) = fg2_k_all(i) 
                         ! The values of constraints in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj) 
                                   fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)             
                               END DO
                           END IF                            
                         ! The value of A in the individual search direction problem (current iteration point)
                           AandB_k_all_sep(1) = AB_obj_k(i) 
                         ! The values of B in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   AandB_k_all_sep(j+1) = AB_const_k(i,j)
                               END DO
                           END IF              
                           
                         ! The value of objective in the individual search direction problem (new iteration point)
                           fg1_new_all_sep(1) = fg1_new_all(i) 
                           fg2_new_all_sep(1) = fg2_new_all(i) 
                         ! The values of constraints in the individual search direction problem (new iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_new_all_sep(j+1) = fg1_new_all(j+number_of_obj) 
                                   fg2_new_all_sep(j+1) = fg2_new_all(j+number_of_obj)             
                               END DO
                           END IF
                       
                        ! The subgradient of the individual objctive in the individual search direction problem (new iteration point) 
                           DO j = 1, user_n
                              grad1_new_all_sep(1,j) = grad1_new_all(i,j) 
                              grad2_new_all_sep(1,j) = grad2_new_all(i,j) 
                           END DO
                               
                           IF (n_const>0) THEN  
                             ! The subradients of constraints in the individual search direction problem (new iteration point) 
                               DO j = 1, n_const  
                                 DO k = 1, user_n
                                    grad1_new_all_sep(j+1,k) = grad1_new_all(j+number_of_obj,k) 
                                    grad2_new_all_sep(j+1,k) = grad2_new_all(j+number_of_obj,k)                    
                                 END DO                
                               END DO                  
                           END IF

                         ! the values of components A and B for H_1 at x_new in the individual search direction problem
                           DO j = 1, 1+n_const                 
                               AandB_new_all_sep(j) = AorB_i_sep(fg1_new_all_sep, fg2_new_all_sep, &
                                           & fg1_k_all_sep, fg2_k_all_sep ,j)
                           END DO                          
                           
                           H1_new_sep(i) = H1_sep(AandB_new_all_sep)              ! the value of the DC component H_1 at x_new in the individual search direction problem
                           H2_new_sep(i) = H2_sep(fg2_new_all_sep)                ! the value of the DC component H_2 at x_new in the individual search direction problem
                           
                           DO j = 1, n_const+1
                                grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, j, user_n)

                                lin_err1 = AandB_k_all_sep(j) - AandB_new_all_sep(j) + &
                                           & DOT_PRODUCT(common_d,grad1)                        ! a new linearization error for A_i or B_l  

                                CALL add_element_b1(B1(i), grad1, lin_err1, j)
                                           
                           END DO

                           grad2 = subgradient_H2_sep(grad2_new_all_sep, user_n)              ! a subgradient of f_2 at x_new  
                           lin_err2 = H2_current_sep(i) - H2_new_sep(i) + DOT_PRODUCT(common_d,grad2)  ! a new linearization error for f_2  
                            
                           CALL add_element_b2(B2(i), grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                                       
                       END DO
                       
                    END IF         
               
                 !------------------------------------------
                 !                POINT 2
                 !------------------------------------------
                 
                   DO i = 1, number_of_obj
                      
                       ind = (user_n)*(i-1) 
                       DO j = 1, user_n
                          common_d(j) = mD(ind+j)
                       END DO
                       
                     ! The point x_new is calculated using the separate direction   
                       x_new = x_current + common_d
                           
                     ! The values of DC components g1_l and g2_l at x_new                 
                       DO j = number_of_obj+1, number_of_func       
                            fg1_new_all(j) = g1(x_new,problem1(j),factors_one(j),user_n) 
                            fg2_new_all(j) = g2(x_new,problem2(j),factors_one(j),user_n)
                            fg_new_all(i) = fg1_new_all(i)-fg2_new_all(i)                                                                   
                       END DO  
                       IF (number_of_const>0 .AND. use_const_in_sep_d ) THEN
                           g_counter = g_counter+1
                       END IF      
                       
                     ! Function values of f_1 and f_2 at the new point x_new
                       fg1_new_all(i) = f1(x_new,problem1(i),factors_one(i),user_n) 
                       fg2_new_all(i) = f2(x_new,problem2(i),factors_one(i),user_n) 
                       fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)
                       f_counter(i) = f_counter(i) + 1
                       
                     ! Subgradients of f_1 and f_2 at the new point x_new
                       grad1 = subgradient_f1(x_new, problem1(i),factors_one(i),user_n)
                       grad2 = subgradient_f2(x_new, problem2(i),factors_one(i),user_n)   
                       subgrad1_counter(i) = subgrad1_counter(i) + 1
                       subgrad2_counter(i) = subgrad2_counter(i) + 1   
                       
                       IF (n_const==0) THEN 
     
                           linerr1 = fg1_k_all(i) - fg1_new_all(i) + DOT_PRODUCT(common_d,grad1)       ! a new linearization error for f_1
                           linerr2 = fg2_k_all(i) - fg2_new_all(i) + DOT_PRODUCT(common_d,grad2)       ! a new linearization error for f_2
                        
                           CALL add_element_b1(B1(i), grad1, linerr1, i)
                           CALL add_element_b2(B2(i), grad2, linerr2)
                       
                       ELSE
                           DO j = 1, user_n
                              grad1_new_all(i,j) = grad1(j)
                              grad2_new_all(i,j) = grad2(j)
                           END DO                  
                       END IF
                       
                      IF (n_const>0) THEN
                      ! The values of subgradients for DC components g1_l and g2_l at x_0
                         DO j = number_of_obj+1, number_of_func
                           grad1 = subgradient_g1(x_new, problem1(j),factors_one(j),user_n)
                           grad2 = subgradient_g2(x_new, problem2(j),factors_one(j),user_n)                      
                           DO k = 1, user_n
                              grad1_new_all(j,k) = grad1(k)
                              grad2_new_all(j,k) = grad2(k)
                           END DO
                         END DO
                         subgrad1_g_counter = subgrad1_g_counter + 1
                         subgrad2_g_counter = subgrad2_g_counter + 1                         
                      END IF               
           
                      IF (n_const>0) THEN
                      ! The subdradients are calculated in the individual search direction problems
                                              
                         ! The value of objective in the individual search direction problem (current iteration point)
                           fg1_k_all_sep(1) = fg1_k_all(i) 
                           fg2_k_all_sep(1) = fg2_k_all(i) 
                         ! The values of constraints in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj) 
                                   fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)             
                               END DO
                           END IF                            
                         ! The value of A in the individual search direction problem (current iteration point)
                           AandB_k_all_sep(1) = AB_obj_k(i) 
                         ! The values of B in the individual search direction problem (current iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   AandB_k_all_sep(j+1) = AB_const_k(i,j)
                               END DO
                           END IF              
                           
                         ! The value of objective in the individual search direction problem (new iteration point)
                           fg1_new_all_sep(1) = fg1_new_all(i) 
                           fg2_new_all_sep(1) = fg2_new_all(i) 
                         ! The values of constraints in the individual search direction problem (new iteration point)
                           IF (n_const>0) THEN  
                               DO j = 1, n_const
                                   fg1_new_all_sep(j+1) = fg1_new_all(j+number_of_obj) 
                                   fg2_new_all_sep(j+1) = fg2_new_all(j+number_of_obj)             
                               END DO
                           END IF
                       
                        ! The subgradient of the individual objctive in the individual search direction problem (new iteration point) 
                           DO j = 1, user_n
                              grad1_new_all_sep(1,j) = grad1_new_all(i,j) 
                              grad2_new_all_sep(1,j) = grad2_new_all(i,j) 
                           END DO
                               
                           IF (n_const>0) THEN  
                             ! The subradients of constraints in the individual search direction problem (new iteration point) 
                               DO j = 1, n_const  
                                 DO k = 1, user_n
                                    grad1_new_all_sep(j+1,k) = grad1_new_all(j+number_of_obj,k) 
                                    grad2_new_all_sep(j+1,k) = grad2_new_all(j+number_of_obj,k)                    
                                 END DO                
                               END DO                  
                           END IF

                         ! the values of components A and B for H_1 at x_new in the individual search direction problem
                           DO j = 1, 1+n_const                 
                               AandB_new_all_sep(j) = AorB_i_sep(fg1_new_all_sep, fg2_new_all_sep, &
                                           & fg1_k_all_sep, fg2_k_all_sep ,j)
                           END DO                          
                           
                           H1_new_sep(i) = H1_sep(AandB_new_all_sep)              ! the value of the DC component H_1 at x_new in the individual search direction problem
                           H2_new_sep(i) = H2_sep(fg2_new_all_sep)                ! the value of the DC component H_2 at x_new in the individual search direction problem
                           
                           DO j = 1, n_const+1
                                grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, j, user_n)

                                lin_err1 = AandB_k_all_sep(j) - AandB_new_all_sep(j) + &
                                           & DOT_PRODUCT(common_d,grad1)                        ! a new linearization error for A_i or B_l  

                                CALL add_element_b1(B1(i), grad1, lin_err1, j)
                                           
                           END DO

                           grad2 = subgradient_H2_sep(grad2_new_all_sep, user_n)              ! a subgradient of f_2 at x_new  
                           lin_err2 = H2_current_sep(i) - H2_new_sep(i) + DOT_PRODUCT(common_d,grad2)  ! a new linearization error for f_2  
                            
                           CALL add_element_b2(B2(i), grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                  
                      END IF     

 
                   END DO
 
              !------------------------------------------------------- 
              !                ** END: NULL STEP  **
              !------------------------------------------------------- 
 
                END IF 
                
              
              ELSE
              
              
              !------------------------------------------------------- 
              !                  ** SERIOUS STEP  **
              !-------------------------------------------------------                
                
                counter_null = 0                              ! The counter for consecutive null step can be initialized

              ! The new better iteration point
                x_new = x_current + tau_smallest * common_d
                
                tot_ratio = 0.0_dp
                
              ! The values of DC components f1_i and f2_i at x_new
                DO i = 1, number_of_obj                       
                   fg1_new_all(i) = f1(x_new,problem1(i),factors_one(i),user_n) 
                   fg2_new_all(i) = f2(x_new,problem2(i),factors_one(i),user_n)
                   fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)                             
                END DO         
               
              ! One function value evaluated for each DC component f1_i      
                f_counter = f_counter+1         
                
                IF (number_of_const>0) THEN
                  ! The values of DC components g1_l and g2_l at x_new                 
                    DO i = number_of_obj+1, number_of_func       
                      fg1_new_all(i) = g1(x_new,problem1(i),factors_one(i),user_n) 
                      fg2_new_all(i) = g2(x_new,problem2(i),factors_one(i),user_n)
                      fg_new_all(i) = fg1_new_all(i)-fg2_new_all(i)                                                           
                    END DO 
                    IF (use_const_in_sep_d) THEN                    
                        g_counter = g_counter+1 
                    END IF                      
                END IF

             ! The values of subgradients for DC components f1_i and f2_i at x_0   
               DO i = 1, number_of_obj             

                   ! Reset of bundles                  
                    IF (use_reset_b1) THEN
                        CALL reset_b1(B1(i))
                    END IF
                   
                    IF (use_reset_b2) THEN
                        CALL reset_b2(B2(i))
                    END IF
               
                   grad1 = subgradient_f1(x_new, problem1(i),factors_one(i),user_n)
                   grad2 = subgradient_f2(x_new, problem2(i),factors_one(i),user_n)
                   
                   IF (n_const==0) THEN 
                       subgrad1_counter(i) = subgrad1_counter(i) + 1
                       change1 = fg1_new_all(i) - fg1_k_all(i)
                    
                       subgrad2_counter(i) = subgrad2_counter(i) + 1
                       change2 = fg2_new_all(i) - fg2_k_all(i)             
                    
                       dd = x_new - x_current                           ! the search direction
                       change = change1 - change2                       ! the change in the objective function value                           
                
                       CALL update_b1(B1(i), grad1, dd, change1, i)     ! bundle update for B_1
                       CALL update_b2(B2(i), grad2, dd, change2)        ! bundle update for B_2
                   ELSE               
                       subgrad1_counter(i) = subgrad1_counter(i) + 1
                       subgrad2_counter(i) = subgrad2_counter(i) + 1

                       DO j = 1, user_n
                          grad1_new_all(i,j) = grad1(j)
                          grad2_new_all(i,j) = grad2(j)
                       END DO    
                   END IF                  
               END DO  
                   
                IF (n_const>0) THEN
                ! The values of subgradients for DC components g1_l and g2_l at x_0
                   DO i = number_of_obj+1, number_of_func
                       grad1 = subgradient_g1(x_new, problem1(i),factors_one(i),user_n)
                       grad2 = subgradient_g2(x_new, problem2(i),factors_one(i),user_n)                       
                       DO j = 1, user_n
                          grad1_new_all(i,j) = grad1(j)
                          grad2_new_all(i,j) = grad2(j)
                       END DO
                   END DO
                   subgrad1_g_counter = subgrad1_g_counter + 1
                   subgrad2_g_counter = subgrad2_g_counter + 1                     
                END IF   
                  
                IF (n_const>0) THEN
                  ! The subdradients are calculated in the individual search direction problems
                   DO i = 1, number_of_obj
                                           
                     ! The value of objective in the individual search direction problem (current iteration point)
                       fg1_k_all_sep(1) = fg1_k_all(i) 
                       fg2_k_all_sep(1) = fg2_k_all(i) 
                     ! The values of constraints in the individual search direction problem (current iteration point)
                       IF (n_const>0) THEN  
                           DO j = 1, n_const
                               fg1_k_all_sep(j+1) = fg1_k_all(j+number_of_obj) 
                               fg2_k_all_sep(j+1) = fg2_k_all(j+number_of_obj)             
                           END DO
                       END IF                            
                     ! The value of A in the individual search direction problem (current iteration point)
                       AandB_k_all_sep(1) = AB_obj_k(i) 
                     ! The values of B in the individual search direction problem (current iteration point)
                       IF (n_const>0) THEN  
                           DO j = 1, n_const
                               AandB_k_all_sep(j+1) = AB_const_k(i,j)
                           END DO
                       END IF              
                       
                     ! The value of objective in the individual search direction problem (new iteration point)
                       fg1_new_all_sep(1) = fg1_new_all(i) 
                       fg2_new_all_sep(1) = fg2_new_all(i) 
                     ! The values of constraints in the individual search direction problem (new iteration point)
                       IF (n_const>0) THEN  
                           DO j = 1, n_const
                               fg1_new_all_sep(j+1) = fg1_new_all(j+number_of_obj) 
                               fg2_new_all_sep(j+1) = fg2_new_all(j+number_of_obj)             
                           END DO
                       END IF
                   
                    ! The subgradient of the individual objctive in the individual search direction problem (new iteration point) 
                       DO j = 1, user_n
                          grad1_new_all_sep(1,j) = grad1_new_all(i,j) 
                          grad2_new_all_sep(1,j) = grad2_new_all(i,j) 
                       END DO
                           
                       IF (n_const>0) THEN  
                         ! The subradients of constraints in the individual search direction problem (new iteration point) 
                           DO j = 1, n_const  
                             DO k = 1, user_n
                                grad1_new_all_sep(j+1,k) = grad1_new_all(j+number_of_obj,k) 
                                grad2_new_all_sep(j+1,k) = grad2_new_all(j+number_of_obj,k)                    
                             END DO                
                           END DO                  
                       END IF

                     ! the values of components A and B for H_1 at x_new in the individual search direction problem
                       DO j = 1, 1+n_const                 
                           AandB_new_all_sep(j) = AorB_i_sep(fg1_new_all_sep, fg2_new_all_sep, &
                                       & fg1_new_all_sep, fg2_new_all_sep ,j)
                       END DO
                       
                       AB_obj_new(i) = AandB_new_all_sep(1)
                       DO j = 1, n_const
                           AB_const_new(i,j) = AandB_new_all_sep(j+1)              
                       END DO                              
                       
                       H1_new_sep(i) = H1_sep(AandB_new_all_sep)              ! the value of the DC component H_1 at x_0 in the individual search direction problem
                       H2_new_sep(i) = H2_sep(fg2_new_all_sep)                ! the value of the DC component H_2 at x_0 in the individual search direction problem
          
                       ind = index_H1_sep(AandB_new_all_sep)                                           ! the index of the component A_i or B_l yielding the subgradient of H_1 in the individual search direction problem
                       grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, ind, user_n) ! the subgradient of H_1 at x_0 in the individual search direction problem
                       change1 = H1_new_sep(i) - H1_current_sep(i)

                       grad2 = subgradient_H2_sep(grad2_new_all_sep, user_n)                           ! the subgradient of H_2 at x_0 in the individual search direction problem
                       change2 = H2_new_sep(i) - H2_current_sep(i)  
                    
                        dd = x_new - x_current                              ! the search direction
                        change = H1_new_sep(i) - H2_new_sep(i) 
                        change = change - (H1_current_sep(i) -H2_current_sep(i)) ! the change in the objective function value    
                
                        CALL update_b1_imp(B1(i), grad1, ind, dd, AandB_new_all_sep, AandB_k_all_sep, &  ! bundle update for B_1
                                    & fg1_k_all_sep, fg2_k_all_sep, fg1_new_all_sep, fg2_new_all_sep, &
                                    & 1, n_const+1 )     

                        CALL update_b2(B2(i), grad2, dd, change2)                            ! bundle update for B_2
                       
                       DO j = 1, n_const+1
                          IF ( j/= ind) THEN 
                             grad1 = subgradient_AorB_i_sep(grad1_new_all_sep, grad2_new_all_sep, j, user_n)        ! the subgradients of A_i and B_l at x_0 in the individual search direction problem
                             alpha = 0.0_dp
                             CALL add_element_b1(B1(i), grad1, alpha, j)
                          END IF
                       END DO                  
                       
                   END DO
                   
                END IF
 
                x_current = x_new                                   ! update of the current iteration point
                fg1_k_all = fg1_new_all                             ! update of the function value f_1
                fg2_k_all = fg2_new_all                             ! update of the function value f_2  
                fg_k_all = fg1_k_all - fg2_k_all                
            

                H1_current_sep = H1_new_sep
                H2_current_sep = H2_new_sep
                AB_obj_k = AB_obj_new
                AB_const_k = AB_const_new               
            
            
                CALL cpu_time(finish_round)                         ! The end time for the current round of the algorithm
                
                
              ! The objective values at the new point without scaling 
                DO i = 1, number_of_func
                   fg1_new_all(i) = fg1_new_all(i)/factors_one(i)
                   fg2_new_all(i) = fg2_new_all(i)/factors_one(i)
                   fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)
                END DO



                IF (iprint > 1) THEN    
                    IF (number_of_const==0) THEN                  
                         WRITE(*,*) 'Round:', iter_counter, 'f(x):', fg_new_all, 'time', finish_round-start_round, &
                            & 'tau', tau_smallest        
                    ELSE
                        WRITE(*,*) 'Round:', iter_counter, 'f(x):', fg_k_all(1:number_of_obj), &
                                       & 'g(x):', fg_k_all(number_of_obj+1:number_of_func), &
                                       & 'time', finish_round-start_round, 'tau', tau_smallest 
                    END IF  
                ELSE IF (iprint==-1) THEN 
                              WRITE(*,*) 'Round:', iter_counter, 'x=', x_current                    
                END IF 

                !----------------------------------------------------------------------------------------
                !----------------------------------------------------------------------------------------


              !------------------------------------------------------- 
              !                ** END: SERIOUS STEP ** 
              !-------------------------------------------------------
              END IF              


               
           END DO
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------         
 
         ! the maximum number of 'main iterations' is executed and this causes the termination 
           IF ( iter_counter >= mit ) THEN  
               termination = 4 
           END IF
           
           
         ! the solution to the minimization problem 
           x_solution = x_current               
           
         ! The objective values at x_solution without scaling
           DO i = 1, number_of_obj
              fg1_new_all(i) = f1(x_solution,problem1(i),1.0_dp,user_n) 
              fg2_new_all(i) = f2(x_solution,problem2(i),1.0_dp,user_n) 
              fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)
           END DO          
 
           DO i = number_of_obj+1, number_of_func
              fg1_new_all(i) = g1(x_solution,problem1(i),1.0_dp,user_n) 
              fg2_new_all(i) = g2(x_solution,problem2(i),1.0_dp,user_n) 
              fg_new_all(i) = fg1_new_all(i) - fg2_new_all(i)
           END DO          
           
         ! the objective function value at the solution  
           f_solution =  fg_new_all                       
           
         ! the values of different counters
           counter(1,1) = iter_counter
           counter(9,1) = escape_counter_H
           counter(10,1) = escape_f_counter_H
           counter(11,1) = escape_sub_counter_H 
           
           DO i = 1, number_of_obj
              counter(2,i) = subprob_counter(i)
              counter(3,i) = f_counter(i)
              counter(4,i) = subgrad1_counter(i)
              counter(5,i) = subgrad2_counter(i) 
              counter(6,i) = escape_counter_d(i)
              counter(7,i) = escape_f_counter_d(i)
              counter(8,i) = escape_sub_counter_d(i)              
              counter(12,i) = f_counter_line(i)
           END DO
           
           IF (number_of_const>0) THEN
               DO i = 1, number_of_const
                  counter_g(1,i) = g_counter(i) 
                  counter_g(2,i) = subgrad1_g_counter(i) 
                  counter_g(3,i) = subgrad2_g_counter(i) 
                  counter_g(4,i) = g_counter_line(i)
                  counter_g(5,i) = escape_g_counter_d(i)
                  counter_g(6,i) = escape_sub_g_counter_d(i) 
                  counter_g(7,i) = escape_f_counter_H 
                  counter_g(8,i) = escape_sub_counter_H                   
               END DO
           END IF 
               
            CALL cpu_time(finish_time)         ! Stop CPU timing
            CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
          ! Calculate the elapsed 'clock' time in seconds:
            elapsed_time=(1.0_dp*clock_end-clock_start)/clock_rate   
            time = elapsed_time         
 
 
            IF (iprint >= 1 .OR. iprint <=-1) THEN
                WRITE(*,*)              
                WRITE(*,*)              
                IF (number_of_const==0) THEN
                   WRITE(*,*) 'Finish:', 'f(x):', fg_new_all, 'time', elapsed_time              
                ELSE
                   WRITE(*,*) 'Finish:', 'f(x):', fg_k_all(1:number_of_obj), &
                                       & 'g(x):', fg_k_all(number_of_obj+1:number_of_func), &
                                       & 'time', elapsed_time
                END IF 
                IF (iprint>2 .OR. iprint<-2) THEN 
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i, ')=',x_solution(i)
                   END DO
                END IF 
                WRITE(*,*) 'time', elapsed_time
                WRITE(*,*) 'Rounds', iter_counter
                WRITE(*,*)
                WRITE(*,*) 'termination', termination
                WRITE(*,*)
                WRITE(*,*) 'f_counter', f_counter
                WRITE(*,*) 'sub1_counter', subgrad1_counter
                WRITE(*,*) 'sub2_counter', subgrad2_counter
                WRITE(*,*)
                WRITE(*,*) 'For improvement function:'
                WRITE(*,*) 'escape_counter_H', escape_counter_H
                WRITE(*,*) 'escape_f_counter_H', escape_f_counter_H
                WRITE(*,*) 'escape_sub_counter_H', escape_sub_counter_H
                WRITE(*,*) 
                WRITE(*,*) 'For each objective:'
                WRITE(*,*) 'escape_counter_d', escape_counter_d
                WRITE(*,*) 'escape_f_counter_d', escape_f_counter_d
                WRITE(*,*) 'escape_sub_counter_d', escape_sub_counter_d
                
                IF (number_of_const>0) THEN 
                    WRITE(*,*) 'g_counter', g_counter
                    WRITE(*,*) 'sub1_g_counter', subgrad1_g_counter
                    WRITE(*,*) 'sub2_g_counter', subgrad2_g_counter             
                END IF
                
            END IF    
            
             ! Deallocation of bundles B_1 and B_2 for each objective function
               DO i = 1, number_of_obj                       
                  CALL deallocation_b1(B1(i))    
                  CALL deallocation_b2(B2(i))            
               END DO 
           
           CLOSE(40)   
           
           END SUBROUTINE new_algorithm      
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        
        
         
       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |               THE DESCENT DIRECTION FOR ONE OBJECTIVE                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE separate_direction(ind_f, scaling, x_k, f1_k, f2_k, f_val, user_n, &
                            & x_new , f1_new, f2_new, f_0, B1, B2, &
                            & crit_tol, eps, m, c, r_dec, r_inc, m_escape, step_tol, mrounds, mrounds_escape, &
                            & reason_for_stop, iter_counter, stop_cond_counter, &
                            & subprob_counter, fi_counter, subgrad1_counter, subgrad2_counter,&
                            & agg_in_use, t, t_min, t_max, t_adjust, i_t, decrease, & 
                            & clarke_f_counter, clarke_sub_counter, number_of_const, AB_k_all, &
                            & fg1_k_all, fg2_k_all, scaling_const )                           
            !
            ! Executes the subroutine finding a search direction at the current iteration point for a fixed objective function. 
            !
            ! INPUT: * 'ind_f'             : The index of the objective function
            !        * 'scaling'           : The scaling factor for the objective
            !        * 'x_k'               : The current iteration point
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !        * 'f_val'             : The value of objective function utilized in the descent condition
            !        * 'f_0'               : The value of the objective function f at a starting point x_0
            !        * 'AB_k_all'          : The values of the components A_i and B_l at the current iteration point
            !
            !        * 'fg1_k_all'         : The values of the DC components f1_i and g1_l at the current iteration point 
            !        * 'fg2_k_all'         : The values of the DC components f2_i and g2_l at the current iteration point
            !
            !        * 'user_n'            : The dimension of the problem 
            !        * 'crit_tol'          : The stopping tolerance for the escape procedure
            !        * 'eps'               : The enlargement parameter
            !        * 'm'                 : The descent parameter when the search direction is determined
            !        * 'c'                 : The decrease parameter
            !        * 'r_dec' and 'r_inc' : The decrease and increase parameters   
            !
            !        * 'm_escape'          : The descent parameter in escape procedure  
            !        * 'step_tol'          : The step-size parameter in escape procedure  
            !
            !        * 'mrounds'           : The maximum number of possible rounds during the search direction determination
            !        * 'mrounds_escape'    : The maximum number of possible rounds during the escape procedure
            !
            !        * 'agg_in_use'        : If .TRUE. then aggregation is used when search direction is determined
            !
            !        * 'number_of_const'   : The number of constraints
            !
            ! OUTPUT: * 'x_new'              : the new iteration point
            !         * 'f1_new' and 'f2_new': the new value of DC components f_1 and f_2 at 'x_new'
            !         * 'reason_for_stop'    : Indicates the cause of the termination during the search direction determination
            !         * 'decrease'           : approximation of descent in objective f          
            !
            !         * 'iter_counter'       : the number of rounds needed during the search direction determination
            !         * 'stop_cond_counter'  : the number of times the escape procedure executed (either 1 or 0)
            !         * 'subprob_counter'    : the number of subproblems solved during the search direction determination
            !         * 'fi_counter'         : the number of function values evaluated for a DC component during the search direction determination (same for f_1 and f_2) 
            !         * 'subgrad1_counter'   : the number of subgradients calculated for f_1 during the search direction determination
            !         * 'subgrad2_counter'   : the number of subgradients calculated for f_2 during the search direction determination 
            !
            ! INOUT: * 'B_1' and B_2'        : the bundles of the DC components f_1 and f_2
            !
            !        * 't'                   : the proximity parameter
            !        * 't_min'               : the lower bound for t
            !        * 't_max'               : the upper bound for t
            !        * 't_adjust'            : if .TRUE. then proximity parameter adjustment is in use
            !
            !        * 'i_t'                 : the index used in update procedure of proximity parameter t (in the future perhaps?)           
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimension of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE separate_direction is used in SUBROUTINE new_algorithm.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER *********************************  
            
               TYPE(kimppu1), INTENT(INOUT) :: B1                   ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                   ! the bundle B_2 for the DC component f_2
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 during the search direction determination
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fg1_k_all  ! the values of DC components f1_i at the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fg2_k_all  ! the values of DC components f2_i at the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: AB_k_all   ! the values of components A_i at the starting point x_k
               
               REAL(KIND=dp), INTENT(IN)  :: scaling_const   ! the values of components A_i at the starting point x_k

               REAL(KIND=dp), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), INTENT(INOUT)  :: f_val           ! the value utilized in descent condition
               REAL(KIND=dp), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' during the search direction determination
               REAL(KIND=dp), INTENT(IN)  :: f_0             ! the value of f at the starting point x_0
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol       ! crit_tol=stopping tolerance 
               REAL(KIND=dp), INTENT(IN) :: eps            ! eps= enlargement parameter
               REAL(KIND=dp), INTENT(IN) :: m              ! m=descent parameter
               REAL(KIND=dp), INTENT(IN) :: c              ! c=decrease parameter 
               REAL(KIND=dp), INTENT(IN) :: r_dec, r_inc   ! r_dec=decrease parameter and R_inc=increase parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol       ! step-size tolerance in escape procedure
               REAL(KIND=dp), INTENT(IN) :: m_escape       ! descent parameter in escape procedure 
               
               REAL(KIND=dp), INTENT(INOUT) :: t             ! the proximity parameter t
               REAL(KIND=dp), INTENT(INOUT) :: t_min, t_max  ! the bounds for proximity parameter t            
               INTEGER, INTENT(INOUT) :: i_t                 ! the index used in update procedure of proximity parameter t   
               
               REAL(KIND=dp), INTENT(IN) :: scaling          ! scaling factor used for the objective f  
               
               REAL(KIND=dp), INTENT(OUT) :: decrease        ! approximation of descent in objective f                 
               
               INTEGER, INTENT(IN) :: user_n                 ! the dimension of the problem
               INTEGER, INTENT(IN) :: ind_f                  ! Determines the index of the objective function
               
               INTEGER, INTENT(IN) :: mrounds                ! the maximum number of possible rounds during the search direction determination
               INTEGER, INTENT(IN) :: mrounds_escape         ! the maximum number of possible rounds during the search direction determination
                                                             ! If mrounds<=0, then DEFAULT value 500 is used (The correct value is determined in the SUBROUTINE new_algorithm()

               INTEGER, INTENT(OUT) :: reason_for_stop       !  0 - a new iteration point found
                                                             !  1  -  the stopping condition satisfied for a separate objective (Clarke stationarity)
                                                             !  2  -  the approximate stopping condition satisfied for a separate objective (i.e. the step-length beta* < step_tol_d)  
                                                             !  3  -  the biggest possible number of rounds executed during the descent direction determination
                                                             !  5  -  the biggest possible number of rounds executed in the 'Escape procedure' when descent direction determination is done
               
               INTEGER, INTENT(OUT) :: iter_counter        ! the number of rounds used during the descent direction determination
               INTEGER, INTENT(OUT) :: stop_cond_counter   ! the number of escape procedures (either 1 or 0)
               INTEGER, INTENT(OUT) :: subprob_counter     ! the number of subproblems solved during the search direction determination
               INTEGER, INTENT(OUT) :: fi_counter          ! the number of function values evaluated for a DC component iduring the search direction determination (same for f_1 and f_2) 
               INTEGER, INTENT(OUT) :: subgrad1_counter    ! the number of subgradients calculated for f_1 during the search direction determination 
               INTEGER, INTENT(OUT) :: subgrad2_counter    ! the number of subgradients calculated for f_2 during the search direction determination 

               INTEGER, INTENT(OUT) :: clarke_f_counter    ! the number of function values evaluated for f in 'Escape procedure'
               INTEGER, INTENT(OUT) :: clarke_sub_counter  ! the number of subgradients evaluated for f in 'Escape procedure'
               
               INTEGER, INTENT(IN) :: number_of_const      ! the number of constraints

               LOGICAL, INTENT(IN) :: agg_in_use           ! If .TRUE. then aggregation is used during the search direction determination
               LOGICAL, INTENT(IN) :: t_adjust             ! If .TRUE. then the proximity parameter is adjusted in the algorithm
              

           !***************************** LOCAL VARIABLES ************************************
 
               REAL(KIND=dp), DIMENSION(number_of_const+1)  :: fg1_new_all   ! the values of DC components f1_i at the new iteration point
               REAL(KIND=dp), DIMENSION(number_of_const+1)  :: fg2_new_all  ! the values of DC components f2_i at the new iteration point
               
               REAL(KIND=dp), DIMENSION(number_of_const+1)  :: fg1_y_all   ! the values of DC components f1_i at the new iteration point
               REAL(KIND=dp), DIMENSION(number_of_const+1)  :: fg2_y_all  ! the values of DC components f2_i at the new iteration point
               REAL(KIND=dp), DIMENSION(number_of_const+1)  :: AB_y_all  ! the values of DC components f2_i at the new iteration point
               
               REAL(KIND=dp), DIMENSION(number_of_const+1,user_n)  :: grad1_y_all  ! the values of DC components f2_i at the new iteration point
               REAL(KIND=dp), DIMENSION(number_of_const+1,user_n)  :: grad2_y_all  ! the values of DC components f2_i at the new iteration point
 
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: d_t                    ! the search direction
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: y                      ! the new auxiliary point
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at y
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: agg_grad               ! the new aggregated subgradient of f_1
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: apu_grad               ! 'help' subgradient
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: vect                   ! 'help' vector               
               
               REAL(KIND=dp), DIMENSION(number_of_const+1) :: factor             ! the scaling factors in the individual search direction problem              
               
               REAL(KIND=dp) :: norm1              ! ||\bxi_1(x_k)||
               REAL(KIND=dp) :: max_norm           ! the value of the maximum subgradient norm ||\bxi_{2,max}||
               REAL(KIND=dp) :: div_t              ! 1.0_dp divided by the proximity parameter t
               
               REAL(KIND=dp) :: d_norm             ! the norm of the direction vector d_t
               REAL(KIND=dp) :: delta1, delta2     ! the predicted changes of f_1 and f_2, respectively
               REAL(KIND=dp) :: f1_y, f2_y         ! the function values of f_1 and f_2 at y
               REAL(KIND=dp) :: f_k                ! the function value of f at x_k
               REAL(KIND=dp) :: real_decrease      ! the real decrease of the objective function f
               REAL(KIND=dp) :: lin_err1, lin_err2 ! the linearization errors of f_1 and f_2 (calculated using y)
               REAL(KIND=dp) :: agg_lin_err        ! the new aggregated linearization error of f_1
               
               REAL(KIND=dp) :: help               ! 'help' variable  
               REAL(KIND=dp) :: f_y               ! 'help' variable  

               REAL(KIND=dp) :: test_dec       ! approximation of decrease  

               INTEGER :: f_counter                ! help counter for the number of function values
               INTEGER :: sub_counter              ! help counter for the number of subgradients 
               INTEGER :: clarke_iter_counter      ! help counter for the number  of iterations in Clarke stationary algorithm 
              
               INTEGER :: s_counter                ! the number of subproblems solved during the current iteration round of the search direction determination
               INTEGER :: i, ind, j                   ! help variables
                           
               LOGICAL :: t_changed                ! .TRUE. if t has been changed during the previous round (.TRUE. also when initialization is done)
               LOGICAL :: was_b1_full              ! .TRUE. if during the previous round something was added to B_1 and B_1 was full before this insertion
               LOGICAL :: agg_used                 ! .TRUE. if aggregation element is in use
               LOGICAL :: stop_main_it             ! .TRUE. if the search direction determination can be stopped
               
              

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  DESCENT DIRECTION DETERMINATION STARTS  <>  <>  <>  <>  <>  <>  
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           
           !__________________________________________________________________________________     
           !************************ STEP 0: STOPPING CONDITION ******************************
                
                IF (number_of_const>0) THEN 
                    factor = 1.0_dp
                    factor(1) = scaling
                    factor(2) = scaling_const
                END IF 
                
               iter_counter = 1                   ! the first round is executed
               subprob_counter = 0                ! the number of subproblems solved
               fi_counter = 0                     ! the number of function values evaluated
               subgrad1_counter = 0               ! the number of subgradients calculated for f_1
               subgrad2_counter = 0               ! the number of subgradients calculated for f_2
               stop_cond_counter = 0              ! The number of escape procedures
               clarke_f_counter = 0               ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               clarke_sub_counter = 0             ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               agg_used = .FALSE.                 ! Aggregation is initialized to be .FALSE.
           
               stop_main_it = .FALSE.             ! We are NOT ready to stop
               
               vect = give_subgrad_b1(B1,0) - give_subgrad_b2(B2,0)         ! vector: ' \bxi_1(x_k) - \bxi_2(x_k) '

               start: IF ( SQRT(DOT_PRODUCT(vect,vect)) < (crit_tol)) THEN  ! ||\bxi_1(x_k) - \bxi_2(x_k)|| <  crit_tol 
               !**>>**>>**>> APPROXIMATE CRITICALITY OBTAINED <<**<<**<<**
 
                   d_t = 1.0_dp
                   
                 IF (number_of_const==0) THEN   
                  ! 'Escape proceure' is executed
                    CALL escape_procedure_d( x_k, f1_k, f2_k, reason_for_stop, &
                            & x_new , f1_new, f2_new, d_t, &
                            & crit_tol, m_escape, step_tol, mrounds_escape, &
                            & clarke_iter_counter, f_counter, sub_counter, &
                            & scaling, ind_f, user_n)
                            
                  ELSE
                    
                    CALL  escape_procedure_d_const(ind_f, x_k, f1_k, f2_k, fg1_k_all, fg2_k_all, &
                            & reason_for_stop, x_new , f1_new, f2_new, fg1_new_all, fg2_new_all,&
                            &  d_t, crit_tol, m_escape, step_tol, 0, mrounds_escape, &
                            & clarke_iter_counter, f_counter, sub_counter, factor, user_n, & 
                            & 1, number_of_const, number_of_const+1)
                  END IF                  
                    
                   ! Counters are updated   
                   clarke_f_counter = clarke_f_counter + f_counter                  
                   clarke_sub_counter = clarke_sub_counter + sub_counter                    
           
                   stop_cond_counter = stop_cond_counter + 1    ! the update of the stopping condition counter
                   stop_main_it = .TRUE.                        ! the 'main iteration' can be stopped
                   i_t = 0
                  
           !__________________________________________________________________________________         
           !************************ STEP 0: END *********************************************         
               ELSE start   
               !_______________________________________________________________________________
               !************************ STEP 1: PARAMETER INITIALIZATION *********************
           
                   f_k = f1_k - f2_k                       ! the function value of f at x_k
                   
                   vect = give_subgrad_b1(B1,0)            ! the vector \bxi_1(x_k)
                   norm1 = SQRT(DOT_PRODUCT(vect,vect))    ! the value of the norm ||\bxi_1(x_k)||
                   max_norm = max_norm_value(B2)           ! the value of the maximum subgradient norm in the bundle B_2                                  
                                    
                   t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm )  ! the lower bound for t
                   t_max = r_inc * t_min                                  ! the upper bound for t
                   
                   IF (.NOT. t_adjust) THEN
                       t = 0.8*(t_min + t_max)     ! the parameter t is selected      
                   END IF                      
            
                   t_changed = .TRUE.                 ! the parameter t was changed
                   was_b1_full = .FALSE.              ! B_1 was not full during the previous round since this is the first round of 'main iteration'
 
                   IF (t < t_min) THEN
                     t = t_min
                   END IF
                   IF (t > t_max) THEN
                     t = t_max 
                   END IF 


               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               END IF start

               
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------
           
               DO WHILE ( (.NOT. stop_main_it ) .AND. (iter_counter <= mrounds ))   ! is repeated until the 'main iteration' can be terminated         
               !______________________________________________________________________________
               !************************ STEP 2: SEARCH DIRECTION ****************************
                   
                   
                   IF (number_of_const==0) THEN
   
                       IF (agg_in_use) THEN             ! 'agg_in_use' tells whether we use aggregation or not. If 'agg_in_use=.FALSE. then also 'agg_used'=.FALSE.
                           agg_used = is_agg_used(B1)   ! tells whether aggregation element is used or not? 
                       END IF   
                      
                     ! IF we have .FALSE. here then the algorithm does NOT use AGGREGATION                
                       IF (agg_used) THEN            
                       
                            CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+2, &  ! subproblems are solved with aggregation
                                              &  B1, B2, t, s_counter)
                       
                       ELSE ! Aggregation
                       
                            CALL subproblem_solver(x_k, give_n_b1(B1), give_size_b1(B1)+1, &  ! subproblems are solved without aggregation
                                              &  B1, B2, t, s_counter)
                       END IF   
                                        
                   ELSE 
                        CALL subproblem_solver_const(x_k, AB_k_all, give_n_b1(B1),&  ! subproblems are solved
                                        & give_size_b1(B1)+1, B1, B2, t, s_counter, number_of_const+1)
                   END IF                  
                                        
                   subprob_counter = subprob_counter + s_counter    ! the number of subproblems solved so far
                                   
                   CALL add_glob_index(B2)              ! calculates the index of the subproblem yielding the global solution        
                   d_t = give_solution(B2)              ! the global solution d_t is selected   
                   d_norm = SQRT(DOT_PRODUCT(d_t,d_t))  ! the norm ||d_t|| of the solution d_t is calculated
                           
                 
                   ! Aggregated element is calculated. 
                   ! NOTICE: The aggregated element (if used) is now always added into bundle B_1 (also in those cases when the bundle B_1 is not full)

                   IF (number_of_const==0) THEN 
                       IF (agg_in_use) THEN                     ! Is done only when 'agg_in_use'=.TRUE.
                           ind = give_solution_ind(B2)
                           apu_grad = give_subgrad_b2(B2,ind)
                           div_t = 1.0_dp / t
                           DO i = 1, give_n_b2(B2)
                              agg_grad(i) = (-div_t) * d_t(i) + apu_grad(i)
                           END DO
                           agg_lin_err = - give_decrease(B2) - (div_t) * DOT_PRODUCT(d_t,d_t)
                           agg_lin_err = agg_lin_err + give_linerr_b2(B2,ind)
                           
                           CALL add_agg_element_b1(B1,agg_grad,agg_lin_err)
                       END IF
                   END IF 
                  
               !______________________________________________________________________________
               !************************ STEP 2: END *****************************************     
                   
                   
               !->->->->->-> EXECUTION OF BRANCH BEGINS (2 POSSIBLE BRANCHES) <-<-<-<-<-<-<-<-
                   branches: IF (d_norm < crit_tol) THEN                
               !->->->->->->->->->->->->->-> BRANCH 1 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-             
                   !__________________________________________________________________________
                   !**************** STEP 3: ESCAPE PROCEDURE ********************************
                         
                        d_t = 1.0_dp
                        
                        IF (number_of_const==0) THEN    
                          ! 'Escape procedure' is executed
                            CALL escape_procedure_d( x_k, f1_k, f2_k, reason_for_stop, &
                                & x_new , f1_new, f2_new, d_t, &
                                & crit_tol, m_escape, step_tol, mrounds_escape, &
                                & clarke_iter_counter, f_counter, sub_counter, &
                                & scaling, ind_f, user_n)        
                        ELSE
 
                            CALL  escape_procedure_d_const(ind_f, x_k, f1_k, f2_k, fg1_k_all, fg2_k_all, &
                                    & reason_for_stop, x_new , f1_new, f2_new, fg1_new_all, fg2_new_all,&
                                    &  d_t, crit_tol, m_escape, step_tol, 0, mrounds_escape, &
                                    & clarke_iter_counter, f_counter, sub_counter, factor, user_n, & 
                                    & 1, number_of_const, number_of_const+1)
                        END IF    
 
                      ! Counters are updated
                        clarke_f_counter = clarke_f_counter + f_counter                 
                        clarke_sub_counter = clarke_sub_counter + sub_counter   
                       
                        stop_cond_counter = stop_cond_counter + 1   ! the update of the stopping condition counter
                        
                        stop_main_it = .TRUE.                       ! Main iteration can be stopped
                        i_t = 0             

                   !__________________________________________________________________________
                   !******************** STEP 3: EXIT, ESCAPE PROCEDURE **********************  
                
               !->->->->->->->->->->->->->-> BRANCH 1 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                        
                   ELSE branches
               !->->->->->->->->->->->->->-> BRANCH 2 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                   !__________________________________________________________________________
                   !******************** STEP 4: DESCENT TEST ********************************
       
                       y = x_k + d_t     ! a new auxiliary point y
                       
                       IF (number_of_const==0) THEN 
                            f1_y = f1(y, problem1(ind_f), scaling, user_n)      ! the value of f_1 at y
                            f2_y = f2(y, problem2(ind_f), scaling, user_n)      ! the value of f_2 at y   
                            f_y = f1_y - f2_y
                       ELSE
                           DO i = 1, 1                  ! the values of DC components f1_i and f2_i at y
                              fg1_y_all(1) = f1(y, problem1(ind_f),factor(1),user_n) 
                              fg2_y_all(1) = f2(y, problem2(ind_f),factor(1),user_n) 
                           END DO 
                           
                           f_y = fg1_y_all(1) - fg2_y_all(1)
                           
                           DO i = 1, number_of_const   ! the values of DC components g1_l and g2_l at y
                              fg1_y_all(i+1) = g1(y,problem1(i+number_of_obj0),factor(i+1),user_n) 
                              fg2_y_all(i+1) = g2(y,problem2(i+number_of_obj0),factor(i+1),user_n) 
                           END DO  

                           DO i = 1, number_of_const+1              ! the values of components A_i and B_l for H_1 at y
                              AB_y_all(i) = AorB_i_sep(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
                           END DO                  
                           
                           f1_y = H1_sep(AB_y_all)    ! the value of H_1 at y
                           f2_y = H2_sep(fg2_y_all)   ! the value of H_2 at y
                                       
                       END IF                      
                       
                       fi_counter = fi_counter + 1 ! one more objective function value calculated for a DC component

                       real_decrease = ( f1_y - f2_y ) -  f_k  ! the real decrease in the value of f    

                       IF (number_of_const==0) THEN 
                           test_dec = give_decrease(B2)
                       ELSE
                           test_dec = give_decrease(B2)-f1_k
                       END IF

                      
                  IF (test_dec>0) THEN
                      WRITE(*,*) 'big problems!!!!  predcit. decrease',  give_decrease(B2)
                      reason_for_stop=10
                      stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                      
                  ELSE        
                      
                      
                       branch2: IF ( real_decrease <= (m * test_dec ) ) THEN   ! the descent condition is dependent on the value 'f_val' 
                       !**>>**>>**>> NEW ITERATION POINT FOUND <<**<<**<<**  
                                                    
                           x_new = y               ! the new iteration point is the current auxiliary point
                           f1_new = f1_y           ! the value of f_1 at x_new
                           f2_new = f2_y           ! the value of f_2 at x_new
                           stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                           reason_for_stop = 0     ! the reason for stopping is the new iteration point
    
                           decrease  = give_decrease(B2)
                           
                           IF (number_of_const>0) THEN
                             decrease = decrease - f1_k
                           END IF
                           
                   !__________________________________________________________________________
                   !******************** STEP 4: END *****************************************         
                   
                       ELSE branch2 
                       !----------- SOME UPDATES IN VARIABLES BEGINS -------------------------  
                       
                   !__________________________________________________________________________     
                   !******************** STEP 5: BUNDLE UPDATE *******************************
                          
                          update: IF ( (( f_y - f_0)>0 ) .AND. & 
                                                      & (d_norm > eps)) THEN 
                                                      
                          !----------- PARAMETER t REDUCTION --------------------------------
                              t = t - r_dec * ( t - t_min )  ! the parameter t is updated
                              t_changed = .TRUE.             ! the parameter t was changed
                              i_t = -1                       ! null step is done in main iteration

                          ELSE update
                          !----------- BUNDLE AND PARAMETER UPDATE --------------------------
 
                              IF (number_of_const==0) THEN 
                                  new_grad1 = subgradient_f1(y, problem1(ind_f), scaling, user_n)   ! a subgradient of f_1 at y
                                  subgrad1_counter = subgrad1_counter + 1                           ! a new subgradient was evaluated for f_1            

                                  lin_err1 = f1_k - f1_y + DOT_PRODUCT(d_t,new_grad1)       ! a new linearization error for f_1
                                  
                                !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                                  bundle1: IF (is_full_b1(B1)) THEN   
                                  ! bundle B_1 is full and a new element is added to the bundle B_1
                                  ! the overwritten element is written/taken down 
                                     
                                      CALL add_element_b1(B1, new_grad1, lin_err1, ind_f)
                                      was_b1_full = .TRUE.
                                    
                                  ELSE bundle1
                                  ! bundle B_1 is NOT full and a new element is added to the bundle B_1

                                      CALL add_element_b1(B1, new_grad1, lin_err1, ind_f)
                                      was_b1_full = .FALSE.

                                  END IF bundle1                             
                                !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-
                                
                              ELSE
 
                                  DO i = 1, 1                       ! the values of subgradients for DC components f1_i and f2_i at y
                                    new_grad1 = subgradient_f1(y, problem1(ind_f),factor(i),user_n)
                                    new_grad2 = subgradient_f2(y, problem2(ind_f),factor(i),user_n)
                                    DO j = 1, user_n
                                       grad1_y_all(i,j) = new_grad1(j)
                                       grad2_y_all(i,j) = new_grad2(j)
                                    END DO
                                  END DO      
                                  DO i = 1, number_of_const        ! the values of subgradients for DC components g1_l and g2_l at y
                                    new_grad1 = subgradient_g1(y, problem1(i+number_of_obj0),factor(i+1),user_n)
                                    new_grad2 = subgradient_g2(y, problem2(i+number_of_obj0),factor(i+1),user_n)
                                    DO j = 1, user_n
                                       grad1_y_all(i+1,j) = new_grad1(j)
                                       grad2_y_all(i+1,j) = new_grad2(j)
                                    END DO
                                  END DO                              
                                  subgrad1_counter = subgrad1_counter + 1          ! the subgradient counter is updated
                                  subgrad2_counter = subgrad2_counter + 1          ! the subgradient counter is updated

                                  DO j = 1, number_of_const+1

                                    new_grad1 = subgradient_AorB_i_sep(grad1_y_all, grad2_y_all, j, user_n)    ! a subgradient of A_i or B_l at y

                                    lin_err1 = AB_k_all(j) - AB_y_all(j) + &
                                               & DOT_PRODUCT(d_t,new_grad1)                        ! a new linearization error for A_i or B_l                       
                                  
                                   !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                                     bundle1const: IF (is_full_b1(B1)) THEN   
                                     ! bundle B_1 is full and a new element is added to the bundle B_1                                
                                         CALL add_element_b1(B1, new_grad1, lin_err1, j)
                                         was_b1_full = .TRUE.
                                    
                                     ELSE bundle1const
                                     ! bundle B_1 is NOT full and a new element is added to the bundle B_1

                                         CALL add_element_b1(B1, new_grad1, lin_err1, j)
                                         was_b1_full = .FALSE.

                                     END IF bundle1const                          
                                   !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-
                                
                                  END DO

                              END IF                              

                              IF (real_decrease > -m * test_dec)  THEN         ! adjustment of the proximity parameter if necessary
                                  t =  t - c * ( t - t_min )
                              END IF                          
                                 !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE BEGINS <<<<<<<<<<<-             

                                  bundle2: IF((give_max_size_b2(B2)>0) ) THEN   ! bundle B_2 is always updated

                                      IF ( number_of_const==0) THEN  ! NO constraints
                                           new_grad2 = subgradient_f2(y, problem2(ind_f), scaling, user_n)   ! a subgradient of f_2 at y  
                                           subgrad2_counter = subgrad2_counter + 1                 ! a new subgradient was evaluated for f_2
                                           lin_err2 = f2_k - f2_y + DOT_PRODUCT(d_t,new_grad2)     ! a new linearization error for f_2  
                                    
                                           CALL add_element_b2(B2, new_grad2, lin_err2)     ! a new element is inserted into the bundle B_2
                                           !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.        
                                      ELSE
                                           new_grad2 = subgradient_H2_sep(grad2_y_all, user_n)              ! a subgradient of f_2 at y  
                                           lin_err2 = f2_k - f2_y + DOT_PRODUCT(d_t,new_grad2)  ! a new linearization error for f_2  
                                        
                                           CALL add_element_b2(B2, new_grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                                           !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.
                                      END IF                          
                                  
                   !__________________________________________________________________________     
                   !******************** STEP 5: END *****************************************
                   
                   
                   !__________________________________________________________________________     
                   !******************** STEP 6: PARAMETER UPDATE ****************************
                   
                                  help = SQRT(DOT_PRODUCT(new_grad2,new_grad2)) ! the norm of the new subgradient of f_2 (subgradient is calculated 
                                                                                ! at the new auxiliary point y)
                                  IF(help > max_norm) THEN 
                                      max_norm = help            ! the updated value of the maximum subgradient norm in the bundle B_2
                                      t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                  END IF            
     
                   !__________________________________________________________________________     
                   !******************** STEP 6: END *****************************************
                   
                            END IF bundle2
                            
                            i_t = -1                       ! null step is done in main iteration
                            
                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE ENDS <<<<<<<<<<<<<-
                               
                              t_changed = .FALSE.   ! the parameter t was NOT changed during this ELSE branch of 'update'
                                  
                          END IF update        
                        !------ PARAMETER t REDUCTION & BUNDLE AND PARAMETER UPDATE ENDS -----
                        
                          iter_counter = iter_counter + 1  ! update of the iteration counter                        
                              
                       END IF branch2   
                    END IF     
                     !--------- NEW ITERATION POINT & SOME UPDATES IN VARIABLES ENDS --------- 
                     
               !->->->->->->->->->->->->->-> BRANCH 2 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-   
                   END IF branches
               !->->->->->->-> EXECUTION OF BRANCH ENDS  (2 BRANCHES) <-<-<-<-<-<-<-<-<-<-<-<-
              END DO 
              
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------           
              
              IF ( (.NOT. stop_main_it ) ) THEN 
                  reason_for_stop = 3            ! the maximum number of rounds have been executed 
                  iter_counter = iter_counter -1 
              END IF        

              
           END SUBROUTINE separate_direction
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             Escape procedure                                   |  |
        !  |                      guaranteeing Clarke stationarity                          |  |
        !  |                           for separate objective                               |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE escape_procedure_d( x_k, f1_k, f2_k, reason_for_stop,&
                            & x_new , f1_new, f2_new, d_cause,&
                            & crit_tol, m, step_tol, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, scaling, ind_f, user_n )
                                                
            !
            ! Executes the 'Escape procedure' for the separate objective 'problem1(ind_f)'. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary, then this algorithm generates a descent direction for the separate objective.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'f1_k' and 'f2_k'   : The values of DC components f_1 and f_2 at the current iteration point
            !
            !        * 'd_cause'           : The search direction causing the execution of the 'Escape procedure'
            !       
            !        * 'user_n'            : The dimension of the problem 
            !        * 'ind_f'             : The index of the objective (the objective is 'problem1(ind_f)')
            !        * 'scaling'           : The scaling factor used in the objective 'problem1(ind_f)'
            !
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance (i.e. proximity measure)
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Escape procedure'
            !
            !
            ! OUTPUT: * 'x_new'              : the new iteration point (obtained if 'reason_for_stop' = 0)
            !         * 'f1_new' and 'f2_new': the new values of DC components f_1 and f_2 (obtained if 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'Escape procedure'
            !
            !         * 'iter_counter'       : the number of rounds needed in 'Escape procedure'
            !         * 'f_counter'          : the number of function values evaluated for f in 'Escape procedure' (same for f_1 and f_2)
            !         * 'subgrad_counter'    : the number of subgradients calculated for f in 'Escape procedure' (same for f_1 and f_2)
            !
            ! NOTICE: The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !         'user_n' when SUBROUTINE escape_procedure_d is used in SUBROUTINE separate_direction.
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Escape procedure'

               REAL(KIND=dp), INTENT(IN)  :: f1_k, f2_k      ! the value of f_1 and f_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: f1_new, f2_new  ! the value of f_1 and f_2 at x_new if 'reason_for_stop=0' in the 'Escape procedure'
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance 
               REAL(KIND=dp), INTENT(IN) :: m             ! 'm'=descent parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! 'step_tol'=step-length tolerance (i.e. proximity measure)
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'Escape procedure'
                                                          
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < step_tol)
                                                          ! 5 - the maximum number of rounds executed in 'Escape procedure'
                                                          
               INTEGER, INTENT(OUT) :: iter_counter       ! the number of rounds used in 'Escape procedure'
               INTEGER, INTENT(OUT) :: f_counter          ! the number of function values evaluated for f in 'Escape procedure' (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for f in 'Escape procedure'
               
               REAL(KIND=dp), INTENT(IN) :: scaling       ! scaling factor used for the objective f 'problem1(ind_f)'
               INTEGER, INTENT(IN) :: ind_f               ! Determines the index of the objective function
               INTEGER, INTENT(IN) :: user_n              ! The dimension of the problem
               
           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                         ! the bundle B of f containing subgradients from '\partial f(x)'
               
               REAL(KIND=dp), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: x_d                    ! an auxiliary point
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad               ! the subgradient of f at x
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of f_1 and f_2 at x_d
               REAL(KIND=dp), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! a new auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: test_y                 ! another auxiliary point            
               
               REAL(KIND=dp) :: f_k              ! the value of the objective function at the current iteration point   
               REAL(KIND=dp) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=dp) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=dp) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=dp) :: div_eps          ! 1.0_dp / eps     
               REAL(KIND=dp) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=dp) :: f1_y, f2_y       ! the function values of f_1 and f_2 at y
               REAL(KIND=dp) :: test_f1, test_f2 ! the function values of f_1 and f_2 at test_y
               REAL(KIND=dp) :: direc_der        ! directional derivative of f 
               REAL(KIND=dp) :: stepsize         ! stepsize 
               REAL(KIND=dp) :: descent_app      ! approximation of descent in the objective function value
                            
               INTEGER :: size_b                 ! the biggest possible bundle size for B
               INTEGER :: N                      ! the current bundle size for B with the current element
                            
               LOGICAL :: stop_alg               ! .TRUE. if ''Escape procedure' can be stopped
               LOGICAL :: stop_step_det          ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************              
                        

            ! the size of the bundle B
               SELECT CASE(user_size)
                  
                  CASE(0)
                      size_b = (user_n + 5) * 2
                      
                  CASE(1)
                      size_b = MIN((user_n+5)*2, 250)

                  CASE(2)
                      size_b = MIN((user_n+5)*2, 100)                 
                      
               END SELECT 

               
               eps = (10.0_dp)**(-4)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_dp / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               
               f_k = f1_k - f2_k                      ! The value of f at the current iteration point
               
               d = -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))    ! The search direction used to determine subgradient for f         
               
            ! The bundles B is initialized
               CALL init_bundle_b1(B, size_b, user_n) 
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.          ! Algorithm cannot be stopped
               

           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************         

               
               IF (iter_counter > 100) THEN  ! If over 100 iterations are executed then stepsize 'eps' is modified
                   eps = (10.0_dp)**(-6)
               END IF
           
               x_d = x_k + (eps) * d                      ! Auxialiary point 'x_d' is calculated

               new_grad1 = subgradient_f1(x_d, problem1(ind_f), scaling, user_n)  ! Subgradient of f_1 at 'x_d'
               new_grad2 = subgradient_f2(x_d, problem2(ind_f), scaling, user_n)  ! Subgradient of f_2 at 'x_d'
               
               ! new subgradient for f at x_k 
               new_grad = new_grad1 - new_grad2
                   
               subgrad_counter = subgrad_counter +1 
               
               ! New subgradien is added into the bundle 'B'
               IF (N == 0) THEN 
                  CALL add_first_element_b1(B, new_grad, 0)    
                  N = N + 1
               ELSE
                  CALL add_element_b1(B, new_grad, 0.0_dp, 0)
               END IF 
           
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CALRKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
                              
               ! Minimum norm element 'u_k' is determined in the bundle 'B' (i.e. aggregated subgradient)
               CALL quadratic_solver(x_k, user_n, N, B, 1.0_dp, u_k, obj_u) 
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
    
               norm_u = SQRT( DOT_PRODUCT(u_k,u_k))
            
               IF (norm_u < crit_tol) THEN  
                              
                  reason_for_stop = 1     ! Approximate Clarke stationarity is achieved
                  stop_alg = .TRUE.       ! Algorihtm can be stopped
                
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SERACH DIRECTION *****************************            
               
               d =  u_k / norm_u          ! New search direction
               y = x_k + eps * d          ! New auxialiary point    
               
               f1_y = f1(y, problem1(ind_f), scaling, user_n)      ! f_1 at 'y'
               f2_y = f2(y, problem2(ind_f), scaling, user_n)      ! f_2 at 'y'
               f_counter = f_counter +1 
               
               real_decrease = f1_y - f2_y - f_k      ! Real decraese in the objective function f
               
               direc_der = div_eps * real_decrease    ! Directional derivative at 'x_k' into the direction 'd'       
               
               IF (direc_der <= (-m * norm_u)) THEN   ! If .TRUE. decrease in the objective f is sufficient
                     reason_for_stop = 10             ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.                ! algorithm can be stopped  
               END IF
               
               iter_counter = iter_counter +1         ! one more iteration is executed
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1) .AND. stop_alg ) THEN 
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  stepsize = 1.0_dp                      ! initial stepsize
                  stop_step_det = .FALSE.                ! stepsize determination cannot be stopped
                  DO WHILE((.NOT. stop_step_det))
                      test_y = x_k + stepsize * d        ! auxiliary test point
                      test_f1 = f1(test_y, problem1(ind_f), scaling, user_n)     ! the value of f_1 at test_y
                      test_f2 = f2(test_y, problem2(ind_f), scaling, user_n)     ! the value of f_2 at test_y   
                      f_counter = f_counter + 1          
                      descent_app = -(10.0_dp)**(-4)     ! sufficient descent 
                      IF ( (test_f1 - test_f2 - f_k) > MIN(-m * stepsize * norm_u,descent_app )) THEN 
                         stepsize = 0.5_dp * stepsize    ! stepsize is decreased
                         IF (stepsize < step_tol) THEN  
                               stop_step_det = .TRUE.    
                         END IF
                      ELSE
                         stop_step_det = .TRUE.      
                         y = test_y
                         f1_y = test_f1
                         f2_y = test_f2
                      END IF
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-
                            

                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                  
                     x_new = y               ! the new iteration point is the auxilary point y
                     f1_new = f1_y           ! the value of f_1 at x_new
                     f2_new = f2_y           ! the value of f_2 at x_new
                     reason_for_stop = 0     ! the reason for stop is the new iteration point             
                     
                  ELSE
                     reason_for_stop = 2     ! the reason for stop is that the approximate stopping condition is satisfied (i.e. the step-length beta* < step_tol)
                     
                  END IF 
               
                  iter_counter = iter_counter -1    ! one extra iteration in the counter needs to be removed 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                               
              END IF           
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         
          
           
           END SUBROUTINE escape_procedure_d
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
  


        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             Escape procedure                                   |  |
        !  |                      guaranteeing Clarke stationarity                          |  |
        !  |                           for separate objective                               |  |
        !  |                              with constraints                                  |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************         
           
            SUBROUTINE escape_procedure_d_const(ind_f, x_k, H1_k, H2_k, fg1_k_all, fg2_k_all, &
                            & reason_for_stop, x_new , H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            & d_cause, crit_tol, m, step_tol, iprint, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, factor, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)
                                                
            !
            ! Executes the 'Clarke stationary' algorithm. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary then the algorithm generates a descent direction.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'H1_k' and 'H2_k'   : The values of DC components H_1 and H_2 at the current iteration point (individual direction determination with constraints)
            !       
            !        * 'ind_f'             : Determines the index of the objective function
            !  
            !        * 'fg1_k_all' and 'fg2_k_all'   : The values of DC components f1_i, f2_i, g1_l and g2_l at the current iteration point
            !
            !        * 'factor'            : The values of scaling factors for f_i and g_l
            !        * 'd_cause'           : The search direction which caused the execution of the clarke stationary algorithm
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance  
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !        * 'iprint'            : Specifies print
            
            !        * 'user_n'            : The dimension of the problem
            !
            !
            ! OUTPUT: * 'x_new'                 : Tthe new iteration point (if it is found, i.e., 'reason_for_stop' = 0)
            !         * 'H1_new' and 'H2_new'   : The values of DC components H_1 and H_2 at the new iteration point
            !         * 'fg1_new_all' and 'fg2_new_all': the new value of DC components f1_i, f2_i, g1_l and g2_l (if new iteration point is found, i.e., 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'       : Indicates the cause of the termination in the 'Clarke stationarity' algorithm  
            !
            !         * 'iter_counter'          : the number of rounds needed in 'Clarke stationary' algorithm
            !         * 'f_counter'             : the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !         * 'subgrad_counter'       : the number of subgradients calculated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !
            ! NOTICE: * The dimensions of vectors 'x_k', 'x_new' and 'd_cause' has to be same. (i.e. the dimension of 'x_k', 'x_new' and 'd_cause' has to be
            !           'user_n' when SUBROUTINE guaranteeing_clarke() is used in SUBROUTINE main_iteration().
            !         * The dimension of vectors 'fg1_k_all', 'fg2_k_all', 'fg1_new_all' and 'fg2_new_all' has to be 'number_of_func'
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_k_all    ! the values of DC components f1_i and g1_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_k_all    ! the values of DC components f2_i and g2_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg1_new_all ! the values of DC components f1_i and g1_l at x_new
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg2_new_all ! the values of DC components f2_i and g2_l at x_new

               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: factor       ! the values of scaling factors for f_i and g_l           
               
               REAL(KIND=dp), INTENT(IN)  :: H1_k, H2_k      ! the value of H_1 and H_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: H1_new, H2_new  ! the value of H_1 and H_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance (needed in stopping condition)
               REAL(KIND=dp), INTENT(IN) :: m             ! a descent parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! a step-length tolerance
               
               INTEGER, INTENT(IN) :: user_n              ! the dimension of the problem
               INTEGER, INTENT(IN) :: ind_f               ! Determines the index of the objective function
               
               INTEGER, INTENT(IN) :: number_of_obj      ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_const    ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_func     ! The dimension of the problem            
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm
 
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (approximate Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < beta_0)
                                                          ! 5 - the maximum number of rounds executed in 'Clarke stationary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: iter_counter       ! the number of rounds used in 'Clarke stationary' algorithm
               INTEGER, INTENT(OUT) :: f_counter          ! the number of function values evaluated for H in 'Clarke stationary' algorithm (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for H in 'Clarke stationary' algorithm 
               
               INTEGER, INTENT(IN) :: iprint    ! the variable that specifies print:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result (nothing is printed in main iteration)
                                                !   iprint = -1: basic print of final result without the solution vector (nothing is printed in main iteration)                                             
                                                !   iprint = 2 : extended print of final result (nothing is printed in main iteration)
                                                !   iprint = -2: extended print of final result without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results (nothing is printed in main iteration)
                                                !   iprint = -3: basic print of intermediate results and extended print of final results without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 4 : extended print of intermediate and final results (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = -4: extended print of intermediate and final results without the solution vector (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = 5 : everything is printed step by step (this is the only option which causes print in the 'main iteration') 
               

           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                ! the bundle B of H containing subgradients from \partial H(x)
               
               REAL(KIND=dp), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad               ! the subgradient of H at x
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of H_1 and H_2 at x_d
               REAL(KIND=dp), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! the new auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: test_y                 ! the test auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: lower_y                ! the 'hepl' auxiliary point            
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_y_all    ! the values of DC components f1_i and g1_l at y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_y_all    ! the values of DC components f2_i and g2_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_y_all     ! the values of components A_i and B_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_test_all     ! the values of DC components f1_i and g1_l at test_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_test_all     ! the values of DC components f2_i and g2_l at test_y   
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_test_all      ! the values of components A_i and B_l at test_y              
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_lower_all    ! the values of DC components f1_i and g1_l at lower_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_lower_all    ! the values of DC components f2_i and g2_l at lower_y  
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_lower_all     ! the values of components A_i and B_l at lower_y              

               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg1_grad_all    ! the subgradients of DC components f1_i and g1_l
               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg2_grad_all    ! the subgradients of DC components f2_i and g2_l
               
               REAL(KIND=dp) :: H_k              ! the value of the objective function H at the current iteration point 
               REAL(KIND=dp) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=dp) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=dp) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=dp) :: div_eps          ! 1.0_dp / eps     
               REAL(KIND=dp) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=dp) :: H1_y, H2_y       ! the function values of H_1 and H_2 at y
               REAL(KIND=dp) :: test_H1, test_H2 ! the function values of H_1 and H_2 at test_y
               REAL(KIND=dp) :: lower_H1, lower_H2 ! the function values of H_1 and H_2 at lower_y
               REAL(KIND=dp) :: direc_der        ! directional derivative of H 
               REAL(KIND=dp) :: stepsize         ! stepsize 
               REAL(KIND=dp) :: upper, lower     ! bounds for stepsize 
               REAL(KIND=dp) :: descent_app      ! approximation of descent in the objective function value
                              
               INTEGER :: size_b                 ! the biggest possible bundle size for B
               INTEGER :: N                      ! the current bundle size for B with the current element
               
               INTEGER :: i, j, l 
               
               LOGICAL :: ylos               ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_alg               ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_step_det          ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             
              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '              CLARKE STATIONARY ALGORITHM STARTS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF             
             
           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************              
                        
            ! the size of the bundle B
               SELECT CASE(user_size)
                  
                  CASE(0)
                      size_b = (user_n + 5) * 2
                      
                  CASE(1)
                      size_b = MIN((user_n+5)*2, 250)

                  CASE(2)
                      size_b = MIN((user_n+5)*2, 100)                 
                      
               END SELECT 
               
               eps = (10.0_dp)**(-4)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_dp / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               l = 0
               
               H_k = H1_k - H2_k
            
               d =  -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))       
               
               
            ! The bundles B is initialized
               CALL init_bundle_b1(B, size_b, user_n) 
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.
               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 0: INITIALIATION '
                   WRITE(*,*)''
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
               
               y = x_k + (eps) * d 
               
               DO i = 1, 1                     ! the values of DC components f1_i and f2_i at y
                   fg1_y_all(i) = f1(y,problem1(ind_f),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(ind_f),factor(i),user_n) 
               END DO  
               
               DO i = 1, number_of_const! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i+1) = g1(y,problem1(i+number_of_obj0),factor(i+1),user_n) 
                   fg2_y_all(i+1) = g2(y,problem2(i+number_of_obj0),factor(i+1),user_n) 
               END DO              
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i_sep(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO     

  
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************          
               
               IF (iter_counter > 100) THEN 
                   eps = (10.0_dp)**(-6)
               END IF
                 
               DO i = 1, 1                       ! the values of subgradients for DC components f1_i and f2_i at x_0
                  new_grad1 = subgradient_f1(y, problem1(ind_f),factor(i),user_n)
                  new_grad2 = subgradient_f2(y, problem2(ind_f),factor(i),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i,j) = new_grad1(j)
                     fg2_grad_all(i,j) = new_grad2(j)
                  END DO
               END DO   
               
               DO i = 1,number_of_const  ! the values of subgradients for DC components g1_l and g2_l at x_0
                  new_grad1 = subgradient_g1(y, problem1(i+number_of_obj0),factor(i+1),user_n)
                  new_grad2 = subgradient_g2(y, problem2(i+number_of_obj0),factor(i+1),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i+1,j) = new_grad1(j)
                     fg2_grad_all(i+1,j) = new_grad2(j)
                  END DO
               END DO                
               subgrad_counter = subgrad_counter+1            ! one more subgradient calculated for DC components 

               new_grad1 = subgradient_H1_sep(fg1_grad_all, fg2_grad_all, AB_y_all, user_n) 
               new_grad2 = subgradient_H2_sep(fg2_grad_all, user_n)
               
               ! new subgradient for H at x 
               new_grad = new_grad1 - new_grad2
                              
               IF (N == 0) THEN 
                  CALL add_first_element_b1(B, new_grad,0)  
                  N = N + 1
               ELSE
                  CALL add_element_b1(B, new_grad, 0.0_dp, 0)
               END IF 
           
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 1: NEW SUBGRADIENT '
                   WRITE(*,*)''
                   WRITE(*,*)' Iteration round: ', iter_counter
                   WRITE(*,*)' New subgardient = [', new_grad, ']' 

               END IF
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CALRKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
                              
               CALL quadratic_solver(x_k, user_n, N, B, 1.0_dp, u_k, obj_u) 
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
 
 
               norm_u = SQRT( DOT_PRODUCT(u_k,u_k))

               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 2: CLARKE STATIONARITY '
                   WRITE(*,*)''
                   WRITE(*,*)' The value of minimum norm = ', norm_u
                   WRITE(*,*)' The solution to quadratic norm minimization problem  = [',u_k, ']'               
               END IF              
               
               IF (norm_u < crit_tol) THEN
                  
                  reason_for_stop = 1 
                  stop_alg = .TRUE.
                  
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' CLARKE STATIONARITY ACHIEVED'
                   WRITE(*,*)''                
               END IF                     
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SERACH DIRECTION *****************************            
               
               d =  u_k / norm_u
               y = x_k + eps * d        

               
               DO i = 1, 1                  ! the values of DC components f1_i anf f2_i at y
                   fg1_y_all(i) = f1(y,problem1(ind_f),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(ind_f),factor(i),user_n) 
               END DO      
               
               DO i = 1,number_of_const    ! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i+1) = g1(y,problem1(i+number_of_obj0),factor(i+1),user_n) 
                   fg2_y_all(i+1) = g2(y,problem2(i+number_of_obj0),factor(i+1),user_n) 
               END DO                  
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_const+1                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i_sep(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO   
                              
               H1_y = H1_sep(AB_y_all) 
               H2_y = H2_sep(fg2_y_all) 
               
               real_decrease = H1_y - H2_y - H_k

               direc_der = div_eps * real_decrease

            
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 3: SEARCH DIRECTION '
                   WRITE(*,*)''               
                   WRITE(*,*)' The search direction d = [', d, ']' 
                   WRITE(*,*)''                
               END IF           
               
               IF (direc_der <= (-m * norm_u)) THEN
                     reason_for_stop = 10                 ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.
                     
                     IF(iprint==5) THEN  !prints everything
                         WRITE(*,*)''
                         WRITE(*,*)' SUFFICIENT DESCENT IS OBTAINED.'
                         WRITE(*,*)''                  
                     END IF                      
                     
               END IF
               
               iter_counter = iter_counter +1 
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1) .AND. stop_alg ) THEN 
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  upper = 1.0_dp
                  lower = step_tol
                  stepsize = step_tol
                  stop_step_det = .FALSE.
                  ylos = .FALSE.
                
                  lower_y = x_k + stepsize * d
                  DO i = 1, 1                  ! the values of DC components f1_i and f2_i at lower_y
                      fg1_lower_all(i) = f1(lower_y,problem1(ind_f),factor(i),user_n) 
                      fg2_lower_all(i) = f2(lower_y,problem2(ind_f),factor(i),user_n) 
                  END DO  
                  DO i = 1,number_of_const    ! the values of DC components g1_l and g2_l at lower_y
                      fg1_lower_all(i+1) = g1(lower_y,problem1(i+number_of_obj0),factor(i+1),user_n) 
                      fg2_lower_all(i+1) = g2(lower_y,problem2(i+number_of_obj0),factor(i+1),user_n) 
                  END DO                      
                  f_counter = f_counter+1             
                               
                  DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at lower_y
                      AB_lower_all(i) = AorB_i_sep(fg1_lower_all, fg2_lower_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                  END DO                          
                      
                 lower_H1 = H1_sep(AB_lower_all)     ! the value of H_1 at test_y
                 lower_H2 = H2_sep(fg2_lower_all)    ! the value of H_2 at test_y 

                 descent_app = -m * stepsize * norm_u       
                 
                 IF ((lower_H1 - lower_H2 - H_k) > descent_app ) THEN 
                    stop_step_det = .TRUE.
                    stepsize = -1.0_dp                
                 END IF
                 l = 1  
     
                  
                 DO WHILE((.NOT. stop_step_det))
                  
                      IF (l == 1) THEN
                        stepsize = 100.0_dp * stepsize
                      ELSE IF (l == 2) THEN 
                        IF (upper < 0.02_dp) THEN 
                          stepsize = 0.1_dp * stepsize
                          !stepsize = 0.5_dp*(upper +lower)
                        ELSE
                          stepsize = 10.0_dp*stepsize
                        END IF
                      ELSE IF (l > 2) THEN
                        stepsize = 0.5_dp*(upper +lower)
                      END IF 
                      
                      test_y = x_k + stepsize * d
                      DO i = 1, 1                  ! the values of DC components f1_i and f2_i at test_y
                         fg1_test_all(i) = f1(test_y,problem1(ind_f),factor(i),user_n) 
                         fg2_test_all(i) = f2(test_y,problem2(ind_f),factor(i),user_n) 
                      END DO  
                      DO i = 1, number_of_const    ! the values of DC components g1_l and g2_l at test_y
                         fg1_test_all(i+1) = g1(test_y,problem1(i+number_of_obj0),factor(i+1),user_n) 
                         fg2_test_all(i+1) = g2(test_y,problem2(i+number_of_obj0),factor(i+1),user_n) 
                      END DO                      
                      f_counter = f_counter + 1             ! one more objective function value calculated for a DC component          
                               
                      DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at test_y
                        AB_test_all(i) = AorB_i_sep(fg1_test_all, fg2_test_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                      END DO                          
                      
                      test_H1 = H1_sep(AB_test_all)     ! the value of H_1 at test_y
                      test_H2 = H2_sep(fg2_test_all)    ! the value of H_2 at test_y 

                      IF (l < 7) THEN
                        descent_app = -m * stepsize * norm_u  
                        IF ( (test_H1 - test_H2 - H_k) < descent_app) THEN
                      
                         fg1_lower_all = fg1_test_all
                         fg2_lower_all = fg2_test_all
                         lower_H1 = test_H1
                         lower_H2 = test_H2
                         lower_y = test_y
                         
                         l = l + 1 
                         lower = stepsize
                        ELSE 
                         upper = stepsize
                         l = l+1
                        END IF
                      ELSE  

                         stop_step_det = .TRUE. 
                         y = lower_y
                         H1_y = lower_H1
                         H2_y = lower_H2
                         fg1_y_all = fg1_lower_all
                         fg2_y_all = fg2_lower_all

                      END IF

                      
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-
               
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'-------------------------------------------------------------------------------'               
                     WRITE(*,*)''
                     WRITE(*,*)' STEP 4: STEP-LENGTH '
                     WRITE(*,*)''                 
                     WRITE(*,*)' The step-lenght t =', stepsize
                  END IF                   
               
       
                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                 
                     x_new = y                 ! the new iteration point is the auxilary point y
                     H1_new = H1_y             ! the value of H_1 at x_new
                     H2_new = H2_y             ! the value of H_2 at x_new
                     fg1_new_all = fg1_y_all   ! the values of f1_i and g1_l at x_new
                     fg2_new_all = fg2_y_all   ! the values of f2_i and g2_l at x_new
                     
                     reason_for_stop = 0     ! the reason for stopping is the new iteration point   
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' New itaration point is found and x_new = [', x_new, ']' 
                        WRITE(*,*)' New objective function value f_new = ', H1_new -H2_new                 
                        WRITE(*,*)''                   
                     END IF                      
                     
                  ELSE
                     reason_for_stop = 2
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' APPROXIMATE STOPPING CONDITION HOLD at the current iteration point.'               
                        WRITE(*,*)''                   
                     END IF
                     
                  END IF 
               
                  iter_counter = iter_counter -1 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                  
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'--------------------------------------&
                     &-----------------------------------------'                       
                     WRITE(*,*)''                                     
                     WRITE(*,*)' THE MAXIMUM NUMBER OF ROUNDS HAS BEEN EXECUTED.'              
                     WRITE(*,*)''                  
                  END IF                  
              END IF           
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         

              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '                    CLARKE STATIONARY ALGORITHM ENDS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF               
           
           END SUBROUTINE escape_procedure_d_const
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|  
  

       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             Escape procedure                                   |  |
        !  |                      guaranteeing Clarke stationarity                          |  |
        !  |                         for the improvement function                           |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE escape_procedure_H( x_k, H1_k, H2_k, fg1_k_all, fg2_k_all, &
                            & reason_for_stop, x_new , H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            & d_cause, crit_tol, m, step_tol, iprint, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, factor, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)
                                                
            !
            ! Executes the 'Escape procedure'. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary then the algorithm generates a descent direction.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'H1_k' and 'H2_k'   : The values of DC components H_1 and H_2 at the current iteration point
            !  
            !        * 'fg1_k_all'         : The values of DC components f1_i and g1_l at the current iteration point
            !        * 'fg2_k_all'         : The values of DC components f2_i and g2_l at the current iteration point
            !
            !        * 'factor'            : The values of scaling factors for f_i and g_l
            !        * 'd_cause'           : The search direction which caused the execution of the 'Escape procedure'
            !        * 'user_n'            : The dimension of the problem
            !
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance  
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Escape procedure'
            !        * 'iprint'            : Specifies print
            !
            !        * 'number_of_obj'     : The number of objective functions        
            !        * 'number_of_const'   : The number of constants functions        
            !        * 'number_of_func'    : The sum of objective and constants functions  
            
            !
            ! OUTPUT: * 'x_new'                 : The new iteration point (if it is found, i.e., 'reason_for_stop' = 0)
            !         * 'H1_new' and 'H2_new'   : The values of DC components H_1 and H_2 at the new iteration point
            !         * 'fg1_new_all'           : The new value of DC components f1_i and g1_l (if new iteration point is found, i.e., 'reason_for_stop' = 0)
            !         * 'fg2_new_all'           : The new value of DC components f2_i and g2_l (if new iteration point is found, i.e., 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'       : Indicates the cause of the termination in the 'Escape procedure'
            !
            !         * 'iter_counter'          : the number of rounds needed in the 'Escape procedure'
            !         * 'f_counter'             : the number of function values evaluated for f in the 'Escape procedure' (same for f_1 and f_2)
            !         * 'subgrad_counter'       : the number of subgradients calculated for f in the 'Escape procedure' (same for f_1 and f_2)
            !
            ! NOTICE: * The dimensions of vectors 'x_k', 'x_new' and 'd_cause' has to be same. (i.e. the dimension of 'x_k', 'x_new' and 'd_cause' has to be
            !           'user_n' when SUBROUTINE escape_procedure_H is used in SUBROUTINE new_algorithm.
            !         * The dimension of vectors 'fg1_k_all', 'fg2_k_all', 'fg1_new_all' and 'fg2_new_all' has to be 'number_of_func'
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_k_all    ! the values of DC components f1_i and g1_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_k_all    ! the values of DC components f2_i and g2_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg1_new_all ! the values of DC components f1_i and g1_l at x_new
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg2_new_all ! the values of DC components f2_i and g2_l at x_new

               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: factor       ! the values of scaling factors for f_i and g_l           
               
               REAL(KIND=dp), INTENT(IN)  :: H1_k, H2_k      ! the value of H_1 and H_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: H1_new, H2_new  ! the value of H_1 and H_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance (needed in stopping condition)
               REAL(KIND=dp), INTENT(IN) :: m             ! a descent parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! a step-length tolerance
               
               INTEGER, INTENT(IN) :: user_n              ! the dimension of the problem
               
               INTEGER, INTENT(IN) :: number_of_obj      ! The number of objectives
               INTEGER, INTENT(IN) :: number_of_const    ! The number of contraints 
               INTEGER, INTENT(IN) :: number_of_func     ! The sum of objectives and constraints          
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm
 
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (approximate Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < beta_0)
                                                          ! 5 - the maximum number of rounds executed in the 'Escape procedure'
                                                          
               INTEGER, INTENT(OUT) :: iter_counter       ! the number of rounds used in the 'Escape procedure'
               INTEGER, INTENT(OUT) :: f_counter          ! the number of function values evaluated for H in the 'Escape procedure' (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for H in the 'Escape procedure'
               
               INTEGER, INTENT(IN) :: iprint    ! the variable that specifies print:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result (nothing is printed in main iteration)
                                                !   iprint = -1: basic print of final result without the solution vector (nothing is printed in main iteration)                                             
                                                !   iprint = 2 : extended print of final result (nothing is printed in main iteration)
                                                !   iprint = -2: extended print of final result without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results (nothing is printed in main iteration)
                                                !   iprint = -3: basic print of intermediate results and extended print of final results without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 4 : extended print of intermediate and final results (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = -4: extended print of intermediate and final results without the solution vector (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = 5 : everything is printed step by step (this is the only option which causes print in the 'main iteration') 
               

           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                         ! the bundle B of H containing subgradients from \partial H(x)
               
               REAL(KIND=dp), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad               ! the subgradient of H at x
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of H_1 and H_2 at x_d
               REAL(KIND=dp), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! the new auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: test_y                 ! the test auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: lower_y                ! the 'hepl' auxiliary point            
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_y_all      ! the values of DC components f1_i and g1_l at y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_y_all      ! the values of DC components f2_i and g2_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_y_all       ! the values of components A_i and B_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_test_all     ! the values of DC components f1_i and g1_l at test_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_test_all     ! the values of DC components f2_i and g2_l at test_y   
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_test_all      ! the values of components A_i and B_l at test_y              
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_lower_all    ! the values of DC components f1_i and g1_l at lower_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_lower_all    ! the values of DC components f2_i and g2_l at lower_y  
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_lower_all     ! the values of components A_i and B_l at lower_y              

               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg1_grad_all    ! the subgradients of DC components f1_i and g1_l
               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg2_grad_all    ! the subgradients of DC components f2_i and g2_l
               
               REAL(KIND=dp) :: H_k              ! the value of the objective function H at the current iteration point 
               REAL(KIND=dp) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=dp) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=dp) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=dp) :: div_eps          ! 1.0_dp / eps     
               REAL(KIND=dp) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=dp) :: H1_y, H2_y       ! the function values of H_1 and H_2 at y
               REAL(KIND=dp) :: test_H1, test_H2 ! the function values of H_1 and H_2 at test_y
               REAL(KIND=dp) :: lower_H1, lower_H2 ! the function values of H_1 and H_2 at lower_y
               REAL(KIND=dp) :: direc_der        ! directional derivative of H 
               REAL(KIND=dp) :: stepsize         ! stepsize 
               REAL(KIND=dp) :: upper, lower     ! bounds for stepsize 
               REAL(KIND=dp) :: descent_app      ! approximation of descent in the objective function value
                              
               INTEGER :: size_b                 ! the biggest possible bundle size for B
               INTEGER :: N                      ! the current bundle size for B with the current element
               
               INTEGER :: i, j, l 
               
               LOGICAL :: ylos                   ! 
               LOGICAL :: stop_alg               ! .TRUE. if 'Escape procedure' algorithm can be stopped
               LOGICAL :: stop_step_det          ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             
              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '              ESCAPE PROCEDURE STARTS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF             
             
           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************              
                        
            ! the size of the bundle B
               SELECT CASE(user_size)
                  
                  CASE(0)
                      size_b = (user_n + 5) * 2
                      
                  CASE(1)
                      size_b = MIN((user_n+5)*2, 250)

                  CASE(2)
                      size_b = MIN((user_n+5)*2, 100)                 
                      
               END SELECT 

               
               eps = (10.0_dp)**(-5)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_dp / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               l = 0
               
               H_k = H1_k - H2_k
            
               d =  -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))             
               
            ! The bundles B is initialized
               CALL init_bundle_b1(B, size_b, user_n) 
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.
               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 0: INITIALIATION '
                   WRITE(*,*)''
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
               
               y = x_k + (eps) * d 
               
               DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at y
                   fg1_y_all(i) = f1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(i),factor(i),user_n) 
               END DO  
               
               DO i = number_of_obj+1, number_of_func   ! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i) = g1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = g2(y,problem2(i),factor(i),user_n) 
               END DO              
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO     

  
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************          
               
               IF (iter_counter > 100) THEN 
                   eps = (10.0_dp)**(-6)
               END IF
              
                   
               DO i = 1, number_of_obj                       ! the values of subgradients for DC components f1_i and f2_i at x_0
                  new_grad1 = subgradient_f1(y, problem1(i),factor(i),user_n)
                  new_grad2 = subgradient_f2(y, problem2(i),factor(i),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i,j) = new_grad1(j)
                     fg2_grad_all(i,j) = new_grad2(j)
                  END DO
               END DO   
               
               DO i = number_of_obj+1,number_of_func          ! the values of subgradients for DC components g1_l and g2_l at x_0
                  new_grad1 = subgradient_g1(y, problem1(i),factor(i),user_n)
                  new_grad2 = subgradient_g2(y, problem2(i),factor(i),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i,j) = new_grad1(j)
                     fg2_grad_all(i,j) = new_grad2(j)
                  END DO
               END DO                
               subgrad_counter = subgrad_counter+1            ! one more subgradient calculated for DC components 

               new_grad1 = subgradient_H1(fg1_grad_all, fg2_grad_all, AB_y_all, user_n) 
               new_grad2 = subgradient_H2(fg2_grad_all, user_n)
               
               ! new subgradient for H at x 
               new_grad = new_grad1 - new_grad2
               
               IF (N == 0) THEN 
                  CALL add_first_element_b1(B, new_grad, 0)  
                  N = N + 1
               ELSE
                  CALL add_element_b1(B, new_grad, 0.0_dp, 0)
               END IF 
           
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 1: NEW SUBGRADIENT '
                   WRITE(*,*)''
                   WRITE(*,*)' Iteration round: ', iter_counter
                   WRITE(*,*)' New subgardient = [', new_grad, ']' 

               END IF
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CALRKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
               
               CALL quadratic_solver(x_k, user_n, N, B, 1.0_dp, u_k, obj_u) 
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
                               
               norm_u = SQRT( DOT_PRODUCT(u_k,u_k))

               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 2: CLARKE STATIONARITY '
                   WRITE(*,*)''
                   WRITE(*,*)' The value of minimum norm = ', norm_u
                   WRITE(*,*)' The solution to quadratic norm minimization problem  = [',u_k, ']'               
               END IF              
               
               IF (norm_u < crit_tol) THEN
                  
                  reason_for_stop = 1 
                  stop_alg = .TRUE.
                  
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' CLARKE STATIONARITY ACHIEVED'
                   WRITE(*,*)''                
               END IF                     
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SERACH DIRECTION *****************************            
               
               d =  u_k / norm_u
               y = x_k + eps * d           
               
               DO i = 1, number_of_obj                  ! the values of DC components f1_i anf f2_i at y
                   fg1_y_all(i) = f1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(i),factor(i),user_n) 
               END DO      
               
               DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i) = g1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = g2(y,problem2(i),factor(i),user_n) 
               END DO                  
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO   
               
               H1_y = H1(AB_y_all) 
               H2_y = H2(fg2_y_all) 
               
               real_decrease = H1_y - H2_y - H_k

               direc_der = div_eps * real_decrease

            
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 3: SEARCH DIRECTION '
                   WRITE(*,*)''               
                   WRITE(*,*)' The search direction d = [', d, ']' 
                   WRITE(*,*)''                
               END IF           
               
               IF (direc_der <= (-m * norm_u)) THEN
                     reason_for_stop = 10                 ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.
                     
                     IF(iprint==5) THEN  !prints everything
                         WRITE(*,*)''
                         WRITE(*,*)' SUFFICIENT DESCENT IS OBTAINED.'
                         WRITE(*,*)''                  
                     END IF                      
                     
               END IF
               
               iter_counter = iter_counter +1 
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
                CALL deallocation_b1(B)
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1) .AND. stop_alg ) THEN 
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  upper = 1.0_dp
                  lower = step_tol
                  stepsize = step_tol
                  stop_step_det = .FALSE.
                  ylos = .FALSE.
                
                  lower_y = x_k + stepsize * d
                  DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at lower_y
                      fg1_lower_all(i) = f1(lower_y,problem1(i),factor(i),user_n) 
                      fg2_lower_all(i) = f2(lower_y,problem2(i),factor(i),user_n) 
                  END DO  
                  DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at lower_y
                      fg1_lower_all(i) = g1(lower_y,problem1(i),factor(i),user_n) 
                      fg2_lower_all(i) = g2(lower_y,problem2(i),factor(i),user_n) 
                  END DO                      
                  f_counter = f_counter+1             
                               
                  DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at lower_y
                      AB_lower_all(i) = AorB_i(fg1_lower_all, fg2_lower_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                  END DO                          
                      
                 lower_H1 = H1(AB_lower_all)     ! the value of H_1 at test_y
                 lower_H2 = H2(fg2_lower_all)    ! the value of H_2 at test_y 

                 descent_app = -m * stepsize * norm_u       
                 
                 IF ((lower_H1 - lower_H2 - H_k) > descent_app ) THEN 
                    stop_step_det = .TRUE.
                    stepsize = -1.0_dp                
                 END IF
                 l = 1  
     
                  
                 DO WHILE((.NOT. stop_step_det))
                  
                      IF (l == 1) THEN
                        stepsize = 100.0_dp * stepsize
                      ELSE IF (l == 2) THEN 
                        IF (upper < 0.02_dp) THEN 
                          stepsize = 0.1_dp * stepsize
                          !stepsize = 0.5_dp*(upper +lower)
                        ELSE
                          stepsize = 10.0_dp*stepsize
                        END IF
                      ELSE IF (l > 2) THEN
                        stepsize = 0.5_dp*(upper +lower)
                      END IF 
                      
                      test_y = x_k + stepsize * d
                      DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at test_y
                         fg1_test_all(i) = f1(test_y,problem1(i),factor(i),user_n) 
                         fg2_test_all(i) = f2(test_y,problem2(i),factor(i),user_n) 
                      END DO  
                      DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at test_y
                         fg1_test_all(i) = g1(test_y,problem1(i),factor(i),user_n) 
                         fg2_test_all(i) = g2(test_y,problem2(i),factor(i),user_n) 
                      END DO                      
                      f_counter = f_counter + 1             ! one more objective function value calculated for a DC component          
                               
                      DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at test_y
                        AB_test_all(i) = AorB_i(fg1_test_all, fg2_test_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                      END DO                          
                      
                      test_H1 = H1(AB_test_all)     ! the value of H_1 at test_y
                      test_H2 = H2(fg2_test_all)    ! the value of H_2 at test_y 

                      IF (l < 7) THEN
                        descent_app = -m * stepsize * norm_u  
                        IF ( (test_H1 - test_H2 - H_k) < descent_app) THEN
                      
                         fg1_lower_all = fg1_test_all
                         fg2_lower_all = fg2_test_all
                         lower_H1 = test_H1
                         lower_H2 = test_H2
                         lower_y = test_y
                         
                         l = l + 1 
                         lower = stepsize
                        ELSE 
                         upper = stepsize
                         l = l+1
                        END IF
                      ELSE  

                         stop_step_det = .TRUE. 
                         y = lower_y
                         H1_y = lower_H1
                         H2_y = lower_H2
                         fg1_y_all = fg1_lower_all
                         fg2_y_all = fg2_lower_all

                      END IF

                      
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-
               
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'-------------------------------------------------------------------------------'               
                     WRITE(*,*)''
                     WRITE(*,*)' STEP 4: STEP-LENGTH '
                     WRITE(*,*)''                 
                     WRITE(*,*)' The step-lenght t =', stepsize
                  END IF                   
               
       
                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                 
                     x_new = y                 ! the new iteration point is the auxilary point y
                     H1_new = H1_y             ! the value of H_1 at x_new
                     H2_new = H2_y             ! the value of H_2 at x_new
                     fg1_new_all = fg1_y_all   ! the values of f1_i and g1_l at x_new
                     fg2_new_all = fg2_y_all   ! the values of f2_i and g2_l at x_new
                     
                     reason_for_stop = 0     ! the reason for stopping is the new iteration point   
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' New itaration point is found and x_new = [', x_new, ']' 
                        WRITE(*,*)' New objective function value f_new = ', H1_new -H2_new                 
                        WRITE(*,*)''                   
                     END IF                      
                     
                  ELSE
                     reason_for_stop = 2
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' APPROXIMATE STOPPING CONDITION HOLD at the current iteration point.'               
                        WRITE(*,*)''                   
                     END IF
                     
                  END IF 
               
                  iter_counter = iter_counter -1 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                  
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'--------------------------------------&
                     &-----------------------------------------'                       
                     WRITE(*,*)''                                     
                     WRITE(*,*)' THE MAXIMUM NUMBER OF ROUNDS HAS BEEN EXECUTED.'              
                     WRITE(*,*)''                  
                  END IF                  
              END IF           
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         

              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '                    ESCAPE PROCEDURE ENDS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF               
           
           END SUBROUTINE escape_procedure_H
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|    


        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                          QUADRATIC NORM SOLVER                                 |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE quadratic_solver(x, NF, NA, B, t, d, obj) 
               ! Solves the quadratic norm problem in the 'Escape procedure'.
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B                   ! the bundle B
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: d       ! the vector giving minimum norm
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! the proximity parameter t
               REAL(KIND=dp), INTENT(OUT) :: obj                   ! the optimal value of the quadratic norm problem
               
               INTEGER, INTENT(IN)  :: NF    ! the number of variables
               INTEGER, INTENT(IN)  :: NA    ! the bundle size of B is 'give_size_b1(B) + 1' 
                
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  

               REAL(KIND=dp) :: u       ! 'help' variable
               
               INTEGER :: j, k, l
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               !INTEGER :: IDECF=10                     ! IDECF=10 diagonal matrix
               !INTEGER :: MFP=2                        ! MFP=2 optimum feasible point
               INTEGER :: KBC=0, KBF=0                  ! KBC=0 no linear constraints; KBF=0 no simple bounds
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b           ! subgradient matrix of B


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               grad_m_b = grad_matrix(B)       ! subgradient matrix of B 
               
               ! NO linearization errors
               AF = 0.0_dp      
  
               ! matrix containing subgradients   
               DO j = 1, NA
                   k = (j-1)*NF
                   DO l = 1, NF
                      AG(k+l) = grad_m_b(k+l) 
                   END DO
               END DO   
                        
               !Calls PLQDF1 by Ladislav Luksan
               CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                           & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                           & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                       
                           & UMAX,GMAX,N,ITERQ)

               ! the vector giving minimum norm               
               DO j = 1, NF
                   d(j) = S(j) 
               END DO
              
              ! the value of minimum norm
               obj = 0.5_dp * DOT_PRODUCT(d, d)         
    
                   
               
           END SUBROUTINE quadratic_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        
        
       ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                       THE COMMON DESECNT DIRECTION                             |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE common_direction(x, NF, NA, mD, t, d, obj) 
               ! Determines the common search direction by solving the quadratic minimum norm problem 
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
               REAL(KIND=dp), DIMENSION(NF*NA), INTENT(IN) :: mD   ! the vector containing separate descent directions for objectives

               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: d       ! the vector giving the minimum norm
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! the proximity parameter t
               REAL(KIND=dp), INTENT(OUT) :: obj                   ! the optimal value of the quadratic minimum norm problem
               
               INTEGER, INTENT(IN)  :: NF    ! the number of variables (i.e. user_n)
               INTEGER, INTENT(IN)  :: NA    ! the number of separate direction used in the quadratic minimum norm problem
                
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  

               REAL(KIND=dp) :: u       ! 'help' variable
               
               INTEGER :: j
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               !INTEGER :: IDECF=10                     ! IDECF=10 diagonal matrix
               !INTEGER :: MFP=2                        ! MFP=2 optimum feasible point
               INTEGER :: KBC=0, KBF=0                  ! KBC=0 no linear constraints; KBF=0 no simple bounds
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
                              
               ! NO linearization errors
               AF = 0.0_dp      
  
               ! matrix containing separate directions  
               AG = mD  
                        
               !Calls PLQDF1 by Ladislav Luksan
               CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                           & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                           & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                       
                           & UMAX,GMAX,N,ITERQ)

               ! the vector giving minimum norm               
               DO j = 1, NF
                   d(j) = S(j) 
               END DO
              
              ! the value of minimum norm
               obj = 0.5_dp * DOT_PRODUCT(d, d)         
     
               
           END SUBROUTINE common_direction
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|            

        
     
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           SUBPROBLEM SOLVER                                    |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE subproblem_solver(x, NF, NA, B1, B2, t, subprob_counter) 
            ! Solves the subproblems in the separate search direction problem 
            ! Calls PLQDF1 by Ladislav Luksan
            !
            ! INPUT: * 'x'               : The current iteration point
            !        * 'NF' and 'NA'     : The number of variables and the size of the bundle 'B1', respectively    
            !        * 'B1' and 'B2'     : The bundle 'B1' and 'B2'
            !        * 't'               : The proximity parameter
            !
            ! OUTPUT: * 'subprob_counter'  : the number of subproblems solved
            !
            ! NOTICE: The dimensions of vectors 'x' has to be ' user_n' when SUBROUTINE subproblem_solver is used in SUBROUTINE main_iteration.
              
              IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B1                  ! bundle B_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                  ! bundle B_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! current iteration point
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! proximity parameter t
               
               INTEGER, INTENT(IN)  :: NF    ! number of variables
               INTEGER, INTENT(IN)  :: NA    ! bundle size of B1         IF ( NA = give_size_b1(B1) + 2 ) THEN aggregation is used
               
               INTEGER, INTENT(OUT) :: subprob_counter    ! number of subproblems solved    
                       
              
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative Lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  
               
               REAL(KIND=dp) :: alpha_b2     ! a linearization error of B_2
               REAL(KIND=dp) :: u, obj, a    ! help variables
               
               INTEGER :: i, j, k, l         ! help variables
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               INTEGER :: IDECF=10, KBC=0, KBF=0, MFP=2 ! IDECF=10 diagonal matrix; KBC=0 no linear constraints; KBF=0 no simple bounds; MFP=2 optimum feasible point
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indexes of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF) :: direction       ! Actual direction vector used (notice dimension)
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b1           ! subgradient matrix of B_1
               REAL(KIND=dp), DIMENSION(NA) :: alpha_m_b1             ! linearization error matrix of B_1
               REAL(KIND=dp), DIMENSION(NF) :: grad_b2                ! a subgradient of B_2


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               IF ( (give_size_b1(B1)+2) == NA) THEN     ! if this is TRUE then aggregation is in use
                   grad_m_b1 = grad_matrix_agg(B1)       ! subgradient matrix of B_1 with aggregation
                   alpha_m_b1 = lin_error_matrix_agg(B1) ! linearization error matrix of B_1 with aggregation
               ELSE
                   grad_m_b1 = grad_matrix(B1)       ! subgradient matrix of B_1
                   alpha_m_b1 = lin_error_matrix(B1) ! linearization error matrix of B_1
               END IF
               
               subprob_counter = give_size_b2(B2) 
                    
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                    
                        
                        grad_b2  = give_subgrad_b2(B2, i)  ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i)   ! linearization error of B_2 in the subproblem i

                        
                        DO j = 1, NA
                           k = (j-1)*NF
                           DO l = 1, NF
                               AG(k+l) = grad_m_b1(k+l) - grad_b2(l)
                           END DO
                           AF(j) = - alpha_m_b1(j) + alpha_b2 
                        END DO
                        
                        !Calls PLQDF1 by Ladislav Luksan
                        CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,MFP,KBF, &
                              & KBC,IDECF,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
                              & UMAX,GMAX,N,ITERQ)

                              
                        DO j = 1, NF
                           direction(j) = S(j) 
                        END DO
              
                        a = DOT_PRODUCT(direction, direction)           

                        obj =   XNORM + (a * u) / 2
                        CALL add_solution(B2, i , direction, XNORM, obj )   
                    
                    END DO subproblems1
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      

               subprob_counter = subprob_counter + 1 
               
               
           END SUBROUTINE subproblem_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        
        
        
          
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |           SUBPROBLEM SOLVER WITH SOPHISTICATED CONSTRAINT HANDLING             |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE subproblem_solver_const(x, AB_k_all, NF, NA, B1, B2, t, subprob_counter, number_of_func) 
               ! Solves the subproblems in the separate search direction problem when sophisticated constraint handling is used.
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B1                  ! bundle B_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                  ! bundle B_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! current iteration point
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AB_k_all  ! the values of components A_i and B_l at the current iteration point
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! proximity parameter t
               
               INTEGER, INTENT(IN)  :: NF    ! number of variables
               INTEGER, INTENT(IN)  :: NA    ! bundle size of B1  
               
               INTEGER, INTENT(IN)  :: number_of_func    ! bundle size of B1         
               
               INTEGER, INTENT(OUT) :: subprob_counter    ! number of subproblems solved    
                       
              
            !*************************** LOCAL VARIABLES ************************************
            
               
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  
               
               REAL(KIND=dp) :: alpha_b2     ! a linearization error of B_2
               REAL(KIND=dp) :: u, obj       ! 'help' variable
               
               REAL(KIND=dp) :: a
               INTEGER :: i, j, k, l
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               !INTEGER :: IDECF=10                     ! IDECF=10 diagonal matrix;
               !INTEGER :: MFP=2                        ! MFP=2 optimum feasible point
               INTEGER :: KBC=0, KBF=0                  ! KBC=0 no linear constraints; KBF=0 no simple bounds
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF) :: direction       ! Actual direction vector used (notice dimension)
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b1           ! subgradient matrix of B_1
               REAL(KIND=dp), DIMENSION(NA) :: alpha_m_b1             ! linearization error matrix of B_1
               REAL(KIND=dp), DIMENSION(NF) :: grad_b2                ! a subgradient of B_2


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               grad_m_b1 = grad_matrix(B1)                     ! subgradient matrix of B_1
               alpha_m_b1 = lin_err_and_f1_matrix(B1,AB_k_all,number_of_func) ! linearization error matrix of B_1 with component values A_i(x_h) and B_l(x_h)
               
               subprob_counter = give_size_b2(B2) 
              

               !$OMP PARALLEL DO PRIVATE(grad_b2,alpha_b2,direction,a,obj) & 
               !$OMP FIRSTPRIVATE(NA,X,IX,XL,XU,AFD,IA,IAA) &
               !$OMP FIRSTPRIVATE(AR,AZ,CF,IC,CL,CU,CG,G,H,S,KBF) &
               !$OMP FIRSTPRIVATE(KBC,XNORM,UMAX,GMAX,N,ITERQ) &
               !$OMP PRIVATE(i, AG, AF ) &
               !$OMP SHARED(grad_m_b1,alpha_m_b1,B2)               
                       
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                    

                        grad_b2  = give_subgrad_b2(B2, i)  ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i)   ! linearization error of B_2 in the subproblem i

                        
                        !$OMP CRITICAL  
                        DO j = 1, NA
                           k = (j-1)*NF
                           DO l = 1, NF
                               AG(k+l) = grad_m_b1(k+l) - grad_b2(l)
                           END DO
                           AF(j) = - alpha_m_b1(j) + alpha_b2 
                        END DO
                        !$OMP END CRITICAL  
                        
                        !Calls PLQDF1 by Ladislav Luksan
                        CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                              & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
                              & UMAX,GMAX,N,ITERQ)

                              
                        DO j = 1, NF
                           direction(j) = S(j) 
                        END DO
              
                        a = DOT_PRODUCT(direction, direction)           

                        obj =   XNORM + (a * u) / 2
                        !$OMP CRITICAL
                        CALL add_solution(B2, i , direction, XNORM, obj )   
                        !$OMP END CRITICAL
                    
                    END DO subproblems1
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      
               !$OPM END PARALLEL DO

               subprob_counter = subprob_counter + 1 
               
               
           END SUBROUTINE subproblem_solver_const
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        
       

      END MODULE imbdc