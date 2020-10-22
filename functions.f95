        MODULE functions    
      
        USE constants, ONLY : dp   ! double precision (i.e. accuracy)    
        IMPLICIT NONE
        
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                 INFORMATION SUPPLIED BY THE USER:                        | | 
        !| |                                                                          | |
        !| |                                                                          | |
        !| |    * PROBLEM specification is done in timbdc.f95 by giving, e.g.         | |
        !| |                                                                          | |
        !| |        - the number of objectives:    'number_of_obj'                    | |       
        !| |        - the number of constraints:   'number_of_const'                  | |
        !| |        - the number of functions:     'number_of_func'                   | |
        !| |        - the number of variables:     'user_n'                           | |       
        !| |                                                                          | |       
        !| |        - The DC components f1 and g1: 'problem10'                        | |
        !| |        - The DC components f2 and g2: 'problem20'                        | |       
        !| |                                                                          | | 
        !| |     AND allocated here using                                             | |
        !| |                                                                          | |
        !| |     allocate_prob_data(number_of_obj, number_of_const, number_of_func, & | |
        !| |                               & problem10, problem20)                    | |
        !| |                                                                          | |
        !| |     DEALLOCATION of problem data is done with                            | |
        !| |                                                                          | |
        !| |     deallocate_prob_data()                                               | |      
        !| |                                                                          | |
        !| |                                                                          | |
        !| |    * Different PARAMETERS:                                               | |
        !| |                                                                          | |           
        !| |        - the number of variables:     'user_n'                           | |       
        !| |        - the stopping tolerance:      'user_crit_tol_H'                  | |
        !| |        - the stopping tolerance:      'user_crit_tol_d'                  | |
        !| |        - the extra descent parameter  'user_m_extra'                     | |               
        !| |                                                                          | |       
        !| |        MAIN ITERATION:                                                   | |  
        !| |                                                                          | |               
        !| |        - the size of bundle B_1:      'user_size_b1'                     | |       
        !| |        - the size of bundle B_2:      'user_size_b2'                     | |       
        !| |        - the descent parameter:       'user_m'                           | |       
        !| |        - the decrease parameter:      'user_r_dec'                       | |       
        !| |        - the increase parameter:      'user_r_inc'                       | |           
        !| |        - the small number:            'user_eps'                         | |       
        !| |        - the decrease parameter:      'user_c1'                          | |       
        !| |        - the decrease parameter:      'user_c2'                          | |       
        !| |        - the decrease parameter:      'user_c3'                          | |       
        !| |                                                                          | |       
        !| |        CLARKE STATIONARY ALGORITHM:                                      | | 
        !| |                                                                          | |               
        !| |        - the size of bundle B:        'user_size'                        | |       
        !| |        - the descent parameter:       'user_m_clarke'                    | |       
        !| |        - the step-length tolerance:   'user_step_tol'                    | |            
        !| |                                                                          | |       
        !| |                                                                          | |       
        !| |--------------------------------------------------------------------------| |      
        !| |                       .**..**..**..**..**..**.                           | |       
        !| |                                                                          | |       
        !| |            NOTE THAT EACH PARAMETER ONLY CONTAINS INTEGER VALUE.         | |       
        !| |            ACTUAL VALUES FOR THESE PARAMTERS ARE SET AT THE              | |       
        !| |            BEGINNING OF imbdc.f95 BASED ON THE INTEGER.                  | |       
        !| |                                                                          | |       
        !| |                       .**..**..**..**..**..**.                           | |       
        !| |--------------------------------------------------------------------------| |      
        !| |                                                                          | |       
        !| |                                                                          | |       
        !| |                       .**..**..**..**..**..**.                           | |       
        !| |                              OBJECTIVES:                                 | |
        !| |                       .**..**..**..**..**..**.                           | |               
        !| |                                                                          | |       
        !| |    * Computation of the values of the DC functions f_1 and f_2:          | |
        !| |      - f1(y, problem, user_n)    value of DC component f_1 at a point y  | |
        !| |      - f2(y, problem, user_n)    value of DC component f_2 at a point y  | |           
        !| |                                                                          | |               
        !| |                                                                          | |               
        !| |    * Computation of the subgradient of the DC components f_1 and f_2:    | |
        !| |      - subgradient_f1(y, problem, user_n)    subgradient of f_1 at y     | |
        !| |      - subgradient_f2(y, problem, user_n)    subgradient of f_2 at y     | |       
        !| |                                                                          | |
        !| |                                                                          | |
        !| |                       .**..**..**..**..**..**.                           | |       
        !| |                              CONSTRAINTS:                                | |
        !| |                       .**..**..**..**..**..**.                           | |               
        !| |                                                                          | |       
        !| |    * Computation of the values of the DC functions g_1 and g_2:          | |
        !| |      - g1(y, problem, user_n)    value of DC component g_1 at a point y  | |
        !| |      - g2(y, problem, user_n)    value of DC component g_2 at a point y  | |           
        !| |                                                                          | |               
        !| |                                                                          | |               
        !| |    * Computation of the subgradient of the DC components g_1 and g_2:    | |
        !| |      - subgradient_g1(y, problem1, user_n)   subgradient of g_1 at y     | |
        !| |      - subgradient_g2(y, problem2, user_n)   subgradient of g_2 at y     | |       
        !| |                                                                          | |
        !| |                                                                          | |       
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                          | |
        !| |                                                                                          | |
        !| |                           OTHER FUNCTIONS:                                               | |
        !| |                                                                                          | |
        !| |    The value of H_1:           H1(AandB_y_all)                                           | |
        !| |    The value of H_2:           H2(fg2_y_all)                                             | |
        !| |                                                                                          | |
        !| |    The subgradient of H_1:     subgradient_H1(grad1_y_all, grad2_y_all, AandB_y_all)     | |
        !| |    The subgradient of H_2:     subgradient_H2(grad2_y_all)                               | |
        !| |                                                                                          | |
        !| |    Index of max comp. in H_1:  index_H1(AandB_y_all)                                     | |
        !| |                                                                                          | |
        !| |    The value of A_i or B_l:    AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all, ind)   | |
        !| |    The subgradient of A_i      subgradient_AorB_i(grad1_y_all, grad2_y_all, ind)         | |
        !| |                                                                                          | |
        !| |                                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
 
        
        !--------------------------------------------------------------------------------
        ! ----------------------------------------------------------------------------- |
        ! |                  INFORMATION ABOUT PARAMETERS:                            | |
        ! ----------------------------------------------------------------------------- |
        !--------------------------------------------------------------------------------       
        
      
        !****************** GLOBAL PARAMETERS *******************************************
      
        INTEGER, SAVE :: number_of_obj0                       ! The number of objectives in the multiobjective problem
                                                                        ! IF 'number_of_obj0'=1 THEN we have a single objective case.
        INTEGER, SAVE :: number_of_const0                     ! The number of constraintss in the multiobjective problem
                                                                        ! IF 'number_of_const0'=0 THEN we have an unconstrained problem    
       
        INTEGER, SAVE :: number_of_func0                      ! the sum of objective and constraint functions         
  
        !--------------------------------------------------------------------------------------  
      
        INTEGER, SAVE :: number_of_obj_sep                       ! The number of objectives in the individual search direction problem
                                                                        ! IF 'number_of_obj0'=1 THEN we have a single objective case.
        INTEGER, SAVE :: number_of_const_sep                     ! The number of constraintss in the individual search direction problem
                                                                        ! IF 'number_of_const0'=0 THEN we have an unconstrained problem    
       
        INTEGER, SAVE :: number_of_func_sep                      ! the sum of objective and constraint functions in the individual search direction problem                                                          


        !---------------------------------------------------------------------------------------------------------------------------------
        ! DC component lists where the first 'number_of_obj0' places tell the DC components f1 for the objectives
        ! and the latter 'number_of_const0' places tell the DC components g1 for the constraints      
        !---------------------------------------------------------------------------------------------------------------------------------      
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem1               ! First is DC components f1 of objective functions used and then DC components g1 of constraint functions used
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem2               ! DC components f2 and g2 of objective and constraint functions used 
       
        !__________________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
        !********* PARAMETRES NEEDED ONLY IN SEPARATE DESCENT DIRECTION DETERMINATIONS ************
        
        INTEGER, SAVE :: user_size_b1                               ! The biggest possible size of the bundle B_1 
                                                                    ! If user_size_b1 = 0 then DEFAULT value MIN((user_n+5),50) is used 
                                                                   
        INTEGER, SAVE :: user_size_b2                               ! The biggest possible size of the bundle B_2 
                                                                    ! If user_size_b2 = 0 then DEFAULT value 3 is used  
              
        INTEGER, SAVE :: user_m                                     ! The descent parameter 'm_1'  If user_m = 0 
                                                                    !                        then DEFAULT value 0.2_dp is used

        INTEGER, SAVE :: user_r_dec                                 ! The decrease parameter 'r_1'    
        
        !If user_r_dec = 0 then DEFAULT value is used.
        !                               
        !   DEFAULT value:                          
        !     If user_n < 10:           user_r_dec = 0.75_dp    
        !     If 10 <= user_n < 300:    user_r_dec = the first two decimals of n/(n+5)
        !     If user_n >= 300:         user_r_dec = 0.99_dp
        !
        !   Some examples of the DEFAULT value of the parameter 'user_r_dec':
        !     If user_n=10:     user_r_dec = 0.66_dp                          
        !     If user_n=20:     user_r_dec = 0.80_dp                         
        !     If user_n=25:     user_r_dec = 0.83_dp                         
        !     If user_n=50:     user_r_dec = 0.90_dp                         
        !     If user_n=100:    user_r_dec = 0.95_dp                         
        !     If user_n=150:    user_r_dec = 0.96_dp     
        !     If user_n=200:    user_r_dec = 0.97_dp                      
        !     If user_n=250:    user_r_dec = 0.98_dp    
        !
               
        INTEGER, SAVE :: user_r_inc                                 ! The increase parameter  'R'   

        !If user_r_inc = 9 then DEFAULT value is used.      
        !
        !   DEFAULT VALUES
        !     if user_n < 10:             r_inc = (10.0_dp)**7    
        !     if 10 <= user_n < 100:      r_inc = (10.0_dp)**5                    
        !     if user_n == 100:           r_inc = (10.0_dp)**8                    
        !     if user_n > 100:            r_inc = (10.0_dp)**6      
                                                                 
        INTEGER, SAVE :: user_c                                     ! The decrease parameter 'r_2' If user_c = 0
                                                                    !                        then DEFAULT value 0.1_dp is used 
                                                                    
        INTEGER, SAVE :: user_tol_tau                               ! The step-size tolerance 'tol_tau'  If user_tol_tau = 0 the DEFAULT value is used 
        !
        !  DEFAULT VALUES
        !     if user_n  < 101:       user_tol_tau = (10.0_dp)**(-2)    
        !     if 101 <= user_n:       user_tol_tau = (10.0_dp)**(-3)        

        !____________________________________________________________________________________________                                                       
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETRES NEEDED ONLY IN ESCAPE PROCEDURE ******************************
        
        
        INTEGER, SAVE :: user_size                                  ! The biggest possible bundle size for the bundle used in 'Escape procedure'
                                                                    ! If user_size = 1 then DEFAULT value MIN((user_n+5)*2, 250) is used 
         
        INTEGER, SAVE :: user_escape_tol_H                          ! The tolerance 'gamma' launching escape procedure. common_obj < user_escape_tol_H 
                                                                    ! If user_escape_tol_H = 5 then DEFAULT value is used 
        !
        !  DEFAULT VALUES
        !     if user_n  <= 10:     user_escape_tol_H = (10.0_dp)**(-7)    
        !     if 10 < user_n:       user_escape_tol_H = (10.0_dp)**(-6) 
         
         
        INTEGER, SAVE :: user_m_escape_H                            ! The descent parameter 'm_2' with improvement function: If user_m_escape_H = 0 
                                                                    !                       then DEFAULT value 0.01_dp is used

        INTEGER, SAVE :: user_m_escape_d                            ! The descent parameter 'm_2' when escape procedure is entered during the individual direction determination: 
                                                                    !                       If user_m_escape_d = 0 then DEFAULT value 0.01_dp is used
                                                                    
          
        INTEGER, SAVE :: user_step_tol_H                            ! The step-length tolerance 'eps2' with improvement function: If user_set_tol_H = 5
                                                                    !                       then DEFAULT value (10.0_dp)**(-4) is used 
        !  DEFAULT VALUES
        !     if user_n  < 10:         user_step_tol_H = (10.0_dp)**(-5)    
        !     if 10 <= user_n < 50:    user_step_tol_H = (10.0_dp)**(-4) 
        !     if 50 <= user_n:         user_step_tol_H = (10.0_dp)**(-3) 
 
 
        INTEGER, SAVE :: user_step_tol_d                            ! The step-length tolerance 'eps1' when escape procedure is entered during the individual direction determination: 
                                                                    !                       If user_step_tol_d = 0 then DEFAULT value (10.0_dp)**(-6) is used 
     
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        !****************** PARAMETRES NEEDED IN STOPPING CONDITIONS ********************
              
        ! ** The stopping tolerances 'delta_2'**     
        INTEGER, SAVE :: user_crit_tol_H                       ! If user_crit_tol_H = 0
                                                               ! then DEFAULT value (10.0_dp)**(-5) is used
                                                               
        ! ** The stopping tolerances 'delta_1'**     
        INTEGER, SAVE :: user_crit_tol_d                       ! If user_crit_tol_d = 0
                                                               ! then DEFAULT value (10.0_dp)**(-5) is used

        ! ** The enlargement parameter 'theta' **
        INTEGER, SAVE :: user_eps                              ! If user_eps = 0 then DEFAULT value 0.00005_dp is used
                
       
        !________________________________________________________________________________
        !>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<
        
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |     To EXECUTE the bundle algorithm PBDC the USER needs to DETERMINE:    | | 
        !| |                                                                          | |       
        !| |        f1(y, problem1)     - value of DC component f_1 at a point y      | |
        !| |        f2(y, problem2)     - value of DC component f_2 at a point y      | |
        !| |                                                                          | |
        !| |        subgradient_f1(y, problem1)   - subgradient of f_1 at y           | |
        !| |        subgradient_f2(y, problem2)   - subgradient of f_2 at y           | |
        !| |                                                                          | |
        !| |        g1(y, problem1)     - value of DC component g_1 at a point y      | |
        !| |        g2(y, problem2)     - value of DC component g_2 at a point y      | |
        !| |                                                                          | |
        !| |        subgradient_g1(y, problem1)   - subgradient of g_1 at y           | |
        !| |        subgradient_g2(y, problem2)   - subgradient of g_2 at y           | |       
        !| |                                                                          | |
        !| |                                                                          | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
    
        
        CONTAINS
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                        SUBROUTINES:                                      | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*


       !********************************************************************************
        !                                                                               |
        !                        ALLOCATING PROBLEM DATA                                |
        !                                                                               |
        !********************************************************************************
           SUBROUTINE parameters(ind, scale_func, agg_used, t_adjust, use_reset_b1, use_reset_b2, max_null_step, &
                                   & use_const_in_sep_d) 
               !
               ! Allocates the paramters based on given 'ind' value  
               ! 
               ! 
               ! NOTICE: * 'ind' > 0
                
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               INTEGER, INTENT(IN) :: ind      ! the value defining the used parameters  
               
               INTEGER, INTENT(OUT) :: max_null_step      ! the maximum number of null steps 

               LOGICAL, INTENT(OUT) :: scale_func         ! Scaling of objectives  
               LOGICAL, INTENT(OUT) :: agg_used           ! Aggregation
               LOGICAL, INTENT(OUT) :: t_adjust           ! Adjustment of paramter t (needs to be .FALSE.)
               LOGICAL, INTENT(OUT) :: use_reset_b1       ! B1 resetation   
               LOGICAL, INTENT(OUT) :: use_reset_b2       ! B2 resetation    
               LOGICAL, INTENT(OUT) :: use_const_in_sep_d ! Constraints in the individual search direction problem         
               !**************************** OTHER VARIABLES **************************************
              
                SELECT CASE (ind)
                    
                    ! The default values
                    CASE(1)
                        
                        scale_func = .TRUE.           ! Scaling of objectives is used
                        !scale_func = .FALSE.         ! Scaling of objectives is NOT used
                         
                        agg_used = .TRUE.             ! Aggregation is used
                        
                        t_adjust = .FALSE.            ! CANNOT BE CHANGED!
                        
                        use_reset_b1 = .FALSE.        ! B1 is NOT reset after serious step
                        use_reset_b2 = .TRUE.         ! B2 is reset after serious step
                        
                        !use_const_in_sep_d = .TRUE.   ! Constraints are used when individual search directions are determined
                        use_const_in_sep_d = .FALSE. ! Constraints are NOT used when individual search directions are determined
                        
                        
                        ! Below are default integer values for different parameters
                        ! Actual values for these parameters are set in the beginning of imbdc.f95
                        ! based on the integer given.

                        max_null_step = 0             ! DEFAUL number of null_steps is 5
                        
                        user_size_b1 = 0              ! size of B1
                        user_size_b2 = 0              ! size of B2
                        user_size = 1                 ! size of B in escape procedure
                        
                        user_m = 0                    ! descent parameter
                        user_r_dec = 0                ! decrease parameter
                        user_r_inc = 9                ! increase parameter
                        user_c = 0                    ! decrease parameter
 
                        user_eps = 0                  ! enlargement parameter
                        
                        user_tol_tau = 0              ! step-size tolerance
                        
                        ! For escape procedure inside the individual descent direction determination
                        user_crit_tol_d = 0      ! criticality i.e. stopping tolerance
                        user_m_escape_d = 0      ! descent parameter
                        user_step_tol_d = 0      ! step-length tolerance
                        
                        ! For escape procedure with the improvement function H
                        user_crit_tol_H = 0      ! criticality i.e. stopping tolerance
                        user_escape_tol_H = 5    ! tolerance launching escape procedure                   
                        user_m_escape_H = 0      ! descent parameter
                        user_step_tol_H = 5      ! step-length tolerance
                       
                
                END SELECT

                
           END SUBROUTINE parameters
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                
               
        !********************************************************************************
        !                                                                               |
        !                        ALLOCATING PROBLEM DATA                                |
        !                                                                               |
        !********************************************************************************
           SUBROUTINE allocate_prob_data(number_of_obj, number_of_const, number_of_func, & 
                                        & problem10, problem20 ) 
               !
               ! Allocates the problem data 
               ! 
               ! 
               ! NOTICE: * 'number_of_obj' > 0
               !         * 'number_of_const' >= 0               
               !         * 'number_of_func' = 'number_of_obj' + 'number_of_const'   
               !         * The length of vectors 'problem10' and 'problem20' is 'number_of_func'
               ! 
               !         * 'problem10' is a list where is DC components f1 and g1 of objective and constraint functions used        
               !         * 'problem20' is a list where is DC components f2 and g2 of objective and constraint functions used                    
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               INTEGER, DIMENSION(number_of_func), INTENT(IN) :: problem10       ! DC components f1 and g1 of objective and constraint functions used  (first is objectives)
               INTEGER, DIMENSION(number_of_func), INTENT(IN) :: problem20       ! DC components f2 and g2 of objective and constraint functions used  (first is objectives)              
               
               INTEGER, INTENT(IN) :: number_of_obj       ! the number of objectives
               INTEGER, INTENT(IN) :: number_of_const     ! the number of constants
               INTEGER, INTENT(IN) :: number_of_func      ! the sum of objective and constraint functions     
               !**************************** OTHER VARIABLES **************************************
              
             ! The number of objcetive funtions              
               IF (number_of_obj>0) THEN 
                  number_of_obj0 = number_of_obj
                  number_of_obj_sep = 1
                  
               ELSE
                  WRITE(*,*) 'There is no objective or the number of objectives is negative!'
               END IF 
                
             ! The number of constraint functions  
               IF (number_of_const >= 0) THEN
                  number_of_const0 = number_of_const
                  number_of_const_sep = number_of_const
               ELSE
                  WRITE(*,*) 'Negative value for number of constraints!'
               END IF 
                
             ! The sum of objective and constraint functions
               number_of_func0 = number_of_func             
               number_of_func_sep = 1+number_of_const             
                
              ! Objective and constraint functions used
                ALLOCATE(problem1(number_of_func),problem2(number_of_func))
                problem1 = problem10
                problem2 = problem20
                            
                         
           END SUBROUTINE allocate_prob_data
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

           
        !********************************************************************************
        !                                                                               |
        !                       DEALLOCATING PROBLEM DATA                               |
        !                                                                               |
        !********************************************************************************
           SUBROUTINE deallocate_prob_data()
               !
               ! Deallocates the data matrices.              
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               !**************************** OTHER VARIABLES **************************************
               
               DEALLOCATE(problem1, problem2)    
                  
           
           END SUBROUTINE deallocate_prob_data      
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
 
        !********************************************************************************
        !                                                                               |
        !                                STARTING POINT                                 |
        !                                                                               |
        !********************************************************************************             
           SUBROUTINE starting_point(x_0,ind,user_n)
               !
               ! - Initializes the starting point used             
               ! - The starting points for the test problems presented in the article [1]
               !
               !  [1] Outi Montonen and Kaisa Joki: "Bundle-based descent method for nonsmooth       
               !      multiobjective DC optimization with inequality constraint". Journal of         
               !      Glonbal Optimization, 72 (2018), 403–429.    
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(user_n), INTENT(OUT) :: x_0    ! Starting point           
                INTEGER, INTENT(IN) :: ind                              ! the parameter defining the starting point used
                INTEGER, INTENT(IN) :: user_n                           ! the dimension of the problem 
               
               !**************************** OTHER VARIABLES **************************************
                INTEGER :: i
        
               SELECT CASE(ind)    

                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                   CASE(1)
                     x_0 = 1.0_dp
                     x_0(1) = -1.20_dp
         
                   !-------------------------------------
                   !           Problem   2
                   !-------------------------------------
                   ! x_0 = (/ -0.5_dp, 1.0_dp /)   
                   CASE(2)
                     x_0 = 1.0_dp
                     x_0(1) = -0.5_dp 
                    
                   !-------------------------------------
                   !           Problem   3
                   !-------------------------------------
                   ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                   CASE(3)
                     x_0 = 1.0_dp
                     x_0(1) = -1.20_dp
                    
                   !-------------------------------------
                   !           Problem   4
                   !-------------------------------------
                   ! x_0 = (/ -2.0_dp, 1.0_dp /)
                   CASE(4)
                     x_0 = 1.0_dp 
                     x_0(1) = -2.0_dp 
                    
                   !-------------------------------------
                   !           Problem   5
                   !-------------------------------------
                   ! x_0= (/ 4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp/)
                   CASE(5)
                     DO i = 1, user_n
                       IF ( (i == 1) .OR. (i == 3)) THEN
                          x_0(i) = 4.0_dp
                       ELSE
                          x_0(i) = 2.0_dp              
                       END IF
                     END DO              
                    
                   !-------------------------------------
                   !           Problem   6
                   !-------------------------------------
                   CASE(6)
                     DO i = 1, user_n
                       x_0(i)  = 0.1_dp * i
                     END DO    

                   !-------------------------------------
                   !           Problem   7
                   !-------------------------------------
                   CASE(7)
                     DO i = 1, user_n
                       x_0(i) = 2.0_dp * i 
                     END DO                  
                    
                   !-------------------------------------
                   !           Problem   8
                   !-------------------------------------
                   CASE(8)
                     DO i = 1, user_n
                       x_0(i)  = 0.1_dp * i
                     END DO              
                    
                   !-------------------------------------
                   !           Problem   9
                   !-------------------------------------
                   CASE(9)
                     DO i = 1, user_n
                       x_0(i) = 2.0_dp * i 
                     END DO          
                        
                   !-------------------------------------
                   !           Problem   10
                   !-------------------------------------
                   CASE(10)
                     x_0 = -1.0_dp
                     DO i = 1,user_n, 2
                      x_0(i) = 1.0_dp
                     END DO         

                   !-------------------------------------
                   !           Problem   11
                   !-------------------------------------
                   ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                   CASE(11)
                     x_0 = 1.0_dp
                     x_0(1) = -1.20_dp   

                   !-------------------------------------
                   !           Problem   12
                   !-------------------------------------
                   ! x_0 = (/ 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp /)
                   CASE(12)
                     x_0 = 3.0_dp
                     x_0(1) = 1.0_dp
                     x_0(user_n) = 1.0_dp

                   !-------------------------------------
                   !           Problem   13
                   !-------------------------------------
                   CASE(13)
                     DO i = 1, user_n
                       x_0(i) = 0.1_dp * i 
                     END DO                  
                    
                   !-------------------------------------
                   !           Problem   14
                   !-------------------------------------
                   CASE(14)
                     DO i = 1, user_n
                       x_0(i) = 0.1_dp * i 
                     END DO                  
                    
                   !-------------------------------------
                   !           Problem   15
                   !-------------------------------------
                   CASE(15)
                     DO i = 1, user_n
                       x_0(i) = 0.1_dp * i 
                     END DO                  
                    
                   !-------------------------------------
                   !           Problem   16
                   !-------------------------------------
                   ! x_0 = (/ -0.5_dp, 1.0_dp /)   
                   CASE(16)
                     x_0 = 1.0_dp
                     x_0(1) = -0.5_dp    

                   !-------------------------------------
                   !           Problem   17
                   !-------------------------------------
                   ! x_0= (/ 4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp/)
                   CASE(17)
                     DO i = 1, user_n
                       IF ( (i == 1) .OR. (i == 3)) THEN
                          x_0(i) = 4.0_dp
                       ELSE
                          x_0(i) = 2.0_dp              
                       END IF
                     END DO              

                   !-------------------------------------
                   !           Problem   18
                   !-------------------------------------
                   CASE(18)
                     DO i = 1, user_n
                       x_0(i) = 2.0_dp * i 
                     END DO    

                   !-------------------------------------
                   !           Problem   19
                   !-------------------------------------
                   ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                   CASE(19)
                     x_0 = 1.0_dp
                     x_0(1) = -1.20_dp               

                   !-------------------------------------
                   !           Problem   20
                   !-------------------------------------
                   ! x_0 = (/ 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp /)
                   CASE(20)
                     x_0 = 3.0_dp
                     x_0(1) = 1.0_dp
                     x_0(user_n) = 1.0_dp
                    
                   !-------------------------------------
                   !           Problem   21
                   !-------------------------------------
                   CASE(21)
                     DO i = 1, user_n
                       x_0(i) = 0.1_dp * i 
                     END DO        

                   !-------------------------------------
                   !           Problem   22
                   !-------------------------------------
                   ! x_0 = (/ 0.0_dp, 2.0_dp /)   
                   CASE(22)
                     x_0(1) = 0.0_dp
                     x_0(2) = 2.0_dp    
                     
                    
                  END SELECT
                  
           END SUBROUTINE starting_point   
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                        OBJECTIVE FUNCTIONS:                              | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       

        
        !********************************************************************************
        !                                                                               |
        !              FUNCTION VALUES OF THE DC COMPONENTS f_1 AND f_2                 |
        !                                                                               |
        !********************************************************************************

           FUNCTION f1(y, problem1, factor, user_n) RESULT(f)        
                !
                ! Calculates the function value of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                 ! the objective function selected for the DC component f2
                REAL(KIND=dp), INTENT(IN) :: factor             ! the scaling factor for objective
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component f_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(4) :: a                ! help variables
                REAL(KIND=dp) :: apu, largest                   ! help variables
                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp
                REAL(KIND=dp) :: a1, a2, a3, a4, a5, a6 
                INTEGER :: i, j, ind                                    ! help variable
                
                REAL(KIND=dp), DIMENSION(user_n) :: f_i    ! only first n-1 places are used            
                
                SELECT CASE(problem1)

                
                   !-------------------------------------
                   !           Problem   0
                   !-------------------------------------
                   CASE(0)
                     f = 0.0_dp

                     DO i = 1, user_n   
                       f = f + MAX( -1.0_dp * y(i), 0.5_dp* y(i), -3+2.0_dp*y(i) )
                     END DO   
                   !-------------------------------------

                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1)  
                     a1 = y(1)**4 + y(2)**2
                     a2 = (2.0_dp - y(1))**2 + (2.0_dp - y(2))**2 
                     a3 = 2.0_dp * EXP(-y(1)+y(2))

                     a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                     a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                     a5 = a5 + 4.0_dp
                     a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp               
                
                     f = MAX(a1,a2,a3) + a4 + a5 + a6                  

                     f = f * factor                  
                     
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   2
                   !-------------------------------------
                   CASE(2)  
                     IF ((ABS(y(1)) - y(2) ) > 0.0_dp) THEN
                        f = ABS(y(1) - 1.0_dp) + 200.0_dp * (ABS(y(1)) - y(2) )
                     ELSE
                        f = ABS(y(1) - 1.0_dp)
                     END IF
                     
                     f = f * factor                  
                   !-------------------------------------    
                  
                   !-------------------------------------
                   !           Problem   3
                   !-------------------------------------
                   CASE(3)
                     f = ABS(y(1) - 1.0_dp)
                     f = f + ABS(y(3) - 1.0_dp )                
                     f = f + 4.95_dp * ( ABS(y(2) + y(4) - 2.0_dp ) )
                     f = f + 10.1_dp * ( ABS(y(2)-1.0_dp) + ABS(y(4)-1.0_dp) )      
                
                     IF ((ABS(y(1)) - y(2) ) > 0.0_dp) THEN
                        f = f + 200.0_dp * (ABS(y(1)) - y(2) )
                     END IF 
                
                     IF ((ABS(y(3)) - y(4) ) > 0.0_dp) THEN
                        f = f + 180.0_dp * (ABS(y(3)) - y(4) )
                     END IF  

                     f = f * factor
                   !-------------------------------------                  
                   
                   !-------------------------------------
                   !           Problem   4
                   !-------------------------------------
                   CASE(4)
                     apu = ABS(y(1))
                     DO i = 2, user_n
                       IF (apu < ABS(y(i)) ) THEN 
                          apu = ABS(y(i))
                       END IF
                     END DO
                     f = user_n * apu 

                     f = f * factor
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   5
                   !-------------------------------------
                   CASE(5)
                     apu = ABS(summa(y,1,user_n)) 
                     DO j = 2, 20
                       IF (apu <= ABS(summa(y,j,user_n)) ) THEN 
                         apu = ABS(summa(y,j,user_n))
                       END IF
                     END DO
                     f = 20.0_dp * apu

                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   6
                   !-------------------------------------
                   CASE(6)
                     f = y(2) + 0.1_dp * (y(1)**2 + Y(2)**2)  

                     IF (-y(2) <= 0.0_dp ) THEN 
                         f = f
                     ELSE
                         f = f + 10.0_dp * (-y(2))
                     END IF
                     
                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   7
                   !-------------------------------------
                   CASE(7)
                     f = ABS(y(1)-1.0_dp) + 200.0_dp * MAX(0.0_dp, ABS(y(1))-y(2) )                 
                
                     a1 = y(1)**2 + y(2)**2 + ABS(y(2))
                     a2 = y(1) + y(1)**2 + y(2)**2 + ABS(y(2)) - 0.5_dp
                     a3 = ABS( y(1) - y(2) ) + ABS(y(2)) - 1.0_dp
                     a4 = y(1) + y(1)**2 + y(2)**2
                
                     f = f + 10.0_dp * MAX(a1,a2,a3,a4)
                     
                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   8
                   !-------------------------------------                    
                   CASE(8)
                     f = 9.0_dp - 8.0_dp * y(1) -6.0_dp * y(2) - 4.0_dp * y(3)  
                     f = f + 2.0_dp * ABS(y(1)) + 2.0_dp * ABS(y(2))+ 2.0_dp * ABS(y(3)) 
                     f = f + 4.0_dp * y(1)**2 + 2.0_dp * y(2)**2 + 2.0_dp * y(3)**2 

                     a(1) = y(1) + y(2) + 2.0_dp * y(3) - 3.0_dp
                     a(2) = -y(1)
                     a(3) = -y(2)
                     a(4) = -y(3)
                
                     f = f + 10.0_dp * MAX(0.0_dp, a(1), a(2), a(3), a(4))
                     
                     f = f * factor                  
                   !------------------------------------
                    
                   !-------------------------------------
                   !           Problem   9
                   !-------------------------------------                    
                   CASE(9)
                     f = y(1)**2 + (y(1)-1.0_dp)**2 + 2.0_dp*(y(1)-2.0_dp)**2 
                     f = f + (y(1)-3.0_dp)**2 + 2.0_dp * y(2)**2 + (y(2)-1.0_dp)**2 
                     f = f + 2.0_dp*(y(2)-2.0_dp)**2 + y(3)**2 + (y(3)-1.0_dp)**2 
                     f = f + 2.0_dp*(y(3)-2.0_dp)**2 + (y(3)-3.0_dp)**2 
                     f = f + 2.0_dp*y(4)**2 + (y(4)-1.0_dp)**2 + 2.0_dp*(y(4)-2.0_dp)**2 

                     f = f * factor                  
                   !------------------------------------
                               
                   !-------------------------------------
                   !           Problem   10
                   !-------------------------------------
                   CASE(10)
                     f = 0.0_dp
                     DO i = 1, user_n
                       f =  f + y(i)**2                             
                     END DO      

                     f = f * factor
                   !-------------------------------------

                   !-------------------------------------
                   !           Problem   11
                   !-------------------------------------
                   CASE(11)
                     f = 0.0_dp
                     f = f + 4.0_dp * ABS(y(1)) 
                     f = f + 2.0_dp * ABS(y(2)) 
                     f = f + 22.0_dp * ABS(y(3)) 
                     f = f - 33.0_dp*y(1) + 16.0_dp*y(2) - 24.0_dp*y(3)
                     f = f + 100.0_dp*MAX(0.0_dp, 2.0_dp*ABS(y(2))-3.0_dp*y(1)-7.0_dp)
                     f = f + 100.0_dp*MAX(0.0_dp, ABS(y(3))-4.0_dp*y(1)-11.0_dp)
                     
                     f = f * factor                  
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !           Problem   12
                   !-------------------------------------
                   CASE(12)                
                     f = 0.0_dp
                     
                     DO i = 1, user_n
                        f = f + ABS(y(i))
                     END DO
                     
                     DO i = 1, user_n
                        apu = 2.0_dp * (y(i)*y(i) -y(i) -1.0_dp)
                        IF (apu > 0.0_dp) THEN
                           f = f + 10.0_dp*apu
                        END IF
                     END DO
                     
                     f = f * factor
                 
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !           Problem   13
                   !-------------------------------------
                   CASE(13)             
                    
                     f = 0.0_dp
                     
                     DO i = 1, user_n-1
                        f = f + ABS(y(i)+y(i+1))
                     END DO
                     
                     DO i = 1, user_n-2
                        f = f + ABS(y(i)+y(i+2))
                     END DO
                     
                     f = f + ABS(y(1)+y(9)) 
                     f = f + ABS(y(1)+y(10)) 
                     f = f + ABS(y(2)+y(10)) 
                     f = f + ABS(y(1)+y(5)) 
                     f = f + ABS(y(4)+y(7)) 
                     
                     apu = 0.0_dp
                     DO i = 1, user_n
                        apu = apu + y(i)
                     END DO
                     
                     IF ((apu-1.0_dp)>0.0_dp) THEN 
                        f = f + 10.0_dp * (apu-1.0_dp)
                     END IF
                     
                     DO i = 1, user_n
                        IF (-y(i) > 0.0_dp) THEN 
                            f = f + 10.0_dp*(-y(i))
                        END IF
                     END DO
                    
                     f = f * factor                 
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   14
                   !-------------------------------------
                   CASE(14)         

                     f = 0.0_dp
                     
                     apu = 0.0_dp
                     DO j =1, user_n
                        apu = apu + y(j)/(1+j-1.0_dp)
                     END DO
                     largest = ABS(apu) 
                     
                     DO i = 2, user_n
                        apu = 0.0_dp
                        DO j = 1, user_n
                          apu = apu + y(j)/(i+j-1.0_dp)
                        END DO
                        IF (ABS(apu) > largest) THEN 
                           largest = ABS(apu)
                        END IF
                     END DO
                     
                     f = user_n * largest   

                     f = f * factor                  
                   
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   15
                   !-------------------------------------
                   CASE(15)     

                     f = 0.0_dp
                     
                     DO i = 1, user_n-1
                         a(1) = y(i)**4 + y(i+1)**2
                         a(2) = (2.0_dp-y(i))**2 + (2.0_dp-y(i+1))**2
                         a(3) = 2.0_dp*EXP(-y(i)+y(i+1))
                         IF (a(1)>a(2)) THEN
                            IF (a(1)>a(3)) THEN
                              f_i(i) = a(1)
                            ELSE
                              f_i(i) = a(3)
                            END IF
                         ELSE
                            IF (a(2)>a(3)) THEN
                              f_i(i) = a(2)
                            ELSE
                              f_i(i) = a(3)
                            END IF      
                         END IF
                     END DO
                     
                     largest = f_i(1)
                     ind = 1
                     DO i =2, user_n-1
                         IF (f_i(i) > largest) THEN 
                            largest = f_i(i)
                            ind = i
                         END IF
                     END DO
                     
                     f = (user_n-1.0_dp)*largest
                     
                     f = f * factor
                   
                   !-------------------------------------                      
                   
               
                   !-------------------------------------
                   !              Problem 16:
                   !          Chained Crescent I
                   !-------------------------------------
                   CASE(16)
                     f = 0.0_dp
                     
                     apu = 0.0_dp
                     DO i = 1, user_n-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) -1.0_dp
                     END DO
                     
                    apu = 2.0_dp *apu
                    IF (apu > 0.0_dp) THEN
                       f = apu
                    END IF
                    
                    f = f * factor
                   !------------------------------------- 
               
                   
                   !-------------------------------------
                   !          BAD critical point
                   !-------------------------------------
                    CASE(30)                   
                      f= MAX(y(1)**2, y(1)) 
                      
                   !-------------------------------------
                
                END SELECT
                
           END FUNCTION f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION f2(y, problem2, factor, user_n) RESULT(f)           
                !
                ! Calculates the function value of DC component f_2 at a point 'y'.
                ! Variable 'problem2' identifies the objective function used.
                !
                ! NOTICE: The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                 ! the objective function selected for the DC component f2   
                REAL(KIND=dp), INTENT(IN) :: factor             ! the scaling factor for objective 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                          ! the function value of the DC component f_2 at a point 'y'
                REAL(KIND=dp) :: apu, largest               
                REAL(KIND=dp) :: a, b              ! help variables
                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp
                REAL(KIND=dp), DIMENSION(3) :: aa
                REAL(KIND=dp), DIMENSION(user_n-1) :: f_i
                REAL(KIND=dp), DIMENSION(user_n) :: term
                REAL(KIND=dp) :: a4, a5, a6
                INTEGER :: i, j, ind                                     ! help variable
                
                
                SELECT CASE(problem2)
                
                                        
                   !-------------------------------------
                   !           Problem   0
                   !-------------------------------------
                   CASE(0)
                     f = 0.0_dp
                
                     DO i = 1, user_n
                       f = f + MAX( -0.5_dp * y(i), y(i) )
                     END DO
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------
                   CASE(1)
                     a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                     a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                     a5 = a5 + 4.0_dp
                     a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp
                
                     f = MAX( (a4 + a5) , (a4 + a6) , (a5 + a6) )                  

                     f = f * factor
                   !-------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   2
                   !-------------------------------------                  
                   CASE(2)
                     f = 100.0_dp * ( ABS(y(1)) - y(2) )
                     
                     f = f * factor                  
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Problem   3
                   !-------------------------------------
                   CASE(3)
                     f = 4.95_dp * ( ABS(y(2) - y(4)) )
                     f = f + 90.0_dp * ( ABS(y(3)) - y(4) )
                     f = f + 100.0_dp * ( ABS(y(1)) - y(2) )

                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   4
                   !-------------------------------------
                   CASE(4)
                     f = 0.0_dp
                     DO i = 1, user_n
                        f = f + ABS(y(i))                                   
                     END DO  

                     f = f * factor
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   5
                   !-------------------------------------
                   CASE(5)
                     f = 0.0_dp
                     DO j = 1, 20
                        f = f + ABS( summa(y,j,user_n) )
                     END DO 

                     f = f * factor                  
                   !-------------------------------------                  
                   
                   !-------------------------------------
                   !           Problem   6
                   !-------------------------------------
                   CASE(6)
                     f = 0.0_dp
                
                     IF (y(1) <= 0.0_dp ) THEN 
                         f = f - y(1)
                     ELSE
                         f = f + y(1) 
                     END IF
                
                     IF (y(2) <= 0.0_dp ) THEN 
                         f = f - y(2)
                     ELSE
                         f = f + y(2)
                     END IF   

                     f = f * factor                  
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Problem   7
                   !-------------------------------------
                   CASE(7)
                     f = 10.0_dp * ( y(1)**2 + y(2)**2 + ABS(y(2)) )
                     f = f + 100.0_dp * ( ABS(y(1)) - y(2) )

                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   8
                   !-------------------------------------                   
                   CASE(8)
                     f = ABS( y(1) -y(2) ) + ABS( y(1) -y(3) )
                     
                     f = f * factor                  
                   !------------------------------------
                    
                   !-------------------------------------
                   !           Problem   9
                   !-------------------------------------                    
                   CASE(9)
                     f = 0.0_dp 

                     a = (y(1) -2.0_dp)**2 + y(2)**2
                     b = (y(3) -2.0_dp)**2 + y(4)**2
                
                     IF ( a >= b) THEN 
                       f = f + a
                     ELSE 
                       f = f + b
                     END IF

                     a = (y(1) -2.0_dp)**2 + (y(2)-1.0_dp)**2
                     b = (y(3) -2.0_dp)**2 + (y(4)-1.0_dp)**2
                
                     IF ( a >= b) THEN 
                       f = f + a
                     ELSE 
                       f = f + b
                     END IF

                     a = (y(1) -3.0_dp)**2 + y(2)**2
                     b = (y(3) -3.0_dp)**2 + y(4)**2
                
                     IF ( a >= b) THEN 
                       f = f + a
                     ELSE 
                       f = f + b
                     END IF

                    a = (y(1))**2 + (y(2)-2.0_dp)**2
                     b = (y(3))**2 + (y(4)-2.0_dp)**2
                
                     IF ( a >= b) THEN 
                       f = f + a
                     ELSE 
                       f = f + b
                     END IF

                     a = (y(1)-1.0_dp)**2 + (y(2)-2.0_dp)**2
                     b = (y(3)-1.0_dp)**2 + (y(4)-2.0_dp)**2
                
                     IF ( a >= b) THEN 
                       f = f + a
                     ELSE 
                       f = f + b
                     END IF            

                     f = f * factor                  
                   !------------------------------------  
                   
                   !-------------------------------------
                   !           Problem   10
                   !-------------------------------------
                   CASE(10)
                     f = 0.0_dp
                     DO i = 2, user_n
                         f = f + ABS( y(i) - y(i-1) ) 
                     END DO
                     
                     f = f * factor
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   11
                   !-------------------------------------
                   CASE(11)
                     apu = - 7.0_dp * y(1) 
                     apu = apu + 2.0_dp * ABS(y(2)) - 18.0_dp
                     f = 20.0_dp * apu
                     
                     f = f * factor                  
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   12
                   !-------------------------------------
                   CASE(12)                
                     f = 0.0_dp

                     DO i = 1, user_n
                        apu = (y(i)*y(i) -y(i) -1.0_dp)
                        f = f + 10.0_dp * apu
                     END DO
                     
                     term = 0.0_dp
                     DO i = 1, user_n
                        DO j = 1, user_n 
                           IF (i /= j) THEN
                               term(i) = term(i) + ABS(y(j))
                           END IF
                        END DO
                     END DO
                     
                     largest = term(1)
                     ind = 1 
                     
                     DO i = 2, user_n
                        IF (term(i) > largest) THEN 
                           largest = term(i)
                           ind = i
                        END IF
                     END DO
                     
                     f = f + largest
                     
                    f = f * factor                   
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !           Problem   13
                   !-------------------------------------
                   CASE(13)             
                     f = 0.0_dp
                     
                     DO i = 1, user_n-1
                        f  = f + ABS(y(i)) + ABS(y(i+1))
                     END DO
                     
                     DO i = 1, user_n-2
                        f  = f + ABS(y(i)) + ABS(y(i+2))
                     END DO                  
                     
                     f = f + 3.0_dp*ABS(y(1))
                     f = f + ABS(y(2))
                     f = f + ABS(y(4))
                     f = f + ABS(y(5))
                     f = f + ABS(y(7))
                     f = f + ABS(y(9))
                     f = f + 2.0_dp * ABS(y(10))    
    
                     f = f * factor 
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   14
                   !-------------------------------------
                   CASE(14)         
                  
                     f = 0.0_dp
                     
                     DO i = 1, user_n
                        apu = 0.0_dp
                        DO j = 1, user_n
                          apu = apu + y(j)/(i+j-1.0_dp)
                        END DO
                        f = f + ABS(apu)
                     END DO
                     
                    f = f * factor                   
                  
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   15
                   !-------------------------------------
                   CASE(15)     
                     f = 0.0_dp
                     
                     DO i = 1, user_n-1
                         aa(1) = y(i)**4 + y(i+1)**2
                         aa(2) = (2.0_dp-y(i))**2 + (2.0_dp-y(i+1))**2
                         aa(3) = 2.0_dp*EXP(-y(i)+y(i+1))
                         IF (aa(1)>aa(2)) THEN
                            IF (aa(1)>aa(3)) THEN
                              f_i(i) = aa(1)
                            ELSE
                              f_i(i) = aa(3)
                            END IF
                         ELSE
                            IF (aa(2)>aa(3)) THEN
                              f_i(i) = aa(2)
                            ELSE
                              f_i(i) = aa(3)
                            END IF      
                         END IF
                     END DO

                     DO i =1, user_n-1
                         f = f + f_i(i)
                     END DO 

                     f = f * factor
                    
                   !-------------------------------------                  
                   
                   !-------------------------------------
                   !            Problem 16:
                   !         Chained Crescent I
                   !-------------------------------------
                   CASE(16)
                     f = 0.0_dp

                     DO i = 1, user_n-1
                        f = f + y(i)**2 + (y(i+1)-1.0_dp)**2
                        f = f + y(i+1) -1.0_dp
                     END DO
                     
                     f = f * factor                  
                   !-------------------------------------                  
  
                   !-------------------------------------
                   !          BAD critical point
                   !-------------------------------------
                    CASE(30)                   
                      f= MAX(0.5_dp* y(1)**2, -y(1)) 
                      
                   !-------------------------------------               
                
                END SELECT              
            


           END FUNCTION f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           

        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS f_1 AND f_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_f1(y, problem1, factor, user_n) RESULT(grad)           
                !
                ! Calculates a subgradient of the DC component f_1 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_1 is calculated
                INTEGER, INTENT(IN) :: problem1                 ! the objective function selected for the DC component f1 
                REAL(KIND=dp), INTENT(IN) :: factor             ! the scaling factor for objective     
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(4) :: a
                REAL(KIND=dp), DIMENSION(5) :: b                ! help variable
                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp
                REAL(KIND=dp) :: a1, a2, a3
                REAL(KIND=dp) :: apu                                ! help varaible
                REAL(KIND=dp) :: largest                            ! help varaible
                REAL(KIND=dp), DIMENSION(user_n) :: abs_sign        ! help varaible
                REAL(KIND=dp), DIMENSION(user_n) :: f_i             ! help varaible (only first n-1 placs are used
                INTEGER :: i, j, ind                                ! help variable
                


                SELECT CASE(problem1)
                
                   !-------------------------------------
                   !           Problem   0
                   !-------------------------------------  
                   CASE(0)              
                     DO i = 1, user_n
                       IF ( y(i) <= 0.0_dp) THEN
                          IF (y(i) == 0.0_dp) THEN
                            grad(i) = 0.5_dp                    
                          ELSE
                            grad(i) = -1.0_dp                   
                          END IF                
                       ELSE
                         IF (y(i) <= 2.0_dp ) THEN
                            grad(i) = 0.5_dp                    
                         ELSE
                            grad(i) = 2.0_dp                    
                         END IF                 
                       END IF
                     END DO
                   !-------------------------------------
 
                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------                    
                   CASE(1)              
                     a1 = y(1)**4 + y(2)**2
                     a2 = (2.0_dp - y(1))**2 + (2.0_dp - y(2))**2 
                     a3 = 2.0_dp * EXP(-y(1)+y(2))              
                
                     IF (a1 >= a2) THEN
                       IF (a1 >= a3) THEN
                         grad(1) = 4.0_dp * y(1)**3 
                         grad(2) = 2.0_dp * y(2)                    
                       ELSE
                         grad(1) = -2.0_dp * EXP(-y(1)+y(2))
                         grad(2) = 2.0_dp * EXP(-y(1)+y(2))                 
                       END IF               
                     ELSE
                       IF (a2 >= a3) THEN
                         grad(1) = 2.0_dp * y(1) - 4.0_dp
                         grad(2) = 2.0_dp * y(2) - 4.0_dp                   
                       ELSE
                         grad(1) = -2.0_dp * EXP(-y(1)+y(2))
                         grad(2) = 2.0_dp * EXP(-y(1)+y(2))                 
                       END IF                   
                     END IF
                
                     grad(1) = grad(1) + 2.0_dp * y(1) - 2.0_dp
                     grad(1) = grad(1) + 4.0_dp * y(1) - 5.0_dp
                     grad(1) = grad(1) + 2.0_dp * y(1)
  
                     grad(2) = grad(2) + 2.0_dp * y(2) - 4.0_dp
                     grad(2) = grad(2) + 2.0_dp * y(2) - 2.0_dp
                     grad(2) = grad(2) + 4.0_dp * y(2) - 4.0_dp                

                     grad = grad * factor
                     
                   !-------------------------------------      
                   
                   !-------------------------------------
                   !           Problem   2
                   !-------------------------------------                    
                   CASE(2)              
                     IF ( y(1) <= 1.0_dp ) THEN
                        grad(1) = -1.0_dp
                     ELSE
                        grad(1) = 1.0_dp
                     END IF
                
                     IF( ( ABS(y(1))-y(2) ) > 0.0_dp  ) THEN
                        IF(y(1) <= 0.0_dp) THEN 
                          grad(1) = grad(1) - 200.0_dp
                        ELSE
                          grad(1) = grad(1) + 200.0_dp  
                        END IF                   
                        grad(2) = -200.0_dp
                     ELSE        
                        grad(1) = grad(1) 
                        grad(2) = 0.0_dp
                     END IF
                     
                     grad = grad * factor                    
                   !-------------------------------------                   

                   !-------------------------------------
                   !           Problem   3
                   !-------------------------------------  
                   CASE(3)
                     grad(1) = 0.0_dp
                     grad(2) = 0.0_dp
                     grad(3) = 0.0_dp
                     grad(4) = 0.0_dp
                
                     IF ( y(1) <= 1.0_dp ) THEN
                       grad(1) = -1.0_dp
                     ELSE
                       grad(1) = 1.0_dp
                     END IF
                
                     IF( ( ABS(y(1))-y(2) ) > 0.0_dp  ) THEN
                       IF(y(1) <= 0.0_dp) THEN 
                          grad(1) = grad(1) - 200.0_dp
                       ELSE
                          grad(1) = grad(1) + 200.0_dp  
                       END IF                    
                       grad(2) = -200.0_dp
                     ELSE        
                       grad(1) = grad(1) 
                       grad(2) = 0.0_dp
                     END IF

                     IF( ( ABS(y(3))-y(4) ) > 0.0_dp  ) THEN
                       IF(y(3) <= 0.0_dp) THEN 
                          grad(3) =  - 180.0_dp
                       ELSE
                          grad(3) =  180.0_dp  
                       END IF                    
                       grad(4) = -180.0_dp
                     ELSE        
                       grad(3) = 0.0_dp 
                       grad(4) = 0.0_dp
                     END IF
                
                     IF ( y(3) <= 1.0_dp ) THEN
                       grad(3) = grad(3) - 1.0_dp
                     ELSE
                       grad(3) = grad(3) + 1.0_dp
                     END IF         
                
                     IF ( y(2) <= 1.0_dp ) THEN
                       grad(2) = grad(2) - 10.1_dp  
                     ELSE
                       grad(2) = grad(2) + 10.1_dp
                     END IF         
                
                     IF ( y(4) <= 1.0_dp ) THEN
                       grad(4) = grad(4) - 10.1_dp  
                     ELSE
                       grad(4) = grad(4) + 10.1_dp
                     END IF 
                
                     IF ( (y(2) + y(4) - 2.0_dp) <= 0.0_dp ) THEN
                       grad(2) = grad(2) - 4.95_dp
                       grad(4) = grad(4) - 4.95_dp  
                     ELSE
                       grad(2) = grad(2) + 4.95_dp
                       grad(4) = grad(4) + 4.95_dp
                     END IF 
                     
                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   4
                   !-------------------------------------  
                   CASE(4)
                     grad = 0.0_dp 
                     apu = ABS(y(1))
                     ind = 1
                
                     DO i = 2, user_n
                       IF (apu < ABS(y(i))) THEN 
                         apu = ABS(y(i))
                         ind = i
                       END IF
                     END DO
                
                     IF (y(ind) <= 0.0_dp) THEN 
                         grad(ind) = - user_n
                     ELSE 
                         grad(ind) = user_n
                    END IF  

                     grad = grad * factor
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   5
                   !-------------------------------------  
                   CASE(5)
                     grad = 0.0_dp 
                     apu = ABS(summa(y,1,user_n)) 
                     ind = 1

                     DO j = 2, 20
                       IF (apu <= ABS(summa(y,j,user_n)) ) THEN 
                         apu = ABS(summa(y,j,user_n))
                         ind = j
                       END IF
                     END DO

                     DO i = 1, user_n
                       grad(i) = (0.05_dp * ind)**(i-1)
                     END DO             
                
                     IF ( summa(y,ind,user_n) <= 0.0_dp ) THEN 
                         grad = -20.0_dp * grad
                     ELSE
                         grad = 20.0_dp * grad
                     END IF 
                     
                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   6
                   !-------------------------------------  
                   CASE(6)
                     grad = 0.0_dp 
                
                     grad(1) = 0.2_dp * y(1)
                
                     grad(2) = 1.0_dp + 0.2_dp * y(2) 

                     IF ( -y(2) <= 0.0_dp ) THEN 
                         grad(2) = grad(2)
                     ELSE
                         grad(2) = grad(2) - 10.0_dp
                     END IF        

                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   7
                   !-------------------------------------  
                   CASE(7)
                     grad = 0.0_dp 
                
                     IF ( (y(1)-1.0_dp) <= 0.0_dp ) THEN 
                       grad(1) = -1.0_dp
                     ELSE
                       grad(1) = 1.0_dp
                     END IF 

                     a(1) = y(1)**2 + y(2)**2 + ABS(y(2))
                     a(2) = y(1) + y(1)**2 + y(2)**2 + ABS(y(2)) - 0.5_dp
                     a(3) = ABS( y(1) - y(2) ) + ABS(y(2)) - 1.0_dp
                     a(4) = y(1) + y(1)**2 + y(2)**2                
                
                     ind = 1
                     DO i = 2, 4
                       IF ( a(ind) < a(i) ) THEN 
                         ind = i
                       END IF 
                     END DO
                
                     IF (ind == 1) THEN
                       grad(1) = grad(1) + 20.0_dp * y(1)
                       grad(2) = grad(2) + 20.0_dp * y(2)
                       IF (y(2) <= 0.0_dp) THEN
                          grad(2) = grad(2) - 10.0_dp
                       ELSE
                          grad(2) = grad(2) + 10.0_dp
                       END IF
                     END IF
                
                     IF (ind == 2) THEN
                       grad(1) = grad(1) + 10.0_dp + 20.0_dp * y(1)
                       grad(2) = grad(2) + 20.0_dp * y(2)
                       IF (y(2) <= 0.0_dp) THEN
                          grad(2) = grad(2) - 10.0_dp
                       ELSE
                          grad(2) = grad(2) + 10.0_dp
                       END IF
                     END IF 
                
                     IF (ind == 3) THEN
                        IF ( ( y(1) - y(2) ) <= 0.0_dp ) THEN
                          grad(1) = grad(1) - 10.0_dp 
                          grad(2) = grad(2) + 10.0_dp
                        ELSE
                          grad(1) = grad(1) + 10.0_dp
                          grad(2) = grad(2) - 10.0_dp                   
                        END IF
                    
                        IF (y(2) <= 0.0_dp) THEN
                          grad(2) = grad(2) - 10.0_dp                   
                        ELSE
                          grad(2) = grad(2) +   10.0_dp                 
                        END IF
                     END IF 
                
                     IF (ind == 4) THEN
                       grad(1) = grad(1) + 10.0_dp + 20.0_dp * y(1)
                       grad(2) = grad(2) + 20.0_dp * y(2)
                     END IF
                
                     IF ( (ABS(y(1)) - y(2) ) >= 0.0_dp ) THEN
                       grad(2) = grad(2) - 200.0_dp
                       IF ( y(1) <= 0.0_dp ) THEN
                         grad(1) = grad(1) - 200.0_dp
                       ELSE
                         grad(1) = grad(1) + 200.0_dp
                       END IF 
                     END IF
                     
                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   8
                   !-------------------------------------                      
                   CASE(8)
                     grad(1) = -8.0_dp + 8.0_dp * y(1)
                     grad(2) = -6.0_dp + 4.0_dp * y(2)
                     grad(3) = -4.0_dp + 4.0_dp * y(3)

                     IF ( y(1) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 2.0_dp
                     ELSE
                         grad(1) = grad(1) + 2.0_dp
                     END IF                 
                
                     IF ( y(2) <= 0.0_dp ) THEN 
                         grad(2) = grad(2) - 2.0_dp
                     ELSE
                         grad(2) = grad(2) + 2.0_dp
                     END IF 

                     IF ( y(3) <= 0.0_dp ) THEN 
                         grad(3) = grad(3) - 2.0_dp
                     ELSE
                         grad(3) = grad(3) + 2.0_dp
                     END IF                 
                
                     b(1) = y(1) + y(2) + 2.0_dp * y(3) - 3.0_dp
                     b(2) = -y(1)
                     b(3) = -y(2)
                     b(4) = -y(3)
                     b(5) = 0.0_dp
                
                     ind = 5
                
                     DO i = 1, 4
                       IF ( b(ind) < b(i) ) THEN 
                         ind = i
                       END IF 
                     END DO
                
                     IF (ind == 1) THEN
                       grad(1) = grad(1) + 10.0_dp 
                       grad(2) = grad(2) + 10.0_dp 
                       grad(3) = grad(3) + 20.0_dp 
                     END IF
                
                     IF (ind == 2) THEN
                       grad(1) = grad(1) - 10.0_dp 
                     END IF 
                
                     IF (ind == 3) THEN
                       grad(2) = grad(2) - 10.0_dp 
                     END IF 
                
                     IF (ind == 4) THEN
                       grad(3) = grad(3) - 10.0_dp 
                     END IF
                     
                     grad = grad * factor                    
                   !-------------------------------------

                   !-------------------------------------
                   !           Problem   9
                   !-------------------------------------                       
                   CASE(9)
                     grad(1) = 10.0_dp * y(1) - 16.0_dp
                     grad(2) = 10.0_dp * y(2) - 10.0_dp 
                     grad(3) = 10.0_dp * y(3) - 16.0_dp
                     grad(4) = 10.0_dp * y(4) - 10.0_dp  

                     grad = grad * factor                    
                   !------------------------------------
                    
                   !-------------------------------------
                   !           Problem   10
                   !-------------------------------------  
                   CASE(10)
                      grad = 2.0_dp * y    
                      grad = grad * factor
                   !------------------------------------- 

                   !-------------------------------------
                   !           Problem   11
                   !-------------------------------------
                   CASE(11)
                     grad = 0.0_dp
                     
                     IF (y(1)> 0.0_dp) THEN
                        grad(1) = grad(1) + 4.0_dp
                     ELSE
                        grad(1) = grad(1) -4.0_dp
                     END IF
                     
                     IF (y(2)> 0.0_dp) THEN
                        grad(2) = grad(2) + 2.0_dp
                     ELSE
                        grad(2) = grad(2) -2.0_dp
                     END IF

                     IF (y(3)> 0.0_dp) THEN
                        grad(3) = grad(3) + 22.0_dp
                     ELSE
                        grad(3) = grad(3) -22.0_dp
                     END IF                  
                     
                     grad(1) = grad(1) -33.0_dp 
                     grad(2) = grad(2) +16.0_dp 
                     grad(3) = grad(3) -24.0_dp 

                     apu = 2.0_dp*ABS(y(2))-3.0_dp*y(1)-7.0_dp
                     IF (apu > 0.0_dp) THEN 
                        grad(1) = grad(1) -300.0_dp
                        IF (y(2)>0.0_dp) THEN 
                           grad(2) = grad(2) +200.0_dp
                        ELSE
                           grad(2) = grad(2) -200.0_dp
                        END IF 
                     END IF
                     
                     apu = ABS(y(3))-4.0_dp*y(1)-11.0_dp
                     IF (apu > 0.0_dp) THEN 
                        grad(1) = grad(1) -400.0_dp
                        IF (y(3)>0.0_dp) THEN 
                           grad(3) = grad(3) +100.0_dp
                        ELSE
                           grad(3) = grad(3) -100.0_dp
                        END IF 
                     END IF  

                     grad = grad * factor                    
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   12
                   !-------------------------------------
                   CASE(12)                
                     grad = 0.0_dp
                     
                     DO i = 1, user_n
                        IF (y(i) > 0.0_dp) THEN 
                           grad(i) = grad(i) +1.0_dp
                        ELSE
                           grad(i) = grad(i) -1.0_dp
                        END IF               
                     END DO                  

                     DO i = 1, user_n
                         apu = 2.0_dp * (y(i)*y(i)- y(i)-1.0_dp) 
                         IF (apu > 0.0_dp) THEN 
                            grad(i) = grad(i) + 40.0_dp * y(i) -20.0_dp
                         END IF
                     END DO
                     
                     grad = grad * factor
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !           Problem   13
                   !-------------------------------------
                   CASE(13)             
                     grad = 0.0_dp
                     
                     DO i = 1, user_n-1
                        IF ( (y(i)+y(i+1)) > 0.0_dp) THEN 
                           grad(i) = grad(i) + 1.0_dp
                           grad(i+1) = grad(i+1) + 1.0_dp
                        ELSE
                           grad(i) = grad(i) - 1.0_dp
                           grad(i+1) = grad(i+1) - 1.0_dp
                        END IF
                     END DO
                     
                     DO i = 1, user_n-2
                        IF ( (y(i)+y(i+2)) > 0.0_dp) THEN 
                           grad(i) = grad(i) + 1.0_dp
                           grad(i+2) = grad(i+2) + 1.0_dp
                        ELSE
                           grad(i) = grad(i) - 1.0_dp
                           grad(i+2) = grad(i+2) - 1.0_dp
                        END IF
                     END DO
                     
                     IF ( (y(1)+y(9)) > 0.0_dp) THEN 
                        grad(1) = grad(1) + 1.0_dp
                        grad(9) = grad(9) + 1.0_dp
                     ELSE
                        grad(1) = grad(1) - 1.0_dp
                        grad(9) = grad(9) - 1.0_dp
                     END IF 
                     
                     IF ( (y(1)+y(10)) > 0.0_dp) THEN 
                        grad(1) = grad(1) + 1.0_dp
                        grad(10) = grad(10) + 1.0_dp
                     ELSE
                        grad(1) = grad(1) - 1.0_dp
                        grad(10) = grad(10) - 1.0_dp
                     END IF 
                     
                     IF ( (y(2)+y(10)) > 0.0_dp) THEN 
                        grad(2) = grad(2) + 1.0_dp
                        grad(10) = grad(10) + 1.0_dp
                     ELSE
                        grad(2) = grad(2) - 1.0_dp
                        grad(10) = grad(10) - 1.0_dp
                     END IF 
                     
                     IF ( (y(1)+y(5)) > 0.0_dp) THEN 
                        grad(1) = grad(1) + 1.0_dp
                        grad(5) = grad(5) + 1.0_dp
                     ELSE
                        grad(1) = grad(1) - 1.0_dp
                        grad(5) = grad(5) - 1.0_dp
                     END IF 
                     
                     IF ( (y(4)+y(7)) > 0.0_dp) THEN 
                        grad(4) = grad(4) + 1.0_dp
                        grad(7) = grad(7) + 1.0_dp
                     ELSE
                        grad(4) = grad(4) - 1.0_dp
                        grad(7) = grad(7) - 1.0_dp
                     END IF                      
                 
                     apu = 0.0_dp
                     DO i = 1, user_n
                        apu = apu + y(i)
                     END DO
                     
                     IF ((apu-1.0_dp)>0.0_dp) THEN 
                        grad = grad + 10.0_dp 
                     END IF
                     
                     DO i = 1, user_n
                        IF (-y(i) > 0.0_dp) THEN 
                            grad(i) = grad(i) - 10.0_dp
                        END IF
                     END DO     

                     grad = grad * factor                    
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   14
                   !-------------------------------------
                   CASE(14)         
                     grad = 0.0_dp
                   
                     apu = 0.0_dp
                     DO j =1, user_n
                        apu = apu + y(j)/(1+j-1.0_dp)
                     END DO
                     IF (apu >= 0.0_dp) THEN 
                           abs_sign(1) = 1.0_dp
                     ELSE   
                           abs_sign(1) = -1.0_dp
                     END IF        
                     largest = ABS(apu) 
                     ind = 1
                     
                     DO i = 2, user_n
                        apu = 0.0_dp
                        DO j = 1, user_n
                          apu = apu + y(j)/(i+j-1.0_dp)
                        END DO
                        IF (apu >= 0.0_dp) THEN 
                           abs_sign(i) = 1.0_dp
                        ELSE   
                           abs_sign(i) = -1.0_dp
                        END IF
                        IF (ABS(apu) > largest) THEN 
                            largest = ABS(apu)
                            ind = i
                        END IF
                     END DO
                     
                     DO j = 1, user_n
                        grad(j) = abs_sign(ind) * user_n / (ind+j-1)
                     END DO                      
    
                     grad = grad * factor   
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   15
                   !-------------------------------------
                   CASE(15)     
                     grad = 0.0_dp
                     
                     DO i = 1, user_n-1
                         a(1) = y(i)**4 + y(i+1)**2
                         a(2) = (2.0_dp-y(i))**2 + (2.0_dp-y(i+1))**2
                         a(3) = 2.0_dp*EXP(-y(i)+y(i+1))
                         
                         IF (a(1)>a(2)) THEN
                            IF (a(1)>a(3)) THEN
                              f_i(i) = a(1)
                            ELSE
                              f_i(i) = a(3)
                            END IF
                         ELSE
                            IF (a(2)>a(3)) THEN
                              f_i(i) = a(2)
                            ELSE
                              f_i(i) = a(3)
                            END IF      
                         END IF
                     END DO
                     
                     largest = f_i(1)
                     ind = 1
                     DO i =2, user_n-1
                         IF (f_i(i) > largest) THEN 
                            largest = f_i(i)
                            ind = i
                         END IF
                     END DO
                     
                     
                         a(1) = y(ind)**4 + y(ind+1)**2
                         a(2) = (2.0_dp-y(ind))**2 + (2.0_dp-y(ind+1))**2
                         a(3) = 2.0_dp*EXP(-y(ind)+y(ind+1))
                         IF (a(1)>a(2)) THEN
                            IF (a(1)>a(3)) THEN
                              grad(ind) = 4.0_dp*y(ind)**3 
                              grad(ind+1) = 2.0_dp*y(ind+1)
                            ELSE
                              grad(ind) = -2.0_dp * EXP(-y(ind)+y(ind+1)) 
                              grad(ind+1) = 2.0_dp * EXP(-y(ind)+y(ind+1)) 
                            END IF
                         ELSE
                            IF (a(2)>a(3)) THEN
                              grad(ind) = -2.0_dp * (2.0_dp -y(ind))
                              grad(ind+1) = -2.0_dp * (2.0_dp -y(ind+1))
                            ELSE
                              grad(ind) = -2.0_dp * EXP(-y(ind)+y(ind+1)) 
                              grad(ind+1) = 2.0_dp * EXP(-y(ind)+y(ind+1)) 
                            END IF      
                         END IF                  
                     
                     grad = (user_n-1.0_dp)*grad

                     grad = grad * factor
                     
                   !-------------------------------------   
                       
                   
                   !-------------------------------------
                   !             Problem 16:
                   !         Chained Crescent I
                   !-------------------------------------
                   CASE(16)
                     grad = 0.0_dp
                     
                     apu = 0.0_dp
                     DO i = 1, user_n-1
                        apu = apu + y(i)**2 + (y(i+1)-1.0_dp)**2
                        apu = apu + y(i+1) -1.0_dp
                     END DO
                     
                     apu = 2.0_dp *apu
                     IF (apu > 0.0_dp) THEN
                       DO i = 1, user_n-1
                          grad(i) = grad(i) + 4.0_dp *y(i)
                          grad(i+1) = grad(i+1) + 4.0_dp *(y(i+1)-1.0_dp) +2.0_dp
                       END DO
                     END IF
                    
                     grad = grad * factor
                     
                   !-------------------------------------
   
                   
                   !-------------------------------------
                   !          BAD critical point
                   !-------------------------------------
                    CASE(30)                   
                      IF ((y(1)**2) >= y(1)) THEN
                         grad = 2*y(1)
                      ELSE
                         grad = 1.0_dp
                      END IF
                   !-------------------------------------       

                   
 
                END SELECT                      
                
           END FUNCTION subgradient_f1      
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION subgradient_f2(y, problem2, factor, user_n) RESULT(grad)                
                !
                ! Calculate a subgradient of the DC component f_2 at a point 'y'.
                ! Variable 'problem1' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component f_2 is calculated
                INTEGER, INTENT(IN) :: problem2                 ! the objective function selected for the DC component f2
                REAL(KIND=dp), INTENT(IN) :: factor             ! the scaling factor for objective      
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component f_2 at a point 'y'
                REAL(KIND=dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
                REAL(KIND=dp), PARAMETER :: two=2.0_dp
                REAL(KIND=dp) :: a,b 
                REAL(KIND=dp) :: a4, a5, a6 
                REAL(KIND=dp) :: apu                                ! help varaible
                REAL(KIND=dp) :: largest, abs_sign                  ! help varaible
                INTEGER :: i, j, ind                                ! help variable
                
                REAL(KIND=dp), DIMENSION(3) :: aa               ! sign of the active partial objective for f1 and for f2    
                REAL(KIND=dp), DIMENSION(user_n) :: term        ! sign of the active partial objective for f1 and for f2    
                
                
                SELECT CASE(problem2)
                
                   !-------------------------------------
                   !           Problem   0
                   !-------------------------------------  
                   CASE(0)              
                     grad = 0.0_dp
                            
                     DO i =1, user_n
                       IF ( y(i) <= 0.0_dp ) THEN
                         IF ( y(i) == 0.0_dp ) THEN
                           grad(i) = 0.5_dp 
                         ELSE
                           grad(i) = -0.5_dp
                         END IF             
                       ELSE
                         grad(i) = 1.0_dp                   
                       END IF
                     END DO
                   !-------------------------------------

                   !-------------------------------------
                   !           Problem   1
                   !-------------------------------------                    
                   CASE(1)              
                     grad = 0.0_dp
                
                     a4 = y(1)**2 - 2.0_dp * y(1) + y(2)**2 - 4.0_dp * y(2) + 4.0_dp
                     a5 = 2.0_dp * y(1)**2 - 5.0_dp * y(1) + y(2)**2 - 2.0_dp * y(2) 
                     a5 = a5 + 4.0_dp
                     a6 = y(1)**2 + 2.0_dp * y(2)**2 - 4.0_dp * y(2) + 1.0_dp               
                
                     IF ( (a4 + a5) >= (a4 + a6) ) THEN
                        IF ( (a4 + a5) >= (a5 + a6)) THEN
                          grad(1) = 6.0_dp * y(1) - 7.0_dp
                          grad(2) = 4.0_dp * y(2) - 6.0_dp
                        ELSE
                          grad(1) = 6.0_dp * y(1) - 5.0_dp
                          grad(2) = 6.0_dp * y(2) - 6.0_dp
                        END IF              
                     ELSE
                        IF ( (a4 + a6) >= (a5 + a6) ) THEN
                          grad(1) = 4.0_dp * y(1) - 2.0_dp
                          grad(2) = 6.0_dp * y(2) - 8.0_dp
                        ELSE
                          grad(1) = 6.0_dp * y(1) - 5.0_dp
                          grad(2) = 6.0_dp * y(2) - 6.0_dp 
                        END IF                  
                     END IF
                     
                     grad = grad * factor                    
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Problem   2
                   !-------------------------------------                    
                   CASE(2)              
                     IF(y(1) <= 0) THEN
                       grad(1) = -100.0_dp 
                     ELSE
                       grad(1) = 100.0_dp
                     END IF
                
                     grad(2) = -100.0_dp
                     
                     grad = grad * factor                    
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Problem   3
                   !-------------------------------------  
                   CASE(3)
                     IF(y(1) <= 0.0_dp) THEN
                       grad(1) = -100.0_dp 
                     ELSE
                       grad(1) = 100.0_dp
                     END IF
                
                     grad(2) = -100.0_dp
                
                     IF(y(3) <= 0.0_dp) THEN
                       grad(3) = -90.0_dp 
                     ELSE
                       grad(3) = 90.0_dp
                     END IF
                
                     grad(4) = -90.0_dp             
                
                     IF( (y(2) - y(4) ) <= 0.0_dp) THEN
                       grad(2) = grad(2) - 4.95_dp 
                       grad(4) = grad(4) + 4.95_dp
                     ELSE
                       grad(2) = grad(2) + 4.95_dp 
                       grad(4) = grad(4) - 4.95_dp                 
                     END IF         

                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   4
                   !-------------------------------------  
                   CASE(4)
                     DO i = 1, user_n
                       IF (y(i) <= 0.0_dp ) THEN 
                         grad(i) = -1.0_dp
                       ELSE
                         grad(i) = 1.0_dp
                       END IF
                     END DO
                     
                     grad = grad * factor
                   !-------------------------------------
                    
                   !-------------------------------------
                   !           Problem   5
                   !-------------------------------------  
                   CASE(5)
                     grad = 0.0_dp 

                     DO j = 1, 20
                      IF (summa(y,j,user_n) <= 0.0_dp) THEN 
                        DO i = 1, user_n
                            grad(i) = grad(i) - (0.05_dp * j)**(i-1)
                        END DO                  
                      ELSE 
                        DO i = 1, user_n
                            grad(i) = grad(i) + (0.05_dp * j)**(i-1)
                        END DO                   
                      END IF
                     END DO   

                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   6
                   !-------------------------------------  
                   CASE(6)
                     DO i = 1, user_n
                       IF (y(i) <= 0.0_dp ) THEN 
                         grad(i) = -1.0_dp
                       ELSE
                         grad(i) = 1.0_dp
                       END IF
                     END DO 

                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   7
                   !-------------------------------------  
                   CASE(7)
                     grad = 0.0_dp
                
                     IF (y(2) <= 0.0_dp ) THEN 
                         grad(2) = 20.0_dp * y(2) - 10.0_dp 
                     ELSE
                         grad(2) = 20.0_dp * y(2) + 10.0_dp 
                     END IF
                
                     grad(2) = grad(2) - 100.0_dp
                     grad(1) = 20.0_dp * y(1)

                     IF (y(1) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 100.0_dp
                     ELSE
                         grad(1) = grad(1) + 100.0_dp
                     END IF
                     
                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   8
                   !-------------------------------------                    
                   CASE(8)
                     grad = 0.0_dp
                
                     IF ( (y(1)-y(2)) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 1.0_dp 
                         grad(2) = grad(2) + 1.0_dp
                     ELSE
                         grad(1) = grad(1) + 1.0_dp
                         grad(2) = grad(2) - 1.0_dp 
                     END IF
                
                     IF ( (y(1)-y(3)) <= 0.0_dp ) THEN 
                         grad(1) = grad(1) - 1.0_dp 
                         grad(3) = grad(3) + 1.0_dp
                     ELSE
                         grad(1) = grad(1) + 1.0_dp
                         grad(3) = grad(3) - 1.0_dp     
                     END IF 
                     
                     grad = grad * factor                    
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   9
                   !-------------------------------------                     
                   CASE(9)
                     grad = 0.0_dp 

                     a = (y(1) -2.0_dp)**2 + y(2)**2
                     b = (y(3) -2.0_dp)**2 + y(4)**2
                
                     IF ( a >= b) THEN 
                       grad(1) = grad(1) + 2.0_dp * (y(1)-2.0_dp)
                       grad(2) = grad(2) + 2.0_dp * y(2)
                     ELSE 
                       grad(3) = grad(3) + 2.0_dp * (y(3)-2.0_dp)
                       grad(4) = grad(4) + 2.0_dp * y(4)
                     END IF

                     a = (y(1) -2.0_dp)**2 + (y(2)-1.0_dp)**2
                     b = (y(3) -2.0_dp)**2 + (y(4)-1.0_dp)**2
                
                     IF ( a >= b) THEN 
                       grad(1) = grad(1) + 2.0_dp * (y(1)-2.0_dp)
                       grad(2) = grad(2) + 2.0_dp * (y(2)-1.0_dp)
                     ELSE 
                       grad(3) = grad(3) + 2.0_dp * (y(3)-2.0_dp)
                       grad(4) = grad(4) + 2.0_dp * (y(4)-1.0_dp)
                     END IF

                     a = (y(1) -3.0_dp)**2 + y(2)**2
                     b = (y(3) -3.0_dp)**2 + y(4)**2
                
                     IF ( a >= b) THEN 
                       grad(1) = grad(1) + 2.0_dp * (y(1)-3.0_dp)
                       grad(2) = grad(2) + 2.0_dp * (y(2))
                     ELSE 
                       grad(3) = grad(3) + 2.0_dp * (y(3)-3.0_dp)
                       grad(4) = grad(4) + 2.0_dp * (y(4)) 
                     END IF

                     a = (y(1))**2 + (y(2)-2.0_dp)**2
                     b = (y(3))**2 + (y(4)-2.0_dp)**2
                
                     IF ( a >= b) THEN 
                       grad(1) = grad(1) + 2.0_dp * (y(1))
                       grad(2) = grad(2) + 2.0_dp * (y(2)-2.0_dp)
                     ELSE 
                       grad(3) = grad(3) + 2.0_dp * (y(3))
                       grad(4) = grad(4) + 2.0_dp * (y(4)-2.0_dp)
                     END IF

                     a = (y(1)-1.0_dp)**2 + (y(2)-2.0_dp)**2
                     b = (y(3)-1.0_dp)**2 + (y(4)-2.0_dp)**2
                
                     IF ( a >= b) THEN 
                       grad(1) = grad(1) + 2.0_dp * (y(1)-1.0_dp)
                       grad(2) = grad(2) + 2.0_dp * (y(2)-2.0_dp)
                     ELSE 
                       grad(3) = grad(3) + 2.0_dp * (y(3)-1.0_dp)
                       grad(4) = grad(4) + 2.0_dp * (y(4)-2.0_dp)
                     END IF           

                     grad = grad * factor                    
                   !------------------------------------
                    
                   !-------------------------------------
                   !           Problem   10
                   !-------------------------------------  
                   CASE(10)
                     grad = 0.0_dp
                
                     DO i = 2, user_n
                       IF ( y(i) - y(i-1) <= 0.0_dp) THEN
                         grad(i-1) = grad(i-1) + 1.0_dp 
                         grad(i) = grad(i) - 1.0_dp 
                       ELSE
                         grad(i-1) = grad(i-1) - 1.0_dp  
                         grad(i) = grad(i) + 1.0_dp     
                       END IF
                     END DO
                     
                     grad = grad * factor
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Problem   11
                   !-------------------------------------  
                   CASE(11)
                     grad = 0.0_dp
                
                     grad(1) = -140.0_dp
                     
                     IF (y(2)>0.0_dp) THEN 
                        grad(2) = 40.0_dp
                     ELSE
                        grad(2) = -40.0_dp
                     END IF    

                     grad = grad * factor                    
                   !-------------------------------------       

                   !-------------------------------------
                   !           Problem   12
                   !-------------------------------------
                   CASE(12)                
                     grad = 0.0_dp

                     DO i = 1, user_n
                        grad(i) = grad(i) +20.0_dp*y(i) - 10.0_dp
                     END DO
                     
                     term = 0.0_dp
                     DO i = 1, user_n
                        DO j = 1, user_n 
                           IF (i /= j) THEN
                               term(i) = term(i) + ABS(y(j))
                           END IF
                        END DO
                     END DO
                     
                     largest = term(1)
                     ind = 1 
                     
                     DO i = 2, user_n
                        IF (term(i) > largest) THEN 
                           largest = term(i)
                           ind = i
                        END IF
                     END DO
                     
                     DO i = 1, user_n
                        IF (i /= ind) THEN
                           IF (y(i) > 0.0_dp) THEN 
                              grad(i) = grad(i) +1.0_dp
                           ELSE
                              grad(i) = grad(i) -1.0_dp
                           END IF
                        END IF                   
                     END DO
                     
                     grad = grad * factor                    
                   !-------------------------------------   
                   
                   !-------------------------------------
                   !           Problem   13
                   !-------------------------------------
                   CASE(13)             
                     grad = 0.0_dp
                     
                     DO i = 1, user_n-1
                        IF ( y(i) > 0.0_dp) THEN
                           grad(i) = grad(i) + 1.0_dp
                        ELSE
                           grad(i) = grad(i) - 1.0_dp
                        END IF
                        IF ( y(i+1) > 0.0_dp) THEN
                           grad(i+1) = grad(i+1) + 1.0_dp
                        ELSE
                           grad(i+1) = grad(i+1) - 1.0_dp
                        END IF                      
                     END DO
                     
                     DO i = 1, user_n-2
                        IF ( y(i) > 0.0_dp) THEN
                           grad(i) = grad(i) + 1.0_dp
                        ELSE
                           grad(i) = grad(i) - 1.0_dp
                        END IF
                        IF ( y(i+2) > 0.0_dp) THEN
                           grad(i+2) = grad(i+2) + 1.0_dp
                        ELSE
                           grad(i+2) = grad(i+2) - 1.0_dp
                        END IF      
                     END DO                  
                     
                    
                     IF ( y(1) > 0.0_dp) THEN
                        grad(1) = grad(1) + 3.0_dp
                     ELSE
                        grad(1) = grad(1) - 3.0_dp
                     END IF

                     IF ( y(2) > 0.0_dp) THEN
                        grad(2) = grad(2) + 1.0_dp
                     ELSE
                        grad(2) = grad(2) - 1.0_dp
                     END IF
                     
                     IF ( y(4) > 0.0_dp) THEN
                        grad(4) = grad(4) + 1.0_dp
                     ELSE
                        grad(4) = grad(4) - 1.0_dp
                     END IF 

                     IF ( y(5) > 0.0_dp) THEN
                        grad(5) = grad(5) + 1.0_dp
                     ELSE
                        grad(5) = grad(5) - 1.0_dp
                     END IF

                     IF ( y(7) > 0.0_dp) THEN
                        grad(7) = grad(7) + 1.0_dp
                     ELSE
                        grad(7) = grad(7) - 1.0_dp
                     END IF

                     IF ( y(9) > 0.0_dp) THEN
                        grad(9) = grad(9) + 1.0_dp
                     ELSE
                        grad(9) = grad(9) - 1.0_dp
                     END IF

                     IF ( y(10) > 0.0_dp) THEN
                        grad(10) = grad(10) + 2.0_dp
                     ELSE
                        grad(10) = grad(10) - 2.0_dp
                     END IF     

                     grad = grad * factor                    
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   14
                   !-------------------------------------
                   CASE(14)         
                     grad = 0.0_dp
                     
                     DO i = 1, user_n
                        apu = 0.0_dp
                        DO j = 1, user_n
                          apu = apu + y(j)/(i+j-1.0_dp)
                        END DO
                        IF (apu >= 0.0_dp) THEN 
                           abs_sign = 1.0_dp
                        ELSE   
                           abs_sign = -1.0_dp
                        END IF
                        DO j = 1, user_n
                          grad(j) = grad(j) + abs_sign / (i+j-1.0_dp)
                        END DO
                     END DO  

                     grad = grad * factor                    
                   !-------------------------------------   

                   !-------------------------------------
                   !           Problem   15
                   !-------------------------------------
                   CASE(15)     
                     grad = 0.0_dp
                                 
                     DO ind = 1, user_n-1
                         aa(1) = y(ind)**4 + y(ind+1)**2
                         aa(2) = (2.0_dp-y(ind))**2 + (2.0_dp-y(ind+1))**2
                         aa(3) = 2.0_dp*EXP(-y(ind)+y(ind+1))
                         IF (aa(1)>aa(2)) THEN
                            IF (aa(1)>aa(3)) THEN
                              grad(ind) = grad(ind) + 4.0_dp*y(ind)**3 
                              grad(ind+1) = grad(ind+1) + 2.0_dp*y(ind+1)
                            ELSE
                              grad(ind) = grad(ind) -2.0_dp * EXP(-y(ind)+y(ind+1)) 
                              grad(ind+1) = grad(ind+1) + 2.0_dp * EXP(-y(ind)+y(ind+1)) 
                            END IF
                         ELSE
                            IF (aa(2)>aa(3)) THEN
                              grad(ind) = grad(ind) - 2.0_dp * (2.0_dp -y(ind))
                              grad(ind+1) = grad(ind+1) - 2.0_dp * (2.0_dp -y(ind+1))
                            ELSE
                              grad(ind) = grad(ind) -2.0_dp * EXP(-y(ind)+y(ind+1)) 
                              grad(ind+1) = grad(ind+1) +2.0_dp * EXP(-y(ind)+y(ind+1)) 
                            END IF      
                         END IF 
                     END DO  
                                      
                     grad = grad * factor
                     
                   !-------------------------------------                      
                   
                   !-------------------------------------
                   !              Problem 16:
                   !         Chained Crescent I
                   !-------------------------------------
                   CASE(16)
                     grad = 0.0_dp

                     DO i = 1, user_n-1
                        grad(i) = grad(i) + 2.0_dp * y(i) 
                        grad(i+1) = grad(i+1) + 2.0_dp * (y(i+1) -1.0_dp) +1.0_dp
                     END DO

                     grad = grad * factor
                     
                   !-------------------------------------
    
                   
                   !-------------------------------------
                   !          BAD critical point
                   !-------------------------------------
                    CASE(30)                   
                      IF ( (0.5_dp*y(1)**2) >= -y(1)) THEN
                         grad = y(1)
                      ELSE
                         grad = -1.0_dp
                      END IF
                   !-------------------------------------                      
                   
                
                END SELECT  

           END FUNCTION subgradient_f2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
           
           
           FUNCTION summa(y, j, user_n) RESULT(f)           
                !
                ! Calculates the sum used in f_1 and f_2 for parameter t_j
                !
                ! NOTICE: The dimension of 'y' has to be 'n'. 'j' needs to be an integer from interval [1,20]
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component f_2 is calculated
                INTEGER, INTENT(IN) ::  j       ! used to determines the parameter t_j
                INTEGER, INTENT(IN) ::  user_n  ! the dimension of the problem 
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f, t, apu                          ! the function value of the sum used in f_1 and parameter t_j
                INTEGER :: i                                    ! help variable
                
                f = 0.0_dp
                apu = 1.0_dp/ user_n
                
                DO i = 1, user_n 
                    t = (0.05_dp * j )**(i-1)
                    f = f + (y(i)-apu)*t            
                END DO 
            
           END FUNCTION summa

           
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                          | |
        !| |                       CONSTRAINT FUNCTIONS:                              | |
        !| |                                                                          | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   

        
           
        !********************************************************************************
        !                                                                               |
        !              FUNCTION VALUES OF THE DC COMPONENTS g_1 AND g_2                 |
        !                                                                               |
        !********************************************************************************

           FUNCTION g1(y, problem, factor, user_n) RESULT(g)      
                !
                ! For the constraint calculates the function value of the DC component g_1 at a point 'y'.
                ! Variable 'problem' identifies the constraint function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component g_1 is calculated
                REAL(KIND=dp), INTENT(IN) :: factor             ! a scaling factor
                INTEGER, INTENT(IN) :: problem                  ! identifies the constraint 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem

                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: g       ! the function value of the DC component g_1 at a point 'y'
                REAL(KIND=dp) :: a1       
               
                g = 0.0_dp
                SELECT CASE(problem)
                   
                    
                   !-------------------------------------
                   !           Constraint   1
                   !-------------------------------------                    
                   CASE(1)              
                       
                      a1 = (y(1)+1.5_dp)**2+y(2)**2-4.0_dp
                      a1 = a1+(y(1)-1.0_dp)**2+(y(2)-1.0_dp)**2-1.0_dp
                      IF(a1>=0)THEN
                        g = a1
                      ELSE 
                        g = 0.0_dp
                      END IF
                      g = g*factor
                    
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Constraint   2
                   !-------------------------------------                    
                   CASE(2)
                      g = 0.0_dp * factor

                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Constraint   3
                   !-------------------------------------  
                   CASE(3)
                   
                       g = user_n * 0.5 * factor

                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   4
                   !-------------------------------------  
                   CASE(4)
                       
                       a1 = (y(1))**2+(y(2))**2-8.0_dp
                    
                       g = factor * a1

                    
                   !-------------------------------------
                    
                   !-------------------------------------
                   !           Constraint   5
                   !-------------------------------------  
                   CASE(5)
                       
                       
                 
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   6
                   !-------------------------------------  
                   CASE(6)
                   
               
                   !-------------------------------------

                
                END SELECT              
            


           END FUNCTION g1
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           
           FUNCTION g2(y, problem, factor, user_n) RESULT(g)      
                !
                ! For the constraint calculates the function value of the DC component g_2 at a point 'y'.
                ! Variable 'problem' identifies the constraint function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the function value of the DC component g_2 is calculated
                REAL(KIND=dp), INTENT(IN) :: factor             ! a scaling factor
                INTEGER, INTENT(IN) :: problem                  ! identifies the constraint 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: g       ! the function value of the DC component g_2 at a point 'y'
                REAL(KIND=dp) :: a1, a2         
                INTEGER :: i
                
                g = 0.0_dp
                
                
                SELECT CASE(problem)

                   
                   !-------------------------------------
                   !           Constraint   1
                   !-------------------------------------                    
                   CASE(1)
                   
                      g = (y(1)-1.0_dp)**2+(y(2)-1.0_dp)**2-1.0_dp
                      g = g*factor
                   
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Constraint   2
                   !-------------------------------------                    
                    CASE(2)
                    
                      a1 = y(1)+y(2)+y(3)+y(4)-5.5_dp
                      a2 = y(1)**2+y(2)**2+y(3)**2+y(4)**2-10.0_dp
                      
                      IF (a1>=a2) THEN
                         g = a1
                      ELSE   
                         g = a2
                      END IF
                      g = g * factor                  
                      
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Constraint   3
                   !-------------------------------------  
                   CASE(3)
                   
                     g = 0.0_dp                 
                     
                     DO i = 2, user_n, 2
                       g = g + (y(i)-0.5_dp)**2
                     END DO
                     DO i = 1, user_n, 2
                       g = g + (y(i)+0.5_dp)**2
                     END DO
                     g = g * factor 
                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   4
                   !-------------------------------------  
                   CASE(4)
                   
                     g = 0.0_dp                 
                     
                   !-------------------------------------
                    
                   !-------------------------------------
                   !           Constraint   5
                   !-------------------------------------  
                   CASE(5)

                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   6
                   !-------------------------------------  
                   CASE(6)
                   

               
                   !-------------------------------------
                              
                
                END SELECT              
            


           END FUNCTION g2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           

        !********************************************************************************
        !                                                                               |
        !                SUBGRADIENTS OF THE DC COMPONENTS g_1 AND g_2                  |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION subgradient_g1(y, problem, factor, user_n) RESULT(grad)                
                !
                ! For the constraint calculates a subgradient of the DC component g_1 at a point 'y'.
                ! Variable 'problem' identifies the constraint function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y    ! a point where the subgradient of the DC component g_1 is calculated
                REAL(KIND=dp), INTENT(IN) :: factor             ! a scaling factor
                INTEGER, INTENT(IN) :: problem                  ! identifies the constraint 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem
            
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component g_1 at a point 'y'
                REAL(KIND=dp) :: a1
                
                grad = 0.0_dp
 
                SELECT CASE(problem)

                   !-------------------------------------
                   !           Constraint   1
                   !-------------------------------------                    
                   CASE(1)              
                       
                       a1 = (y(1)+1.5_dp)**2+y(2)**2-4.0_dp
                       a1 = a1+(y(1)-1.0_dp)**2+(y(2)-1.0_dp)**2-1.0_dp
                       IF(a1 >= 0)THEN
                          grad(1) = 2.0_dp*(y(1)+1.5_dp)+2.0_dp*(y(1)-1.0_dp)
                          grad(2) = 2.0_dp*y(2)+2.0_dp*(y(2)-1.0_dp)
                       ELSE 
                          grad = 0.0_dp
                       END IF
                       grad = grad*factor                     
                    
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Constraint   2
                   !-------------------------------------                    
                   CASE(2)
                       
                       grad = 0.0_dp               
                       
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Constraint   3
                   !-------------------------------------  
                   CASE(3)
                       
                       grad = 0.0_dp
                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   4
                   !-------------------------------------  
                   CASE(4)

                       grad(1) = 2.0_dp*(y(1))
                       grad(2) = 2.0_dp*(y(2))
                         
                       grad = grad*factor                     
                  
                   !-------------------------------------
                    
                   !-------------------------------------
                   !           Constraint   5
                   !-------------------------------------  
                   CASE(5)
                       
                   
                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   6
                   !-------------------------------------  
                   CASE(6)
                   
               
                   !-------------------------------------
                

                
                END SELECT  

           END FUNCTION subgradient_g1
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
           
           FUNCTION subgradient_g2(y, problem, factor, user_n) RESULT(grad)                
                !
                ! For the constraint calculates a subgradient of the DC component g_2 at a point 'y'.
                ! Variable 'problem' identifies the constraint function used.
                !
                ! NOTICE: * The dimension of 'y' has to be 'user_n'.
                !         * The dimension of 'grad' is also 'user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: y         ! a point where the subgradient of the DC component g_2 is calculated
                REAL(KIND=dp), INTENT(IN) :: factor                  ! a scaling factor
                INTEGER, INTENT(IN) :: problem                       ! identifies the constraint 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem              
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(SIZE(y)) :: grad       ! the subgradient of the DC component g_2 at a point 'y'
                REAL(KIND=dp) :: a1, a2       
                INTEGER :: i
                
                grad = 0.0_dp
               
                SELECT CASE(problem)

                   !-------------------------------------
                   !           Constraint   1
                   !-------------------------------------                    
                   CASE(1) 
                   
                       grad(1) = 2.0_dp*(y(1)-1.0_dp)
                       grad(2) = 2.0_dp*(y(2)-1.0_dp)
                       grad = factor * grad   
                       
                   !------------------------------------- 
                   
                   !-------------------------------------
                   !           Constraint   2
                   !-------------------------------------                    
                   CASE(2)
                   
                       a1 = y(1)+y(2)+y(3)+y(4)-5.5_dp
                       a2 = y(1)**2+y(2)**2+y(3)**2+y(4)**2-10.0_dp
                       
                       IF (a1 >= a2) THEN
                          grad = 1.0_dp
                       ELSE   
                          DO i = 1, user_n
                             grad(i) = 2.0_dp * y(i)
                          END DO
                       END IF
                       
                       grad = factor * grad                  
                   !-------------------------------------                   
                   
                   !-------------------------------------
                   !           Constraint   3
                   !-------------------------------------  
                   CASE(3)
                       grad = 0.0_dp
                     
                       DO i = 2, user_n, 2
                         grad(i) = grad(i) + 2.0_dp*(y(i)-0.5_dp)
                       END DO
                       DO i = 1, user_n, 2
                         grad(i) = grad(i) + 2.0_dp*(y(i)+0.5_dp)
                       END DO
                       grad = grad * factor 
                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   4
                   !-------------------------------------  
                   CASE(4)
                       grad = 0.0_dp
                  
                   
                   !-------------------------------------
                    
                   !-------------------------------------
                   !           Constraint   5
                   !-------------------------------------  
                   CASE(5)

                   
                   !-------------------------------------
                   
                   !-------------------------------------
                   !           Constraint   6
                   !-------------------------------------  
                   CASE(6)
                   
        
                   !-------------------------------------
                    

                
                END SELECT  

           END FUNCTION subgradient_g2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -              
           
           
        !----------------------------------------------------------------------------------   
        ! ------------------------------------------------------------------------------- |  
        ! |                                                                             | |
        ! |                            FUCTIONS NEEDED IN                               | |
        ! |                         MULTIOBJECTIVE OPTIMIZATION                         | |
        ! |                                                                             | |
        ! ------------------------------------------------------------------------------- |
        !----------------------------------------------------------------------------------

       

        !********************************************************************************
        !                                                                               |
        !   FUNCTION VALUES OF THE DC COMPONENTS H_1 AND H_2 IN MULTIOBJECTIVE MODEL    |
        !                                                                               |
        !********************************************************************************


           FUNCTION H1(AandB_y_all) RESULT(f)        
                !
                ! Calculates the function value of the DC component H_1=max{A_i,B_l} in multiobjective model when the componets used are 'AandB_y_all' 
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func0'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all ! the values of the components in H_1
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component H_1 at a point 'y'

                INTEGER :: i, ind                               ! help variable
           
                   !-------------------------------------
                     ind = 1 
                     DO i = 2, number_of_func0
                        IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                           ind = i
                        END IF 
                     END DO
     
                     f = AandB_y_all(ind)
                   !-------------------------------------                  

           END FUNCTION H1   
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION H2(fg2_y_all) RESULT(f)          
                !
                ! Calculates the function value of the DC component H_2 when the DC components f2_i and g2_l have values 'fg2_y_all'.
                !
                ! NOTICE: * The dimension of 'fg2_y_all' has to be 'number_of_func0'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_y_all    ! function values for all DC components f2_i and g2_l at point y
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                                     ! the function value of the DC component H_2 at a point 'y'
                
                INTEGER :: i
               
                    !-------------------------------------             
                      f = 0.0_dp
                      DO i = 1, number_of_func0
                          f = f + fg2_y_all(i)
                      END DO
                    !-------------------------------------

           END FUNCTION H2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           
           
        !********************************************************************************
        !                                                                               |
        !   SUBGRADIENTS OF THE DC COMPONENTS H_1 AND H_2 IN MULTIOBJECTIVE MODEL       |
        !                                                                               |
        !********************************************************************************          
           
            FUNCTION subgradient_H1(grad1_y_all, grad2_y_all, AandB_y_all, user_n) RESULT(grad)
                !
                ! Calculates the subgradient of the DC component H_1=max{A_i,B_l} at 'y' in multiobjective model when 
                ! the componets used at 'y' are 'AandB_y_all'
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func0'.
                !         * The dimension of 'grad1_y_all' and 'grad2_y_all' has to be 'number_of_func0:user_n'.             
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all       ! the values of the components A_i in H_1 at y
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad1_y_all     ! all subgradients for the DC components f1_i and g1_l at y             
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all     ! all subgradients for the DC components f2_i and g2_l at y     
                INTEGER, INTENT(IN) :: user_n                          ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad        ! the subgradient of the DC component H_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_f2     ! the sum of the subgradients of the DC components f2_i at 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_g2     ! the sum of the subgradients of the DC components g2_i at 'y'  
                
                INTEGER :: i, j, ind                              ! help variable

                   !-------------------------------------
                     ind = 1 
                     DO i = 2, number_of_func0
                        IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                           ind = i
                        END IF 
                     END DO
                   !-------------------------------------
                    
                   IF (ind < number_of_obj0+1) THEN  ! subgradient of A_i is calculated                  
                   !-------------------------------------
                     grad_g2 = 0.0_dp
                     DO i = number_of_obj0+1, number_of_func0
                       DO j = 1, user_n
                          grad_g2(j) = grad_g2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = 1, number_of_obj0
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_g2
                   !-------------------------------------
                   ELSE ! subgradient of B_l is calculated
                   !-------------------------------------
                     grad_f2 = 0.0_dp
                     DO i = 1, number_of_obj0
                       DO j = 1, user_n
                          grad_f2(j) = grad_f2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = number_of_obj0+1, number_of_func0
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_f2
                   !-------------------------------------
                   END IF                    
    
                
           END FUNCTION subgradient_H1  
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION subgradient_H2(grad2_y_all, user_n) RESULT(grad)              
                !
                ! Calculates a subgradient of the DC component H_2 at a point 'y' in multiobjective model when 
                ! subgradients of DC componets f2_i and g2_l at y are 'grad2_y_all'.
                !
                ! NOTICE: * The dimension of 'grad_fg2_y_all' has to be 'number_of_func0:user_n'.             
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all   ! all subgradients for the DC components f2_i and g2_l at y 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem

                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad         ! the subgradient of the DC component H_2 at a point 'y'

                INTEGER :: i , j

                   !-------------------------------------
                      grad = 0.0_dp
                      DO i = 1, number_of_func0
                        DO j = 1, user_n
                          grad(j) = grad(j) + grad2_y_all(i,j)
                        END DO  
                      END DO                  
                   !------------------------------------- 

           END FUNCTION subgradient_H2
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
           
 
           
        !********************************************************************************
        !                                                                               |
        !           GIVES THE INDEX OF THE COMPONENT YIELDING THE VALUE OF H_1          |
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION index_H1(AandB_y_all) RESULT(ind)
                !
                ! Calculates the index of the component yielding the value of H_1=max{A_i,B_l} when
                ! components used are 'AandB_y_all'
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func0'
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all     ! the values of components used in H_1

                !**************************** OTHER VARIABLES **************************************
                INTEGER :: ind   ! the index of the component yielding the value of H_1                
                INTEGER :: i     ! help variable
              
                     ind = 1 
                     DO i = 2, number_of_func0
                        IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                           ind = i
                        END IF 
                     END DO
                                
                
           END FUNCTION index_H1                   
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           
        
        !********************************************************************************
        !                                                                               |
        !    FUNCTION VALUE AND SUBGRADIENT OF COMPONENT A_i IN THE DC COMPONENT H_1    |
        !                                                                               |
        !********************************************************************************


           FUNCTION AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all, ind) RESULT(f)        
                !
                ! Calculates the value of the component A_i or B_l in the DC component H_1 at a point 'y',
                ! when the model is formed at 'x_k'. Variable 'ind' will identify the component function used.
                !
                ! NOTICE: * The dimension of 'fg1_y_all' and 'fg2_y_all' has to be 'number_of_func0'.
                !         * The dimension of 'fg1_k_all' and 'fg2_k_all' has to be 'number_of_func0'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_y_all   ! the values of the DC components f1_i and g1_l at point y
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_y_all   ! the values of the DC components f2_i and g2_l at point y
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_k_all   ! the values of the DC components f1_i and g1_l at point x_k
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_k_all   ! the values of the DC components f2_i and g2_l at point x_k
                INTEGER, INTENT(IN) :: ind                             ! index of the component used
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                                 ! the value of the component A_'ind' in the DC component H_1 at a point 'y'
                REAL(KIND=dp) :: sum_f2                            ! the sum of values of the all the DC components f2_i at a point 'y'
                REAL(KIND=dp) :: sum_g2                            ! the sum of values of the all the DC components g2_i at a point 'y'

                INTEGER :: i, j                                    ! help variable
                    
                    IF (ind < number_of_obj0+1) THEN  ! A_i is calculated 
                   !-------------------------------------
                   
                       sum_g2 = 0.0_dp
                       DO i = number_of_obj0+1, number_of_func0
                         sum_g2 = sum_g2 + fg2_y_all(i)
                       END DO
                    
                       f = fg1_y_all(ind) 
                       f = f -(fg1_k_all(ind)-fg2_k_all(ind))
     
                       DO j = 1, number_of_obj0
                           IF (ind /= j) THEN 
                              f = f + fg2_y_all(j)
                           END IF 
                       END DO    
                       f = f + sum_g2
                   !-------------------------------------  
                    ELSE  ! B_l is calculated
                   !-------------------------------------
                   
                       sum_f2 = 0.0_dp
                       DO i = 1, number_of_obj0
                         sum_f2 = sum_f2 + fg2_y_all(i)
                       END DO                  
                       
                       f =  fg1_y_all(ind) 
     
                       DO j = number_of_obj0+1, number_of_func0
                           IF (ind /= j) THEN 
                              f = f + fg2_y_all(j)
                           END IF 
                       END DO    
                       f = f + sum_f2                   
                   
                   !-------------------------------------
                    END IF                 
                
           END FUNCTION AorB_i 
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
           
           FUNCTION subgradient_AorB_i(grad1_y_all, grad2_y_all, ind, user_n) RESULT(grad)
                !
                ! Calculates a subgradient of the component A_i or B_l in the DC component H_1 at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'grad1_y_all' and 'grad2_y_all' has to be 'number_of_func0:user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad1_y_all    ! subgradients of DC components f1_i and g1_l at point y
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all    ! subgradients of DC components f2_i and g2_l at point y
                INTEGER, INTENT(IN) :: ind                                  ! the index of the component used
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem         
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad        ! the subgradient of the component 'ind' in the DC component H_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_f2     ! the sum of the subgradients of the DC components f2_i at 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_g2     ! the sum of the subgradients of the DC components g2_i at 'y'

                INTEGER :: i, j                                 ! help variable

                    
                   IF (ind < number_of_obj0+1) THEN  ! subgradient of A_i is calculated                  
                   !-------------------------------------
                     grad_g2 = 0.0_dp
                     DO i = number_of_obj0+1, number_of_func0
                       DO j = 1, user_n
                          grad_g2(j) = grad_g2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = 1, number_of_obj0
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_g2
                   !-------------------------------------
                   ELSE ! subgradient of B_l is calculated
                   !-------------------------------------
                     grad_f2 = 0.0_dp
                     DO i = 1, number_of_obj0
                       DO j = 1, user_n
                          grad_f2(j) = grad_f2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = number_of_obj0+1, number_of_func0
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_f2
                   !-------------------------------------
                   END IF

           END FUNCTION subgradient_AorB_i     

   
           
        !----------------------------------------------------------------------------------   
        ! ------------------------------------------------------------------------------- |  
        ! |                                                                             | |
        ! |                            FUCTIONS NEEDED IN                               | |
        ! |                         CONSTRAINED OPTIMIZATION                            | |
        ! |                                  WHEN                                       | |
        ! |                       SOPHISTICATED STRATEGY IS USED                        | |
        ! |                                                                             | |
        ! ------------------------------------------------------------------------------- |
        !----------------------------------------------------------------------------------

       

        !********************************************************************************************
        !                                                                                           |
        !   FUNCTION VALUES OF THE DC COMPONENTS H_1 AND H_2 IN CONSTRAINED MODEL WITH OBJECTIVE k  |
        !                                                                                           |
        !********************************************************************************************


           FUNCTION H1_sep(AandB_y_all) RESULT(f)        
                !
                ! Calculates the function value of the DC component H_1=max{A_i,B_l} in the constraint model with objective k when the componets used are 'AandB_y_all' 
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func0'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all ! the values of the components in H_1
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                              ! the function value of the DC component H_1 at a point 'y'

                INTEGER :: i, ind                               ! help variable
           
                   !-------------------------------------
                     ind = 1 
                     IF ( number_of_func_sep > 1 ) THEN
                         DO i = 2, number_of_func_sep
                            IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                               ind = i
                            END IF 
                         END DO
                     END IF
                     f = AandB_y_all(ind)
                   !-------------------------------------                  

           END FUNCTION H1_sep   
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION H2_sep(fg2_y_all) RESULT(f)          
                !
                ! Calculates the function value of the DC component H_2 in the constraint model with objective k when the DC components f2_i and g2_l have values 'fg2_y_all'.
                !
                ! NOTICE: * The dimension of 'fg2_y_all' has to be 'number_of_func0'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_y_all    ! function values for all DC components f2_i and g2_l at point y
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                                     ! the function value of the DC component H_2 at a point 'y'
                
                INTEGER :: i
               
                    !-------------------------------------             
                      f = 0.0_dp
                      DO i = 1, number_of_func_sep
                          f = f + fg2_y_all(i)
                      END DO
                    !-------------------------------------

           END FUNCTION H2_sep
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           
           
        !********************************************************************************************
        !                                                                                           |
        !   SUBGRADIENTS OF THE DC COMPONENTS H_1 AND H_2 IN MULTIOBJECTIVE MODEL WITH OBJECTIVE k  |      
        !                                                                                           |
        !********************************************************************************************          
           
            FUNCTION subgradient_H1_sep(grad1_y_all, grad2_y_all, AandB_y_all, user_n) RESULT(grad)
                !
                ! Calculates the subgradient of the DC component H_1=max{A_i,B_l} in the constraint model with objective k 
                ! at 'y' in multiobjective model when the componets used at 'y' are 'AandB_y_all'
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func_sep'.
                !         * The dimension of 'grad1_y_all' and 'grad2_y_all' has to be 'number_of_func_sep:user_n'.             
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all       ! the values of the components A_i in H_1 at y
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad1_y_all     ! all subgradients for the DC components f1_i and g1_l at y             
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all     ! all subgradients for the DC components f2_i and g2_l at y     
                INTEGER, INTENT(IN) :: user_n                          ! The dimension of the problem
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad        ! the subgradient of the DC component H_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_f2     ! the sum of the subgradients of the DC components f2_i at 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_g2     ! the sum of the subgradients of the DC components g2_i at 'y'  
                
                INTEGER :: i, j, ind                              ! help variable

                   !-------------------------------------
                     ind = 1 
                     IF ( number_of_func_sep > 1 ) THEN                  
                         DO i = 2, number_of_func_sep
                            IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                               ind = i
                            END IF 
                         END DO
                     END IF
                   !-------------------------------------
                    
                   IF (ind < number_of_obj_sep+1) THEN  ! subgradient of A_i is calculated                  
                   !-------------------------------------
                     grad_g2 = 0.0_dp
                     DO i = number_of_obj_sep+1, number_of_func_sep
                       DO j = 1, user_n
                          grad_g2(j) = grad_g2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = 1, number_of_obj_sep
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_g2
                   !-------------------------------------
                   ELSE ! subgradient of B_l is calculated
                   !-------------------------------------
                     grad_f2 = 0.0_dp
                     DO i = 1, number_of_obj_sep
                       DO j = 1, user_n
                          grad_f2(j) = grad_f2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = number_of_obj_sep+1, number_of_func_sep
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_f2
                   !-------------------------------------
                   END IF                    
    
                
           END FUNCTION subgradient_H1_sep  
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           FUNCTION subgradient_H2_sep(grad2_y_all, user_n) RESULT(grad)              
                !
                ! Calculates a subgradient of the DC component H_2 at a point 'y' in the constraint model with objective k when 
                ! subgradients of DC componets f2_i and g2_l at y are 'grad2_y_all'.
                !
                ! NOTICE: * The dimension of 'grad_fg2_y_all' has to be 'number_of_func0:user_n'.             
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all   ! all subgradients for the DC components f2_i and g2_l at y 
                INTEGER, INTENT(IN) :: user_n                   ! The dimension of the problem

                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad         ! the subgradient of the DC component H_2 at a point 'y'

                INTEGER :: i , j

                   !-------------------------------------
                      grad = 0.0_dp
                      DO i = 1, number_of_func_sep
                        DO j = 1, user_n
                          grad(j) = grad(j) + grad2_y_all(i,j)
                        END DO  
                      END DO                  
                   !------------------------------------- 

           END FUNCTION subgradient_H2_sep
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -          
           
 
           
        !********************************************************************************
        !                                                                               |
        !  GIVES THE INDEX OF THE COMPONENT YIELDING THE VALUE OF H_1 WITH OBJECTIVE k  |        
        !                                                                               |       
        !********************************************************************************       
        
           FUNCTION index_H1_sep(AandB_y_all) RESULT(ind)
                !
                ! Calculates the index of the component yielding the value of H_1=max{A_i,B_l} when
                ! components used are 'AandB_y_all'
                !
                ! NOTICE: * The dimension of 'AandB_y_all' has to be 'number_of_func0'
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AandB_y_all     ! the values of components used in H_1

                !**************************** OTHER VARIABLES **************************************
                INTEGER :: ind   ! the index of the component yielding the value of H_1                
                INTEGER :: i     ! help variable
              
                     ind = 1 
                     DO i = 2, number_of_func_sep
                        IF (AandB_y_all(i) > AandB_y_all(ind)) THEN 
                           ind = i
                        END IF 
                     END DO
                                
                
           END FUNCTION index_H1_sep                   
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           
           
        
        !***********************************************************************************************
        !                                                                                              |
        !    FUNCTION VALUE AND SUBGRADIENT OF COMPONENT A_i IN THE DC COMPONENT H_1 WITH OBJECTIVE k  |   
        !                                                                                              |
        !***********************************************************************************************


           FUNCTION AorB_i_sep(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all, ind) RESULT(f)        
                !
                ! Calculates the value of the component A_i or B_l in the DC component H_1 at a point 'y' in the constraint model with objective k, 
                ! when the model is formed at 'x_k'. Variable 'ind' will identify the component function used.
                !
                ! NOTICE: * The dimension of 'fg1_y_all' and 'fg2_y_all' has to be 'number_of_func_sep'.
                !         * The dimension of 'fg1_k_all' and 'fg2_k_all' has to be 'number_of_func_sep'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_y_all   ! the values of the DC components f1_i and g1_l at point y
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_y_all   ! the values of the DC components f2_i and g2_l at point y
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_k_all   ! the values of the DC components f1_i and g1_l at point x_k
                REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_k_all   ! the values of the DC components f2_i and g2_l at point x_k
                INTEGER, INTENT(IN) :: ind                             ! index of the component used
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp) :: f                                 ! the value of the component A_'ind' in the DC component H_1 at a point 'y'
                REAL(KIND=dp) :: sum_f2                            ! the sum of values of the all the DC components f2_i at a point 'y'
                REAL(KIND=dp) :: sum_g2                            ! the sum of values of the all the DC components g2_i at a point 'y'

                INTEGER :: i, j                                    ! help variable
                    
                    IF (ind < number_of_obj_sep+1) THEN  ! A_i is calculated 
                   !-------------------------------------
                   
                       sum_g2 = 0.0_dp
                       DO i = number_of_obj_sep+1, number_of_func_sep
                         sum_g2 = sum_g2 + fg2_y_all(i)
                       END DO
                    
                       f = fg1_y_all(ind) 
                       f = f -(fg1_k_all(ind)-fg2_k_all(ind))
     
                       DO j = 1, number_of_obj_sep
                           IF (ind /= j) THEN 
                              f = f + fg2_y_all(j)
                           END IF 
                       END DO    
                       f = f + sum_g2
                   !-------------------------------------  
                    ELSE  ! B_l is calculated
                   !-------------------------------------
                   
                       sum_f2 = 0.0_dp
                       DO i = 1, number_of_obj_sep
                         sum_f2 = sum_f2 + fg2_y_all(i)
                       END DO                  
                       
                       f =  fg1_y_all(ind) 
     
                       DO j = number_of_obj_sep+1, number_of_func_sep
                           IF (ind /= j) THEN 
                              f = f + fg2_y_all(j)
                           END IF 
                       END DO    
                       f = f + sum_f2                   
                   
                   !-------------------------------------
                    END IF                 
                
           END FUNCTION AorB_i_sep 
           
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
           
           FUNCTION subgradient_AorB_i_sep(grad1_y_all, grad2_y_all, ind, user_n) RESULT(grad)
                !
                ! Calculates a subgradient of the component A_i or B_l in the DC component H_1 for objective ind at a point 'y'.
                ! Variable 'ind' identifies the objective function used.
                !
                ! NOTICE: * The dimension of 'grad1_y_all' and 'grad2_y_all' has to be 'number_of_func0:user_n'.
                !
                IMPLICIT NONE
                !**************************** NEEDED FROM USER *************************************
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad1_y_all    ! subgradients of DC components f1_i and g1_l at point y
                REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: grad2_y_all    ! subgradients of DC components f2_i and g2_l at point y
                INTEGER, INTENT(IN) :: ind                                  ! the index of the component used
                INTEGER, INTENT(IN) :: user_n                               ! The dimension of the problem         
                
                !**************************** OTHER VARIABLES **************************************
                REAL(KIND=dp), DIMENSION(user_n) :: grad        ! the subgradient of the component 'ind' in the DC component H_1 at a point 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_f2     ! the sum of the subgradients of the DC components f2_i at 'y'
                REAL(KIND=dp), DIMENSION(user_n) :: grad_g2     ! the sum of the subgradients of the DC components g2_i at 'y'

                INTEGER :: i, j                                 ! help variable

                    
                   IF (ind < number_of_obj_sep+1) THEN  ! subgradient of A_i is calculated                  
                   !-------------------------------------
                     grad_g2 = 0.0_dp
                     DO i = number_of_obj_sep+1, number_of_func_sep
                       DO j = 1, user_n
                          grad_g2(j) = grad_g2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = 1, number_of_obj_sep
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_g2
                   !-------------------------------------
                   ELSE ! subgradient of B_l is calculated
                   !-------------------------------------
                     grad_f2 = 0.0_dp
                     DO i = 1, number_of_obj_sep
                       DO j = 1, user_n
                          grad_f2(j) = grad_f2(j) + grad2_y_all(i,j)
                       END DO
                     END DO
                       
                     grad = 0.0_dp
                     DO i = 1, user_n
                        grad(i) = grad1_y_all(ind,i)
                     END DO     

                     DO j = number_of_obj_sep+1, number_of_func_sep
                        IF ( ind /= j ) THEN 
                           DO i = 1, user_n
                             grad(i) = grad(i) + grad2_y_all(j,i)
                           END DO   
                        END IF
                     END DO
                     
                     grad = grad + grad_f2
                   !-------------------------------------
                   END IF

           END FUNCTION subgradient_AorB_i_sep                        
           
      END MODULE functions     