        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                IMB-DC - THE INTERACTIVE MULTIBUNDLE METHOD FOR                   | |
        !| |                     CONSTRAINED NONSMOOTH DC OPTIMIZATION                        | | 
        !| |                                 (version 1)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                       by Kaisa Joki (last modified September 2020)               | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |     The software is free for academic teaching and research purposes but I       | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    Utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver.          |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !|                                                                                      |
        !|                                                                                      |
        !|   Code include:                                                                      |
        !|                                                                                      |
        !|   timbdc.f95         - Main program for IMB-DC (this file)                           |
        !|   constants.f95      - Double precision (also some parameters)                       |
        !|   bundle1.f95        - Bundle of DC component f_1                                    |
        !|   bundle2.f95        - Bundle of DC component f_2                                    |
        !|   functions.f95      - User-specified DC components f_1 and f_2 together with        |
        !|                        subgradients of DC components. Contains also user-specified   |
        !|                        initial values for parameters                                 |
        !|   imbdc.f95          - IMB-DC method with actual values for parameters               |
        !|                                                                                      |
        !|   plqdf1.f           - Quadratic solver by Ladislav Luksan                           |
        !|                                                                                      |
        !|   Makefile           - Makefile                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software MODIFY   timbdc.f95, imbdc.f95 and functions.f95   as needed   |
        !|                                                                                      |
        !|                                                                                      |
        !|   Reference:                                                                         |
        !|                                                                                      |
        !|   [1] Outi Montonen and Kaisa Joki: "Interactive multibundle method for constrained  |
        !|       multiobjective DC optimization." (manuscript)                                  |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       


    PROGRAM timbdc
      
         USE constants, ONLY : dp   ! double precision (i.e. accuracy)
         USE functions              ! INFORMATION from the USER
         USE bundle1                ! The BUNDLE of the DC component f_1
         USE bundle2                ! The BUNDLE of the DC component f_2
         USE imbdc                  ! The new method for nonsmooth multiobjective DC optimization
         
        IMPLICIT NONE   
        
        ! 'user_n' is the number of variables in the problem (USER specifies this in MODULE functions.f95)
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: x_solution         ! The solution obtained to the problem
        
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f_solution         ! The objective function values at the solution 'x_solution'
        
        REAL(KIND=dp) :: scaling_const                                 ! The constant used in the scaling of the objective
                
        INTEGER, DIMENSION(:,:), ALLOCATABLE  :: counter        ! Contains the values of different counteres:
                                                                !   counter(1,1) = iter_counter:         the number of 'main iterations' executed during the new method 
                                                                !   counter(2,j) = subprob_counter:      the number of subproblems solved during all separate descent direction determinations for f_j 
                                                                !   counter(3,j) = f_counter:            the number of function values evaluated for the objective f_j (the same value holds for the DC components also)
                                                                !   counter(4,j) = subgrad1_counter:     the number of subgradients calculated for the first DC component of f_j (i.e. f_1)
                                                                !   counter(5,j) = subgrad2_counter:     the number of subgradients calculated for the second DC component of f_j (i.e f_2)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(6,1) = stop_cond_counter:    the number of times 'Escape procedure' is executed for the separate decsent direction determination
                                                                !   counter(7,1) = escape_f_counter:     the number of function values evaluated for each objective f_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                !   counter(8,1) = escape_sub_counter:   the number of subgradients caluculated for each objective f_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(9,1) = stop_cond_counter:    the number of times 'Escape procedure' is executed for the improvement function
                                                                !   counter(10,1) = escape_f_counter:    the number of function values evaluated for each objective f_i (in 'Escape procedure' for the improvement function)
                                                                !   counter(11,1) = escape_sub_counter:  the number of subgradients caluculated for each objective f_i (in 'Escape procedure' for the improvement function)
                                                                !--------------------------------------------------------------------------------------------------------------------------                 
                                                                !   counter(12,1) = f_counter_line:      the number of function values evaluated for the objective f_j during line search (the same value holds for the DC components also)
                                                                !--------------------------------------------------------------------------------------------------------------------------  
 

        INTEGER, DIMENSION(:,:), ALLOCATABLE  :: counter_g        ! Contains the values of different counteres:
                                                                  !   counter_g(1,j) = g_counter:              the number of function values evaluated for the constraint g_j (the same value holds for the DC components also)
                                                                  !   counter_g(2,j) = subgrad1_g_counter:     the number of subgradients calculated for the first DC component of g_j (i.e. g_1)
                                                                  !   counter_g(3,j) = subgrad2_g_counter:     the number of subgradients calculated for the second DC component of g_j (i.e g_2)
                                                                  !   counter_g(4,1) = g_counter_line:         the number of function values evaluated for the constraint g_j during line search (the same value holds for the DC components also)
                                                                  !--------------------------------------------------------------------------------------------------------------------------     
                                                                  !   counter_g(5,i) = escape_g_counter_d:       the number of function values evaluated for each constraint g_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                  !   counter_g(6,i) = escape_sub_g_counter_d:   the number of subgradients caluculated for each constraint g_i (in 'Escape procedure' for the separate decsent direction determination)
                                                                  !--------------------------------------------------------------------------------------------------------------------------     
                                                                  !   counter_g(7,1) = escape_f_counter_H:    the number of function values evaluated for each constraint g_i (in 'Escape procedure' for the improvement function)
                                                                  !   counter_g(8,1) = escape_sub_counter_H:  the number of subgradients caluculated for each constraint g_i (in 'Escape procedure' for the improvement function)
                                                                  !-----------------------------------------------------------------------------------------------------------------------
                                                                  
        INTEGER :: iprint                   ! Variable that specifies print option (specified by USER): 
                                            !   iprint = 0 : print is suppressed
                                            !   iprint = 1 : basic print of final result 
                                            !   iprint = -1: basic print of final result (without the solution vector)
                                            !   iprint = 2 : extended print of final result 
                                            !   iprint = -2: extended print of final result (without the solution vector)
                                            !   iprint = 3 : basic print of intermediate results and extended print of final results
                                            !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                            !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                            !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                            !   iprint = 5 : prints each step of the IMB-DC algorithm (i.e. everything is printed step by step) (NOT recommended)
                                            !
                                            ! If 'iprint' <= -5 .OR. 'iprint' >= 6 then DEFAULT value 'iprint'=1 is used    
                                            
        
        INTEGER :: mit                      ! The maximum number of 'main iterations' (specified by USER).
                                            ! If 'mit' <=0 then DEFAULT value 'mit'=1000 is used

        INTEGER :: mrounds                  ! The maximum number of rounds during descent direction determination
                                            ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used

        INTEGER :: mrounds_escape           ! The maximum number of rounds during one 'Escape procedure' (specified by USER).
                                            ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_escape'=5000 is used
                                            
        INTEGER :: termination              ! The reason for termination in the new method:
                                            !   1 - the stopping condition is satisfied (i.e. approximate Clarke stationarity)
                                            !   2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < epsilon)
                                            !   3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                            !   4 - the maximum number of 'main iterations' is executed 
                                            !   5 - the maximum number 'mrounds_clarke' of rounds is executed in one 'Clarke stationary' algorithm

        LOGICAL :: scale_func               ! If .TRUE. then the scaling procedure is performed at the beginning of the new method
        
        LOGICAL :: agg_used                 ! If .TRUE. then aggregation is used during separate descent direction determination for the objective f_i
             
        LOGICAL :: t_adjust                 ! NEEDS TO BE ALWAYS .FALSE. (There is no adjustment for t available)                       
        LOGICAL :: use_reset_b1             ! If .TRUE. then bundle B1 is reset after serious step and escape procedure. Otherwise .FALSE.        
        LOGICAL :: use_reset_b2             ! If .TRUE. then bundle B2 is reset after serious step and escape procedure. Otherwise .FALSE.        
        
        LOGICAL :: use_const_in_sep_d       ! If .TRUE. then constraints are used in the individual search direction problem     
        
        
        !--------------------------------------------------------------------------------------------------------------------------------
        !      *****   WHAT IS NEEDED TO DEFINE PROBLEM !    *****
        !--------------------------------------------------------------------------------------------------------------------------------
        INTEGER :: number_of_obj                           ! The number of objectives in the multiobjective problem
                                                           ! IF 'number_of_obj'=1 THEN we have a single objective case.
        INTEGER :: number_of_const                         ! The number of constraintss in the multiobjective problem
                                                           ! IF 'number_of_const'=0 THEN we have an unconstrained problem    
        
        INTEGER :: number_of_func                          ! The sum of objective and constraint functions                                                               

        !---------------------------------------------------------------------------------------------------------------------------------
        ! DC component lists where the first 'number_of_obj' places tell the DC components f1 for the objectives
        ! and the latter 'number_of_const' places tell the DC components g1 for the constraints         
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem10      ! First is DC components f1 of objective functions used and then DC components g1 of constraint functions used
         
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem20      ! DC components f2 and g2 of objective and constraint functions used 
        
        INTEGER :: startpoint                                ! Starting point used

        !--------------------------------------------------------------------------------------------------------------------------------
      
        INTEGER :: dim_loop                      ! Defines the dimension(s) in the considered problem 
                                                 !    1 - the dimension is 2
                                                 !    2 - the dimension is 4
                                                 !    3 - the dimension is 10
                                                 !    4 - the dimensions are 10, 50, 100, 250 and 500
        
        INTEGER :: ind_start
        INTEGER :: ind_finish
        
        INTEGER :: tot_num                       ! The total number of problems        
        
        INTEGER :: user_n                        ! The number of variables in the problem
        INTEGER :: max_null_step                 ! The maximum number of consecutive null steps
        
        REAL(KIND=dp) :: time                    ! The cpu time
        CHARACTER(LEN=80) :: outfile1 
               
        INTEGER :: i_user             ! The maximum number of times the user can steer the algorithm
        
        INTEGER :: teht_dim           ! The dimensions of in varying dimension test problems
                                      !  0  -  10, 50, 100, 250, 500
                                      !  1  -  10, 50, 100, 250 
                                      !  2  -  10, 50, 100 
 
        !--------------------------------------------------------------------------------------------------------------------------------
               
        INTEGER :: i, j, k, ind                       ! help variables

        outfile1 = 'results_IMB-DC.txt'
        
        OPEN(46,file=outfile1)    
        
        ! Values from scale_func, agg_used, t_adjust, use_reset_b1, use_reset_b2, max_null_step and use_const_in_sep_d are set
        CALL parameters(1, scale_func, agg_used, t_adjust, use_reset_b1, use_reset_b2, max_null_step, & 
                          & use_const_in_sep_d )   
        
      ! Some parameters of the method 
        mrounds = 1000              ! maximum number of rounds during separate descent direction determination
        mit = 1000                  ! maximum number of 'main iterations'
        mrounds_escape = 1000       ! maximum number of rounds in Clarke stationary algorithm 
          
        iprint = 2                  ! basic print of intermediate results and extended print of final results
        
        i_user = 0                  ! The maximum number of times the user can steer the algorithm
                                    ! If i_user=0 then the user cannot steer the algorithm 
    
        teht_dim = 0
                       ! If teht_dim=0: values for n are 10, 50, 100, 250, 500s
                       ! If teht_dim=1: values for n are 10, 50, 100, 250
                       ! If teht_dim=2: values for n are 10, 50, 100
                       ! If teht_dim=3: values for n are 10, 50
                       
        tot_num = 0   ! The total number of solved problems so far
        
        
        WRITE(46,*)
        WRITE(46,*) 'Parameters ind:', ind
        WRITE(46,*)
        WRITE(46,*) 'mrounds', mrounds
        WRITE(46,*) 'mit', mit
        WRITE(46,*) 'mrounds_escape', mrounds_escape
        WRITE(46,*)
        WRITE(46,*) 'max_null_step', max_null_step
        WRITE(46,*)
        WRITE(46,*) 'scale_func', scale_func
        WRITE(46,*) 'agg_used', agg_used
        WRITE(46,*) 't_adjust', t_adjust
        WRITE(46,*) 'use_reset_b1', use_reset_b1
        WRITE(46,*) 'use_reset_b2', use_reset_b2
        WRITE(46,*)
        WRITE(46,*) 'user_size_b1=B1', user_size_b1
        WRITE(46,*) 'user_size_b2=B2', user_size_b2
        WRITE(46,*) 'user_escape_b=B_escape', user_size
        WRITE(46,*) 
        WRITE(46,*) 'user_m=m1', user_m
        WRITE(46,*) 'user_r_dec=r', user_r_dec
        WRITE(46,*) 'user_r_inc=R', user_r_inc
        WRITE(46,*) 'user_c=c', user_c
        WRITE(46,*) 
        WRITE(46,*) 'user_eps=theta', user_eps        
        WRITE(46,*) 
        WRITE(46,*) 'For escape procedure inside the descent direction determination:'
        WRITE(46,*) 'user_crit_tol_d=delta1', user_crit_tol_d
        WRITE(46,*) 'user_m_escape=m2', user_m_escape_d
        WRITE(46,*) 'user_step_tol=eps1', user_step_tol_d        
        WRITE(46,*) 
        WRITE(46,*) 'For escape procedure with the improvement function:'
        WRITE(46,*) 'user_crit_tol=delta2', user_crit_tol_H        
        WRITE(46,*) 'user_escape_tol=gamma', user_escape_tol_H
        WRITE(46,*) 'user_m_escape=m2', user_m_escape_H
        WRITE(46,*) 'user_step_tol=eps', user_step_tol_H
        WRITE(46,*) 
        WRITE(46,*) 'In constrained problems constraints are used in the individual search direction problem', &
                     & use_const_in_sep_d
        WRITE(46,*)
            
           
        WRITE(46,*)  ' Problem ', ' user_n ', ' | ', ' obj1 ', ' obj2 ', ' obj3 ', 'g ', ' | ', &
                      & ' total_f1 ', ' total_f2 ', ' total_f3 ', 'total_g', &
                      & ' total_sub_f1 ', ' total_sub_f2 ', ' total_sub_f3 ', 'total_sub_g ',' | ', &
                      & ' n_f_obj1 ', ' n_f_obj2 ', ' n_f_obj3 ', ' n_g ', ' | ',& 
                      & ' n_f_line_obj1 ', ' n_f_line_obj2 ', ' n_f_line_obj3 ', ' n_g_line ',' | ', & 
                      & ' n_sub1_obj1 ', ' n_sub1_obj2 ', ' n_sub1_obj3 ', ' n_sub1_g ',' | ', &
                      & ' n_sub2_obj1 ', ' n_sub2_obj2 ', ' n_sub2_obj3 ', ' n_sub2_g ', ' | ',&                
                      & ' n_f_escape_d1 ', ' n_sub_escape_d1 ', ' times_escape_d1 ', ' | ',&  
                      & ' n_f_escape_d2 ', ' n_sub_escape_d2 ', ' times_escape_d2 ', ' | ',&  
                      & ' n_f_escape_d3 ', ' n_sub_escape_d3 ', ' times_escape_d3 ', ' | ',&  
                      & ' n_f_escape_H ', ' n_sub_escape_H ', ' times_escape_H ', ' | ',' time ', ' n_iter ', ' termination '                     

                      
        
      ! Different test problems are looked through
        DO i = 1,21
        
            SELECT CASE(i)

              !-------------------------------------
              !           Problem   1
              !-------------------------------------           
               CASE(1)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                  
                ! The functions used 
                  problem10 = (/ 2,6 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 1
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)

              !-------------------------------------
              !           Problem   2
              !-------------------------------------           
               CASE(2)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 2
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   3
              !-------------------------------------           
               CASE(3)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 3
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   4
              !-------------------------------------           
               CASE(4)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 6,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 4
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   5
              !-------------------------------------           
               CASE(5)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,9 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 5
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   6
              !-------------------------------------           
               CASE(6)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,13 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 6
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 3
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   7
              !-------------------------------------           
               CASE(7)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 12,13 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 7
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 3
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   8
              !-------------------------------------           
               CASE(8)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 8
                  
                ! The dimension(s) of the problem looked through 
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8
                  END SELECT
                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   9
              !-------------------------------------           
               CASE(9)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,12 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 9
                  
                ! The dimension(s) of the problem looked through 
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8
                  END SELECT
                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   10
              !-------------------------------------           
               CASE(10)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 14,15 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 10
                  
                ! The dimension(s) of the problem looked through 
                  !SELECT CASE(teht_dim)
                  !teht_dim = 1
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8                    
                  END SELECT
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   11
              !-------------------------------------           
               CASE(11)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 11
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)   


              !-------------------------------------
              !           Problem   12
              !-------------------------------------           
               CASE(12)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,4,9 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 12
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   13
              !-------------------------------------           
               CASE(13)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,12 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 13
                  
                ! The dimension(s) of the problem looked through 
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8                    
                  END SELECT
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   14
              !-------------------------------------           
               CASE(14)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,16 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 14
                  
                ! The dimension(s) of the problem looked through 
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8                    
                  END SELECT
          
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)
                          
              !-------------------------------------
              !           Problem   15
              !-------------------------------------           
               CASE(15)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,14,15 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 15
                  
                ! The dimension(s) of the problem looked through 
                  SELECT CASE(teht_dim)
                    CASE(0)
                      dim_loop = 4
                    CASE(1)
                      dim_loop = 6
                    CASE(2) 
                      dim_loop = 7
                    CASE(3) 
                      dim_loop = 8                    
                  END SELECT
                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   16
              !-------------------------------------           
               CASE(16)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7,1 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 16
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   17
              !-------------------------------------           
               CASE(17)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,9,2 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 17
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   18
              !-------------------------------------           
               CASE(18)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,12,3 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 18
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   19
              !-------------------------------------           
               CASE(19)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,7,1 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 19
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   20
              !-------------------------------------           
               CASE(20)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,4,9,2 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 20
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)


              !-------------------------------------
              !           Problem   21
              !-------------------------------------           
               CASE(21)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,16,3 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 21
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)
                  
 
              !-------------------------------------
              !           Problem   22
              !-------------------------------------           
               CASE(22)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,4 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 22
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)

                  use_const_in_sep_d = .FALSE.
          

              !-------------------------------------
              !           Problem   23
              !-------------------------------------           
               CASE(23)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,4 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 22
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)

                 use_const_in_sep_d = .TRUE.

              !-------------------------------------
              !           Problem   24
              !-------------------------------------           
               CASE(24)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                  
                ! The functions used 
                  problem10 = (/ 2,6 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 22
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                                  
                 CALL allocate_prob_data(number_of_obj,number_of_const,number_of_func, &
                          & problem10, problem20)
     

            END SELECT      

            ALLOCATE(f_solution(number_of_func),counter(12,number_of_obj))
            ALLOCATE(counter_g(8,number_of_const))

            SELECT CASE(dim_loop)   
 
                CASE(1) ! (n=2)
                    ind_start = 1
                    ind_finish = 1
                    
                CASE(2) ! (n=4)
                    ind_start = 2
                    ind_finish = 2
                    
                CASE(3) ! (n=10)
                    ind_start = 3
                    ind_finish = 3
                    
                CASE(4) ! (n= 10, 50, 100, 250, 500)
                    ind_start = 3
                    ind_finish = 7 
                    
                CASE(5) ! (n= 250)
                    ind_start = 3
                    ind_finish = 6

                CASE(6) ! (n= 10, 50, 100, 250)
                    ind_start = 3
                    ind_finish = 6

                CASE(7) ! (n= 10, 50, 100)
                    ind_start = 3
                    ind_finish = 5
                    
                CASE(8) ! (n= 10, 50)
                    ind_start = 3
                    ind_finish = 4

            END SELECT          
            
            DO j = ind_start, ind_finish
            
                SELECT CASE(j)
                   
                    CASE(1)
                       user_n = 2
                    
                    CASE (2)
                       user_n = 4
                    
                    CASE (3)
                       user_n = 10
                    
                    CASE (4)
                       user_n = 50
                    
                    CASE (5)
                       user_n = 100
                    
                    CASE (6)
                       user_n = 250
                    
                    CASE (7)
                       user_n = 500
                           
                END SELECT
                
                SELECT CASE(i)
                
                   CASE(16)
                      scaling_const = (10.0_dp)**(1)

                   CASE(17)
                      scaling_const = (10.0_dp)**(1)
                   
                   CASE(18)
                      SELECT CASE(user_n)
                        CASE(10)
                            scaling_const = (10.0_dp)**(1)
                        CASE(50)
                            scaling_const = (10.0_dp)**(1)
                        CASE(100)
                            scaling_const = (10.0_dp)**(1)
                        CASE(250) 
                            scaling_const = (10.0_dp)**(1)
                        CASE(500)
                            scaling_const = (10.0_dp)**(1)                        
                        
                      END SELECT
                   
                   CASE(19)
                       scaling_const = (10.0_dp)**(1)
                   
                   CASE(20)
                       scaling_const = (10.0_dp)**(1)

                   CASE(21)
                      SELECT CASE(user_n)
                        CASE(10)
                            scaling_const = (10.0_dp)**(1)
                        CASE(50)
                            scaling_const = (10.0_dp)**(0)
                        CASE(100)
                            scaling_const = (10.0_dp)**(0)
                        CASE(250) 
                            scaling_const = (10.0_dp)**(1)
                        CASE(500)
                            scaling_const = (10.0_dp)**(0)                      
                      END SELECT              
                
                    CASE(22)
                       scaling_const = (10.0_dp)**(1)    
                    
                    CASE(23)
                       scaling_const = (10.0_dp)**(1)    
                       
                END SELECT
                
              
                ALLOCATE(x_solution(user_n))
               

                IF (iprint >= 1 .OR. iprint<=-1) THEN 
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) '** START ** START ** START ** START ** START ** START ** START **'  
                WRITE(*,*) '------------------------------------------------------------------'         
                WRITE(*,*) ' '
                WRITE(*,*) ' ', 'Problem', i, 'user_n', user_n
                WRITE(*,*) ' '
                WRITE(*,*) ' ', 'Objectives:', problem10
                WRITE(*,*) ' '
                END IF 

                CALL new_algorithm( x_solution, f_solution, mit, mrounds, &
                                    & mrounds_escape, termination, counter,  &
                                    & iprint, agg_used, scale_func, t_adjust, use_reset_b1, use_reset_b2, &
                                    & max_null_step, & 
                                    & user_n, startpoint, number_of_obj, number_of_const, number_of_func, &
                                    & time, i_user, counter_g, use_const_in_sep_d, scaling_const)
                  
                  
                IF (iprint >= 1 .OR. iprint<=-1) THEN 
                WRITE(*,*) ' '
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) '** END ** END ** END ** END ** END ** END ** END ** END ** END **'  
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) ' '      
                END IF 
                
               
                tot_num = tot_num + 1     ! One new problem solved
                
                IF (number_of_const==0) THEN
                    IF (number_of_obj == 1) THEN 
                        WRITE(46,*)   i, user_n , ' | ', f_solution(1), ' . ', ' . ', ' . ', ' | ',  &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1) , &
                              & ' . ', &                              
                              & ' . ', &
                              & ' . ', &
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & ' . ', &
                              & ' . ', &
                              & ' . ', ' | ', &  
                              & counter(3,1), ' . ', ' . ', ' . ',' | ', & 
                              & counter(12,1), ' . ', ' . ', ' . ',' | ', & 
                              & counter(4,1), ' . ', ' . ', ' . ',' | ', & 
                              & counter(5,1), ' . ', ' . ', ' . ',' | ', &           
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &   
                              ' . ', ' . ', ' . ', ' | ', &                     
                              ' . ', ' . ', ' . ', ' | ', &                                                   
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination
                       
                    ELSE IF (number_of_obj == 2) THEN
                        WRITE(46,*)   i, user_n ,' | ', f_solution(1), f_solution(2), ' . ', ' . ',   ' | ', &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1), &
                              & counter(3,2) + counter(12,2) + counter(7,2) + counter(10,1), &
                              & ' . ', &                              
                              & ' . ', &                              
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & Max(counter(4,2), counter(5,2)) + counter(8,2) + counter(11,1), &
                              & ' . ' ,  &
                              & ' . ' , ' | ', &
                              & counter(3,1), counter(3,2), ' . ', ' . ',' | ', & 
                              & counter(12,1), counter(12,2), ' . ', ' . ',' | ', & 
                              & counter(4,1), counter(4,2), ' . ', ' . ',' | ', & 
                              & counter(5,1), counter(5,2), ' . ', ' . ',' | ', &            
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &                     
                              & counter(7,2), counter(8,2), counter(6,2), ' | ', &                     
                              ' . ', ' . ', ' . ', ' | ', &                     
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination
                  
                    ELSE IF (number_of_obj == 3) THEN
                        WRITE(46,*)   i, user_n ,' | ', f_solution(1), f_solution(2), f_solution(3), ' . ',  ' | ',  &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1), &
                              & counter(3,2) + counter(12,2) + counter(7,2) + counter(10,1), &
                              & counter(3,3) + counter(12,3) + counter(7,3) + counter(10,1), ' . ',  &
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & Max(counter(4,2), counter(5,2)) + counter(8,2) + counter(11,1), &
                              & Max(counter(4,3), counter(5,3)) + counter(8,3) + counter(11,1), ' . ', ' | ', & 
                              & counter(3,1), counter(3,2), counter(3,3), ' . ',' | ', & 
                              & counter(12,1), counter(12,2), counter(12,3), ' . ',' | ', & 
                              & counter(4,1), counter(4,2), counter(4,3), ' . ', ' | ', & 
                              & counter(5,1), counter(5,2), counter(5,3), ' . ', ' | ', &             
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &                     
                              & counter(7,2), counter(8,2), counter(6,2), ' | ', &                     
                              & counter(7,3), counter(8,3), counter(6,3), ' | ', &                     
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination 

                              
                    END IF
                ELSE                    
                    IF (number_of_obj == 1) THEN 
                        WRITE(46,*)   i, user_n , ' | ', f_solution(1), ' . ', ' . ', f_solution(number_of_func), ' | ',  &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1) , &
                              & ' . ', &                              
                              & ' . ', &
                              & counter_g(1,1)+counter_g(4,1)+counter_g(5,1)+counter_g(7,1),   &
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & ' . ', &
                              & ' . ', &
                              & Max(counter_g(2,1),counter_g(3,1))+counter_g(6,1)+counter_g(8,1),   &
                              & ' | ', &  
                              & counter(3,1), ' . ', ' . ', counter_g(1,1) ,' | ', & 
                              & counter(12,1), ' . ', ' . ', counter_g(4,1) ,' | ', & 
                              & counter(4,1), ' . ', ' . ', counter_g(2,1) ,' | ', & 
                              & counter(5,1), ' . ', ' . ', counter_g(3,1) ,' | ', &           
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &   
                              ' . ', ' . ', ' . ', ' | ', &                     
                              ' . ', ' . ', ' . ', ' | ', &                                                   
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination
                       
                    ELSE IF (number_of_obj == 2) THEN
                        WRITE(46,*)   i, user_n ,' | ', f_solution(1), f_solution(2), ' . ', &
                              & f_solution(number_of_func),  ' | ',  &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1), &
                              & counter(3,2) + counter(12,2) + counter(7,2) + counter(10,1), &
                              & ' . ', & 
                              & counter_g(1,1)+counter_g(4,1)+counter_g(5,1)+counter_g(7,1),   &                              
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & Max(counter(4,2), counter(5,2)) + counter(8,2) + counter(11,1), &
                              & ' . ' , &
                              & Max(counter_g(2,1),counter_g(3,1))+counter_g(6,1)+counter_g(8,1),   &
                              & ' | ', &  
                              & counter(3,1), counter(3,2), ' . ', counter_g(1,1) ,' | ', & 
                              & counter(12,1), counter(12,2), ' . ',counter_g(4,1) ,' | ', & 
                              & counter(4,1), counter(4,2), ' . ',counter_g(2,1) ,' | ', & 
                              & counter(5,1), counter(5,2), ' . ',counter_g(3,1) ,' | ', &            
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &                     
                              & counter(7,2), counter(8,2), counter(6,2), ' | ', &                     
                              ' . ', ' . ', ' . ', ' | ', &                     
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination
                  
                    ELSE IF (number_of_obj == 3) THEN
                        WRITE(46,*)   i, user_n ,' | ', f_solution(1), f_solution(2), f_solution(3),  &
                              & f_solution(number_of_func),  ' | ',  &
                              & counter(3,1) + counter(12,1) + counter(7,1) + counter(10,1), &
                              & counter(3,2) + counter(12,2) + counter(7,2) + counter(10,1), &
                              & counter(3,3) + counter(12,3) + counter(7,3) + counter(10,1), &
                              & counter_g(1,1)+counter_g(4,1)+counter_g(5,1)+counter_g(7,1),   &                              
                              & Max(counter(4,1), counter(5,1)) + counter(8,1) + counter(11,1), &
                              & Max(counter(4,2), counter(5,2)) + counter(8,2) + counter(11,1), &
                              & Max(counter(4,3), counter(5,3)) + counter(8,3) + counter(11,1), & 
                              & Max(counter_g(2,1),counter_g(3,1))+counter_g(6,1)+counter_g(8,1),   &
                              & ' | ', &
                              & counter(3,1), counter(3,2), counter(3,3),counter_g(1,1),' | ', & 
                              & counter(12,1), counter(12,2), counter(12,3),counter_g(4,1),' | ', & 
                              & counter(4,1), counter(4,2), counter(4,3),counter_g(2,1),' | ', & 
                              & counter(5,1), counter(5,2), counter(5,3),counter_g(3,1),' | ', &             
                              & counter(7,1), counter(8,1), counter(6,1), ' | ', &                     
                              & counter(7,2), counter(8,2), counter(6,2), ' | ', &                     
                              & counter(7,3), counter(8,3), counter(6,3), ' | ', &                     
                              & counter(10,1), counter(11,1), counter(9,1),' | ', time, counter(1,1), termination 

                              
                    END IF

                END IF
            
                
                DEALLOCATE(x_solution)
                
            END DO   ! END: Different dimensions
      
            DEALLOCATE(problem10,problem20,f_solution,counter,counter_g)
            CALL deallocate_prob_data()
             
        END DO   ! END: Different test problems
        
        WRITE(46,*)
        WRITE(46,*) 'Total number of problems', tot_num
        WRITE(46,*)
        WRITE(46,*) '------------------------------------------------------------------------------', &
                     & '------------------------------------------------------------------------------', &
                     & '------------------------------------------------------------------------------'
        
        
        CLOSE(46)
 
    END PROGRAM timbdc






















