      MODULE bundle1  
      
        USE constants, ONLY : dp   ! double precision (i.e. accuracy)
        USE functions              ! double precision (i.e. accuracy)

        IMPLICIT NONE 
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |           THE BUNDLE ELEMENT AND THE BUNDLE OF THE DC COMPONENT F_1              | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   

        TYPE bundle_element1 ! bundle element of F_1
           PRIVATE
           REAL (KIND=dp), DIMENSION(:), POINTER  :: subgrad    ! subgradient of the bundle element
           REAL (KIND=dp) :: lin_error                          ! linearization error of the bundle element
           INTEGER :: component                                 ! Index of the objective function for which the element is calculated
        END TYPE bundle_element1        

        TYPE kimppu1 ! bundle of F_1
           PRIVATE
           TYPE(bundle_element1), DIMENSION(:), POINTER :: b_elements   ! bundle elements (does NOT contain the 'current element' and 'agg_element')
           TYPE(bundle_element1) :: current_element ! bundle element calculated at the current iteration point ('current element') 
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           TYPE(bundle_element1) :: agg_element ! the aggregated bundle element ('agg_element')
           
           INTEGER :: n         ! number of variables (also the length of subgradients)
           INTEGER :: b_maxsize ! 'maximum size of the bundle' - 1, (i.e. b_maxsize=size(b_elements) NOTICE: the 'current element' and 'agg_element' are stored separately)        
           INTEGER :: b_size    ! the current size of the bundle without the 'current element' and 'agg_element' (the actual size of the bundle is 'b_size+1' and 'agg_element is NOT taken into account in this value)
           INTEGER :: indeksi   ! the place where the next bundle element is tried to be added in the bundle element table 'b_elements'  
          
           LOGICAL :: full      ! tells whether the bundle is full or not             
           ! NOTICE: if the aggregated bundle element 'agg_element' is used, then the actual size of the bundle is b_size+2, since the 'agg_element' is also stored separately.
           LOGICAL :: agg       ! tells whether the aggregated bundle element was inserted into the bundle during the previous round           
        END TYPE kimppu1 


        CONTAINS
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                           CONTAINS SUBROUTINES:                                  | | 
        !| |                                                                                  | |
        !| |    INITIALIZATION         : init_bundle_b1(set, set_size, grad_length)           | |
        !| |    ADD ELEMENT            : add_element_b1(set, grad, alpha, num)                | |
        !| |    ADD AGGREGATED ELEMENT : add_agg_element_b1(set, grad, alpha)                 | |
        !| |    ADD 1. CURRENT ELEMENT : add_first_element_b1(set, grad, num)                 | |
        !| |    UPDATE BUNDLE          : update_b1(set, new_grad, new_ind, d, &               | |
        !| |                                        & value_change, x_old, x_new,             | |
        !| |                                        &  number_of_obj, number_of_func)         | |
        !| |    UPDATE FACTORS IN BUNDLE: factor_update_b1(set, factor_old, factor_new)       | |
        !| |    UPDATE BUNDLE (imp.func.): update_b1_imp(set, new_grad, new_ind, d,  &        | |
        !| |                                             & AB_new, AB_old, fg1_old, &         | |
        !| |                                             & fg2_old, fg1_new, fg2_new, &       | |
        !| |                                             & number_of_obj, number_of_func)     | |
        !| |    SCALING OF THE BUNDLE  : scale_elements(set, ind, increase)                   | |
        !| |    DELETE FROM BUNDLE     : delete_b1(set, eps)                                  | |
        !| |    RESET BUNDLE           : reset_b1(set)                                        | |
        !| |    RESET AGG. ELEMENT     : reset_agg_element(set)                               | |
        !| |    DEALLOCATION           : deallocate_b1(set)                                   | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |              CONTAINS FUNCTIONS GIVING DIFFERENT VALUES:                         | |   
        !| |                                                                                  | |
        !| |    MATRIX OF SUBGRADIENTS      : grad_matrix(set)                                | |
        !| |    MATRIX OF LIN.ERR.-AB_k     : lin_err_and_f1_matrix(set,AB_multi)             | |
        !| |    INDEX OF LAST ELEMENT       : give_last_element_ind_b1(set)                   | |
        !| |    BUNDLE SIZE                 : give_size_b1(set)                               | |
        !| |    NUMBER OF VARIABLES         : give_n_b1(set)                                  | |
        !| |    IS BUNDLE FULL?             : is_full_b1(set)                                 | |
        !| |    IS AGGREGATION USED?        : is_agg_used(set)                                | |
        !| |    SUBGRADIENT OF ELEMENT i    : give_subgrad_b1(set, i)                         | | 
        !| |    LIN. ERROR OF ELEMENT i     : give_linerr_b1(set, i)                          | | 
        !| |    PREDICTED DECREASE          : give_decrease_b1(set,d)                         | |       
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| |                                                                                  | |       
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        
        
        
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                SUBROUTINES                                       | | 
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       
        
        
        !**************************************************************************************
        !                                                                                     |
        !                               INITIALIZATION                                        |
        !                                                                                     |
        !************************************************************************************** 
           
           SUBROUTINE init_bundle_b1(set, set_size, grad_length) 
               !
               ! Initializes the bundle 'set'. The size of the bundle is 'set_size' and the length of subgradients is 'grad_size'.
               ! 
               ! 
               ! NOTICE: * 'grad_length' >= 1
               !         * IF (set_size < 2 ) THEN the size of the bundle is set to be 1 and only the 'current element' is stored (If aggregation is used, then also the aggregated element 'agg_element' is stored)               
               !         * 'set_size' does NOT include the 'aggregated element'. So if aggregation is used, then the actual size of the bundle is 'set_size+1'. 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set          ! bundle
               INTEGER, INTENT(IN):: set_size, grad_length  ! bundle size and the length of subgardients
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i, allocstat
               
                           
               IF (set_size < 2) THEN
                    set%b_maxsize = 0                       ! only the 'current element' and the 'agg_element' (if used) are stored 
                    set%full = .TRUE.
                
               ELSE     
                    set%b_maxsize = set_size - 1            ! the biggest possible size of the bundle without the 'current element' (the 'agg_element' is not taken into account here)
                    set%indeksi = 1
                    set%full = .FALSE.
               END IF
               
               set%b_size = 0                    ! the number of stored bundle elements in the table 'b_elements' ( ! without the 'current element' and 'agg_element' ! )  
               set%n = grad_length               ! the number of variables (this is also the length of subgradients)
               set%agg = .FALSE.
               
               ALLOCATE(set%b_elements(set%b_maxsize), STAT=allocstat)  ! initializes the maximum size of the bundle table 'b_elements' 
               ALLOCATE(set%current_element%subgrad(grad_length), &     ! initializes the length of the subgradient in the 'current element'
                         & STAT=allocstat)             
               ALLOCATE(set%agg_element%subgrad(grad_length), &         ! initializes the length of the subgradient in the 'aggregated element'
                         & STAT=allocstat)  

               DO i=1, set%b_maxsize
                    ALLOCATE(set%b_elements(i)%subgrad(grad_length), &  ! initializes the length of subgradients in the table 'b_elements'
                        & STAT=allocstat)     
               END DO 
               
           END SUBROUTINE init_bundle_b1
           
           

        !**************************************************************************************
        !                                                                                     |
        !                     ADD ELEMENT INTO TO THE BUNDLE                                  | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_element_b1(set, grad, alpha, num)
               !
               ! Adds the element '(grad, alpha)' into the bundle 'set' (i.e. into the bundle element table 'b_elements'). 
               ! 'num' is the index of the objective for which the element is calculated.
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * 'num' is the index of the objective for which the element is calculated
               !         * the dimension of the vector 'grad' has to be 'set%n'.
               !         * IF the size of the bundle is 1, THEN nothing is added to the bundle element table 'b_elements'.         
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad   ! the subgradient                     
               REAL(KIND=dp), INTENT(IN) :: alpha                   ! the corresponding linearization error
               INTEGER, INTENT(IN) :: num                           ! the index of the objective
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i
               
               IF (set%b_maxsize > 0 ) THEN ! executed if bundle is larger than 1 (i.e. something can be stored into the table 'b_elements')
                   IF ( set%indeksi > set%b_maxsize ) THEN
                       set%indeksi = 1         
                   END IF
               
                   i = set%indeksi
                   set%b_elements(i)%subgrad = grad     ! adds the new subgradient into position i
                   set%b_elements(i)%lin_error = alpha  ! adds the new linearization error into position i
                   set%b_elements(i)%component = num    ! adds the index of the objective into position i
                   set%indeksi = i + 1                  ! the position where the next element is tried to be added
                   
                   IF ( .NOT. set%full ) THEN           ! if the bundle was not full during the previous round, then the size of the bundle is increased with 1
                       set%b_size = set%b_size + 1
                   END IF                  
                   
                   IF(set%b_size == set%b_maxsize) THEN  ! we test: Is the bundle full ?
                       set%full = .TRUE.
                   ELSE
                       set%full = .FALSE.
                   END IF 
                   
               END IF
               
           END SUBROUTINE add_element_b1
            
           
           
        !**************************************************************************************
        !                                                                                     |
        !                     ADD AGGREGATED ELEMENT INTO TO THE BUNDLE                       | 
        !                                                                                     |
        !**************************************************************************************        

           SUBROUTINE add_agg_element_b1(set, grad, alpha)
               !
               ! Adds the aggregated element '(grad, alpha)' into the bundle 'set'.
               !
               ! NOTICE: * 'grad' is the subgradient and 'alpha' is the corresponding linearizatio error. 
               !         * the dimension of the vector 'grad' has to be 'set%n'.           
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                  ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad   ! the aggregated subgradient                      
               REAL(KIND=dp), INTENT(IN) :: alpha                   ! the corresponding linearization error
               !**************************** OTHER VARIABLES **************************************
                           
               set%agg_element%subgrad = grad
               set%agg_element%lin_error = alpha 
               set%agg = .TRUE.
               
           END SUBROUTINE add_agg_element_b1           

           
        !**************************************************************************************
        !                                                                                     |     
        !                  INITIALIZE/ADD THE FIRST CURRENT ELEMENT                           |
        !                                                                                     |     
        !**************************************************************************************            
           
           SUBROUTINE add_first_element_b1(set, grad, num)
               !
               ! Adds the element '(grad, 0)' calculated at the first iteration point x_0 into the bundle 'set'.
               ! 'num' is the index of the objective function for which the element is calculated.
               !
               ! NOTICE: * the dimension of the 'grad' has to be 'set%n'.
               !         * 'num' is the index of the objective function for which the element is calculated.
               !         * the linearization error of the first current element is always zero.            
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                   ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: grad    ! the subgradient    
               INTEGER, INTENT(IN):: num                             ! the index of the objective  
               !**************************** OTHER VARIABLES **************************************            
               
               set%current_element%subgrad = grad
               set%current_element%lin_error = 0.0_dp ! the linearization error is zero at the iteration point x_0
               set%current_element%component = num    ! the linearization error is zero at the iteration point x_0
               
           END SUBROUTINE add_first_element_b1  

           
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE update_b1(set, new_grad, d, value_change, ind)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this updates also all the other linearization errors in the bundle 'set' 
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * f1(x_{k+1}) - f1(x_k) is the 'value_change'             
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: new_grad   ! the subgradient calculated at the new iteration point               
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=dp), INTENT(IN) :: value_change                ! f1(x_{k+1}) - f1(x_k), i.e. the value change in the objective function
               INTEGER, INTENT(IN) :: ind                               ! The index of the objective function
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i            
               
               ! The old 'current element' is added into the bundle set 'b_elements' and after that 
               ! the 'current element' can be updated with the new element
               ! (the linearization error of the current element is always zero, thus it is not changed).                  
               CALL add_element_b1(set, set%current_element%subgrad, 0.0_dp, ind)
               set%current_element%subgrad = new_grad
               set%current_element%component = ind
           
               ! Linearization error update in the bundle set 'b_elements'
                
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error + &
                          & value_change - DOT_PRODUCT(set%b_elements(i)%subgrad, d)
               END DO
               
               IF (set%agg) THEN
                   set%agg_element%lin_error = set%agg_element%lin_error + &
                          & value_change - DOT_PRODUCT(set%agg_element%subgrad, d)           ! update of aggregated element  
               END IF
               
           
           END SUBROUTINE update_b1 
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                        UPDATE THE FACTORS IN THE BUNDLE                             |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE factor_update_b1(set, factor_old, factor_new) 
               !
               ! Updates the bundle elements based on the changed scaling factor in the objective function
               !          
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set          ! The bundle
               REAL(KIND=dp), INTENT(IN):: factor_old       ! the old scaling factor in the objective function         
               REAL(KIND=dp), INTENT(IN) :: factor_new      ! the new scaling factor in the objective function
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: multiplier
               INTEGER :: i            
               
               multiplier = factor_new/factor_old
               
               ! the 'current element' is updated 
               ! (the linearization error of the current element is always zero, thus it is not changed).                  
               set%current_element%subgrad = set%current_element%subgrad*multiplier
           
               ! Updates in the bundle set 'b_elements'                
               DO i = 1, set%b_size
                   set%b_elements(i)%lin_error = set%b_elements(i)%lin_error * multiplier
                   set%b_elements(i)%subgrad = set%b_elements(i)%subgrad * multiplier
               END DO
               
              ! update of aggregated element    
               IF (set%agg) THEN
                   set%agg_element%lin_error = set%agg_element%lin_error * multiplier          
                   set%agg_element%subgrad = set%agg_element%subgrad * multiplier          
               END IF
               
           
           END SUBROUTINE factor_update_b1 
           
           
        !**************************************************************************************
        !                                                                                     |     
        !                                UPDATE THE BUNDLE                                    |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE update_b1_imp(set, new_grad, new_ind, d, AB_new, AB_old, &
                                & fg1_old, fg2_old, fg1_new, fg2_new, & 
                                & number_of_obj, number_of_func)
               !
               ! Updates the 'current element' with the bundle element calculated at the new iteration point x_(k+1)
               ! and due to this also the linearization errors are updated in the bundle 'set' (utilize when improvement function is used)
               !
               ! NOTICE: * the dimension of vectors 'new_grad' and 'd' has to be 'set%n'
               !         * the vector 'd' is the new search direction d^k = x_{k+1} - x_k
               !         * the dimension of vectors 'AB_new', 'AB_old', fg1_old, fg2_old, fg1_new, and fg2_new has to be 'number_of_func'         
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN):: new_grad   ! the subgradient calculated at the new iteration point                               
               REAL(KIND=dp), DIMENSION(set%n), INTENT(IN) :: d         ! d^k = x_{k+1} - x_k, i.e. the search direction
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: AB_new     ! the values of new A_i and B_l components 
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: AB_old     ! the values of old A_i and B_l components
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: fg1_old    ! the old values of the DC components f1_i and g1_l
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: fg2_old    ! the old values of the DC components f2_i and g2_l
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: fg1_new    ! the new values of the DC components f1_i and g1_l
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN) :: fg2_new    ! the new values of the DC components f2_i and g2_l
               INTEGER, INTENT(IN) :: new_ind                                     ! the index of the objective at the new bundle element
               INTEGER, INTENT(IN) :: number_of_obj                               ! the index of the objective at the new bundle element
               INTEGER, INTENT(IN) :: number_of_func                              ! the index of the objective at the new bundle element
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, ind               
               
               ! The old 'current element' is added into the bundle set 'b_elements' and after that 
               ! the 'current element' can be updated with the new element
               ! (the linearization error of the current element is always zero, thus it is not changed).                  
               CALL add_element_b1(set, set%current_element%subgrad, 0.0_dp, &
                                    & set%current_element%component)
               set%current_element%subgrad = new_grad
               set%current_element%component = new_ind
           
               ! Linearization error update in the bundle set 'b_elements'
                
               DO i = 1, set%b_size
                   ind = set%b_elements(i)%component
                   IF (ind < (number_of_obj_sep+1)) THEN ! update for A_i
                        set%b_elements(i)%lin_error = set%b_elements(i)%lin_error &
                           & + AB_new(ind) - AB_old(ind) &  
                           & - DOT_PRODUCT(set%b_elements(i)%subgrad, d)&                                                                              
                           & - fg1_old(ind) + fg2_old(ind) &  
                           & + fg1_new(ind) - fg2_new(ind)   
                    ELSE ! update for B_l
                        set%b_elements(i)%lin_error = set%b_elements(i)%lin_error &
                           & + AB_new(ind) - AB_old(ind) &  
                           & - DOT_PRODUCT(set%b_elements(i)%subgrad, d)                    
                    END IF
      
               END DO
               
               IF (set%agg) THEN
                 !  No update for aggregated element  
               END IF
               
           
           END SUBROUTINE update_b1_imp 
           
  
        !**************************************************************************************
        !                                                                                     |
        !                    CHANGE THE SCALING OF ELEMENTS IN THE BUNDLE                     |
        !                                                                                     |
        !**************************************************************************************    

           SUBROUTINE scale_elements(set, ind, increase, user_n)
               !
               ! Changes the scaling of the elements in the bundle 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set  ! bundle
               INTEGER, INTENT(IN):: ind            ! the index of objective or constraint which elements are scaled    
               INTEGER, INTENT(IN):: increase       ! the coefficient for scaling
               INTEGER, INTENT(IN):: user_n         ! The dimension of the problem
                           
               !**************************** OTHER VARIABLES **************************************  
               REAL(KIND=dp), DIMENSION(user_n) :: subgrad       ! the subgradient of objective or constraint which is scaled    
               REAL(KIND=dp) :: linerr                           ! the linearization error of objective or constraint which is scaled                  
               INTEGER :: i, comp
              
               DO i = 1, set%b_size
                   comp = set%b_elements(i)%component
                   IF (ind == comp) THEN 
                       subgrad = ((10.0_dp)**increase)*set%b_elements(i)%subgrad
                       linerr = ((10.0_dp)**increase)*set%b_elements(i)%lin_error
                       
                       set%b_elements(i)%subgrad = subgrad
                       set%b_elements(i)%lin_error = linerr
                   END IF
               END DO
               
               comp = set%current_element%component
               IF (ind == comp) THEN 
                  subgrad = ((10.0_dp)**increase)*set%current_element%subgrad
                  linerr = ((10.0_dp)**increase)*set%current_element%lin_error
                       
                  set%current_element%subgrad = subgrad
                  set%current_element%lin_error = linerr
               END IF
               
           END SUBROUTINE scale_elements       
           

        !**************************************************************************************
        !                                                                                     |
        !                            DELETE FROM THE BUNDLE                                   |
        !                                                                                     |
        !**************************************************************************************    

           SUBROUTINE delete_b1(set, eps)
               !
               ! Deletes from the bundle 'set' elements for which the linearization error 'alpha' > 'eps'.
               !
               ! NOTICE: the linearization error for the' current element' is always zero so it is never deleted/removed.              
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set  ! bundle
               REAL(KIND=dp), INTENT(IN):: eps      ! epsilon    
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: j, k, up
              
               k = 1 
               up = set%b_size
               
               DO j = 1, up 
                   IF( set%b_elements(j)%lin_error + 0.000000001_dp <= eps ) THEN
                      IF (j /= k ) THEN 
                          set%b_elements(k)%subgrad = set%b_elements(j)%subgrad     
                          set%b_elements(k)%lin_error = set%b_elements(j)%lin_error
                      END IF
                      k = k + 1 
                   END IF
               END DO
               
               IF (set%b_size /= k - 1) THEN                ! If .TRUE. then something is deleted/removed from the bundle 
                   set%b_size = k - 1                       ! the new number of stored subgradients
                   set%indeksi = k                          ! the place where the next bundle element will be added
                   IF (set%b_size /= set%b_maxsize) THEN    ! we test whether the bundle is full or not
                      set%full = .FALSE.
                   ELSE 
                      set%full = .TRUE.
                   END IF
               END IF
               
           END SUBROUTINE delete_b1

        !**************************************************************************************
        !                                                                                     |     
        !                                DEALLOCATION                                         |
        !                                                                                     |
        !**************************************************************************************

           SUBROUTINE deallocation_b1(set)
               !
               ! Deallocates arrays     
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set                      ! bundle
               !**************************** OTHER VARIABLES **************************************
               INTEGER :: i, allocstat             
               
               DEALLOCATE(set%current_element%subgrad, STAT=allocstat)     ! deallocates the subgradient in the 'current element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"              
               DEALLOCATE(set%agg_element%subgrad, STAT=allocstat)         ! deallocates the subgradient in the 'aggregated element'
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"                          

               DO i=1, set%b_maxsize
                    DEALLOCATE(set%b_elements(i)%subgrad, STAT=allocstat)  ! deallocates subgradients in the table 'b_elements'
                    IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"             
                    
               END DO            
               DEALLOCATE(set%b_elements, STAT=allocstat)  ! deallocates the bundle table 'b_elements' 
               IF (allocstat /= 0) STOP "*** Could Not release memory B1 ***"              

               
           END SUBROUTINE deallocation_b1      
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                               RESET THE BUNDLE                                      |
        !                                                                                     |
        !************************************************************************************** 
        
           SUBROUTINE reset_b1(set)
               !
               ! deletes all the elements from the bundle 'set' except the 'current element' (Also the 'aggregated element' is deleted)
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set ! bundle

               IF (set%b_maxsize > 0) THEN      ! Reset is executed if it is possible that we have something in the bundle element table 'b_elements'
                   set%b_size = 0               ! the number of stored subgradient in 'b_elements' after reset
                   set%indeksi = 1              ! the place where the next element is placed
                   set%full = .FALSE.           ! the bundle is not full because elements from 'b_elements' are removed
               END IF

               set%agg = .FALSE.            ! the aggregated element is deleted
               
           END SUBROUTINE reset_b1      
           
           
        !**************************************************************************************
        !                                                                                     |
        !                               RESET THE AGGREGATED ELEMENT                          |
        !                                                                                     |
        !************************************************************************************** 
        
           SUBROUTINE reset_agg_element(set)
               !
               ! Resets the 'aggregated element' 
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(INOUT) :: set ! bundle

               set%agg = .FALSE.            ! the aggregated element is deleted
               
           END SUBROUTINE reset_agg_element                
           
           
           
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                        FUNCTIONS GIVING DIFFERENT VALUES                         | |   
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   


        !**************************************************************************************
        !                                                                                     |
        !                             MATRIX OF SUBGRADIENTS                                  |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+1)'.
               !         * The subgradient corresponding to the current element is in the last position.
               !         * The aggregated element is NOT taken into account.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=dp), DIMENSION(set%n*(set%b_size+1)) :: m  ! the subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, j, length, start
               
               length = set%n               ! the length of subgradients
                       
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! the place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO
           END FUNCTION grad_matrix
           
           
        !**************************************************************************************
        !                                                                                     |
        !                           VECTOR OF LINEARIZATION ERRORS                            |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+1'.
               !         * The linearization error corresponding to the current element is in the last position.    
               !         * The aggregated element is NOT taken into account.               
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=dp), DIMENSION(set%b_size+1) :: m  ! the linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the last position

           END FUNCTION lin_error_matrix

        !**************************************************************************************
        !                                                                                     |
        !                         MATRIX OF SUBGRADIENTS WITH AGGREGATION                     |
        !                                                                                     |
        !**************************************************************************************
        
           PURE FUNCTION grad_matrix_agg(set) RESULT(m)
               !
               ! Returns the subgradient matrix 'm' (with the aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the matrix 'm' is 'set%n*(set%b_size+2)'.
               !         * The subgradient corresponding to the current element is in the 'last but one' position.
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set                     ! bundle
               REAL(KIND=dp), DIMENSION(set%n*(set%b_size+2)) :: m  ! subgradient matrix formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: i, j, length, start
               
               length = set%n               ! the length of subgradients
               
               DO i = 1, set%b_size         ! each subgradient is copied from the bundle element table 'b_elements'
                  start = (i-1)*length      ! place where the subgradient is copied
                  DO j = 1, length 
                       m(start+j) = set%b_elements(i)%subgrad(j)
                  END DO           
               END DO
               
               start = set%b_size * length  ! the place where the 'current element' is copied
               DO j = 1, length 
                   m(start + j) = set%current_element%subgrad(j)
               END DO

               start = (set%b_size + 1) * length  ! the place where the 'aggregated element' is copied             
               DO j = 1, length 
                   m(start + j) = set%agg_element%subgrad(j)
               END DO
               
           END FUNCTION grad_matrix_agg
        
        
        
        !**************************************************************************************
        !                                                                                     |
        !                  VECTOR OF LINEARIZATION ERRORS WITH AGGREGATION                    |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION lin_error_matrix_agg(set) RESULT(m)
               !
               ! Returns the linearization error vector 'm' (with aggregation) formed from the bundle 'set' of the DC component f_1.
               ! Needed when we calculate the search direction. 
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+2'.
               !         * The linearization error corresponding to the current element is in the 'last but one' position.  
               !         * The subgradient corresponding to the aggregated element is in the last position.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=dp), DIMENSION(set%b_size+2) :: m  ! linearization error vector formed from the bundle 'set'
               !**************************** OTHER VARIABLES **************************************            
               INTEGER ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       m(j) = set%b_elements(j)%lin_error          
               END DO

               m(set%b_size+1) = set%current_element%lin_error  ! the lin. error corresponding to the 'current element' is copied into the 'last but one' position
               m(set%b_size+2) = set%agg_element%lin_error      ! the lin. error corresponding to the 'aggregated element' is copied into the last position

           END FUNCTION lin_error_matrix_agg               
        
          
        
        !**************************************************************************************
        !                                                                                     |
        !                  VECTOR OF LINEARIZATION ERRORS + A_i(x_k,x_k)/B_l(x_k)             |
        !                                                                                     |
        !**************************************************************************************     

            PURE FUNCTION lin_err_and_f1_matrix(set,AB_multi,number_of_func) RESULT(m)
               !
               ! Substracts vector 'AB_multi' from linearization errors vector formed from the bundle 'set' of the DC component f_1 and returns the corresponding vector 'm'.
               ! Needed when we calculate the search direction.
               !
               ! NOTICE: * The size of the vector 'm' is 'set%b_size+1'.
               !         * The linearization error corresponding to the current element is in the last position.    
               !         * The aggregated element is NOT taken into account.  
               !         * The size of 'AB_multi' is 'number_of_func'.             
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set             ! bundle
               REAL(KIND=dp), DIMENSION(number_of_func), INTENT(IN):: AB_multi   ! the values of A_i and B_l at the current iteration point                 
               INTEGER, INTENT(IN) :: number_of_func
               
               
               REAL(KIND=dp), DIMENSION(set%b_size+1) :: m  ! the "linearization error" vector formed from the bundle 'set'
               INTEGER :: ind
               
               !**************************** OTHER VARIABLES **************************************            
               INTEGER ::  j
               
               DO j = 1, set%b_size                             ! linearization errors corresponding to the table 'b_elements' are copied
                       ind = set%b_elements(j)%component
                       m(j) = set%b_elements(j)%lin_error - AB_multi(ind)                  
               END DO
               
               ind = set%current_element%component
               m(set%b_size+1) = set%current_element%lin_error - AB_multi(ind) ! the lin. error corresponding to the 'current element' is copied into the last position

           END FUNCTION lin_err_and_f1_matrix
           
                  
        !**************************************************************************************
        !                                                                                     |
        !                             INDEX OF THE LAST ELEMENT                               |
        !                                                                                     |
        !**************************************************************************************     

           PURE FUNCTION give_last_element_ind_b1(set) RESULT(ind)
               !
               ! Gives the index of the place where the last bundle element was added in the bundle element table 'b_elements'.
               !
               ! NOTICE: the index 'ind' is zero if there is nothing in the bundle element table 'b_elements'.
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: ind                   ! index of the place where the last element was added
               
               IF (set%b_size /= 0) THEN        ! there is something in the bundle element table 'b_elements'
                   ind = set%indeksi - 1
               ELSE                             ! there in nothing in the bundle element table 'b_elements'
                   ind = 0
               END IF 
               
           END FUNCTION give_last_element_ind_b1        


           
        !**************************************************************************************
        !                                                                                     |
        !                                    BUNDLE SIZE                                      |
        !                                                                                     |
        !**************************************************************************************
       
           PURE FUNCTION give_size_b1(set) RESULT(bundle_size)
               !
               ! Gives the current size of the bundle 'set' (NOTICE: ! without current element and the aggregated element ! )
               !
               IMPLICIT NONE 
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               !**************************** OTHER VARIABLES **************************************            
               INTEGER :: bundle_size           ! size of the bundle
               
               bundle_size = set%b_size
               
           END FUNCTION give_size_b1
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                              NUMBER OF VARIABLES                                    |
        !                                                                                     |
        !**************************************************************************************        
               
           PURE FUNCTION give_n_b1(set) RESULT(variable_number)
               !
               ! Gives the number of varibles in the minimization problem (this is also the length of subgradients)
               ! Number of varibles is (i.e. has to be) same as in kimppu2 when used in the algorithm.
               !
               IMPLICIT NONE 
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               INTEGER :: variable_number        ! the number of variables
               
               variable_number = set%n
               
           END FUNCTION give_n_b1          
           
           
           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS BUNDLE FULL?                                     |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_full_b1(set) RESULT(isfull)
               !
               ! Returns .TRUE. if bundle 'set' (i.e. the bundle element table 'b_element') is full otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isfull                 ! tells whether the bundle is full or not
 
               isfull = set%full 
               
            END FUNCTION is_full_b1         

           
        !**************************************************************************************
        !                                                                                     |
        !                                 IS AGGREGATION USED?                                |
        !                                                                                     |
        !**************************************************************************************        

           PURE FUNCTION is_agg_used(set) RESULT(isUsed)
               !
               ! Returns .TRUE. if aggregation is used in the bundle 'set' otherwise retuns .FALSE.
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set  ! bundle
               !**************************** OTHER VARIABLES **************************************            
               LOGICAL :: isUsed                 ! tells whether the aggregation is used or not
 
               isUsed = set%agg
               
            END FUNCTION is_agg_used                    
        
        
        !**************************************************************************************
        !                                                                                     |
        !                            SUBGRADIENT OF ELEMENT i                                 |
        !                                                                                     |
        !**************************************************************************************
        
           FUNCTION give_subgrad_b1(set, i) RESULT(grad)
               !
               ! Gives the subgradient of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the subgradient of the 'current element'.    
               !         * If 'i=-1' then gives the subgradient of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)

               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set         ! bundle
               INTEGER, INTENT(IN) :: i                 ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp), DIMENSION(set%n) :: grad  ! subgradient at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! subgradient is from the bundle element table 'b_elements'
                   grad = set%b_elements(i)%subgrad
               ELSE IF (i == 0) THEN                        ! subgradient is from the 'current element'
                   grad = set%current_element%subgrad
               ELSE IF (i == -1) THEN                       ! subgradient is from the 'aggregated element'
                   grad = set%agg_element%subgrad                  
               ELSE                                         ! otherwise we are outside the bundle     
                   WRITE(*,*) 'CANNOT RETURN SUBGRADIENT! index ' &
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_subgrad_b1
           
           
        !**************************************************************************************
        !                                                                                     |
        !                          LINEARIZATION ERROR OF ELEMENT i                           |
        !                                                                                     |
        !**************************************************************************************        
    
           FUNCTION give_linerr_b1(set, i) RESULT(error)
               !
               ! Gives the linearization error of the bundle element at position 'i'.
               !
               ! NOTICE: * -1 <= 'i' <= 'set%b_size' (otherwise we are surely outside the bundle).
               !         * If 'i'=0 then gives the linearization error of the 'current element'.    
               !         * If 'i=-1' then gives the linearization error of the aggregated element (NOTICE: Does not take into account whether aggregation is really used or not.)
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               INTEGER, INTENT(IN) :: i         ! the position
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: error           ! the linearization error at the position 'i'
         
               IF ( (i <= set%b_size) .AND. (i > 0) ) THEN  ! the linearization error is from the bundle element table 'b_elements'
                   error = set%b_elements(i)%lin_error
               ELSE IF (i==0) THEN                          ! the linearization error is from the 'current element'
                   error = set%current_element%lin_error
               ELSE IF (i==-1) THEN                         ! the linearization error is from the 'aggregated element'      
                   error = set%agg_element%lin_error               
               ELSE                                         ! otherwise we are outside the bundle
                   WRITE(*,*) 'CANNOT RETURN LINEARIZATION ERROR! index '&
                         & , i , 'outside the bundle (size without the current element:',& 
                         & set%b_size, ')'
               END IF              
           END FUNCTION give_linerr_b1    

        !**************************************************************************************
        !                                                                                     |
        !                          PREDICTED DECREASE FOR DIRECTION d                         |
        !                                                                                     |
        !**************************************************************************************        
    
           FUNCTION give_decrease_b1(set, d) RESULT(decrease)
               !
               ! Gives the predicted decraese fro direction d
               !
               ! NOTICE: * the dimension of d shoud be user_n
               !
               IMPLICIT NONE
               !**************************** NEEDED FROM USER *************************************
               TYPE(kimppu1), INTENT(IN) :: set ! bundle
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: d         ! the direction
               !**************************** OTHER VARIABLES **************************************            
               REAL(KIND=dp) :: decrease   ! The predicted decrease
               REAL(KIND=dp) :: max_piece  ! The predicted decrease at the moment
               REAL(KIND=dp) :: piece      ! 
               INTEGER :: i, j, current_size, length 
               
               length = set%n               ! the length of subgradients
               current_size = set%b_size    ! the size of the bundle 
               
               max_piece = 0.0_dp
               DO j = 1, length
                  max_piece = max_piece + set%current_element%subgrad(j) * d(j)
               END DO
               max_piece = max_piece - set%current_element%lin_error
               
               DO i = 1, current_size
                   piece = 0.0_dp
                   DO j = 1, length
                      piece = piece + set%b_elements(i)%subgrad(j) * d(j)
                   END DO
                   piece = piece - set%b_elements(i)%lin_error             
               
                   IF (piece > max_piece) THEN 
                       max_piece = piece 
                   END IF
               END DO      
               
               decrease = max_piece
                           
           END FUNCTION give_decrease_b1         
        
        

      END MODULE bundle1