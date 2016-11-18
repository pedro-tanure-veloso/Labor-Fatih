!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!DEFINE THE UTILITY FUNCTION AS LOG!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function utility(c)
    IMPLICIT NONE
    double precision ans, c

    ans=log(c)

    utility=ans

END function utility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!DEFINE THE PRODUCTION FUNCTION!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!COBB DOUGLAS LIKE

double precision function production_fun(k,delta,a,alpha)
    IMPLICIT NONE
    double precision k,delta,a,alpha,ans

    ans=a*k**alpha + (1-delta)*k

    production_fun=ans

END function production_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!FIND THE STEADY STATE OF CAPITAL!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


double precision function k_ss(delta,a,alpha, beta)
    IMPLICIT NONE
    double precision k,delta,a,alpha,ans, beta

    ans=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))

    k_ss=ans

END function k_ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!JUST THE ITERATION OF VALUE FUNCTION !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I DEFINE THE ITERATION SEPARTLY FROM THE REST OF THE SUBRUTINE 
!TWO LOOPS TO FIND THE OPTIMAL AND THEN THE FIRST LOOP IS TO ITERATE
!
SUBROUTINE iterate(beta,delta,a,alpha,num_points,initial_guess, policy_funtion,grid_k,value_function)

    IMPLICIT NONE
    integer num_points,i,j
    double precision beta,delta,a,alpha,grid_k(num_points), produc,value_function(num_points)
    double precision initial_guess(num_points),tolerance,error,max_val(num_points)
    double precision policy_funtion(num_points),val_funtion(num_points),utils(num_points)
    double precision, external :: production_fun,k_ss,utility
    !double precision, dimension(num_points) ::iterate

    tolerance=.0000001
    error=100

    val_funtion=initial_guess

    max_val=-100


do while (error>tolerance)


    do  i=1,num_points


        do  j=1,num_points

            produc= production_fun(grid_k(i),delta,a,alpha)

            if (produc-grid_k(j)>0) Then
                utils(j)= utility(produc-grid_k(j))+beta*val_funtion(j)
            else
                utils(j)=-1000000
            END IF

            if (utils(j)>max_val(i)) Then
                max_val(i)=utils(j)
                policy_funtion(i)=grid_k(j)
            END IF

         end  do

    end  do

    error=MAXVAL(abs(val_funtion-max_val))
    !print*, "error is ", error
    val_funtion=max_val
    max_val=-100

 end  do

value_function=val_funtion

END SUBROUTINE iterate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!THIS IS WHAT DOES THE VFI AND IS CALLED IN MAIN!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FIND THE SS K AND CREATES A GRID AROUND IT. THEN IT CALLS THE ITERATE SUBROUTINE 
!I WANTED TO MAKE IT FUNCTION THAT RETUNRS THE VALUE FUNCTION ARRAY BUT I COULDNT DO FIND A WAY, SO ITS A SUBROUTINE THAT YOU !NEED TO IT A BLANK VALUE FUNCTION AND POLICY FUNCTION


SUBROUTINE value_interation(beta,delta,a,alpha,num_points,initial_guess,value_function,policy_function)!,grid_k)
IMPLICIT NONE
integer num_points,mid,i
double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), initial_guess(num_points),error
double precision equally_spaced,value_function(num_points),policy_function(num_points)
double precision, external :: production_fun,k_ss,utility

!double precision, dimension(num_points) ::value_interation


mid=num_points/2+1             !FIND THE MID POINT SO K_SS IS MIDDLE VALUE IN THE GRID
k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))


equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))

grid_k(1)=0.0001               !STARTING IT AT 0 CAN BE PROBLEMATIC

do  i=2,num_points

    if (i==mid) Then
        grid_k(i)=k_at_ss
    else
        grid_k(i)=grid_k(i-1)+ equally_spaced
    END IF
end  do




call iterate(beta,delta,a,alpha,num_points,initial_guess, policy_function,grid_k,value_function,)

end SUBROUTINE value_interation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!WILL GIVE ANALYTICAL SOLUTION FOR DELTA =1     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!USE THIS AS A CHECK, TO SEE IF CONVERGANCE IS CORRECT


SUBROUTINE analytical(a,alpha,beta,delta,VF,PF,num_points)
IMPLICIT NONE
integer num_points,mid,i
double precision equally_spaced, a, alpha, beta, delta, grid_k(num_points), ans, a0, a1, VF(num_points), PF(num_points), k_at_ss
double precision, external :: k_ss
  
mid=num_points/2+1
k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))
equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))

grid_k(1)=0.0001

do  i=2,num_points

    if (i==mid) Then
        grid_k(i)=k_at_ss
    else
        grid_k(i)=grid_k(i-1)+ equally_spaced
    END IF
    end  do
    
a0 = (1/(1-beta))*((alpha*beta)/(1-alpha*beta)*log(alpha*beta)+log(1-alpha*beta))

a1 = alpha/(1-alpha*beta)


VF(1) = 0.0
PF(1) = 0.0

do i=2,num_points  
    VF(i) = a0 + a1*log(grid_k(i))
    PF(i) = alpha*beta*grid_k(i)**alpha
    end do
    
end SUBROUTINE analytical


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!THE ITERATION OF VALUE FUNCTION USING THE MCQUEEN BOUNDS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I DEFINE THE ITERATION SEPARTLY FROM THE REST OF THE SUBRUTINE 
!TWO LOOPS TO FIND THE OPTIMAL AND THEN THE FIRST LOOP IS TO ITERATE
! HOWEVER, NOW IT UPDATES USING THE MIDPOINT OF THE ERRORS
!IT ALSO ALLOWS THE USER TO DICTATE HOW OFTEN TO APPLY THOSE BOUNDS USING NUM
!




SUBROUTINE iterate_c(beta,delta,a,alpha,num_points,initial_guess, policy_funtion,grid_k,value_function, num)
    
    IMPLICIT NONE
    integer num_points,i,j,k,num
    double precision beta,delta,a,alpha,grid_k(num_points), produc,value_function(num_points)
    double precision initial_guess(num_points),tolerance,error,max_val(num_points)
    double precision beta_fraction, b_min, b_max, error_max(num_points), error_min(num_points),error_max_e, error_min_e
    double precision policy_funtion(num_points),val_funtion(num_points),utils(num_points)  
    double precision, external :: production_fun,k_ss,utility
    !double precision, dimension(num_points) ::iterate

    tolerance=.001
    error=100
    
    beta_fraction= beta/(1.0-beta)
   
    k=1                    !THIS KEEPS TRACK OF THE ITERATION TO SEE IF WE SHOULD APPLY THE BOUNDS
    
    val_funtion=initial_guess
    
    max_val=-100
    
    do while (error>tolerance)
        
        
        do  i=1,num_points
        
        
            do  j=1,num_points
                produc= production_fun(grid_k(i),delta,a,alpha)
                if (produc-grid_k(j)>0) Then
                    utils(j)= utility(produc-grid_k(j))+beta*val_funtion(j)
                else  
                    utils(j)=-1000000
                END IF                
        
                if (utils(j)>max_val(i)) Then
                    max_val(i)=utils(j)
                    policy_funtion(i)=grid_k(j)
                    error_min(i)=max_val(i)-val_funtion(i)
                    error_max(i)=max_val(i)-val_funtion(i)
                END IF
 
             end  do
            
        end  do
  
        !error=MAXVAL(abs(val_funtion-max_val))
        
        !print*, "max_val is", max_val
        
        
        !error_min_e=MINVAL(error_min)-MINVAL(max_val-val_funtion)  
        !error_max_e=MAXVAL(error_max)-MAXVAL(max_val-val_funtion) 
        
        !print*, "error_min_e is ", error_min_e
        !print*, "error_max_e is ", error_max_e        
        
        b_max= beta_fraction*MAXVAL(error_max)
        b_min= beta_fraction*MINVAL(error_min)
        
        
        if (k>num-1.) then
        
            !print*, "error_min is ", error_min
            !print*, "error_max is ", error_max
            error=(abs(b_max-b_min))
            !print*,"=============================="
            !print*, "error is ", error  
            !print*, "error_min is ", b_min
            !print*, "error_max is ", b_max
            !print*,"=============================="

            val_funtion=max_val+(b_max*.5+(b_min*.5))
            max_val=-100
        else 
            error=MAXVAL(abs(val_funtion-max_val))
            !print*,"=============================="
            !print*, "error is ", error  
            !print*, "error_min is ", b_min
            !print*, "error_max is ", b_max
            !print*,"=============================="

            val_funtion=max_val!+(b_max*.5+(b_min*.5))
            max_val=-100
        end if
    
     end  do
    
    value_function=val_funtion
    
END SUBROUTINE iterate_c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!THIS IS WHAT DOES THE MCQUEEN BOUNDS AND IS CALLED IN MAIN!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FIND THE SS K AND CREATES A GRID AROUND IT. THEN IT CALLS THE ITERATE SUBROUTINE 
!WE COULD PROBABLY MERGE THIS WITH THE ONE FOR PURE VFI
!IT HAS THE EXTRA PARAMETER(NUM) THAT DICTATES HOW OFTEN TO APPLY THE BOUNDS




SUBROUTINE value_interation_c(beta,delta,a,alpha,num_points,initial_guess,value_function,policy_function,num)!,grid_k_c)
    IMPLICIT NONE
    integer num_points,mid,i,num
    double precision beta,delta,a,alpha, k_at_ss
    double precision grid_k_c(num_points), initial_guess(num_points),error
    double precision equally_spaced,value_function(num_points),policy_function(num_points)
    double precision, external :: production_fun,k_ss,utility

    mid=num_points/2+1
    k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1.0-delta)))**(1.0/(alpha-1.0))
    error=10.0
    
    equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))
    

    grid_k_c(1)=0.0001

    do  i=2,num_points
        
        if (i==mid) Then
            grid_k_c(i)=k_at_ss
        else 
            grid_k_c(i)=grid_k_c(i-1)+ equally_spaced
        END IF
    end  do
   
    
    call iterate_c(beta,delta,a,alpha,num_points,initial_guess, policy_function,grid_k_c,value_function,num)
    
    

end SUBROUTINE value_interation_c



SUBROUTINE howard(it,beta,delta,a,alpha,num_points,initial_guess,policy_function,grid_k,value_function)

    IMPLICIT NONE
    integer num_points,i,j,n,m,it !"it" is the number of iterations that we will have in Howard's improvement algorithm
    double precision beta,delta,a,alpha,grid_k(num_points), produc,value_function(num_points),kp_index(num_points)
    double precision initial_guess(num_points),tolerance,error,v0(num_points),v1(num_points)
    double precision policy_function(num_points),val_funtion(num_points),utils(num_points)
    double precision, external :: production_fun,k_ss,utility
    
    tolerance = .01
    error = 100

    val_funtion = initial_guess

    v0 = -100
    n = 0


    do while (error>tolerance)


        do  i=1,num_points
            produc= production_fun(grid_k(i),delta,a,alpha)
            do  j=1,num_points
                if (produc-grid_k(j)>0) Then
                    utils(j) = utility(produc-grid_k(j))+beta*val_funtion(j) ! getting the "value function"
                else
                    utils(j) = -1000000
                END IF
                if (utils(j)>v0(i)) Then
                    v0(i) = utils(j)
                    policy_function(i) = grid_k(j) ! gets the policy function
                    kp_index(i) = j
                END IF
            end  do
        end  do
    
        do m=1,it
           do i=1,num_points
              produc= production_fun(grid_k(i),delta,a,alpha)
              v1(i) = utility(produc-policy_function(i)) + beta*v0(int(kp_index(i))) ! does the policy iteration
            end do
           v0 = v1
        end do
    
        error = MAXVAL(abs(val_funtion - v0))
        !print*, "error is ", error
        val_funtion = v0
        v0 = -100
        n = n+1
    end  do
    !print*, "number of iterations", n
    value_function=val_funtion

END SUBROUTINE howard


SUBROUTINE policy_interation_howard(it,beta,delta,a,alpha,num_points,initial_guess,value_function,policy_function)
    IMPLICIT NONE
    integer num_points,mid,i,it
    double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), initial_guess(num_points),error
    double precision equally_spaced,value_function(num_points),policy_function(num_points)
    double precision, external :: production_fun,k_ss,utility

    mid=num_points/2+1
    k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))
    error=10.0

    equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))

    grid_k(1)=0.0001

    do  i=2,num_points

        if (i==mid) Then
            grid_k(i)=k_at_ss
        else
            grid_k(i)=grid_k(i-1)+ equally_spaced
        END IF
    end  do

    call howard(it,beta,delta,a,alpha,num_points,initial_guess, policy_function,grid_k,value_function)

end SUBROUTINE policy_interation_howard










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!MAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!THIS PROGRAM WILL FIND THE VALUE FUNCTION AND POLICY FUNCTION 3 DIFFERENT WAYS:
!1.) USING THE REGULAR VFI
!2.) USING THE MCQUEEN PORTEUS BOUNDS
!3.) USING THE MCQUEEN PORTEUS BOUNDS EVERY OFTEN MODUL
!
!THEN IT WILL CALCULATE HOW MUCH TIME EACH TOOK AND HOW FAR AWAY THEY WERE FROM SIMPLE VFI





program main

IMPLICIT NONE
    integer points,modul,howard_iterations
    double precision beta,delta,a,alpha!,grid_k_c(1001),grid_k(1001)
    
    
    !EVERTHING IS ALLOCATABLE SO THAT WE JUST HAVE TO CHANGE A SINGLE PARAMETER IN ORDER DO MORE GRID POINTS
    
    double precision, ALLOCATABLE :: initial_guess(:)
    double precision, ALLOCATABLE :: value_function(:)
    double precision, ALLOCATABLE :: policy_function(:)
    double precision, ALLOCATABLE :: value_function_c(:)
    double precision, ALLOCATABLE :: policy_function_c(:)
    double precision, ALLOCATABLE :: value_function_fiv(:)
    double precision, ALLOCATABLE :: policy_function_fiv(:)    
    double precision, ALLOCATABLE :: anal_value_function(:)
    double precision, ALLOCATABLE :: anal_policy_function(:)
    double precision, ALLOCATABLE :: value_function_howard(:)
    double precision, ALLOCATABLE :: policy_function_howard(:)
        
    
    real :: start_vfi, finish_vfi
    real :: start_MPB_one, finish_MPB_one
    real :: start_MPB_five, finish_MPB_five
    real :: start_how, finish_how




    points=1001             !NUMBER OF POINTS USED FOR THE GRID
    beta=.96               !FUTURE DISCOUNTING 
    delta=1.0              !DEPRECIATION
    a=1.0                  !PRODUCTIVITY 
    alpha=1.0/3.0          !OUTPUT ELASTICITIES  
    initial_guess=0.0      !INITIAL GUESS FOR VALUE FUNCTION
    modul=5                !HOW OFTEN TO APPLY MCQUEEN BOUNDS
    howard_iterations=5    !NUMBER OF HOWARD STEPS
    
    
    
    !ALLOCATES MEMORY FOR THE ARRAYS OF SIZE (POINTS)


    allocate ( initial_guess(points))
    
    allocate ( value_function(points))
    allocate ( policy_function(points))
    
    allocate ( value_function_c(points))
    allocate ( policy_function_c(points))
    
    allocate ( value_function_fiv(points))
    allocate ( policy_function_fiv(points))
    
    allocate ( anal_value_function(points))
    allocate ( anal_policy_function(points))

    allocate ( value_function_howard(points))
    allocate ( policy_function_howard(points))


    call cpu_time(start_vfi)

    call value_interation(beta,delta,a,alpha,points,initial_guess,value_function,policy_function)!,grid_k)

    call cpu_time(finish_vfi)
    print '("Time it took for VFI was = ",f6.3," seconds.")',finish_vfi-start_vfi



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!PRINT THE RESULTS DEVIATIONS AND TIME IT TOOK !!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !call analytical(a,alpha,beta,delta,anal_value_function,anal_policy_function,points)


    !print*, "the max difference in policy functions is ",MAXVAL(abs(anal_policy_function-policy_function))
    !print*, "the max difference in value functions is ",MAXVAL(abs(anal_value_function-value_function))
    
    call cpu_time(start_MPB_one)
    call value_interation_c(beta,delta,a,alpha,points,initial_guess,value_function_c,policy_function_c,1)
    call cpu_time(finish_MPB_one)
    print '("Time it took for MacQueen-Porteus every 1 was = ",f6.3," seconds.")',finish_MPB_one-start_MPB_one
    
    call cpu_time(start_MPB_five)
    call value_interation_c(beta,delta,a,alpha,points,initial_guess,value_function_fiv,policy_function_fiv,modul)
    call cpu_time(finish_MPB_five)
    print '("Time it took for MacQueen-Porteus every 5 was = ",f6.3," seconds.")',finish_MPB_five-start_MPB_five
    
    

    call cpu_time(start_how)

    call policy_interation_howard(howard_iterations,beta,delta,a,alpha,points,initial_guess,value_function_howard,policy_function_howard)


    call cpu_time(finish_how)
    print '("Time it took for Howard step every 5 was =  = ",f6.3," seconds.")',finish_how-start_how
    
    
    
    value_function(1)=0.0
    value_function_c(1)=0.0
    value_function_fiv(1)=0.0
    value_function_howard(1)=0.0
    
    
    !print*, "policy_function is ",policy_function
    !print*, "initial guess is ",initial_guess
    !print*, "value_function is ",value_function_c
    
    
    
    
    print*, "the max difference in policy functions between MacQueen-Porteus every 1 and VFI is ",MAXVAL(abs(policy_function_c-policy_function))
    print*, "the max difference in value functions  between MacQueen-Porteus every 1 and VFI is ",MAXVAL(abs(value_function_c-value_function))


    print*, "the max difference in policy functions between MacQueen-Porteus every 5 and VFI is ",MAXVAL(abs(policy_function_fiv-policy_function))
   print*, "the max difference in value functions  between MacQueen-Porteus every 5 and VFI is ",MAXVAL(abs(value_function_fiv-value_function))
    
    
        print*, "the max difference in policy functions between Howard steps and VFI is ",MAXVAL(abs(policy_function_howard-policy_function))
   print*, "the max difference in value functions  between Howard steps and VFI is ",MAXVAL(abs(value_function_howard-value_function))

end program main