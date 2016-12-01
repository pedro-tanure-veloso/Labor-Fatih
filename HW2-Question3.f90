!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!DEFINE THE G FUNCTION !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Change this to find for a different function


double precision function func_g(p)
    IMPLICIT NONE
    double precision ans, p

    ans=.5*(p**(-.5)+p**(-.2))

    func_g=ans

END function func_g



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!DEFINE THE F FUNCTION !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Makes the problem about finding a root 
!!Not necessary but useful 


double precision function func_f(p,target_value)
    IMPLICIT NONE
    double precision ans, p,target_value
    double precision, external :: func_g
    ans=func_g(p)-target_value                

    func_f=ans

END function func_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Bisection  Method!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Will take mid point of the two extremes and until it gets close enough to the root



SUBROUTINE bisection(tol,start,ending,iterations,x_hat,target_value)

    IMPLICIT NONE
    
    integer iterations
    double precision tol,start,ending,x_hat,mid_point_func,mid_point,n_start,n_end,n_start_fun,target_value
    double precision, external :: func_f
    
        
    n_start=start                                    !Will hold the two end points of the Interval 
                                                     !the loop is working with are keeping track of
    n_end=ending
    
    
    iterations=1                                    !keeps track of iterations
    mid_point=(n_start+n_end)/2.0                   !holds the midpoint between 
    mid_point_func=func_f(mid_point,target_value)   ! holds f evaluated at the midpoint

    
    n_start_fun=func_f(n_start,target_value)

    
    if((func_f(start,target_value)*func_f(ending,target_value))>0.0) then !checks inital conditions are well definied
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*, "!!! INVALID INITIAL CONDITIONS !!!!!"
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        return
    end if

    do while(abs(mid_point_func)>tol)            !while its not close enough to zer0
        
       
        iterations=iterations+1                  !update the iteration
        
        
    
        if(n_start_fun*mid_point_func>0.0) then  !if f(start) and f(midpoint) have same sign, take it to be the new n_start
            n_start=mid_point                    !creates new interval
        else 
            n_end=mid_point                      !otherwise make it the other endpoint in the interval
        end if
        
        
        mid_point=(n_start+n_end)/2.0            !find new midpoint 
        mid_point_func=func_f(mid_point,target_value)  !find new f(midpoint )
        n_start_fun=func_f(n_start,target_value) !find new comparison point
        
    
    end do  

    x_hat=mid_point                              !return the answer 


END SUBROUTINE bisection




program main

    IMPLICIT NONE
    integer iterations
    double precision tol,start,ending,x_hat,target_value
    double precision, external :: func_f
    
    
    
    target_value=.75    !target value we want for g(p)
    
    tol=10.0**(-6)      !tolerance value for g(p_estimate)-target
    
    
    start=10            !inital guess
    
    ending=1   !inital guess.... make sure g(start) and g(end) are both in different sides of target
    
    call  bisection(tol,start,ending,iterations,x_hat,target_value)
    
    
    print*,"p_hat is ",x_hat
    print*,"f of p_hat is ",func_f(x_hat,target_value)
    print*,"total num of iterations is",iterations

end program main