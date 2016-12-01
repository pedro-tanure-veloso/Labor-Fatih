double precision function g(p)
    implicit none
    double precision p
    
    g = .5*p**(-0.5)+.5*p**(-0.2)
    
end function g

double precision function diff(p,target_value)
    implicit none
    double precision p,target_value
    
    diff = .5*p**(-0.5)+.5*p**(-0.2) - target_value
    
end function diff

subroutine bissection(func,ftarget,x1,x2,tol,x3,it)
    implicit none
    external func
    integer it
    double precision func,x1,x2,err,ftarget,tol,x3,fint,mid
    
    it = 0
    err = 10.
        
    do while (err > tol)
        mid = (x1+x2)/2
        fint = func(mid,ftarget) ! "Intermediary" function
        err = abs(fint)
        !print *, "error is",err
        if (fint*func(x2,ftarget) > 0.) then
            x2 = mid
        else 
            x1 = mid
        end if
        it = it + 1
    end do
    
    x3 = mid
    print *, "Closest estimate for p is:", x3
    print *, "Number of iterations required to reach goal is:", it
end subroutine bissection

subroutine secant(func,ftarget,x1,x2,x3,tol,it)
    implicit none
    external func
    integer it
    double precision func,ftarget,tol,x1,x2,x3,x,err,beta,alpha
    
    it = 0
    err = 10.
    
    do while (err > tol)
        beta = (func(x1,ftarget)-func(x2,ftarget))/(x1-x2)
        alpha = func(x1,ftarget) - beta*x1
        x = -alpha/beta
        print *, "x", x
        err = abs(func(x,ftarget))
        if (err > tol) then
            x1 = x2
            x2 = x
        end if
        it = it + 1
    end do
    
    x3 = x
    print *, "Closest estimate for p is:", x3
    print *, "Number of iterations required to reach goal is:", it
    
    
end subroutine secant

! Starting the program
program Q3_4

implicit none
integer it
double precision, external:: g,diff
double precision x1,x2,x3,x4,ftarget,tol

tol = 10**(-6)
! Initial guesses: make sure x1<x2!
x1 = 0.1
x2 = 3.
! Targets
ftarget = 0.75

!! Applying the methods
! Bisection
call bissection(diff,ftarget,x1,x2,tol,x3,it)

! Secant
x1 = 0.1 ! We must reinitialize the guesses
x2 = 3.
call secant(diff,ftarget,x1,x2,x4,tol,it)


end program Q3_4