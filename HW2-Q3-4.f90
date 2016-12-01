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
    double precision func,x1,x2,err,ftarget,tol,x3,fint,mid,a,b
    
    it = 0
    err = 10.
    a = x1
    b = x2
        
    do while (err > tol)
        mid = (a+b)/2
        fint = func(mid,ftarget) ! "Intermediary" function
        err = abs(fint)
        !print *, "error is",err
        if (fint*func(b,ftarget) > 0.) then
            b = mid
        else 
            a = mid
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
    double precision func,ftarget,tol,x1,x2,x3,x,err,beta,alpha,a,b
    
    it = 0
    err = 10.
    a = x1
    b = x2
    
    do while (err > tol)
        beta = (func(a,ftarget)-func(b,ftarget))/(a-b) ! beta = delta(y)/delta(x)
        alpha = func(a,ftarget) - beta*a                 ! finding intercept
        x = -alpha/beta                                    ! finding x such that y=0
        err = abs(func(x,ftarget))
        if (err > tol) then
            a = b                ! stores last value of previous iteration as first one in the next
            b = x                 ! last value is now the previous x such that y=0
        end if
        it = it + 1
    end do
    
    x3 = x
    print *, "Closest estimate for p is:", x3
    print *, "Number of iterations required to reach goal is:", it
    
    
end subroutine secant

subroutine brent(func,ftarget,x1,x2,x3,tol,it)

    implicit none
    integer itmax,it,iter
    double precision func,ftarget,x1,x2,x3,tol,eps
    parameter (itmax=100,eps=3e-8)
    double precision a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    
    a = x1
    b = x2
    fa = func(a,ftarget)
    fb = func(b,ftarget)
    
    if ((fa > 0 .and. fb > 0) .or. (fa < 0 .and. fb < 0)) then
        print *, "root must be bracketed for Brent's method"
        stop
    end if
    
    c = b
    fc = fb
    
    do iter=1,itmax
        if (((fb > 0) .and. fc > 0) .or. (fb < 0 .and. fc < 0)) then
            c = a
            fc = fa
            d = b-a
            e = d
        end if
        if (abs(fc) < abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end if
        tol1 = 2.*eps*abs(b) + .5*tol
        xm = .5*(c-b)
        if (abs(xm) <= tol1 .or. fb == 0.) then
            x3 = b
            it = iter
            print *, "Closest estimate for p is:", x3
            print *, "Number of iterations required to reach goal is:", it
            return
        end if
        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s = fb/fa
            if (a == c) then
                p = 2.*xm*s
                q = 1.-s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                q = (q-1.)*(r-1.)*(s-1.)
            end if
            if (p > 0) then
                q = -q
            end if
            p = abs(p)
            if (2.*p < min(3.*xm*q-abs(tol1*q),abs(e*q))) then
                e = d
                d = p/q
            else
                d = xm
                e = d
            end if
        else
            d = xm
            e = d
        end if
        a = b
        fa = fb
        if (abs(d) > tol1) then
            b = b+d
        else
            b = b+sign(tol1,xm)
        end if
        fb = func(b,ftarget)
    end do
    print *, "Brent's method exceeding maximum iterations"
    x3 = b
    it = itmax
    print *, "Closest estimate for p is:", x3
    print *, "Number of iterations required to reach goal is:", it
        
end subroutine brent

! Starting the program
program Q3_4

implicit none
integer it
double precision, external:: g,diff
double precision x1,x2,x3,x4,x5,ftarget,tol

tol = 10**(-6)
! Initial guesses: make sure x1<x2!
x1 = 0.1
x2 = 3.
! Targets
ftarget = 0.75

!! Applying the methods
! Bisection
print *, "Bissection method"
call bissection(diff,ftarget,x1,x2,tol,x3,it)

! Secant
print *, "Secant method"
call secant(diff,ftarget,x1,x2,x4,tol,it)

! Brent's method
print *, "Brent's method"
call brent(diff,ftarget,x1,x2,x5,tol,it)


end program Q3_4