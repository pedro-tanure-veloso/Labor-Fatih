! Question 4a - Homework 1
! Here, we will interpolate three types of utility functions using 3 methods: linear, Chebyshev polynomials and cubic splines.
! The functions are: log, square root and standard CRRA. We will vary the relative risk aversion coefficient as well.

! Step 1: define the utility functions

double precision function util1(c)
    IMPLICIT NONE
    double precision ans, c

    ans = log(c) ! log utility

    util1 = ans

END function util1

double precision function util2(c)
    IMPLICIT NONE
    double precision ans, c

    ans = c**.5 ! square root utility

    util2 = ans

END function util2

double precision function util3(c,sigma)
    IMPLICIT NONE
    double precision ans, c, sigma

    ans = (c**(1-sigma))/(1-sigma) ! CRRA utility

    util3 = ans

END function util3

! Step 2: interpolaion subroutines (these ones come from Numerical Recipes in Fortran 77)

SUBROUTINE polint(xa,ya,n,x,y,dy)
    INTEGER n, NMAX
    real dy,x,y,xa(n),ya(n)
    PARAMETER (NMAX=10) ! Largest anticipated value of n
        ! Given arrays xa and ya, each of length n, and given a value x, routine returns value y, and error estimate dy
        ! If P(x) is a polynomial of degree N-1 such that P(xa_i)=ya_i, returned value is y=P(x)
    INTEGER i,m,ns
    REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

    ns = 1
    dif = abs(x-xa(1))
    
    do i=1,n              ! Here we find the index ns of the closest table entry
        dift - abs(x-xa(i))
        if (dift < dif) then
            ns = i
            dif = dift
        endif
        c(i) = ya(i)     ! and initialize the tableau of c's and d's
        d(i) = ya(i)
    end do
    
    y = ya(ns)           ! This is the initial approximation to y
    ns = ns-1
    do m=1,n-1           ! For each column of the tableau, we loop over the current c's and d's and update them
        do i=1,n-1
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1)-d(i)
            den = ho-hp
            if (den == 0) pause 'failure point'  ! this error can occur only if two input xa's are (to within roundoff) identical
                den = w/den
                d(i) = hp*den   ! Here, ths c's and d's are updated
                c(i) = ho*den
        end do
        if (2*ns < n-m) then
            dy = c(ns+1)
        else
            dy = d(ns)
            ns = ns-1
        end if
        y = y+dy
    end do
    return
END SUBROUTINE polint

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    INTEGER n,NMAX
    REAL yp1,ypn,x(n),y(n),y2(n)
    PARAMETER (NMAX=500)
        ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i=f(x_i), with
        ! x_1<x_2<...<x_N, and given values yp1 and ypn for the first derivative of the interpolating function
        ! at points 1 and n, respectively, this routine returns an array y2(1:n) of length n which contains the
        ! second derivatives of the interpolating function at the tabulated points x_i. If yp1 and/or ypn are
        ! equal to 10^3 or larger, the routine is signaled to set the corresponding boundary condition for a natural
        ! spline, with zero second derivative on that boundary.
        ! Parameter NMAX is the largest anticipated value of n.
    INTEGER i,k
    REAL p,qn,sig,un,u(NMAX)
    
    
    if (yp1 > 99e30) then
        y2(1) = 0.        ! The lower boundary condition is the set either to be "natural"
        u(1) = 0.         ! That is, 2nd derivatives at the boundaries are set to zero
    else                  ! ...or else to have a specified first derivative
        y2(1) = -0.5
        u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1)-yp1)
    endif
    
    do i=2,n-1                                ! This is the decomposition loop of the tridiagonal algorithm.
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))   ! y2 and u are used for temporary storage of the decomposed factors
        p = sig*y2*i-1)+2
        y2(i) = (sig-1)/p
        u(i) = (6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
*                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn > 99e30) then                   ! The upper boundary condition is set either to be "natural"
        qn = 0.
        un = 0.
    else                                    ! or else to have a specified first derivative
        qn = 0.5
        un = (3./(x(n)-x(n-1)))*((ypn-(y(n)-y(n-1))/(x(n)-x(n-1))
    endif
    
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
    
    do k=n-1,1,-1                         ! This is the backsubstitution loop of the tridiagonal algorithm
        y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    
    return
END

SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER n
    REAL x,y,xa(n),y2a(n),ya(n)
        ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's in order),
        ! and given the array y2a(1:n), which is the output from spline above., and given a value of x, this routine
        ! returns a cubic-spline interpolated value y
    INTEGER k,khi,klo
    REAL a,b,h           
    
    klo = 1           ! We will find the right place in the table by means of bisection. This is optimal if sequential
    khi=n             ! calls to this routine are at random values of x. If sequential calls are in order, and closely
    if (khi-klo > 1) then     ! spaced, one would do better to store previous values of klo and khi and test if they
        k = (khi+klo)/2       ! remain appropriate on the next call
        if (xa(k) > x) then
            khi = k
        else
            klo = k
        endif
    endif     ! klo and khi now bracket the inout value of x
    h = xa(khi)-xa(klo)
    if (h == 0) pause 'bad xa input in splint'   ! the xa's must be distinct
    a = (xa(khi)-x)/h            ! Cubic spline in polynomial is now evaluated
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
    
    return
END

SUBROUTINE chebtf(a,b,c,n,func)

    INTEGER n,NMAX
    REAL a,b,c(n),func
    DOUBLE PRECISION PI
    EXTERNAL func
    PARAMETER (NMAX=50,PI=3.141592653589793d0)
        ! Chebyshev fit: given a function func, lower and upper limits of the interval [a,b], and
        ! a maximum degree n, this routine computes the n coefficients c_k such that func(x)~~
        ! [Sum_(k=1)^n c_kT_(k-1)(y)]-c_1/2, where y and x are related by y=(x-(b+a)/2)/((b-a)/2).
        ! This routine is to be used with moderately large n (e.g., 30 or 50), the array of c's
        ! subsequently to be truncated at the smaller value m such that c_(m+1) and subsequent
        ! elements are negligible.
        ! Parameters: Maximum expected value of n, and pi
    INTEGER j,k
    REAL bma,bpa,fac,y,f(NMAX)
    DOUBLE PRECISION sum
    
    bma = 0.5*(b-a)
    bpa = 0.5*(b+a)
    do k=1,n            ! We evaluate the function at the n points required by eq. 5.8.7 in the book
        y = cos(PI*(k-0.5d0)/n
        f(k) = func(y*bma+bpa)
    end do
    fac = 2./n
    do j=1,n            ! We will accumulate the sum in double precision, a nicety that you can ignore
        sum = 0.d0
        do k=1,n
            sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
        end do
        c(j) = fac*sum
    end do
    return
END

FUNCTION chebev(a,b,c,m,x)
    INTEGER m
    REAL chebev,a,b,x,c(m)
        ! CHebyshev evaluation: All arguments are input. c(1:m) is an array of Chebyshev coefficients,
        ! the first m elements of c output from chebtf (which must have been called with the same a 
        ! and b). The Chebyshev polynomial Sum_(k=1)^m c_kT_(k-1)(y)-c_1/2 is evaluated at a point
        ! y=(x-(b+a)/2)/((b-a)/2), and the result is returned as the function value.
    INTEGER j
    REAL d,dd,sv,y,y2
    
    if ((x-a)*(x-b) > 0) pause 'x not in range in chebev'
    d = 0.
    dd = 0.
    y = (2.*x-a-b)/(b-a)       ! Change of variable
    y2 = 2.*y
    do j=m,2,-1                ! Clenshaw's recurrence
        sv = d
        d = y2*d-dd+c(j)
        dd = sv
    end do
    chebev = y*d-dd+0.5*c(1)  ! Last step is different
    return
END

program main

IMPLICIT none
integer n
real alpha,x(n),y(n),grid(n)

grid(1) = 0.05
grid(n) = 2

! equally spaced grid
do i=2,n-1
    grid(i) = 

end program main