! Question 4a - Homework 1
! Here, we will interpolate three types of utility functions using 3 methods: linear, Chebyshev polynomials and cubic splines.
! The functions are: log, square root and standard CRRA. We will vary the relative risk aversion coefficient as well.

!INCLUDE 'Question4-Routines.f90'

program main

!use Question4-Routines

IMPLICIT none
integer n,i,factor
double precision equally_spaced, equally_spacedl, alpha, dif
double precision, ALLOCATABLE :: x(:),y1(:),y2(:),y3(:),xl(:),yl1(:),yl2(:),yl3(:),int1(:),int2(:),int3(:),dy1(:),yx1(:),yx2(:),yx3(:)
double precision, external :: util1, util2, util3

! Risk aversion parameter
alpha = 2.

! size of the grid
n = 100
factor = 3
allocate (x(factor*n))
allocate (y1(factor*n))
allocate (y2(factor*n))
allocate (y3(factor*n))
allocate (xl(n))   ! this is the x grid with less points
allocate (yl1(n))   ! this is the y grid with less points
allocate (yl2(n))   
allocate (yl3(n))   ! this is the y grid with less points
allocate (yx1(n))   
allocate (yx2(n))   
allocate (yx3(n))   
allocate (int1(factor*n))
allocate (int2(factor*n))
allocate (int3(factor*n))
allocate (dy1(factor*n))

! equally spaced grid
x(1) = 0.05
xl(1) = 0.05
x(factor*n) = 2.
xl(n) = 2.
equally_spaced = (x(factor*n)-x(1))/(factor*n)
equally_spacedl = (xl(n)-xl(1))/n

do  i=2,factor*n-1
    x(i) = x(i-1)+ equally_spaced
end  do

do i=2,n-1
    xl(i) = xl(i-1)+ equally_spacedl
end  do

! calculating utilities
do i=1,factor*n
    y1(i) = util1(x(i))
    !print *, 'utility in consumption', x(i), 'is', y(i)
end do

do i=1,n
    yl1(i) = util1(xl(i))
    !print *, 'utility in consumption', xl(i), 'is', yl1(i)
end do

do i=1,factor*n
    y2(i) = util2(x(i))
    !print *, 'utility in consumption', x(i), 'is', y(i)
end do

do i=1,n
    yl2(i) = util2(xl(i))
    !print *, 'utility in consumption', xl(i), 'is', yl1(i)
end do

! doing the interpolation
! log
call spline(xl,yl1,n,1/xl(1),1/xl(n),yx1)
call splint(xl,yl1,yx1,n,x,int1)

print *, 'log'
do i=1,factor*n
    dif = ABS(int1(i)-y1(i))
    print *, 'difference between interpolation and real utility is', dif
end do

! sqrt
call spline(xl,yl2,n,0.5*xl(1)**(-0.5),0.5*xl(n)**(-0.5),yx2)
call splint(xl,yl2,yx2,n,x,int2)

print *, 'sqrt'
do i=1,factor*n
    dif = ABS(int2(i)-y2(i))
    print *, 'difference between interpolation and real utility is', dif
end do

! crra
call spline(xl,yl3,n,xl(1)**(-alpha),xl(n)**(-alpha),yx3)
call splint(xl,yl3,yx3,n,x,int3)

print *, 'crra'
do i=1,factor*n
    dif = ABS(int3(i)-y3(i))
    print *, 'difference between interpolation and real utility is', dif
end do

end program main