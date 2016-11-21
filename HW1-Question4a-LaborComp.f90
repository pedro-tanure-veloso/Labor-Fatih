! Question 4a - Homework 1
! Here, we will interpolate three types of utility functions using 3 methods: linear, Chebyshev polynomials and cubic splines.
! The functions are: log, square root and standard CRRA. We will vary the relative risk aversion coefficient as well.

program main

IMPLICIT none
integer n,i,factor,natural,err
double precision equally_spaced, equally_spacedl, alpha, dif,theta,nat,maxd
double precision, ALLOCATABLE :: x(:),y1(:),y2(:),y3(:),xl(:),yl1(:),yl2(:),yl3(:)
double precision, ALLOCATABLE :: int1(:),int2(:),int3(:),dy1(:),yx1(:),yx2(:),yx3(:)
double precision, external :: util1,util2,util3
integer, parameter :: out_unit=20

maxd = 0.

! Use natural = 1 if you want natural spline
natural = 1
nat = 99e31

! Risk aversion parameter
alpha = 2.

! size of the grid
n = 900
factor = 10
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

! equally-spaced grid
x(1) = 0.05
x(factor*n) = 2.
equally_spaced = (x(factor*n)-x(1))/(factor*n)

do  i=2,factor*n-1
    x(i) = x(i-1)+ equally_spaced
end  do

! exponential grid
theta = 2.3
call expo_grid(n,xl,x(1),x(factor*n),theta)

! calculating utilities
do i=1,factor*n
    y1(i) = util1(x(i))
    y2(i) = util2(x(i))
    y3(i) = util3(x(i),alpha)
end do
do i=1,n
    yl1(i) = util1(xl(i))
    yl2(i) = util2(xl(i))
    yl3(i) = util3(xl(i),alpha)
end do

! doing the interpolation
! log
do i=1,factor*n
    if (natural == 1) then
        call spline(xl,yl1,n,nat,nat,yx1)
        else
        call spline(xl,yl1,n,1/xl(1),1/xl(n),yx1)
    endif
    call splint(xl,yl1,yx1,n,x(i),int1(i))
    dif = ABS(ABS(int1(i)-y1(i))/y1(i))*100
    if (dif > maxd) then
        maxd = dif
        err = i
    endif
end do

print *, 'log'
print *, 'maximum difference between interpolation and real utility is', maxd
print *, 'in x = ',x(err),'interpolated is', int1(err),' and actual is', y1(err)
print *, 'interpolated in the end is', int1(factor*n),' and actual is', y1(factor*n)
maxd = 0.

! sqrt
do i=1,factor*n
    if (natural == 1) then
        call spline(xl,yl2,n,nat,nat,yx2)
        else
        call spline(xl,yl2,n,0.5*xl(1)**(-0.5),0.5*xl(n)**(-0.5),yx2)
    endif
    call splint(xl,yl2,yx2,n,x(i),int2(i))
    dif = ABS(ABS(int2(i)-y2(i))/y2(i))*100
    if (dif > maxd) then
        maxd = dif
        err = i
    endif
end do

print *, 'sqrt'
print *, 'maximum difference between interpolation and real utility is', maxd
print *, 'in x = ',x(err),'interpolated is', int2(err),' and actual is', y2(err)
print *, 'interpolated in the end is', int2(factor*n),' and actual is', y2(factor*n)
maxd = 0.

! crra
do i=1,factor*n
    if (natural == 1) then
        call spline(xl,yl3,n,nat,nat,yx3)
        else
        call spline(xl,yl3,n,xl(1)**(-alpha),xl(n)**(-alpha),yx3)
    endif
    call splint(xl,yl3,yx3,n,x(i),int3(i))
    dif = ABS(ABS(int3(i)-y3(i))/y3(i))*100
    if (dif > maxd) then
        maxd = dif
        err = i
    endif
end do

print *, 'crra'
print *, 'maximum difference between interpolation and real utility is', maxd
print *, 'in x = ',x(err),'interpolated is', int3(err),' and actual is', y3(err)
print *, 'interpolated in the end is', int3(factor*n),' and actual is', y3(factor*n)

open (unit=out_unit,file="log.csv",action="write",status="replace")
write (out_unit,*) y1
write (out_unit,*) int1
close (out_unit)

open (unit=out_unit,file="sqrt.csv",action="write",status="replace")
write (out_unit,*) y2
write (out_unit,*) int2
close (out_unit)

open (unit=out_unit,file="crra.csv",action="write",status="replace")
write (out_unit,*) y3
write (out_unit,*) int3
close (out_unit)

end program main