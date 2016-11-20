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




double precision function linear_inter(x1,y1,x2,y2,xval) result(yval)
    double precision ,intent(in) :: x1,y1,x2,y2,xval
    double precision :: frac
    frac = ( xval - x1 ) / ( x2 - x1 )
    yval = y1 + frac * ( y2 - y1 )
end function linear_inter


program main


    double precision start, finish, equally_spaced
    integer num_points, factor, track, i, last, next
    double precision, ALLOCATABLE :: grid_real_y(:)
    double precision, ALLOCATABLE :: grid_x(:)
    double precision, ALLOCATABLE :: grid_interpol_y(:)
    
    
    
    start=0.05
    finish=2.0
    num_points=1000
    factor=5
    track=1
    last=1
    next=last+factor-1
    
    equally_spaced=((start-finish))/(REAL(num_points))
    
    
    allocate ( grid_real_y(num_points))
    allocate ( grid_interpol_y(num_points))
    allocate ( grid_x(num_points))
    

    grid_x(1)=start
   

    do  i=2,num_points
        grid_x(i)=grid_x(i-1)+ equally_spaced
    end  do

    grid_real_y(num_points)=util2(grid_x(num_points))
    grid_interpol_y(num_points)=util2(grid_x(num_points))


    do  i=1,num_points-1
        grid_real_y(i)=util2(grid_x(i))
        if (track>factor-1) then 
            grid_interpol_y(i)=grid_real_y(i)
            track=1
            last=i
            if( last+factor-1<num_points) then
                next=last+factor-1
                grid_real_y(next)=util2(grid_x(next))
             else 
                 next=num_points
             end if 
        else 
           grid_interpol_y(i)=linear_inter(grid_x(last),grid_real_y(last),grid_x(next),grid_real_y(next),grid_x(i))
           track=track+1
        end if 
            
    end  do
    




end program main