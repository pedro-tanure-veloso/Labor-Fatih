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

    ans = (c**(1.0-sigma))/(1.0-sigma) ! CRRA utility

    util3 = ans

END function util3



SUBROUTINE expo_grid(points, grid_x, start, finish, theta) 
    implicit none
    integer points, i
    double precision :: grid_x(points),start, finish,theta, equally_spaced
    double precision :: grid_x_hat(points)
    
    
    equally_spaced=(1.0)/(REAL(points))

    grid_x_hat(1)=0.0
   

    do  i=2,points
        grid_x_hat(i)=grid_x(i-1)+ equally_spaced
    end  do
    
    do  i=2,points
        grid_x(i)=start+(finish-start)*(grid_x_hat(i)**theta)
    end  do  
    
    
    
end SUBROUTINE expo_grid






SUBROUTINE linear_inter(x1,y1,x2,y2,xval,yval) 
    double precision :: x1,y1,x2,y2,xval,yval
    double precision :: frac
    
    
    
    frac = ( xval - x1 ) / ( x2 - x1 )
    
    
    yval = y1 + frac * ( y2 - y1 )
    
    
end SUBROUTINE linear_inter




SUBROUTINE linear_interpolation(start,finish,num_points,factor,grid_real_y,grid_x,grid_interpol_y, which)


    double precision start, finish, equally_spaced,error,max_error
    integer num_points, factor, track, i, last, next,which
    double precision :: grid_real_y(num_points)
    double precision :: grid_x(num_points)
    double precision :: grid_interpol_y(num_points)
    double precision sigma
    
    
    

    track=1
    last=1
    next=last+factor-1
    
    equally_spaced=((finish-start))/(REAL(num_points))

    grid_x(1)=start
   

    do  i=2,num_points
        grid_x(i)=grid_x(i-1)+ equally_spaced
    end  do
    if(which==1) then 

        print*,"=================================="
        print*,"=======log utility================"
        print*,"=================================="

        grid_real_y(num_points)=util1(grid_x(num_points))
        grid_interpol_y(num_points)=util1(grid_x(num_points))
        grid_real_y(1)=util1(grid_x(1))
        grid_interpol_y(1)=util1(grid_x(1))

        grid_real_y(next)=util1(grid_x(next))








        do  i=2,num_points-1
            grid_real_y(i)=util1(grid_x(i))
            if (track>factor-1) then 
                grid_interpol_y(i)=grid_real_y(i)
                track=1
                last=i
                if( last+factor-1<num_points) then
                    next=last+factor-1
                    grid_real_y( next)=util1(grid_x( next))
                 else 
                     next=num_points
                 end if 
            else 
               call linear_inter(grid_x(last),grid_real_y(last),grid_x( next),grid_real_y( next),grid_x(i),grid_interpol_y(i))



               track=track+1
            end if 

        end  do




    else if(which==2) then 

        print*,"=================================="
        print*,"=======square root utility========"
        print*,"=================================="

        grid_real_y(num_points)=util2(grid_x(num_points))
        grid_interpol_y(num_points)=util2(grid_x(num_points))
        grid_real_y(1)=util2(grid_x(1))
        grid_interpol_y(1)=util2(grid_x(1))

        grid_real_y(next)=util2(grid_x(next))








        do  i=2,num_points-1
            grid_real_y(i)=util2(grid_x(i))
            if (track>factor-1) then 
                grid_interpol_y(i)=grid_real_y(i)
                track=1
                last=i
                if( last+factor-1<num_points) then
                    next=last+factor-1
                    grid_real_y( next)=util2(grid_x( next))
                 else 
                     next=num_points
                 end if 
            else 
               call linear_inter(grid_x(last),grid_real_y(last),grid_x( next),grid_real_y( next),grid_x(i),grid_interpol_y(i))



               track=track+1
            end if 

        end  do
    
    else if(which==3) then 

        print*,"=================================="
        print*,"=======CRRA utility==============="
        print*,"=================================="


        sigma=2.0
        
        
        grid_real_y(num_points)=util3(grid_x(num_points),sigma)
        grid_interpol_y(num_points)=util3(grid_x(num_points),sigma)
        grid_real_y(1)=util3(grid_x(1),sigma)
        grid_interpol_y(1)=util3(grid_x(1),sigma)

        grid_real_y(next)=util3(grid_x(next),sigma)





        do  i=2,num_points-1
            grid_real_y(i)=util3(grid_x(i),sigma)
            if (track>factor-1) then 
                grid_interpol_y(i)=grid_real_y(i)
                track=1
                last=i
                if( last+factor-1<num_points) then
                    next=last+factor-1
                    grid_real_y( next)=util3(grid_x(next),sigma)
                 else 
                     next=num_points
                 end if 
            else 
               call linear_inter(grid_x(last),grid_real_y(last),grid_x( next),grid_real_y( next),grid_x(i),grid_interpol_y(i))



               track=track+1
            end if 

        end  do
        
    end if
        
    
    max_error=-1
        
    do  i=1,num_points
        error=abs((grid_real_y(i)-grid_interpol_y(i))/(grid_real_y(i)))
        if(error>max_error) then
            max_error=error
        end if
    end  do
    
    print*, "max error is ",max_error

    
    
    
    
    !do  i=1,num_points
    !    print*,grid_real_y(i)
    !end  do
end SUBROUTINE linear_interpolation

program main


    double precision start, finish, equally_spaced,error,max_error
    integer num_points, factor, track, i, last, next, which
    double precision, ALLOCATABLE :: grid_real_y(:)
    double precision, ALLOCATABLE :: grid_x(:)
    double precision, ALLOCATABLE :: grid_interpol_y(:)
    
    
    which=3                            !!!! 1 is log utility, 2  is square root utility, and CRRA utility
    start=0.05
    finish=2.0
    num_points=10000
    factor=222
    track=1
    last=1
    next=last+factor-1
    
    equally_spaced=((finish-start))/(REAL(num_points))
    
    
    allocate ( grid_real_y(num_points))
    allocate ( grid_interpol_y(num_points))
    allocate ( grid_x(num_points))
    
    grid_x=0
    grid_interpol_y=0
    grid_real_y=0
    
    
    call linear_interpolation(start,finish,num_points,factor,grid_real_y,grid_x,grid_interpol_y, which)
    

end program main
