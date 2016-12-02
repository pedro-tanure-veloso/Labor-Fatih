INCLUDE 'AR1_discret.f90'

PROGRAM Q5
    USE AR1_discret
    IMPLICIT NONE
    INTEGER n(3),i
    INTEGER, PARAMETER :: out_unit = 20
    DOUBLE PRECISION rho,sigma_y,sigma_eta
    DOUBLE PRECISION, ALLOCATABLE :: y1(:),y2(:),y3(:),P1(:,:),P2(:,:),P3(:,:),s1(:),s2(:),s3(:)
    
    n(1) = 5
    n(2) = 11
    n(3) = 21
        
    ! Preallocating stuff
    allocate (y1(n(1)))
    allocate (y2(n(2)))
    allocate (y3(n(3)))
    allocate (P1(n(1),n(1)))
    allocate (P2(n(2),n(2)))
    allocate (P3(n(3),n(3)))
    allocate (s1(n(1)))
    allocate (s2(n(2)))
    allocate (s3(n(3)))
    
    ! Parameters, part 1
    rho = .9
    sigma_eta = .2
    sigma_y = sqrt(sigma_eta**2/(1-rho**2))
    
    ! 5-state Markov Chain
    open (unit=out_unit,file="5stMC9.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P1,y1,s1)
    do i=1,n(1)
        write(out_unit,*) P1(i,1),',',P1(i,2),',',P1(i,3),',',P1(i,4),',',P1(i,5),&
                ',', y1(i), ',', s1(i)
    end do
    close(out_unit)
    
    
    ! 11-state Markov Chain
    open (unit=out_unit,file="11stMC9.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P2,y2,s2)
    do i=1,n(2)
        write(out_unit,*) P2(i,1),',',P2(i,2),',',P2(i,3),',',P2(i,4),',',P2(i,5),&
                ',',P2(i,6),',',P2(i,7),',',P2(i,8),',',P2(i,9),',',P2(i,10),',',P2(i,11),&
                ',', y2(i), ',', s2(i)
    end do
    close(out_unit)
    
    ! 21-state Markov Chain
    open (unit=out_unit,file="21stMC9.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P3,y3,s3)
    do i=1,n(3)
        write(out_unit,*) P3(i,1),',',P3(i,2),',',P3(i,3),',',P3(i,4),',',P3(i,5),&
                ',',P3(i,6),',',P3(i,7),',',P3(i,8),',',P3(i,9),',',P3(i,10),',',P3(i,11),&
                ',',P3(i,12),',',P3(i,13),',',P3(i,14),',',P3(i,15),',',P3(i,16),',',P3(i,17),&
                ',',P3(i,18),',',P3(i,19),',',P3(i,20),',',P3(i,21),&
                ',', y3(i), ',', s3(i)
    end do
    close(out_unit)
    
    ! Parameters, part 2
    rho = .999
    sigma_eta = .2
    sigma_y = sqrt(sigma_eta**2/(1-rho**2))
    
    ! 5-state Markov Chain
    open (unit=out_unit,file="5stMC999.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P1,y1,s1)
    do i=1,n(1)
        write(out_unit,*) P1(i,1),',',P1(i,2),',',P1(i,3),',',P1(i,4),',',P1(i,5),&
                ',', y1(i), ',', s1(i)
    end do
    close(out_unit)
    
    
    ! 11-state Markov Chain
    open (unit=out_unit,file="11stMC999.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P2,y2,s2)
    do i=1,n(2)
        write(out_unit,*) P2(i,1),',',P2(i,2),',',P2(i,3),',',P2(i,4),',',P2(i,5),&
                ',',P2(i,6),',',P2(i,7),',',P2(i,8),',',P2(i,9),',',P2(i,10),',',P2(i,11),&
                ',', y2(i), ',', s2(i)
    end do
    close(out_unit)
    
    ! 21-state Markov Chain
    open (unit=out_unit,file="21stMC999.csv",action="write",status="replace")
    call rouwenhorst(rho,sigma_y,P3,y3,s3)
    do i=1,n(3)
        write(out_unit,*) P3(i,1),',',P3(i,2),',',P3(i,3),',',P3(i,4),',',P3(i,5),&
                ',',P3(i,6),',',P3(i,7),',',P3(i,8),',',P3(i,9),',',P3(i,10),',',P3(i,11),&
                ',',P3(i,12),',',P3(i,13),',',P3(i,14),',',P3(i,15),',',P3(i,16),',',P3(i,17),&
                ',',P3(i,18),',',P3(i,19),',',P3(i,20),',',P3(i,21),&
                ',', y3(i), ',', s3(i)
    end do
    close(out_unit)
    
    ! From here on, simulations will be done in Python

END PROGRAM