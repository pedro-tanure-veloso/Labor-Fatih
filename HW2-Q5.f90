PROGRAM Q5
    USE AR1_discret
    IMPLICIT NONE
    INTEGER n(3),i,j
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
    
    call rouwenhorst(rho,sigma_y,P1,y1,s1)
    do i=1,n(1)
        do j=1,n(1)
            write (*,*) P1(i,j)
        end do
    end do
            

END PROGRAM