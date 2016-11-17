double precision function utility(c)
    IMPLICIT NONE
    double precision ans, c

    ans=log(c)

    utility=ans

END function utility

double precision function production_fun(k,delta,a,alpha)
    IMPLICIT NONE
    double precision k,delta,a,alpha,ans

    ans=a*k**alpha + (1-delta)*k

    production_fun=ans

END function production_fun

double precision function k_ss(delta,a,alpha, beta)
    IMPLICIT NONE
    double precision k,delta,a,alpha,ans, beta

    ans=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))

    k_ss=ans

END function k_ss


SUBROUTINE iterate(beta,delta,a,alpha,num_points,initial_guess, policy_funtion,grid_k,value_function)

    IMPLICIT NONE
    integer num_points,i,j,n
    double precision beta,delta,a,alpha,grid_k(num_points), produc,value_function(num_points)
    double precision initial_guess(num_points),tolerance,error,max_val(num_points)
    double precision policy_funtion(num_points),val_funtion(num_points),utils(num_points)
    double precision, external :: production_fun,k_ss,utility
    !double precision, dimension(num_points) ::iterate

    tolerance=.01
    error=100

    val_funtion=initial_guess

    max_val=-100
    n=0


do while (error>tolerance)


    do  i=1,num_points


        do  j=1,num_points

            produc= production_fun(grid_k(i),delta,a,alpha)

            if (produc-grid_k(j)>0) Then
                utils(j)= utility(produc-grid_k(j))+beta*val_funtion(j)
            else
                utils(j)=-1000000
            END IF

            if (utils(j)>max_val(i)) Then
                max_val(i)=utils(j)
                policy_funtion(i)=grid_k(j)
            END IF

         end  do

    end  do

    error=MAXVAL(abs(val_funtion-max_val))
    print*, "error is ", error
    val_funtion=max_val
    max_val=-100
    n=n+1
 
 end  do
print*, "number of iterations", n
value_function=val_funtion

END SUBROUTINE iterate


SUBROUTINE value_interation(beta,delta,a,alpha,num_points,initial_guess,value_function,policy_function)
IMPLICIT NONE
integer num_points,mid,i
double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), initial_guess(num_points),error
double precision equally_spaced,value_function(num_points),policy_function(num_points)
double precision, external :: production_fun,k_ss,utility

!double precision, dimension(num_points) ::value_interation
mid=num_points/2+1
k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))
error=10.0

equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))

grid_k(1)=0.001

do  i=2,num_points

    if (i==mid) Then
        grid_k(i)=k_at_ss
    else
        grid_k(i)=grid_k(i-1)+ equally_spaced
    END IF
end  do

call iterate(beta,delta,a,alpha,num_points,initial_guess, policy_function,grid_k,value_function)

end SUBROUTINE value_interation

SUBROUTINE analytical(a,alpha,beta,delta,VF,PF,num_points)
IMPLICIT NONE
integer num_points,mid,i
double precision equally_spaced, a, alpha, beta, delta, grid_k(num_points), ans, a0, a1, VF(num_points), PF(num_points), k_at_ss
double precision, external :: k_ss
  
mid=num_points/2+1
k_at_ss=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))
print*, "steady state level of capital is", k_at_ss
equally_spaced=(REAL(2.0*k_at_ss))/(REAL(num_points))

grid_k(1)=0.001

do  i=2,num_points

    if (i==mid) Then
        grid_k(i)=k_at_ss
    else
        grid_k(i)=grid_k(i-1)+ equally_spaced
    END IF
    end  do
    
a0 = (1/(1-beta))*((alpha*beta)/(1-alpha*beta)*log(alpha*beta)+log(1-alpha*beta))
print*, "a0 is",a0
a1 = alpha/(1-alpha*beta)
print*, "a1 is",a1

VF(1) = 0.0
PF(1) = 0.0

do i=2,num_points  
    VF(i) = a0 + a1*log(grid_k(i))
    PF(i) = alpha*beta*grid_k(i)**alpha
    end do
    
end SUBROUTINE analytical


program main

IMPLICIT NONE
integer points
real :: start, finish
double precision beta,delta,a,alpha
double precision, ALLOCATABLE :: initial_guess(:)
double precision, ALLOCATABLE :: value_function(:)
double precision, ALLOCATABLE :: policy_function(:)
double precision, ALLOCATABLE :: anal_value_function(:)
double precision, ALLOCATABLE :: anal_policy_function(:)

points=11
beta=.96
delta=.1
a=1.0
alpha=1.0/3.0

allocate ( initial_guess(points))
allocate ( value_function(points))
allocate ( policy_function(points))
allocate ( anal_value_function(points))
allocate ( anal_policy_function(points))

initial_guess=0.0

! Let's starts timing
call cpu_time(start)
              
! Now we do the iterations
call value_interation(beta,delta,a,alpha,points,initial_guess,value_function,policy_function)

if (delta==1) then
call analytical(a,alpha,beta,delta,anal_value_function,anal_policy_function,points)
end if

print*, "initial guess is ",initial_guess
print*, "value_function is ",value_function
print*, "policy_function is ",policy_function
if (delta==1) then
print*, "analytical_value_function is ",anal_value_function
print*, "analytical_policy_function is ",anal_policy_function
end if

call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start


end program main
