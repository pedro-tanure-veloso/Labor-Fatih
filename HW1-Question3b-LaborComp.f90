! Question 3b - Homework 1

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


SUBROUTINE howard(it,beta,delta,a,alpha,num_points,initial_guess,policy_function,grid_k,value_function)

    IMPLICIT NONE
    integer num_points,i,j,n,m,it !"it" is the number of iterations that we will have in Howard's improvement algorithm
    double precision beta,delta,a,alpha,grid_k(num_points), produc,value_function(num_points),kp_index(num_points)
    double precision initial_guess(num_points),tolerance,error,v0(num_points),v1(num_points)
    double precision policy_function(num_points),val_funtion(num_points),utils(num_points)
    double precision, external :: production_fun,k_ss,utility
    
    tolerance = .01
    error = 100

    val_funtion = initial_guess

    v0 = -100
    n = 0


do while (error>tolerance)


    do  i=1,num_points
        produc= production_fun(grid_k(i),delta,a,alpha)
        do  j=1,num_points
            if (produc-grid_k(j)>0) Then
                utils(j) = utility(produc-grid_k(j))+beta*val_funtion(j) ! getting the "value function"
            else
                utils(j) = -1000000
            END IF
            if (utils(j)>v0(i)) Then
                v0(i) = utils(j)
                policy_function(i) = grid_k(j) ! gets the policy function
                kp_index(i) = j
            END IF
        end  do
    end  do
    
    do m=1,it
           do i=1,num_points
               produc= production_fun(grid_k(i),delta,a,alpha)
               v1(i) = utility(produc-policy_function(i)) + beta*v0(int(kp_index(i))) ! does the policy iteration
           end do
           v0 = v1
    end do
    
    error = MAXVAL(abs(val_funtion - v0))
    print*, "error is ", error
    val_funtion = v0
    v0 = -100
    n = n+1
end  do
print*, "number of iterations", n
value_function=val_funtion

END SUBROUTINE howard


SUBROUTINE policy_interation(it,beta,delta,a,alpha,num_points,initial_guess,value_function,policy_function)
IMPLICIT NONE
integer num_points,mid,i,it
double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), initial_guess(num_points),error
double precision equally_spaced,value_function(num_points),policy_function(num_points)
double precision, external :: production_fun,k_ss,utility

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

call howard(it,beta,delta,a,alpha,num_points,initial_guess, policy_function,grid_k,value_function)

end SUBROUTINE policy_interation


program main

IMPLICIT NONE
integer points, howard_iterations
real :: start, finish
double precision beta,delta,a,alpha
double precision, ALLOCATABLE :: initial_guess(:)
double precision, ALLOCATABLE :: value_function(:)
double precision, ALLOCATABLE :: policy_function(:)

howard_iterations = 5
points = 11
beta = .96
delta = .1
a=1.0
alpha=1.0/3.0

allocate ( initial_guess(points))
allocate ( value_function(points))
allocate ( policy_function(points))

initial_guess=0.0

! Let's starts timing
call cpu_time(start)
              
! Now we do the iterations
call policy_interation(howard_iterations,beta,delta,a,alpha,points,initial_guess,value_function,policy_function)

print*, "initial guess is ",initial_guess
print*, "value_function is ",value_function
print*, "policy_function is ",policy_function

call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start

end program main