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
    double precision k,delta,a,alpha,ans
    
    ans=((1.0/(alpha*a))*(1.0/beta-(1-delta)))**(1.0/(alpha-1))
    
    k_ss=ans
    
END function k_ss



double precision function iterate(beta,delta,a,alpha,num_points,initial_guess, policy_funtion)
    
    IMPLICIT NONE
    integer num_points,mid, equally_spaced,i
    double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), produc
    double precision initial_guess(num_points),tolerance,error,
    double precision policy_funtion(num_points),val_funtion(num_points),utils(num_points)  
    double precision, external :: production_fun,k_ss,utility

    tolerance=.0001
    error=100
    
    val_funtion=initial_guess
    
    do while (error>tolerance)
    
        do  i=1,num_points
            do  j=1,num_points
               
                produc= production_fun(grid_k(i),delta,a,alpha)
                utils= utility(produc-grid_k(j))+beta*val_funtion(j)
 
            
            
        end  do
    
    
    end  do
    
END function iterate






double precision function value_interation(beta,delta,a,alpha,num_points,initial_guess)
    IMPLICIT NONE
    integer num_points,mid, equally_spaced,i
    double precision beta,delta,a,alpha,k_at_ss,grid_k(num_points), initial_guess  
    double precision, external :: production_fun,k_ss,utility
    
    


    mid=num_points/2+1
    k_at_ss=k_ss(delta,a,alpha, beta)
    error=10
    
    
    
    equally_spaced=(2*k_at_ss)/(REAL(num_points))
    
    grid_k(1)=0
    
    do  i=2,num_points
        
        if (i==mid) Then
            grid_k(i)=k_at_ss
        else 
            grid_k(i)=grid_k(i-1)+ equally_spaced
        END IF
    end  do
    

    
    
    


end function value_interation



program main

    double precision k

    k=testing(5)

end program main