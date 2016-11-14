double precision function phi_power(n)
    IMPLICIT NONE
    double precision phi,ans
    integer n
     
    phi=((5.0**(.5)-1.0)/2.0)
    
    ans=phi**(n)
    
    phi_power=ans
    
END function phi_power


double precision function How_Close( a, b)
    IMPLICIT NONE
    double precision a, b ,close_
    
    
    close_=abs(a-b)
    
    print *,"power method is " ,a
    print *,"diff method is " ,b
    print *," and they are " ,close_, "close "
    
    How_Close=close_

end function How_Close




program main
    IMPLICIT NONE    
    double precision Method_diff(100), Method_pow(100), closeness
    INTEGER n 
    double precision, external :: phi_power, How_Close
    
   

    do  n=0,99
    
        Method_pow(n+1)= phi_power(n)
    
        if (n==0) Then
            Method_diff(n+1)=1
        else if (n==1) Then
            Method_diff(n+1)= ((5.0**(.5)-1.0)/2.0)
        else 
            Method_diff(n+1)=Method_diff(n-1)-Method_diff(n)
        END IF
        
        closeness=How_Close( Method_pow(n+1), Method_diff(n+1))

        
    end  do

end program main