SUBROUTINE qdratic_formula (a, b, c, ans_1, ans_2)
    IMPLICIT NONE
    double precision a, b ,c, ans_1,ans_2
    
    
    
    ans_1=((-b)+(b**(2.0)-(4.0*a*c))**.5)/(2.0*a)
    ans_2=((-b)-(b**(2.0)-(4.0*a*c))**.5)/(2.0*a)  
    
    
END SUBROUTINE qdratic_formula


SUBROUTINE new_way (a, b, c, ans_1, ans_2)
    IMPLICIT NONE
    double precision a, b ,c, ans_1,ans_2,q
    
    if (b<0) Then
        q=(-.5)*(b-(b**(2.0)-4.0*a*c)**.5)
    else 
        q=(-.5)*(b+(b**(2.0)-4.0*a*c)**.5)
    END IF
    
    
    ans_2=q/a
    ans_1=c/q  
    
    
END SUBROUTINE new_way

double precision function How_Close( a, b)
    IMPLICIT NONE
    double precision a, b ,close_
    
    
    close_=abs(a-b)
    
    print *,"quadratic formula is " ,a
    print *,"the new way is " ,b
    print *," and they are " ,close_, "close "
    
    How_Close=close_

end function How_Close


program main

    double precision a, b, c(8),d, Method_1(8,2), Method_2(8,2), closeness
    
    
    d = 10.0
    a = 1.0
    b = 100000.0
    
    do  i=1,8
        c(i) = d**(-i)
        call qdratic_formula (a,b,c(i),Method_1(i,0),Method_1(i,1))
        call new_way (a,b,c(i),Method_2(i,0),Method_2(i,1))
        
        print *," for n equals ", i, "the roots are::"

        closeness=How_Close(Method_2(i,0), Method_1(i,0))

        
        
    end  do



end program main