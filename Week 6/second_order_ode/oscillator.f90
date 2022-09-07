program oscillator
    implicit none

    real*8:: x_0, v_0, t_0, x, v, t, t_n, dt
    real*8:: f, g

    real*8:: f0, f1, f2, f3, g0, g1, g2, g3
    real*8:: x_temp1, x_temp2, x_temp3, v_temp1, v_temp2, v_temp3

    integer:: i, n_iter

    x_0 = 0.100 ; v_0 = 1.00 ; t_0 = 0.00 ; t_n = 50
    dt = 0.01
    

    x = x_0; v = v_0; t = t_0
    n_iter = int(t_n/dt)
    open(11, file = "output/oscillator.plt")
    write(11,*) t, ',', x, ',', v


    do i = 1, n_iter
        f0 = f(x, v)
        x_temp1 = x + (dt/2)*f0
        g0 = g(x, v)
        v_temp1 = v + (dt/2)*g0

        f1 = f(x_temp1, v_temp1)
        x_temp2 = x + (dt/2)*f1
        g1 = g(x_temp1, v_temp1)
        v_temp2 = v + (dt/2)*g1

        f2 = f(x_temp2, v_temp2)
        x_temp3 = x + dt*f2
        g2 = g(x_temp2, v_temp2)
        v_temp3 = v + dt*g2

        f3 = f(x_temp3, v_temp3)
        g3 = g(x_temp3, v_temp3)

        ! PRINT*, g0, g1, g2, g3

        t = t + dt
        x = x + (dt/6)*(f0 + 2*f1 + 2*f2 + f3)
        v = v + (dt/6)*(g0 + 2*g1 + 2*g2 + g3)

        write(11,*) t, ',', x, ',', v
        PRINT*, "t = ",t, ", x = ", x, ", v = ", v

    end do

end program oscillator

!---------------------------------------------
!   FUNCTIONS
!---------------------------------------------
function f(x, v)
    implicit none 
    real*8:: x, v
    real*8:: f 
    
    f = v

end function f


function g(x, v)
    implicit none 
    real*8:: x, v
    real*8:: g 

    g = -1*sin(x)

    ! PRINT*, "sin(x) = ", g

end function g