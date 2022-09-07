program DE_solver
    implicit none

    real*8:: x_0, y_0, x, y, y_temp, dx, x_n
    real*8:: fxy, f_0, f_1
    
    integer i, n_iter

    real*8:: dy
    real*8:: y_a = 92.620 
    
    ! Functions
    real*8:: func, y_update_rk4

    character (len = 100) :: solver
    print *, "Enter solver type"
    read *, solver

    x_0 = 0; y_0 = 0; x_n = 1.56
    x = x_0; y = y_0

    if (solver == 'euler') then
        
        dx = 0.001
        n_iter = int (x_n/dx) + 1

        open(11, file = "output/euler_soln.plt")
        write(11,*) x , ',', y

        do i = 1,n_iter
            fxy = func(x, y)
            y = y + dx*fxy
            x = x + dx
            write(11,*) x, ',', y
            ! PRINT*, x, ',', y
            
        end do

        dy = y_a - y
        PRINT*, "x = ", x, ",  dy = y_a - y_e = ", dy    
        close(11)

    endif


    if (solver == 'modified_euler') then
        
        dx = 0.001
        n_iter = int (x_n/dx) + 1
        
        open(12, file = "output/modified_euler_soln.plt")
        write(12,*) x, ',', y
    
        do i = 1,n_iter
            y_temp = y + (dx/2)*func(x, y)
            fxy = func(x, y_temp)
            y = y + dx*fxy
            x = x + dx
            write(12,*) x, ',', y
            ! PRINT*, x, ',', y

        end do

        dy = y_a - y
        PRINT*, "x = ", x, ",  dy = y_a - y_e = ", dy
        close(12)

    endif

    if (solver == 'improved_euler') then
        
        dx = 0.001
        n_iter = int (x_n/dx) + 1

        open(13, file = "output/improved_euler_soln.plt")
        write(13,*) x, ',', y
    
        do i = 1,n_iter
            f_0 = func(x,y)
            y_temp = y + dx*f_0

            f_1 = func(x, y_temp)

            fxy = (f_0 + f_1)/2
            y = y + dx*fxy
            x = x + dx
            write(13,*) x, ',', y
            ! PRINT*, x, ',', y

        end do

        dy = y_a - y
        PRINT*, "x = ", x, ",  dy = y_a - y_e = ", dy
        
        close(13)

    endif

    if (solver == 'rk4') then
    
        dx = 0.01
        n_iter = int (x_n/dx) + 1

        open(14, file = "output/rk4_soln.plt")
        write(14,*) x, ',', y
    
        do i = 1,n_iter
            y = y_update_rk4(dx, x, y)
            x = x + dx
            write(14,*) x, ',', y
            ! PRINT*, x, ',', y

        end do

        dy = y_a - y
        PRINT*, "x = ", x, ",  dy = y_a - y_e = ", dy
        
        close(14)

    endif


end program DE_solver

function func(x, y)
    implicit none 
    real*8:: x,y
    real*8:: func 

    func = 1 + y**2

end function func


function y_update_rk4(dx, x, y)
    implicit none
    real*8:: dx, x, y, y_update_rk4, func
    real*8:: f0, f1, f2, f3, y_temp1, y_temp2, y_temp3

    f0 = func(x, y)
    y_temp1 = y + (dx/2)*f0

    f1 = func(x + dx/2, y_temp1)
    y_temp2 = y + (dx/2)*f1

    f2 = func(x + dx/2, y_temp2)
    y_temp3 = y + dx*f2

    f3 = func(x + dx, y_temp3)

    y_update_rk4 = y + (dx/6)*(f0 + 2*f1 + 2*f2 + f3)

end function y_update_rk4