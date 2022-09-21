program laplace
    implicit none

    integer, parameter :: lx =34, ly = 34 ! Analog of DEF in cpp. integer, parameter defines lx and ly as aliases.
    real*8 :: old_temp(1:lx, 1:ly), temp(1:lx, 1:ly)
    real*8 :: bound_temp, increment_temp, dx, dy, prefactor, tolerance
    integer :: i, j, ii, jj
    
    character(len=30) :: filename
    integer :: iunit, ci, test, counter      
    ! counter to count no. of iterations, test to check convergence test, iunit for writing pipeline number.

    ! PARAMATERS
    dx = 0.01d0; dy = 0.01d0
    tolerance = 0.0001d0
    bound_temp = 0.0d0         ! boundary temperature is the offset in the temperaure at (1,1) boundary
    increment_temp = 0.2d0     ! incremental steps of increasing T along the boundary

    
    old_temp = 0.0d0 ! can we set the value of the entire array to 0 just like this?

    ! Setting up the Boundary Conditions.

    do j = 1, ly        
        old_temp(1,j) = 7.4
        old_temp(lx,j) = 0.8
        ! old_temp(1,j) = bound_temp + dfloat(j)*increment_temp
        ! old_temp(lx,j) = dfloat(ly)*increment_temp - dfloat(j)*increment_temp 
    end do

    do i = 2, lx-1
        old_temp(i, 1) = old_temp(1, 1) - dfloat(i)*increment_temp
        old_temp(i,ly) = old_temp(1,ly) - dfloat(i)*increment_temp
    end do

    temp = old_temp

    iunit = 71; ci = 0
    write(filename, '("initialize_",i0,".dat")' ) ci
    
    open(newunit = iunit, file = filename)   ! Writing initial lattice config at 0th iteration

    do ii = 1,lx
        do jj = 1,ly
            write(iunit,*) ii, jj, old_temp(ii, jj)
        end do
    end do
    
    close(iunit)

    test = 0; counter = 0; prefactor = (0.50d0*dx*dx*dy*dy)/(dx*dx + dy*dy)

    do   ! THE MAIN DO LOOP
        counter = counter + 1
        test = 0

        ! Update loop
        do jj = 2, ly - 1
            do ii = 2, lx - 1
                temp(ii, jj) = prefactor*( (old_temp(ii + 1, jj) + old_temp(ii - 1, jj))/(dx*dx) & 
                + (old_temp(ii, jj + 1) + old_temp(ii, jj - 1))/(dy*dy) ) 
            end do
        end do


        ! Convergence Check loop
        do jj = 2, ly - 1
            do ii = 2, lx - 1
                if ( (abs(temp(ii, jj) - old_temp(ii, jj))).gt.tolerance) test = 1
            end do
        end do

        
        if (test.eq.0) exit   
        old_temp = temp

    end do

    PRINT*, "number of iterations to achieve convergence = ", counter 


    ! Writing down the converged result

    iunit = 71; ci = 10000

    write( filename, '("initialize_",i0,".dat")' ) ci

    open(newunit = iunit, file = filename)

    do ii = 1, lx
        do jj = 1, ly
            write(iunit, *) ii, jj, temp(ii,jj)
        end do
    end do

    close(iunit)

end program laplace