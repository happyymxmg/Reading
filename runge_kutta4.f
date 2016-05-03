!	This subroutine does the evolution 
!	Input is f_ini at x=x_ini, and output is f_ini at x=x_ini+h, where h is given in the main program
!	It calls the subroutine derivs(x, f, k) to calculate the dirivative
!	So all the physics information is contained in the subroutine derivs
	subroutine rk4(x_ini, f_ini, h, f_fin)

	implicit none
    	
	include "common.inc"
               
	double precision x_ini, f_ini(n_dim, n_z + 1), h, f_fin(n_dim, n_z + 1)
                                                   
 	double precision x, f(n_dim, n_z + 1)
	double precision k_1(n_dim, n_z + 1), k_2(n_dim, n_z + 1), k_3(n_dim, n_z + 1), k_4(n_dim, n_z + 1) 

	integer i_z, i_dim
        
	x = x_ini
	do i_z = 2, n_z
		do i_dim = 1, n_dim
			f(i_dim, i_z) = f_ini(i_dim, i_z)
		end do
	end do
	call derivs(x, f, k_1)
       
	x = x_ini + h/2.D0
 	do i_z = 2, n_z
		do i_dim = 1, n_dim 
			f(i_dim, i_z) = f_ini(i_dim, i_z) + h/2.D0 * k_1(i_dim, i_z)
		end do
	end do
	call derivs(x, f, k_2)
	
	x = x_ini + h/2.D0
 	do i_z = 2, n_z
 		do i_dim = 1, n_dim 
			f(i_dim, i_z) = f_ini(i_dim, i_z) + h/2.D0 * k_2(i_dim, i_z)
		end do
	end do
	call derivs(x, f, k_3)        
      	
	x = x_ini + h 
 	do i_z = 2, n_z
 		do i_dim = 1, n_dim  
			f(i_dim, i_z) = f_ini(i_dim, i_z) + h * k_3(i_dim, i_z)
		end do
	end do
	call derivs(x, f, k_4)
	
	do i_z = 2, n_z
 		do i_dim = 1, n_dim

!			2nd order
!  			f_fin(i_dim, i_z) = f_ini(i_dim, i_z) + h * k_2(i_dim, i_z)

!			4th order
			f_fin(i_dim, i_z) = f_ini(i_dim, i_z) + h/6.D0 * (k_1(i_dim, i_z) + 2.D0 * k_2(i_dim, i_z) + 2.D0 * k_3(i_dim, i_z) + k_4(i_dim, i_z))
		end do
	end do

	return

	end
