!	This subroutine gives the derivative dfdx at a given x, where the function/array f(x) is given as input	
!	The physics is contained here and its subroutines
!	The derivative has two contributions: vacuum and medium
!	Both contribution are calculated by calling integration subroutines: int_vacuum, int_medium
!	The common factor alpha_s/(2PI) has been extracted here, and does not appear integration subroutine anymore
	subroutine derivs(x, f, dfdx)

	implicit none

	include "common.inc"

!	These are the arrays to store the values of FF's

        integer nmax
        parameter (nmax = 1000)

	double precision dt_Prob1 
	common /dt_Prob1/ dt_Prob1(nmax+1)

	double precision dt_Prob2
	common /dt_Prob2/ dt_Prob2(nmax+1)

	double precision dt_Prob3 
	common /dt_Prob3/ dt_Prob3(nmax+1)

	double precision dt_Prob4 
	common /dt_Prob4/ dt_Prob4(nmax+1)
                                                    
	double precision x, f(n_dim, n_z + 1), dfdx(n_dim, n_z + 1)
	integer i_z
	
	double precision Q
               
	double precision z, int_1, int_2, int_3, int_4
	double precision int_m1, int_m2, int_m3, int_m4

        double precision L_medium
        double precision Q_EL
 
        if (nmax .lt. n_z) then
                write (*,*) "nmax < n_z!"
                stop
        end if

!	Prob1: for gluons
!	Prob2: for u and ubar quarks
!	Prob3: for d and dbar quarks 
!	Prob4: for s and sbar quarks 
	do i_z = 2, n_z
  		
		dt_Prob1(i_z) = f(1, i_z)
		dt_Prob2(i_z) = f(2, i_z)  
		dt_Prob3(i_z) = f(3, i_z)    
		dt_Prob4(i_z) = f(4, i_z)    

	end do
 	
!	Get the values for ending points (i_z = 0, i_z = n_z + 1)
                          
	dt_Prob1(1) = 2.D0*dt_Prob1(2) - dt_Prob1(3)
	dt_Prob2(1) = 2.D0*dt_Prob2(2) - dt_Prob2(3) 
	dt_Prob3(1) = 2.D0*dt_Prob3(2) - dt_Prob3(3)          
	dt_Prob4(1) = 2.D0*dt_Prob4(2) - dt_Prob4(3)
                                 
!	dt_Prob1(1) = 3.D0*dt_Prob1(2) - 3.D0*dt_Prob1(3) + dt_Prob1(4)           
!	dt_Prob2(1) = 3.D0*dt_Prob2(2) - 3.D0*dt_Prob2(3) + dt_Prob2(4)            
!	dt_Prob3(1) = 3.D0*dt_Prob3(2) - 3.D0*dt_Prob3(3) + dt_Prob3(4)                     
!	dt_Prob4(1) = 3.D0*dt_Prob4(2) - 3.D0*dt_Prob4(3) + dt_Prob4(4)           
 	                       
	dt_Prob1(n_z+1) = 0.D0
	dt_Prob2(n_z+1) = 0.D0  
	dt_Prob3(n_z+1) = 0.D0     
	dt_Prob4(n_z+1) = 0.D0
                                  
	Q = sqrt(exp(x))
	Q_com = Q
	alpha_s =  2.D0 * PI / ((11.D0 - 2.D0/3.D0*N_f) * log(Q/Lambda_QCD))
 
	do i_z = 2, n_z
		
		z = z_min + delta_z * (i_z - 1)
             
!		Vacuum  
		int_1 = 0.D0
		int_2 = 0.D0
		int_3 = 0.D0
		int_4 = 0.D0

!		Medium
		int_m1 = 0.D0
		int_m2 = 0.D0
		int_m3 = 0.D0
		int_m4 = 0.D0
		 
		if (flag_vacuum .EQ. 1) call vacuum_int(z, int_1, int_2, int_3, int_4)
                L_medium = length_jim
                Q_EL = sqrt(Ejet_com/L_medium*hbar_c)
		if (Q_com .gt. Q_EL) then
		        if (flag_medium .EQ. 1) call medium_int(z, int_m1, int_m2, int_m3, int_m4)
                end if
                    
		dfdx(1, i_z) = alpha_s/(2.D0*PI) * (int_1 + int_m1)
		dfdx(2, i_z) = alpha_s/(2.D0*PI) * (int_2 + int_m2)
		dfdx(3, i_z) = alpha_s/(2.D0*PI) * (int_3 + int_m3)
		dfdx(4, i_z) = alpha_s/(2.D0*PI) * (int_4 + int_m4)
                     
!		write (*,*) "compare:", z, alpha_s/(2.D0*PI) * int_2, alpha_s/(2.D0*PI) * int_m2, int_m2/int_2
	end do 

	return

	end

