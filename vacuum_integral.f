!	This subroutine does the integration for the vacuum part
!	Give a value of z, it gives the contribution from gluon, u, d, s quarks
!	Prob1: for gluons
!	Prob2: for u and ubar quarks
!	Prob3: for d and dbar quarks 
!	Prob4: for s and sbar quarks 

	subroutine vacuum_int(z, int_1, int_2, int_3, int_4)

	implicit none

	include "common.inc"

!	These are the values to be shared with its subroutines                  
	double precision z_com
	common /z_com/ z_com                 

	double precision zcut_min, zcut_max
	common /z_cut_min_max/ zcut_min, zcut_max

	double precision Q, mu
           
	double precision z, int_1, int_2, int_3, int_4
                                              
	double precision Prob1, Prob2, Prob3, Prob4
        
	double precision dfdx_1, dfdx_2, dfdx_3, dfdx_4
	external dfdx_1, dfdx_2, dfdx_3, dfdx_4 
	double precision dfdx_v1, dfdx_vv1, dfdx_v2, dfdx_v3, dfdx_v4
	external dfdx_v1, dfdx_vv1, dfdx_v2, dfdx_v3, dfdx_v4
	double precision int_min, int_max

	DOUBLE PRECISION ABSERR,EPSABS,EPSREL,RESULT,WORK
	DIMENSION IWORK(100),WORK(400)
	INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
 
	EPSABS = 0.0D0
	EPSREL = 1.0D0
	KEY = 6
	LIMIT = 100
	LENW = LIMIT*4
         
	z_com = z

	Q = Q_com  

!	Any cut to be put for its subroutines dfdx(x)
!	zcut_min = 0.5D0-0.5D0*sqrt(1.D0-Q**2/Ejet_com**2)
!	zcut_max = 0.5D0+0.5D0*sqrt(1.D0-Q**2/Ejet_com**2)

	zcut_min = 0.d0
	zcut_max = 1.d0
       
	int_1 = 0.D0
	int_2 = 0.D0
	int_3 = 0.D0
	int_4 = 0.D0

!	The integration range is always from 0-1
!	Some hard range are implemented in the integral function dfdx(x)
	int_min = 0.D0
	int_max = 1.D0

!	dfdx_1: gluon to gluon/quark
!	dfdx_2: u/ubar quark to gluon/quark      
!	dfdx_3: d/dbar quark to gluon/quark      
!	dfdx_4: s/sbar quark to gluon/quark      

	call trapezoidal(dfdx_1,int_min,int_max,n_z,result)
!	call simpson38(dfdx_1,int_min,int_max,n_z,result)
!	CALL DQAG(dfdx_1,int_min,int_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
	int_1 = int_1 + result

	call trapezoidal(dfdx_2,int_min,int_max,n_z,result)                        
!	call simpson38(dfdx_2,int_min,int_max,n_z,result)
!	CALL DQAG(dfdx_2,int_min,int_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
	int_2 = int_2 + result
             
	call trapezoidal(dfdx_3,int_min,int_max,n_z,result)  
!	call simpson38(dfdx_3,int_min,int_max,n_z,result)
!	CALL DQAG(dfdx_3,int_min,int_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
	int_3 = int_3 + result
              
	call trapezoidal(dfdx_4,int_min,int_max,n_z,result)  
!	call simpson38(dfdx_4,int_min,int_max,n_z,result)     
!	CALL DQAG(dfdx_4,int_min,int_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
	int_4 = int_4 + result

	end

!!!!!!!!
             
!	dfdx_1: gluon to gluon/quark
!	So we need Prob1
	function dfdx_1(x)
        
	implicit none

	include "common.inc"

	double precision dfdx_1, x, z
                
	double precision z_com
	common /z_com/ z_com

	double precision zcut_min, zcut_max
	common /z_cut_min_max/ zcut_min, zcut_max
            
	double precision Prob1, Prob2, Prob3, Prob4 
                                 
	double precision q2q, q2g, g2q, g2g

	double precision P1, P2 
 
	z = z_com

	dfdx_1 = 0.D0
                         
	if (x .LT. zcut_min) return
	if (x .GT. zcut_max) return

	P1 = - Prob1(z)                     
	if (x .GT. z) P1 = P1 + Prob1(z/x)/x
	if (x .LT. 1.D0 - z) P1 = P1 + Prob1(z/(1.D0 - x))/(1.D0 - x)
	g2g = (1.D0/2.D0 * 6.D0 * (x/(1.D0-x) + x*(1.D0-x) + (1.D0-x)/x)) * P1

	P2 = - Prob1(z) * N_f                                            
	if (x .GT. z) P2 = P2 + (Prob2(z/x) + Prob3(z/x) + Prob4(z/x))/x  * 2.D0
!	if (x .GT. z) P2 = P2 + (Prob2(z/x) + Prob3(z/x) + Prob4(z/x))/x 
!	if (x .LT. 1.D0 - z) P2 = P2 + (Prob2(z/(1.D0-x)) + Prob3(z/(1.D0-x)) + Prob4(z/(1.D0-x)))/(1.D0-x) 
	g2q = 1.D0/2.D0 * (x**2+(1.D0-x)**2) * P2

	dfdx_1 = g2g + g2q

     	return
 
	end

!	dfdx_2: u/ubar quark to gluon/quark      
!	So we need Prob2
	function dfdx_2(x)

	implicit none

	include "common.inc"

	double precision dfdx_2, x, z
                
	double precision z_com
	common /z_com/ z_com

	double precision zcut_min, zcut_max
	common /z_cut_min_max/ zcut_min, zcut_max
                                    
	double precision Prob1, Prob2, Prob3, Prob4 
           
	double precision q2q, q2g, g2q, g2g

	double precision P1, P2
  
	z = z_com

	dfdx_2 = 0.D0

	if (x .LT. zcut_min) return
	if (x .GT. zcut_max) return

	P1 = - Prob2(z)
	if (x .GT. z) P1 = P1 + Prob2(z/x)/x
	if (x .LT. 1.D0 - z) P1 = P1 + Prob1(z/(1.D0-x))/(1.D0-x)
	q2q = 4.D0/3.D0 * (1.D0+x**2)/(1.D0-x) * P1

!	Since we are evolving qqhbar together, this is zero
!	P2 = Prob1(z/x)/x  
!	q2g = 4.D0/3.D0 * (1.D0+(1.D0-x)**2)/x * P2
	q2g = 0.D0

	dfdx_2 = q2q + q2g

     	return

	end

!	dfdx_3: d/dbar quark to gluon/quark      
!	So we need Prob3
	function dfdx_3(x)

	implicit none

	double precision dfdx_3, x, z
 
	double precision z_com
	common /z_com/ z_com

	double precision zcut_min, zcut_max
	common /z_cut_min_max/ zcut_min, zcut_max

	double precision Prob1, Prob2, Prob3, Prob4 
       
	double precision q2q, q2g, g2q, g2g

	double precision P1, P2
              
	z = z_com

	dfdx_3 = 0.D0
        
	if (x .LT. zcut_min) return
	if (x .GT. zcut_max) return

	P1 = - Prob3(z)
	if (x .GT. z) P1 = P1 + Prob3(z/x)/x 
	if (x .LT. 1.D0 - z) P1 = P1 + Prob1(z/(1.D0-x))/(1.D0-x) 
	q2q = 4.D0/3.D0 * (1.D0+x**2)/(1.D0-x) * P1

!	Since we are evolving qqhbar together, this is zero
!	P2 = Prob1(z/x)/x 
!	q2g = 4.D0/3.D0 * (1.D0+(1.D0-x)**2)/x * P2
	q2g = 0.D0

	dfdx_3 = q2q + q2g

     	return

	end

!	dfdx_4: s/sbar quark to gluon/quark      
!	So we need Prob4
	function dfdx_4(x)
        
	implicit none

	include "common.inc"

	double precision dfdx_4, x, z
                
	double precision z_com
	common /z_com/ z_com

	double precision zcut_min, zcut_max
	common /z_cut_min_max/ zcut_min, zcut_max

	double precision Prob1, Prob2, Prob3, Prob4
                                 
	double precision q2q, q2g, g2q, g2g

	double precision P1, P2
 
	z = z_com

	dfdx_4 = 0.D0

	if (x .LT. zcut_min) return
	if (x .GT. zcut_max) return

	P1 = - Prob4(z)
	if (x .GT. z) P1 = P1 + Prob4(z/x)/x 
	if (x .LT. 1.D0 - z) P1 = P1 + Prob1(z/(1.D0-x))/(1.D0-x) 
	q2q = 4.D0/3.D0 * (1.D0+x**2)/(1.D0-x) * P1

!	Since we are evolving qqhbar together, this is zero
!	P2 = Prob1(z/x)/x 
!	q2g = 4.D0/3.D0 * (1.D0+(1.D0-x)**2)/x * P2
	q2g = 0.D0

	dfdx_4 = q2q + q2g

     	return

	end
 
