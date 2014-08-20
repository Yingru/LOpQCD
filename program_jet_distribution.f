	program jet_intial_distribution

	implicit none
	
	DOUBLE PRECISION	ABSERR,EPSABS,EPSREL,RESULT,WORK
	DIMENSION IWORK(100),WORK(400)
	INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL

	double precision funjet, xa_min, xa_max
	external funjet

	INTEGER init,itmx,ncall,ndim,nprn
	double precision tgral,chi2a,sd,regionpp(4)

	double precision PI
	common /constants/ PI
	
	double precision alpha_em, alpha_s
	common /parameters/ alpha_em, alpha_s
		
	double precision atom_A, atom_B, s_NN, n_coll, sigma_in 
	common /experiments/ atom_A, atom_B, s_NN, n_coll, sigma_in 
	
	double precision proton_A, proton_B
	common /proton/ proton_A, proton_B
	
	double precision k_jet, k_gamma, k_brem
	common /NLO/ k_jet, k_gamma, k_brem 
	
	logical shadowing, dA, eloss, owens
	common /shadowing/ shadowing, dA,  eloss, owens

	integer flag
	common /flag/ flag

	double precision E_com, Rperd, Rperd2, p_T, x_T, y
	common /variables/ E_com, Rperd, Rperd2, p_T, x_T, y
	
	double precision z_c, pT_c, xT_c, Q, factor, const, temp

     	double precision jet_u, jet_d, jet_s, jet_charm, jet_g 	
	double precision jet_ubar, jet_dbar, jet_sbar, jet_charmbar
	double precision jet_avr_q, jet_avr_qbar, jet_avr_qqbar, jet_sum_qqbar, jet_sum_ccbar, sum_parton

	double precision ratio_u, ratio_d, ratio_s, ratio_ubar, ratio_dbar, ratio_sbar

	integer jet
	common /jet/ jet
	
	double precision zc_com
	common /zc_com/ zc_com
	 
	double precision M_charm
	common /M_charm/ M_charm

	double precision rate_nn, rate_da, rate_aa

	double precision parton_pt
	double precision M_T

        double precision dsigmaOVERdy

        dsigmaOVERdy = 0d0

	z_c = 1.d0
	zc_com = z_c
	
	PI = 3.141592653589793d0
	alpha_em = 1.d0/137.d0

!	M_charm = 1.27d0
	M_charm = 4.19d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Experiments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Au+Au @ 200GeV
	atom_A = 197
	atom_B = 197
	proton_A = 79
	proton_B = 79
	s_NN = 62.4**2
	k_jet = 1.7
 
!	Pb+Pb @ 5500GeV
!	atom_A = 208
!	atom_B = 208
!	proton_A = 82
!	proton_B = 82
!	s_NN = 2760**2
!	k_jet = 1.7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	E_com = sqrt(s_NN) / 2.d0
	Rperd = 1.18d0 * atom_A**(1.d0/3.d0)	
!	write (*,*) "please give the rapidity y:"
!	read (*,*) y
	y = 0.d0

	parton_pt = 0.5D0
	do while (parton_pt .le. 70D0)	! the parton energy

		p_T = parton_pt		
		M_T = sqrt(p_T**2 + M_charm**2)   
		x_T = M_T / E_com
		
 		pT_c = p_T / z_c	! the momentum of the c parton
 		xT_c = pT_c / E_com

! 		shadowing = .FALSE.
!		dA = .FALSE.
	 	 
!		shadowing = .TRUE.
!		dA = .FALSE.
	 	
		shadowing = .TRUE.
		dA = .FALSE.
		eloss = .FALSE.
	 	      
 		EPSABS = 0.0E0
 		EPSREL = 1.0E0
 		KEY = 6
 		LIMIT = 100
 		LENW = LIMIT*4
 	
 		xa_min = x_T * exp(y) / (2.d0 - x_T * exp(-y))
 		xa_max = 1.d0
		
		jet = 0
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_g = result * 0.388d0	 

 		jet = 1
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_u = result * 0.388d0	  

		jet = -1
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_ubar = result * 0.388d0	  
		
		jet = 2
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_d = result * 0.388d0	  
		
		jet = -2
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_dbar = result * 0.388d0	     
   		
		jet = 3
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_s = result * 0.388d0	        
  		
		jet = -3
 		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     c		NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
 		jet_sbar = result * 0.388d0	     
	      		
		jet = 4
		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
		jet_charm= result * 0.388d0 
		
		jet = -4
		CALL DQAG(funjet,xa_min,xa_max,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
		jet_charmbar = result * 0.388d0 

		open (20,file='distribution_g.dat',status="unknown")
		open (21,file='distribution_u.dat',status="unknown")	
		open (22,file='distribution_ubar.dat',status="unknown")	
		open (23,file='distribution_d.dat',status="unknown")	
		open (24,file='distribution_dbar.dat',status="unknown")	
		open (25,file='distribution_s.dat',status="unknown")	
		open (26,file='distribution_sbar.dat',status="unknown")	
		open (27,file='distribution_c.dat',status="unknown")	
		open (28,file='distribution_cbar.dat',status="unknown")	
		
		rate_aa = jet_g 
!		write (20,*) p_T, rate_aa 
	
		rate_aa = jet_u 
!		write (21,*) p_T, rate_aa	
	
		rate_aa = jet_ubar
!		write (22,*) p_T, rate_aa	
	
	     	rate_aa = jet_d 
!		write (23,*) p_T, rate_aa	
	
		rate_aa = jet_dbar 
!		write (24,*) p_T, rate_aa	
	
		rate_aa = jet_s 
!		write (25,*) p_T, rate_aa 
	
	      	rate_aa = jet_sbar 
!		write (26,*) p_T, rate_aa	
	
		rate_aa = jet_charm
!		write (27,*) p_T, rate_aa 
	
	      	rate_aa = jet_charmbar 
!		write (28,*) p_T, rate_aa	

		open (40,file='sum_g.dat',status="unknown")	
		open (50,file='sum_qqbar.dat',status="unknown")	
		open (60,file='sum_ccbar.dat',status="unknown")	
	
		rate_aa = jet_g     
!		write (40,*) p_T, rate_aa 
	
		jet_sum_qqbar = jet_u + jet_d + jet_s + jet_ubar + jet_dbar + jet_sbar
		rate_aa = jet_sum_qqbar    ! yield, dN/d^2p_Tdy = GeV^{-2}
!		write (*,*) p_T, rate_aa
!		write (50,*) p_T, rate_aa

		jet_sum_ccbar = jet_charm + jet_charmbar
		rate_aa = jet_sum_ccbar    ! yield, dN/d^2p_Tdy = GeV^{-2}
		write (*,*) p_T, rate_aa, jet_charm, jet_charmbar
		write (60,*) p_T, rate_aa

                dsigmaOVERdy = dsigmaOVERdy + rate_aa*2d0*PI*p_T*0.5d0

		parton_pt = parton_pt + 0.5D0

	end do

        write(6,*) "dsigma/dy: ", dsigmaOVERdy

	end program

! 	(hc)^2 = (0.197 GeV fm)^2 = 0.197^2 GeV^2 fm^2 = 0.0388d0 GeV^2 fm^2 = 0.388d0 GeV^2 mb = 0.388 * 10^9  GeV^2 pb

