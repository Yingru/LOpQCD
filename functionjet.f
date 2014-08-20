	function funjet(x_a)

	implicit none

	double precision funjet, x_a
	
	DOUBLE PRECISION	ABSERR,EPSABS,EPSREL,RESULT,WORK
      DIMENSION IWORK(100),WORK(400)
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL

	double precision PI
	common /constants/ PI
	
	double precision alpha_em, alpha_s
	common /parameters/ alpha_em, alpha_s
				
	double precision atom_A, atom_B, s_NN, n_coll, sigma_in 
	common /experiments/ atom_A, atom_B, s_NN, n_coll, sigma_in 
	
	double precision proton_A, proton_B
	common /proton/ proton_A, proton_B

	double precision T_i, tau_i, T_c, tau_f, r_d, tau_H
	common /intial/ T_i, tau_i, T_c, tau_f, r_d, tau_H
	
	double precision k_jet, k_gamma, k_brem
	common /NLO/ k_jet, k_gamma, k_brem 
	
	logical shadowing, dA, eloss, owens
	common /shadowing/ shadowing, dA,  eloss, owens

	integer flag
	common /flag/ flag

	double precision E_com, Rperd, Rperd2, p_T, x_T, y
	common /variables/ E_com, Rperd, Rperd2, p_T, x_T, y
	
	double precision x_b, z_c, pT_c, xT_c, Q, factor, const, tmp
	double precision s, t, u
	double precision ratio

	integer Iset
	double precision Ctq5Pdf, eks98r

	double precision ruv, rdv, ru, rd, rs, rc, rb, rt, rg
        double precision QQQ, xxx
        integer order, pset, AAA

	double precision fpa_u, fpa_ubar, fpa_d, fpa_dbar, fpa_s, fpa_sbar, fpa_g, fpa_c, fpa_cbar
	double precision fpb_u, fpb_ubar, fpb_d, fpb_dbar, fpb_s, fpb_sbar, fpb_g, fpb_c, fpb_cbar
	
	double precision fa_u, fa_ubar, fa_d, fa_dbar, fa_s, fa_sbar, fa_g, fa_c, fa_cbar
	double precision fb_u, fb_ubar, fb_d, fb_dbar, fb_s, fb_sbar, fb_g, fb_c, fb_cbar

	double precision sigma1, sigma2, sigma3, sigma4
	double precision sigma5, sigma6, sigma7, sigma8 

	double precision sigma_u1, sigma_u2, sigma_u3, sigma_u4
	double precision sigma_u5, sigma_u6, sigma_u7, sigma_u8

	double precision sum_u, sum_d, sum_s, sum_g
	double precision sum_ubar, sum_dbar, sum_sbar, sum_parton

	integer jet
	common /jet/ jet

	double precision zc_com
	common /zc_com/ zc_com
        
	double precision M_charm
	common /M_charm/ M_charm
	 
	double precision sigma_qq2cc, sigma_gg2cc
	double precision sigma_qc2qc, sigma_gc2gc
	double precision sigma_cq2cq, sigma_cg2cg
	double precision sigma_qc2cq, sigma_gc2cg
	double precision sigma_cq2qc, sigma_cg2gc
	double precision MM
        double precision M_ss, M_tt, M_uu, M_st, M_su, M_tu
	double precision sum_charm, sum_charm_bar
	double precision MT_c, MT_h
	
        integer id_hv

!       for charm
        if (abs(M_charm - 1.3d0) .lt. 0.5d0) then
                id_hv = 4
!       for bottom
        else if (abs(M_charm - 4.2d0) .lt. 0.5d0) then
                id_hv = 5
        else
                write (*, *) "Peterson: not charm nor bottom!"
                stop
        end if

	funjet = 0

	z_c = zc_com 

!        MT_h = sqrt(pT_h**2 + M_meson**2)

!        pT_c = (MT_h+pT_h)/(2.D0*z_c)-(z_c*M_charm**2)/(2.D0*(MT_h+pT_h))
	pT_c = p_T / z_c

	MT_c = sqrt(pT_c**2 + M_charm**2) 
	xT_c = MT_c / E_com
	
	x_b = x_a * xT_c * exp(-y) / ( 2.0 * x_a - xT_c * exp(y) )

	if (x_b .LT. 0.00001) return
	if (x_b .GT. 0.99999) return
	if (x_a .LT. 0.00001) return
	if (x_a .GT. 0.99999) return

!	Q = MT_h
	Q = MT_c
	If (Q .LT. 1.d0) Q = 1.d0
	alpha_s = 12.d0 * PI / ((11.d0*3.d0 - 2.d0*3.d0) * log(Q**2 / 0.2d0**2)) 
                                
	s = s_NN * x_a * x_b
	t = -sqrt(s_NN) * x_a * MT_c * exp(-y)
	u = -sqrt(s_NN) * x_b * MT_c * exp(y)

!	The function Ctq5Pdf (Iparton, X, Q)
!	returns the parton distribution inside the proton for parton [Iparton] 
!	at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
!	Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
!					 for (b, c, s, d, u, g, u_bar, ..., b_bar),
!	Details are:
!	---------------------------------------------------------------------------
!  	Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
!	---------------------------------------------------------------------------
!	1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
!	2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
!	3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
!	4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
!	5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
!	6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
!	7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
!		--------------------------------------------------------
!	8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
!	9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
!	---------------------------------------------------------------------------
! 
!	The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
!	Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
!	The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
!	which is defined as the bottom quark mass, whenever it can be applied.
!
!	The Table_Files are assumed to be in the working directory.
!
!	Before using the PDF, it is necessary to do the initialization by
!		Call SetCtq5(Iset) 
!	where Iset is the desired PDF specified in the above table.
!
!	These programs, as provided, are in double precision.  By removing the
!	"Implicit Double Precision" lines, they can also be run in single precision.

	Iset=3
	Call SetCtq5(Iset)

	fpa_g 	= Ctq5Pdf(0,x_a,Q)
	fpa_u 	=  Ctq5Pdf(1,x_a,Q)
	fpa_d 	=  Ctq5Pdf(2,x_a,Q)
	fpa_s 	=  Ctq5Pdf(3,x_a,Q)
	fpa_c    = Ctq5Pdf(id_hv,x_a,Q)
	fpa_ubar 	=  Ctq5Pdf(-1,x_a,Q) 
	fpa_dbar 	=  Ctq5Pdf(-2,x_a,Q) 
	fpa_sbar 	=  Ctq5Pdf(-3,x_a,Q)  
	fpa_cbar = Ctq5Pdf(-id_hv,x_a,Q)  
		
	fpb_g 	= Ctq5Pdf(0,x_b,Q)
	fpb_u 	=  Ctq5Pdf(1,x_b,Q)
	fpb_d 	=  Ctq5Pdf(2,x_b,Q)
	fpb_s 	=  Ctq5Pdf(3,x_b,Q)
	fpb_c    = Ctq5Pdf(id_hv,x_b,Q)
	fpb_ubar 	=  Ctq5Pdf(-1,x_b,Q)
	fpb_dbar 	=  Ctq5Pdf(-2,x_b,Q) 
	fpb_sbar 	=  Ctq5Pdf(-3,x_b,Q)  
	fpb_cbar = Ctq5Pdf(-id_hv,x_b,Q)  

	if(shadowing .EQV. .FALSE.) then

		ratio = proton_A / Atom_A

		fa_g 		= fpa_g 
		fa_u 		= ratio * fpa_u + (1.d0-ratio) * fpa_d
		fa_d 		= ratio * fpa_d + (1.d0-ratio) * fpa_u
		fa_s 		= fpa_s 
		fa_c 		= fpa_c
		fa_ubar 	= ratio * fpa_ubar + (1.d0-ratio) * fpa_dbar
		fa_dbar 	= ratio * fpa_dbar + (1.d0-ratio) * fpa_ubar
		fa_sbar 	= fpa_sbar  
		fa_cbar 	= fpa_cbar  
	
		ratio = proton_B / Atom_B

		fb_g 		= fpb_g 
		fb_u 		= ratio * fpb_u  + (1.d0-ratio) * fpb_d
		fb_d 		= ratio * fpb_d  + (1.d0-ratio) * fpb_u
		fb_s 		= fpb_s 
		fb_c 		= fpb_c
		fb_ubar 	= ratio * fpb_ubar + (1.d0-ratio) * fpb_dbar
		fb_dbar 	= ratio * fpb_dbar + (1.d0-ratio) * fpb_ubar
		fb_sbar 	= fpb_sbar  
		fb_cbar 	= fpb_cbar  
	
    	end if

C ***************************************************************************
C Instructions:
C
C For given input values of
C
C     order: 1=LO, 2=NLO   ; integer
C     pset : 1...31        ; integer
C            1     = central fit
C            2,3   = error sets S{+1}, S{-1}
C            4,5   = error sets S{+2}, S{-2}
C            ...   ...
C            30,31 = error sets {S+15}, {S-15}
C     A    : atomic number ; integer
C     x    : Bjorken-x     ; double precision
C     Q    : scale in GeV  ; double precision
C
C the command 
C
C   Call EPS09(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, rb, rg)
C
C returns the bound proton nuclear corrections R_f^A(x,Q)
C (in double precision) for
C	
C	ruv = up valence
C	rdv = down valence
C	ru  = up sea
C	rd  = down sea
C	rs  = strange
C	rc  = charm
C	rb  = bottom
C	rg  = gluons
C
C The nuclear corrections for bound neutrons can be obtained
C by the isospin symmetry, e.g. the total up quark distribution
C per nucleon in a nucleus A with Z protons is
C
C  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
C            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]
C
C Note that the parametrization should only be applied at the
C kinematical domain
C
C             1e-6 <= x <= 1
C              1.3 <= Q <= 1000 GeV.
C
C No warning message is displayed if these limits are
C exceeded, and outside these boundaries the modifications
C are frozen to the boundary values, i.e
C
C   for Q > 1000, the modifications at Q=1000 are returned,
C   for Q < 1.3,  the modifications at Q=1.3 are returned,
C   for x > 1,    the modifications at x=1 are returned
C   for x < 1e-6, the modifications at x=1e-6 are returned,
C
C The data used by the program for required order
C and atomic number A, are stored in separate files
C
C   LO : EPS09LOR_A
C   NLO: EPS09NLOR_A
C
C which must be located in the working directory.

	if(shadowing .EQV. .TRUE.) then
		
!		Call eps08(x_a,Q,atom_A,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)

                order = 1
                pset = 1
                AAA = nint(atom_A)
                xxx = x_a
                QQQ = Q
                call EPS09(order, pset, AAA, xxx, QQQ, ruv, rdv, ru, rd, rs, rc, rb, rg)

		ratio = proton_A / Atom_A

		fa_g 		= rg*fpa_g
		fa_u 		= ratio * (ruv*(fpa_u-fpa_ubar)+ru*fpa_ubar) + (1-ratio) * (rdv*(fpa_d-fpa_dbar)+rd*fpa_dbar)
		fa_d 		= ratio * (rdv*(fpa_d-fpa_dbar)+rd*fpa_dbar) + (1-ratio) * (ruv*(fpa_u-fpa_ubar)+ru*fpa_ubar)
		fa_s 		= rs * fpa_s
		if (id_hv .eq. 4) fa_c 		= rc * fpa_c
		if (id_hv .eq. 5) fa_c 		= rb * fpa_c
		fa_ubar 	= ratio * ru * fpa_ubar + (1.d0-ratio) * rd * fpa_dbar
		fa_dbar 	= ratio * rd * fpa_dbar + (1.d0-ratio) * ru * fpa_ubar
		fa_sbar 	= rs * fpa_sbar
		if (id_hv .eq. 4) fa_cbar 	= rc * fpa_cbar
		if (id_hv .eq. 5) fa_cbar 	= rb * fpa_cbar

!		write (*,*) fa_u, fa_d, fa_ubar, fa_dbar, fa_s, fa_sbar

!		Call eps08(x_b,Q,atom_B,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)

                order = 1
                pset = 1
                AAA = nint(atom_B)
                xxx = x_b
                QQQ = Q
                call EPS09(order, pset, AAA, xxx, QQQ, ruv, rdv, ru, rd, rs, rc, rb, rg)

		ratio = proton_B / Atom_B
	
		fb_g 		= rg * fpb_g
		fb_u 		= ratio * (ruv*(fpb_u-fpb_ubar)+ru*fpb_ubar) + (1-ratio) * (rdv*(fpb_d-fpb_dbar)+rd*fpb_dbar)
		fb_d 		= ratio * (rdv*(fpb_d-fpb_dbar)+rd*fpb_dbar) + (1-ratio) * (ruv*(fpb_u-fpb_ubar)+ru*fpb_ubar)
		fb_s 		= rs * fpb_s
		if (id_hv .eq. 4) fb_c 		= rc * fpb_c
		if (id_hv .eq. 5) fb_c 		= rb * fpb_c
		fb_ubar 	= ratio * ru * fpb_ubar + (1.d0-ratio) * rd * fpb_dbar
		fb_dbar 	= ratio * rd * fpb_dbar + (1.d0-ratio) * ru * fpb_ubar
		fb_sbar 	= rs * fpb_sbar
		if (id_hv .eq. 4) fb_cbar 	= rc * fpb_cbar
		if (id_hv .eq. 5) fb_cbar 	= rb * fpb_cbar
	
	end if

!	parton scattering cross section

	factor = PI * alpha_s**2 / s**2
 
	sigma1 = factor * ( 4.0/9.0 * (s**2+u**2)/t**2 )					! qq'->qq'			q'q->q'q
	sigma2 = factor * ( 4.0/9.0 * ((s**2+u**2)/t**2+(s**2 + t**2)/u**2) 
     c	- 8.0/27.0 * s**2/(t*u) )								! qq->qq			
     	sigma3 = factor * ( 4.0/9.0 * (t**2+u**2)/s**2 )					! qq_bar->q'q'_bar	q_bar q->q'_barq'
	sigma4 = factor * ( 4.0/9.0 * ((s**2+u**2)/t**2+(u**2+t**2)/s**2)
     c	- 8.0/27.0 * u**2/(s*u) )								! qq_bar->qq_bar		q_bar q->q_bar q
     	sigma5 = factor * ( -4.0/9.0 * (s/u+u/s)+(s**2 + u**2)/t**2 )			! gq->gq			qg->qg
     	sigma6 = factor * ( 32.0/27.0 * (t/u+u/t) - 8.0/3.0 * (t**2+u**2)/s**2 )	! qq_bar -> gg		q_bar q->gg
	sigma7 = factor * ( 1.0/6.0 * (t/u+u/t) - 3.0/8.0 * (t**2+u**2)/s**2 )		! gg->qq_bar		gg->q_bar q
	sigma8 = factor * ( 9.0/2.0 * (3-t*u/s**2-s*u/t**2-s*t/u**2) )			! gg->gg

!	now exchange 3&4 or u&t, because both 3&4 can be c particle.
!	we can also say we exchange 3&4 because we fix c particle to be a jet.

	tmp = u
	u = t
	t = tmp
		
	sigma_u1 = factor * ( 4.0/9.0 * (s**2+u**2)/t**2 )					! qq'->q'q, 		q'q->qq'
	sigma_u2 = factor * ( 4.0/9.0 * ((s**2+u**2)/t**2+(s**2 + t**2)/u**2) 
     c	- 8.0/27.0 * s**2/(t*u) )								! qq->qq
     	sigma_u3 = factor * ( 4.0/9.0 * (t**2+u**2)/s**2 )					! qq_bar->q'_bar q', 	q_bar q->q'q'_bar
	sigma_u4 = factor * ( 4.0/9.0 * ((s**2+u**2)/t**2+(u**2+t**2)/s**2)
     c	- 8.0/27.0 * u**2/(s*u) )								! qq_bar->q_bar q, 	q_bar q->qq_bar
     	sigma_u5 = factor * ( -4.0/9.0 * (s/u+u/s)+(s**2 + u**2)/t**2 )			! gq->qg, 			qg->qg
     	sigma_u6 = factor * ( 32.0/27.0 * (t/u+u/t) - 8.0/3.0 * (t**2+u**2)/s**2 )	! qq_bar -> gg, 		q_bar q->gg
	sigma_u7 = factor * ( 1.0/6.0 * (t/u+u/t) - 3.0/8.0 * (t**2+u**2)/s**2 )	! gg->q_bar q, 		gg->qq_bar
	sigma_u8 = factor * ( 9.0/2.0 * (3-t*u/s**2-s*u/t**2-s*t/u**2) )			! gg->gg

!       change back to the original 
        tmp = u
        u = t
        t = tmp

!	For different parton species

!	c=u
	if (jet .EQ. 1) then
		sum_u = sigma1 * fa_u * (fb_d + fb_s + fb_dbar + fb_sbar)		!ud->ud, us->us, ud_bar->ud_bar, us_bar->us_bar
     c		+ sigma2 * fa_u * fb_u 							!uu->uu
     c		+ sigma3 * (fa_d * fb_dbar + fa_s * fb_sbar)			!dd_bar->uu_bar, ss_bar->uu_bar
     c		+ sigma4 * fa_u * fb_ubar						!uu_bar->uu_bar
     c		+ sigma5 * fa_u * fb_g							!ug->ug
     c		+ sigma7 * fa_g * fb_g							!gg->uu_bar
     c		+ sigma_u1 * (fa_d + fa_s + fa_dbar + fa_sbar) * fb_u		!du->ud, su->us, d_bar u->ud_bar, s_bar u->us_bar
     c		+ sigma_u3 * (fa_dbar * fb_d + fa_sbar *fb_s)			!d_bar d->uu_bar, s_bar s->uu_bar
     c		+ sigma_u4 * fa_ubar * fb_u 						!u_bar u->uu_bar
     c		+ sigma_u5 * fa_g * fb_u						!gu->ug
	end if

!	c=u_bar
	if (jet .EQ.-1) then
		sum_ubar = sigma1 * fa_ubar * (fb_d + fb_s 
     c		+ fb_dbar + fb_sbar)		!u_bar d->u_bar d , u_bar s->u_bar s, u_bar d_bar->u_bar d_bar, u_bar s_bar->u_bar s_bar
     c		+ sigma2 * fa_ubar * fb_ubar 				!u_bar u_bar->u_bar u_bar
     c		+ sigma3 * (fa_dbar * fb_d + fa_sbar * fb_s)	!d_bar d->u_bar u, s_bar s->u_bar u
     c		+ sigma4 * fa_ubar * fb_u				!u_bar u->u_bar u
     c		+ sigma5 * fa_ubar * fb_g				!u_bar g->u_bar g
     c		+ sigma7 * fa_g * fb_g					!gg->u_bar u
     c		+ sigma_u1 * (fa_d + fa_s 
     c		+ fa_dbar + fa_sbar) * fb_ubar!du_bar->u_bar d, su_bar ->u_bar s, d_bar u_bar->u_bar d_bar, s_bar u_bar ->u_bar s_bar
     c		+ sigma_u3 * (fa_d * fb_dbar + fa_s * fb_sbar)	!dd_bar->u_bar u, ss_bar->u_bar u
     c		+ sigma_u4 * fa_u * fb_ubar 				!uu_bar->u_bar u
     c		+ sigma_u5 * fa_g * fb_ubar				!gu_bar->u_bar u
	end if

!	c=d
	if (jet .EQ.2) then
		sum_d = sigma1 * fa_d * (fb_u + fb_s + fb_ubar + fb_sbar)		!du->du, ds->ds, du_bar->du_bar, ds_bar->ds_bar
     c		+ sigma2 * fa_d * fb_d 							!dd->dd
     c		+ sigma3 * (fa_u * fb_ubar + fa_s * fb_sbar)			!uu_bar->dd_bar, ss_bar->dd_bar
     c		+ sigma4 * fa_d * fb_dbar						!dd_bar->dd_bar
     c		+ sigma5 * fa_d * fb_g							!dg->dg
     c		+ sigma7 * fa_g * fb_g							!gg->dd_bar
     c		+ sigma_u1 * (fa_u + fa_s + fa_ubar + fa_sbar) * fb_d		!ud->du, sd->ds, u_bar d->du_bar, s_bar d->ds_bar
     c		+ sigma_u3 * (fa_ubar * fb_u + fa_sbar *fb_s)			!u_bar u->dd_bar, s_bar s->dd_bar
     c		+ sigma_u4 * fa_dbar * fb_d 						!d_bar d->dd_bar
     c		+ sigma_u5 * fa_g * fb_d						!gd->dg
	end if

!	c=d_bar
	if (jet .EQ. -2) then
		sum_dbar = sigma1 * fa_dbar * (fb_u + fb_s 
     c		+ fb_ubar + fb_sbar)		!d_bar u->d_bar u , d_bar s->d_bar s, d_bar u_bar->d_bar u_bar, d_bar s_bar->d_bar s_bar
     c		+ sigma2 * fa_dbar * fb_dbar 				!d_bar d_bar->d_bar d_bar
     c		+ sigma3 * (fa_ubar * fb_u + fa_sbar * fb_s)	!u_bar u->d_bar d, s_bar s->d_bar d
     c		+ sigma4 * fa_dbar * fb_d				!d_bar d->d_bar d
     c		+ sigma5 * fa_dbar * fb_g				!d_bar g->d_bar g
     c		+ sigma7 * fa_g * fb_g					!gg->d_bar d
     c		+ sigma_u1 * (fa_u + fa_s 
     c		+ fa_ubar + fa_sbar) * fb_dbar!ud_bar->d_bar u, sd_bar ->d_bar s, u_bar d_bar->d_bar u_bar, s_bar d_bar ->d_bar s_bar
     c		+ sigma_u3 * (fa_u * fb_ubar + fa_s * fb_sbar)	!uu_bar->d_bar d, ss_bar->d_bar d
     c		+ sigma_u4 * fa_d * fb_dbar 				!dd_bar->d_bar d
     c		+ sigma_u5 * fa_g * fb_dbar				!gd_bar->d_bar d
	end if

!	c=s
	if (jet .EQ. 3) then
		sum_s = sigma1 * fa_s * (fb_d + fb_u + fb_dbar + fb_ubar)		!sd->sd, su->su, sd_bar->sd_bar, su_bar->su_bar
     c		+ sigma2 * fa_s * fb_s 							!ss->ss
     c		+ sigma3 * (fa_d * fb_dbar + fa_u * fb_ubar)			!dd_bar->ss_bar, uu_bar->ss_bar
     c		+ sigma4 * fa_s * fb_sbar						!ss_bar->ss_bar
     c		+ sigma5 * fa_s * fb_g							!sg->sg
     c		+ sigma7 * fa_g * fb_g							!gg->ss_bar
     c		+ sigma_u1 * (fa_d + fa_u + fa_dbar + fa_ubar) * fb_s		!ds->sd, us->su, d_bar s->sd_bar, u_bar s->su_bar
     c		+ sigma_u3 * (fa_dbar * fb_d + fa_ubar *fb_u)			!d_bar d->ss_bar, u_bar u->ss_bar
     c		+ sigma_u4 * fa_sbar * fb_s 						!s_bar s->ss_bar
     c		+ sigma_u5 * fa_g * fb_s						!gs->sg
	end if

!	c=s_bar
	if (jet .EQ. -3) then
		sum_sbar = sigma1 * fa_sbar * (fb_d + fb_u 
     c		+ fb_dbar + fb_ubar)		!s_bar d->s_bar d , s_bar u->s_bar u, s_bar d_bar->s_bar d_bar, s_bar u_bar->s_bar u_bar
     c		+ sigma2 * fa_sbar * fb_sbar 				!s_bar s_bar->s_bar s_bar
     c		+ sigma3 * (fa_dbar * fb_d + fa_ubar * fb_u)	!d_bar d->s_bar s, u_bar u->s_bar s
     c		+ sigma4 * fa_sbar * fb_s				!s_bar s->s_bar s
     c		+ sigma5 * fa_sbar * fb_g				!s_bar g->s_bar g
     c		+ sigma7 * fa_g * fb_g					!gg->s_bar s
     c		+ sigma_u1 * (fa_d + fa_u 
     c		+ fa_dbar + fa_ubar) * fb_sbar!ds_bar->s_bar d, us_bar ->s_bar u, d_bar s_bar->s_bar d_bar, u_bar s_bar ->s_bar u_bar
     c		+ sigma_u3 * (fa_d * fb_dbar + fa_u * fb_ubar)	!dd_bar->s_bar s, uu_bar->s_bar s
     c		+ sigma_u4 * fa_s * fb_sbar 				!ss_bar->s_bar s
     c		+ sigma_u5 * fa_g * fb_sbar				!gs_bar->s_bar s
	end if
	
!	c=g
	if (jet .EQ. 0) then
		sum_g = sigma5 * fa_g * (fb_u + fb_d + fb_s 
     c		+ fb_ubar + fb_dbar + fb_sbar)		!gu->gu, gd->gd, gs->gs, gu_bar->gu_bar, gd_bar->gd_bar, gs_bar->gs_bar
     c		+ sigma6 * (fa_u * fb_ubar + fa_d * fb_dbar + fa_s * fb_sbar)	!uu_bar->gg, dd_bar->gg, ss_bar->gg
     c		+ sigma8 * fa_g * fb_g 								!gg->gg
     c		+ sigma_u5 * (fa_u + fa_d + fa_s 
     c		+ fa_ubar + fa_dbar + fa_sbar) * fb_g	!ug->gu, dg->gd, sg->gs, u_bar g->gu_bar, d_bar g->gd_bar, s_bar g->gs_bar
     c		+ sigma_u6 * (fa_ubar * fb_u + fa_dbar * fb_d + fa_sbar * fb_s)	!u_bar u->gg, d_bar d->gg, s_bar s->gg
     	end if

!	write (*,*) sum_u, sum_d, sum_ubar, sum_dbar, sum_s, sum_sbar

	if (jet .EQ. 4 .OR. jet .EQ. -4) then

        !       now consider the ccbar production
 
                factor = PI * alpha_s**2 / s**2
               
                MM = 64.0/9.0 / s**2 * ((u-M_charm**2)**2 + (t-M_charm**2)**2 + 2.d0*s*M_charm**2)
                                
!                M_ss = 12.d0 / s**2 * (t-M_charm**2) * (u-M_charm**2) 
!                M_tt = -8.d0/3.d0 / (t-M_charm**2)**2 * (4*M_charm**4 - (t-M_charm**2)*(u-M_charm**2) + 2.D0*M_charm**2*(t-M_charm**2))
!                M_uu = -8.d0/3.d0 / (u-M_charm**2)**2 * (4*M_charm**4 - (u-M_charm**2)*(t-M_charm**2) + 2.D0*M_charm**2*(u-M_charm**2))
!                M_st = 6.d0 / (s*(t-M_charm**2)) * (M_charm**4 - t*(s+t))
!                M_su = 6.d0 / (s*(u-M_charm**2)) * (M_charm**4 - u*(s+u))
!                M_tu = 2.d0/3.d0 * M_charm**2/((t-M_charm**2)*(u-M_charm**2)) * (4.D0*M_charm**2 + (t-M_charm**2) + (u-M_charm**2))
                
                M_ss = 12.d0 / s**2 * (t-M_charm**2) * (u-M_charm**2) 
                M_tt = 8.d0/3.d0 / (t-M_charm**2)**2 * ((t-M_charm**2)*(u-M_charm**2) - 2.D0*M_charm**2*(t+M_charm**2))
                M_uu = 8.d0/3.d0 / (u-M_charm**2)**2 * ((u-M_charm**2)*(t-M_charm**2) - 2.D0*M_charm**2*(u+M_charm**2))
                M_st = 6.d0 / (s*(t-M_charm**2)) * (M_charm**4 - t*(s+t))
                M_su = 6.d0 / (s*(u-M_charm**2)) * (M_charm**4 - u*(s+u))
                M_tu = -2.d0/3.d0 * M_charm**2/((t-M_charm**2)*(u-M_charm**2)) * (s - 4.D0*M_charm**2)

                sigma_qq2cc = factor * (1.D0/16.D0) * MM

                sigma_gg2cc = factor * (1.D0/16.D0) * (M_ss + M_tt + M_uu + M_st + M_su + M_tu) 
 
        !       now consider single c(cbar) production

                factor = PI * alpha_s**2 / (s-M_charm**2)**2
                 
                MM = 64.0/9.0 / t**2 * ((u-M_charm**2)**2 + (s-M_charm**2)**2 + 2.d0*t*M_charm**2)
                                
                M_tt = -12.d0 / t**2 * (s-M_charm**2) * (u-M_charm**2) 
                M_ss = -8.d0/3.d0 / (s-M_charm**2)**2 * ((s-M_charm**2)*(u-M_charm**2) - 2.D0*M_charm**2*(s+M_charm**2))
                M_uu = -8.d0/3.d0 / (u-M_charm**2)**2 * ((u-M_charm**2)*(s-M_charm**2) - 2.D0*M_charm**2*(u+M_charm**2))
                M_st = -6.d0 / (t*(s-M_charm**2)) * (M_charm**4 - s*(t+s))
                M_tu = -6.d0 / (t*(u-M_charm**2)) * (M_charm**4 - u*(t+u))
                M_su = 2.d0/3.d0 * M_charm**2/((s-M_charm**2)*(u-M_charm**2)) * (t - 4.D0*M_charm**2)

                sigma_qc2qc = factor * (1.D0/16.D0) * MM
                sigma_cq2cq = sigma_qc2qc

                sigma_gc2gc = factor * (1.D0/16.D0) * 8.d0/3.d0 * (M_ss + M_tt + M_uu + M_st + M_su + M_tu) 
                sigma_cg2cg = sigma_gc2gc

        !       then change 3 & 4 for qc->qc & gc->gc to get qc->cq & gc->cg (t & u)
                tmp = u
                u = t 
                t = tmp
                  
                MM = 64.0/9.0 / t**2 * ((u-M_charm**2)**2 + (s-M_charm**2)**2 + 2.d0*t*M_charm**2)
                  
                MM = 64.0/9.0 / t**2 * ((u-M_charm**2)**2 + (s-M_charm**2)**2 + 2.d0*t*M_charm**2)
                                
                M_tt = -12.d0 / t**2 * (s-M_charm**2) * (u-M_charm**2) 
                M_ss = -8.d0/3.d0 / (s-M_charm**2)**2 * ((s-M_charm**2)*(u-M_charm**2) - 2.D0*M_charm**2*(s+M_charm**2))
                M_uu = -8.d0/3.d0 / (u-M_charm**2)**2 * ((u-M_charm**2)*(s-M_charm**2) - 2.D0*M_charm**2*(u+M_charm**2))
                M_st = -6.d0 / (t*(s-M_charm**2)) * (M_charm**4 - s*(t+s))
                M_tu = -6.d0 / (t*(u-M_charm**2)) * (M_charm**4 - u*(t+u))
                M_su = 2.d0/3.d0 * M_charm**2/((s-M_charm**2)*(u-M_charm**2)) * (t - 4.D0*M_charm**2)
              
                sigma_qc2cq = factor * (1.D0/16.D0) * MM
                sigma_cq2qc = sigma_qc2cq

                sigma_gc2cg = factor * (1.D0/16.D0) * (M_ss + M_tt + M_uu + M_st + M_su + M_tu) 
                sigma_cg2gc = sigma_gc2cg
              
        !        write (*,*) sigma_qq2cc, sigma_gg2cc, sigma_qc2qc, sigma_gc2gc, sigma_qc2cq, sigma_gc2cg
                 
        !       change back to the original 
                tmp = u
                u = t
                t = tmp
         
        end if

!	c=charm

	if (jet .EQ. 4) then
		sum_charm = sigma_qq2cc * (fa_u * fb_ubar + fa_d * fb_dbar + fa_s * fb_sbar)	!uu_bar->cc_bar, dd_bar->cc_bar, ss_bar->cc_bar
     c  	+ sigma_qq2cc * (fa_ubar * fb_u + fa_dbar * fb_d + fa_sbar * fb_s) 		!u_baru->cc_bar, d_baru->cc_bar, s_bars->cc_bar 
     c		+ sigma_gg2cc * fa_g * fb_g								!gg->cc_bar
     c		+ sigma_cq2cq * fa_c * (fb_u + fb_d + fb_s)					!cq->cq
     c		+ sigma_cg2cg * fa_c * fb_g							!cg->cg
     c		+ sigma_qc2cq * (fa_u + fa_d + fa_s) * fb_c					!qc->cq
     c		+ sigma_gc2cg * fa_g * fb_c							!gc->cg
	end if
        
!	c=charm_bar

	if (jet .EQ. -4) then
		sum_charm_bar = sigma_qq2cc * (fa_u * fb_ubar + fa_d * fb_dbar + fa_s * fb_sbar)
     c  	+ sigma_qq2cc * (fa_ubar * fb_u + fa_dbar * fb_d + fa_sbar * fb_s)
     c		+ sigma_gg2cc * fa_g * fb_g
     c		+ sigma_cq2cq * fa_cbar * (fb_u + fb_d + fb_s)					!cq->cq
     c		+ sigma_cg2cg * fa_cbar * fb_g							!cg->cg
     c		+ sigma_qc2cq * (fa_u + fa_d + fa_s) * fb_cbar					!qc->cq
     c		+ sigma_gc2cg * fa_g * fb_cbar							!gc->cg
	end if

	if (jet .EQ. 0) sum_parton = sum_g
     	if (jet .EQ. 1) sum_parton = sum_u
	if (jet .EQ. 2) sum_parton = sum_d
	if (jet .EQ. 3) sum_parton = sum_s
	if (jet .EQ. -1) sum_parton = sum_ubar
	if (jet .EQ. -2) sum_parton = sum_dbar
	if (jet .EQ. -3) sum_parton = sum_sbar

	if (jet .EQ. 4) sum_parton = sum_charm
	if (jet .EQ. -4) sum_parton = sum_charm_bar

	funjet = sum_parton * 2.0 * x_a * x_b / PI / (2.0 * x_a - xT_c * exp(y))

	funjet = funjet * k_jet	! jet initial distribution(cross section), including NLO effect.

	end function

!	We can classify according to target "c" particle, if c=u, we have:
!	ud->ud,us->us,ud_bar->ud_bar,us_bar->us_bar, uu->uu(*), 
!	dd_bar->uu_bar,ss_bar->uu_bar, uu_bar->uu_bar, ug->ug, gg->uu_bar(*), total 10.	
!	since we fixed "c" parton, we have freedom to specify "a", "b" in the initial states, 
!	we need consider exchanging 1,2, where (*) means we don't need to exchang 1,2.
