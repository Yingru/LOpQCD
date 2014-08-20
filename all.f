c	This is the program to get parton distribution functions.

C============================================================================= 
C============================================================================= 
C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.

C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------

C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C

C   The Table_Files are assumed to be in the working directory.
C   

C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   

C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.

C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
!          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine SetCtq5 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=9)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5m.tbl', 'cteq5d.tbl', 'cteq5l.tbl', 'cteq5hj.tbl'
     > , 'cteq5hq.tbl', 'cteq5f3.tbl', 'cteq5f4.tbl'
     > , 'cteq5m1.tbl', 'ctq5hq1.tbl'  /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >     //'in SetCtq5!!'
      Stop
C                             ********************
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      SUBROUTINE DQAG (F, A, B, EPSABS, EPSREL, KEY, RESULT, ABSERR,
     +   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAG
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (QAG-S, DQAG-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
C             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
C             QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     Declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
C
C            KEY    - Integer
C                     Key for choice of local integration rule
C                     A GAUSS-KRONROD PAIR is used with
C                       7 - 15 POINTS If KEY.LT.2,
C                      10 - 21 POINTS If KEY = 2,
C                      15 - 31 POINTS If KEY = 3,
C                      20 - 41 POINTS If KEY = 4,
C                      25 - 51 POINTS If KEY = 5,
C                      30 - 61 POINTS If KEY.GT.5.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     Which should EQUAL or EXCEED ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for RESULT and ERROR are
C                             Less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C                      ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). HOWEVER, If
C                             this yield no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (I.E. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL) One
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             INTEGRATOR on the SUBRANGES. If possible,
C                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
C                             should be used which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 AND
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set
C                             to zero.
C                             EXCEPT when LENW is invalid, IWORK(1),
C                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
C                             set to zero, WORK(1) is set to A and
C                             WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    Dimensioning parameter for IWORK
C                    Limit determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    If LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for work
C                    LENW must be at least LIMIT*4.
C                    IF LENW.LT.LIMIT*4, the routine will end with
C                    IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least limit, the first K
C                    elements of which contain pointers to the error
C                    estimates over the subintervals, such that
C                    WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K))
C                    form a decreasing sequence with K = LAST If
C                    LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left end
C                    points of the subintervals in the partition of
C                     (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
C                     right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
C                     the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAG
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C***FIRST EXECUTABLE STATEMENT  DQAG
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (LIMIT.GE.1 .AND. LENW.GE.LIMIT*4) THEN
C
C        PREPARE CALL FOR DQAGE.
C
         L1 = LIMIT+1
         L2 = LIMIT+L1
         L3 = LIMIT+L2
C
         CALL DQAGE(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR,NEVAL,
     1     IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C        CALL ERROR HANDLER IF NECESSARY.
C
         LVL = 0
      ENDIF
C
      IF (IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAG',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
*DECK DQAGE
      SUBROUTINE DQAGE (F, A, B, EPSABS, EPSREL, KEY, LIMIT, RESULT,
     +   ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
C***BEGIN PROLOGUE  DQAGE
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral   I = Integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESLT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (QAGE-S, DQAGE-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
C             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
C             QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            KEY    - Integer
C                     Key for choice of local integration rule
C                     A Gauss-Kronrod pair is used with
C                          7 - 15 points if KEY.LT.2,
C                         10 - 21 points if KEY = 2,
C                         15 - 31 points if KEY = 3,
C                         20 - 41 points if KEY = 4,
C                         25 - 51 points if KEY = 5,
C                         30 - 61 points if KEY.GT.5.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.1.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for result and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value
C                             of LIMIT.
C                             However, if this yields no improvement it
C                             is rather advised to analyze the integrand
C                             in order to determine the integration
C                             difficulties. If the position of a local
C                             difficulty can be determined(e.g.
C                             SINGULARITY, DISCONTINUITY within the
C                             interval) one will probably gain from
C                             splitting up the interval at this point
C                             and calling the integrator on the
C                             subranges. If possible, an appropriate
C                             special-purpose integrator should be used
C                             which is designed for handling the type of
C                             difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1) ,
C                             ELIST(1) and IORD(1) are set to zero.
C                             ALIST(1) and BLIST(1) are set to A and B
C                             respectively.
C
C            ALIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the left
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            BLIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the right
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            RLIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the
C                      integral approximations on the subintervals
C
C            ELIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the moduli of the
C                      absolute error estimates on the subintervals
C
C            IORD    - Integer
C                      Vector of dimension at least LIMIT, the first K
C                      elements of which are pointers to the
C                      error estimates over the subintervals,
C                      such that ELIST(IORD(1)), ...,
C                      ELIST(IORD(K)) form a decreasing sequence,
C                      with K = LAST if LAST.LE.(LIMIT/2+2), and
C                      K = LIMIT+1-LAST otherwise
C
C            LAST    - Integer
C                      Number of subintervals actually produced in the
C                      subdivision process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DQK15, DQK21, DQK31, DQK41, DQK51, DQK61,
C                    DQPSRT
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQAGE
C
      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,
     1  BLIST,B1,B2,DEFABS,DEFAB1,DEFAB2,D1MACH,ELIST,EPMACH,
     2  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,
     3  RESABS,RESULT,RLIST,UFLOW
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
     1  NRMAX
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  RLIST(*)
C
      EXTERNAL F
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                      (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGE
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      IF(EPSABS.LE.0.0D+00.AND.
     1  EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      NEVAL = 0
      IF(KEYF.EQ.1) CALL DQK15(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.2) CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.3) CALL DQK31(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.4) CALL DQK41(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.5) CALL DQK51(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.6) CALL DQK61(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF(ABSERR.LE.0.5D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     1  .OR.ABSERR.EQ.0.0D+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        IF(KEYF.EQ.1) CALL DQK15(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.2) CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.3) CALL DQK31(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.4) CALL DQK41(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.5) CALL DQK51(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.6) CALL DQK61(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.1) CALL DQK15(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.2) CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.3) CALL DQK31(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.4) CALL DQK41(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.5) CALL DQK51(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.6) CALL DQK61(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+1
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
        IF(ABS(RLIST(MAXERR)-AREA12).LE.0.1D-04*ABS(AREA12)
     1  .AND.ERRO12.GE.0.99D+00*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
C           EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*
     1  EPMACH)*(ABS(A2)+0.1D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   20   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   40 RESULT = 0.0D+00
      DO 50 K=1,LAST
        RESULT = RESULT+RLIST(K)
   50 CONTINUE
      ABSERR = ERRSUM
   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
  999 RETURN
      END
*DECK DQK15
      SUBROUTINE DQK15 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK15
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK15-S, DQK15-D)
C***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       Result is computed by applying the 15-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 7-POINT GAUSS RULE(RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK15
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.1294849661 6886969327 0611432679 082 D0 /
      DATA WG  (  2) / 0.2797053914 8927666790 1467771423 780 D0 /
      DATA WG  (  3) / 0.3818300505 0511894495 0369775488 975 D0 /
      DATA WG  (  4) / 0.4179591836 7346938775 5102040816 327 D0 /
C
      DATA XGK (  1) / 0.9914553711 2081263920 6854697526 329 D0 /
      DATA XGK (  2) / 0.9491079123 4275852452 6189684047 851 D0 /
      DATA XGK (  3) / 0.8648644233 5976907278 9712788640 926 D0 /
      DATA XGK (  4) / 0.7415311855 9939443986 3864773280 788 D0 /
      DATA XGK (  5) / 0.5860872354 6769113029 4144838258 730 D0 /
      DATA XGK (  6) / 0.4058451513 7739716690 6606412076 961 D0 /
      DATA XGK (  7) / 0.2077849550 0789846760 0689403773 245 D0 /
      DATA XGK (  8) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0229353220 1052922496 3732008058 970 D0 /
      DATA WGK (  2) / 0.0630920926 2997855329 0700663189 204 D0 /
      DATA WGK (  3) / 0.1047900103 2225018383 9876322541 518 D0 /
      DATA WGK (  4) / 0.1406532597 1552591874 5189590510 238 D0 /
      DATA WGK (  5) / 0.1690047266 3926790282 6583426598 550 D0 /
      DATA WGK (  6) / 0.1903505780 6478540991 3256402421 014 D0 /
      DATA WGK (  7) / 0.2044329400 7529889241 4161999234 649 D0 /
      DATA WGK (  8) / 0.2094821410 8472782801 2999174891 714 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK15
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = ABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQK21
      SUBROUTINE DQK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK21
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK21-S, DQK21-D)
C***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the driver program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 21-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 10-POINT GAUSS RULE (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK21
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0666713443 0868813759 3568809893 332 D0 /
      DATA WG  (  2) / 0.1494513491 5058059314 5776339657 697 D0 /
      DATA WG  (  3) / 0.2190863625 1598204399 5534934228 163 D0 /
      DATA WG  (  4) / 0.2692667193 0999635509 1226921569 469 D0 /
      DATA WG  (  5) / 0.2955242247 1475287017 3892994651 338 D0 /
C
      DATA XGK (  1) / 0.9956571630 2580808073 5527280689 003 D0 /
      DATA XGK (  2) / 0.9739065285 1717172007 7964012084 452 D0 /
      DATA XGK (  3) / 0.9301574913 5570822600 1207180059 508 D0 /
      DATA XGK (  4) / 0.8650633666 8898451073 2096688423 493 D0 /
      DATA XGK (  5) / 0.7808177265 8641689706 3717578345 042 D0 /
      DATA XGK (  6) / 0.6794095682 9902440623 4327365114 874 D0 /
      DATA XGK (  7) / 0.5627571346 6860468333 9000099272 694 D0 /
      DATA XGK (  8) / 0.4333953941 2924719079 9265943165 784 D0 /
      DATA XGK (  9) / 0.2943928627 0146019813 1126603103 866 D0 /
      DATA XGK ( 10) / 0.1488743389 8163121088 4826001129 720 D0 /
      DATA XGK ( 11) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0116946388 6737187427 8064396062 192 D0 /
      DATA WGK (  2) / 0.0325581623 0796472747 8818972459 390 D0 /
      DATA WGK (  3) / 0.0547558965 7435199603 1381300244 580 D0 /
      DATA WGK (  4) / 0.0750396748 1091995276 7043140916 190 D0 /
      DATA WGK (  5) / 0.0931254545 8369760553 5065465083 366 D0 /
      DATA WGK (  6) / 0.1093871588 0229764189 9210590325 805 D0 /
      DATA WGK (  7) / 0.1234919762 6206585107 7958109831 074 D0 /
      DATA WGK (  8) / 0.1347092173 1147332592 8054001771 707 D0 /
      DATA WGK (  9) / 0.1427759385 7706008079 7094273138 717 D0 /
      DATA WGK ( 10) / 0.1477391049 0133849137 4841515972 068 D0 /
      DATA WGK ( 11) / 0.1494455540 0291690566 4936468389 821 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
      FC = F(CENTR)
      RESK = WGK(11)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(11)*ABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQK31
      SUBROUTINE DQK31 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK31
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK31-S, DQK31-D)
C***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 31-POINT
C                       GAUSS-KRONROD RULE (RESK), obtained by optimal
C                       addition of abscissae to the 15-POINT GAUSS
C                       RULE (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the modulus,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK31
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 15-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0307532419 9611726835 4628393577 204 D0 /
      DATA WG  (  2) / 0.0703660474 8810812470 9267416450 667 D0 /
      DATA WG  (  3) / 0.1071592204 6717193501 1869546685 869 D0 /
      DATA WG  (  4) / 0.1395706779 2615431444 7804794511 028 D0 /
      DATA WG  (  5) / 0.1662692058 1699393355 3200860481 209 D0 /
      DATA WG  (  6) / 0.1861610000 1556221102 6800561866 423 D0 /
      DATA WG  (  7) / 0.1984314853 2711157645 6118326443 839 D0 /
      DATA WG  (  8) / 0.2025782419 2556127288 0620199967 519 D0 /
C
      DATA XGK (  1) / 0.9980022986 9339706028 5172840152 271 D0 /
      DATA XGK (  2) / 0.9879925180 2048542848 9565718586 613 D0 /
      DATA XGK (  3) / 0.9677390756 7913913425 7347978784 337 D0 /
      DATA XGK (  4) / 0.9372733924 0070590430 7758947710 209 D0 /
      DATA XGK (  5) / 0.8972645323 4408190088 2509656454 496 D0 /
      DATA XGK (  6) / 0.8482065834 1042721620 0648320774 217 D0 /
      DATA XGK (  7) / 0.7904185014 4246593296 7649294817 947 D0 /
      DATA XGK (  8) / 0.7244177313 6017004741 6186054613 938 D0 /
      DATA XGK (  9) / 0.6509967412 9741697053 3735895313 275 D0 /
      DATA XGK ( 10) / 0.5709721726 0853884753 7226737253 911 D0 /
      DATA XGK ( 11) / 0.4850818636 4023968069 3655740232 351 D0 /
      DATA XGK ( 12) / 0.3941513470 7756336989 7207370981 045 D0 /
      DATA XGK ( 13) / 0.2991800071 5316881216 6780024266 389 D0 /
      DATA XGK ( 14) / 0.2011940939 9743452230 0628303394 596 D0 /
      DATA XGK ( 15) / 0.1011420669 1871749902 7074231447 392 D0 /
      DATA XGK ( 16) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0053774798 7292334898 7792051430 128 D0 /
      DATA WGK (  2) / 0.0150079473 2931612253 8374763075 807 D0 /
      DATA WGK (  3) / 0.0254608473 2671532018 6874001019 653 D0 /
      DATA WGK (  4) / 0.0353463607 9137584622 2037948478 360 D0 /
      DATA WGK (  5) / 0.0445897513 2476487660 8227299373 280 D0 /
      DATA WGK (  6) / 0.0534815246 9092808726 5343147239 430 D0 /
      DATA WGK (  7) / 0.0620095678 0067064028 5139230960 803 D0 /
      DATA WGK (  8) / 0.0698541213 1872825870 9520077099 147 D0 /
      DATA WGK (  9) / 0.0768496807 5772037889 4432777482 659 D0 /
      DATA WGK ( 10) / 0.0830805028 2313302103 8289247286 104 D0 /
      DATA WGK ( 11) / 0.0885644430 5621177064 7275443693 774 D0 /
      DATA WGK ( 12) / 0.0931265981 7082532122 5486872747 346 D0 /
      DATA WGK ( 13) / 0.0966427269 8362367850 5179907627 589 D0 /
      DATA WGK ( 14) / 0.0991735987 2179195933 2393173484 603 D0 /
      DATA WGK ( 15) / 0.1007698455 2387559504 4946662617 570 D0 /
      DATA WGK ( 16) / 0.1013300070 1479154901 7374792767 493 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C***FIRST EXECUTABLE STATEMENT  DQK31
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(8)*FC
      RESK = WGK(16)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,8
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(16)*ABS(FC-RESKH)
      DO 20 J=1,15
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQK41
      SUBROUTINE DQK41 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK41
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK41-S, DQK41-D)
C***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the calling program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 41-POINT
C                       GAUSS-KRONROD RULE (RESK) obtained by optimal
C                       addition of abscissae to the 20-POINT GAUSS
C                       RULE (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK41
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(20),FV2(20),XGK(21),WGK(21),WG(10)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 20-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
C
C           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0176140071 3915211831 1861962351 853 D0 /
      DATA WG  (  2) / 0.0406014298 0038694133 1039952274 932 D0 /
      DATA WG  (  3) / 0.0626720483 3410906356 9506535187 042 D0 /
      DATA WG  (  4) / 0.0832767415 7670474872 4758143222 046 D0 /
      DATA WG  (  5) / 0.1019301198 1724043503 6750135480 350 D0 /
      DATA WG  (  6) / 0.1181945319 6151841731 2377377711 382 D0 /
      DATA WG  (  7) / 0.1316886384 4917662689 8494499748 163 D0 /
      DATA WG  (  8) / 0.1420961093 1838205132 9298325067 165 D0 /
      DATA WG  (  9) / 0.1491729864 7260374678 7828737001 969 D0 /
      DATA WG  ( 10) / 0.1527533871 3072585069 8084331955 098 D0 /
C
      DATA XGK (  1) / 0.9988590315 8827766383 8315576545 863 D0 /
      DATA XGK (  2) / 0.9931285991 8509492478 6122388471 320 D0 /
      DATA XGK (  3) / 0.9815078774 5025025919 3342994720 217 D0 /
      DATA XGK (  4) / 0.9639719272 7791379126 7666131197 277 D0 /
      DATA XGK (  5) / 0.9408226338 3175475351 9982722212 443 D0 /
      DATA XGK (  6) / 0.9122344282 5132590586 7752441203 298 D0 /
      DATA XGK (  7) / 0.8782768112 5228197607 7442995113 078 D0 /
      DATA XGK (  8) / 0.8391169718 2221882339 4529061701 521 D0 /
      DATA XGK (  9) / 0.7950414288 3755119835 0638833272 788 D0 /
      DATA XGK ( 10) / 0.7463319064 6015079261 4305070355 642 D0 /
      DATA XGK ( 11) / 0.6932376563 3475138480 5490711845 932 D0 /
      DATA XGK ( 12) / 0.6360536807 2651502545 2836696226 286 D0 /
      DATA XGK ( 13) / 0.5751404468 1971031534 2946036586 425 D0 /
      DATA XGK ( 14) / 0.5108670019 5082709800 4364050955 251 D0 /
      DATA XGK ( 15) / 0.4435931752 3872510319 9992213492 640 D0 /
      DATA XGK ( 16) / 0.3737060887 1541956067 2548177024 927 D0 /
      DATA XGK ( 17) / 0.3016278681 1491300432 0555356858 592 D0 /
      DATA XGK ( 18) / 0.2277858511 4164507808 0496195368 575 D0 /
      DATA XGK ( 19) / 0.1526054652 4092267550 5220241022 678 D0 /
      DATA XGK ( 20) / 0.0765265211 3349733375 4640409398 838 D0 /
      DATA XGK ( 21) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0030735837 1852053150 1218293246 031 D0 /
      DATA WGK (  2) / 0.0086002698 5564294219 8661787950 102 D0 /
      DATA WGK (  3) / 0.0146261692 5697125298 3787960308 868 D0 /
      DATA WGK (  4) / 0.0203883734 6126652359 8010231432 755 D0 /
      DATA WGK (  5) / 0.0258821336 0495115883 4505067096 153 D0 /
      DATA WGK (  6) / 0.0312873067 7703279895 8543119323 801 D0 /
      DATA WGK (  7) / 0.0366001697 5820079803 0557240707 211 D0 /
      DATA WGK (  8) / 0.0416688733 2797368626 3788305936 895 D0 /
      DATA WGK (  9) / 0.0464348218 6749767472 0231880926 108 D0 /
      DATA WGK ( 10) / 0.0509445739 2372869193 2707670050 345 D0 /
      DATA WGK ( 11) / 0.0551951053 4828599474 4832372419 777 D0 /
      DATA WGK ( 12) / 0.0591114008 8063957237 4967220648 594 D0 /
      DATA WGK ( 13) / 0.0626532375 5478116802 5870122174 255 D0 /
      DATA WGK ( 14) / 0.0658345971 3361842211 1563556969 398 D0 /
      DATA WGK ( 15) / 0.0686486729 2852161934 5623411885 368 D0 /
      DATA WGK ( 16) / 0.0710544235 5344406830 5790361723 210 D0 /
      DATA WGK ( 17) / 0.0730306903 3278666749 5189417658 913 D0 /
      DATA WGK ( 18) / 0.0745828754 0049918898 6581418362 488 D0 /
      DATA WGK ( 19) / 0.0757044976 8455667465 9542775376 617 D0 /
      DATA WGK ( 20) / 0.0763778676 7208073670 5502835038 061 D0 /
      DATA WGK ( 21) / 0.0766007119 1799965644 5049901530 102 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
C                    TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK41
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
      FC = F(CENTR)
      RESK = WGK(21)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,10
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,10
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(21)*ABS(FC-RESKH)
      DO 20 J=1,20
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQK51
      SUBROUTINE DQK51 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK51
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK51-S, DQK51-D)
C***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subroutine defining the integrand
C                       function F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the calling program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 51-point
C                       Kronrod rule (RESK) obtained by optimal addition
C                       of abscissae to the 25-point Gauss rule (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   910819  Added WGK(26) to code.  (WRB)
C***END PROLOGUE  DQK51
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(25),FV2(25),XGK(26),WGK(26),WG(13)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 25-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0113937985 0102628794 7902964113 235 D0 /
      DATA WG  (  2) / 0.0263549866 1503213726 1901815295 299 D0 /
      DATA WG  (  3) / 0.0409391567 0130631265 5623487711 646 D0 /
      DATA WG  (  4) / 0.0549046959 7583519192 5936891540 473 D0 /
      DATA WG  (  5) / 0.0680383338 1235691720 7187185656 708 D0 /
      DATA WG  (  6) / 0.0801407003 3500101801 3234959669 111 D0 /
      DATA WG  (  7) / 0.0910282619 8296364981 1497220702 892 D0 /
      DATA WG  (  8) / 0.1005359490 6705064420 2206890392 686 D0 /
      DATA WG  (  9) / 0.1085196244 7426365311 6093957050 117 D0 /
      DATA WG  ( 10) / 0.1148582591 4571164833 9325545869 556 D0 /
      DATA WG  ( 11) / 0.1194557635 3578477222 8178126512 901 D0 /
      DATA WG  ( 12) / 0.1222424429 9031004168 8959518945 852 D0 /
      DATA WG  ( 13) / 0.1231760537 2671545120 3902873079 050 D0 /
C
      DATA XGK (  1) / 0.9992621049 9260983419 3457486540 341 D0 /
      DATA XGK (  2) / 0.9955569697 9049809790 8784946893 902 D0 /
      DATA XGK (  3) / 0.9880357945 3407724763 7331014577 406 D0 /
      DATA XGK (  4) / 0.9766639214 5951751149 8315386479 594 D0 /
      DATA XGK (  5) / 0.9616149864 2584251241 8130033660 167 D0 /
      DATA XGK (  6) / 0.9429745712 2897433941 4011169658 471 D0 /
      DATA XGK (  7) / 0.9207471152 8170156174 6346084546 331 D0 /
      DATA XGK (  8) / 0.8949919978 7827536885 1042006782 805 D0 /
      DATA XGK (  9) / 0.8658470652 9327559544 8996969588 340 D0 /
      DATA XGK ( 10) / 0.8334426287 6083400142 1021108693 570 D0 /
      DATA XGK ( 11) / 0.7978737979 9850005941 0410904994 307 D0 /
      DATA XGK ( 12) / 0.7592592630 3735763057 7282865204 361 D0 /
      DATA XGK ( 13) / 0.7177664068 1308438818 6654079773 298 D0 /
      DATA XGK ( 14) / 0.6735663684 7346836448 5120633247 622 D0 /
      DATA XGK ( 15) / 0.6268100990 1031741278 8122681624 518 D0 /
      DATA XGK ( 16) / 0.5776629302 4122296772 3689841612 654 D0 /
      DATA XGK ( 17) / 0.5263252843 3471918259 9623778158 010 D0 /
      DATA XGK ( 18) / 0.4730027314 4571496052 2182115009 192 D0 /
      DATA XGK ( 19) / 0.4178853821 9303774885 1814394594 572 D0 /
      DATA XGK ( 20) / 0.3611723058 0938783773 5821730127 641 D0 /
      DATA XGK ( 21) / 0.3030895389 3110783016 7478909980 339 D0 /
      DATA XGK ( 22) / 0.2438668837 2098843204 5190362797 452 D0 /
      DATA XGK ( 23) / 0.1837189394 2104889201 5969888759 528 D0 /
      DATA XGK ( 24) / 0.1228646926 1071039638 7359818808 037 D0 /
      DATA XGK ( 25) / 0.0615444830 0568507888 6546392366 797 D0 /
      DATA XGK ( 26) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0019873838 9233031592 6507851882 843 D0 /
      DATA WGK (  2) / 0.0055619321 3535671375 8040236901 066 D0 /
      DATA WGK (  3) / 0.0094739733 8617415160 7207710523 655 D0 /
      DATA WGK (  4) / 0.0132362291 9557167481 3656405846 976 D0 /
      DATA WGK (  5) / 0.0168478177 0912829823 1516667536 336 D0 /
      DATA WGK (  6) / 0.0204353711 4588283545 6568292235 939 D0 /
      DATA WGK (  7) / 0.0240099456 0695321622 0092489164 881 D0 /
      DATA WGK (  8) / 0.0274753175 8785173780 2948455517 811 D0 /
      DATA WGK (  9) / 0.0307923001 6738748889 1109020215 229 D0 /
      DATA WGK ( 10) / 0.0340021302 7432933783 6748795229 551 D0 /
      DATA WGK ( 11) / 0.0371162714 8341554356 0330625367 620 D0 /
      DATA WGK ( 12) / 0.0400838255 0403238207 4839284467 076 D0 /
      DATA WGK ( 13) / 0.0428728450 2017004947 6895792439 495 D0 /
      DATA WGK ( 14) / 0.0455029130 4992178890 9870584752 660 D0 /
      DATA WGK ( 15) / 0.0479825371 3883671390 6392255756 915 D0 /
      DATA WGK ( 16) / 0.0502776790 8071567196 3325259433 440 D0 /
      DATA WGK ( 17) / 0.0523628858 0640747586 4366712137 873 D0 /
      DATA WGK ( 18) / 0.0542511298 8854549014 4543370459 876 D0 /
      DATA WGK ( 19) / 0.0559508112 2041231730 8240686382 747 D0 /
      DATA WGK ( 20) / 0.0574371163 6156783285 3582693939 506 D0 /
      DATA WGK ( 21) / 0.0586896800 2239420796 1974175856 788 D0 /
      DATA WGK ( 22) / 0.0597203403 2417405997 9099291932 562 D0 /
      DATA WGK ( 23) / 0.0605394553 7604586294 5360267517 565 D0 /
      DATA WGK ( 24) / 0.0611285097 1705304830 5859030416 293 D0 /
      DATA WGK ( 25) / 0.0614711898 7142531666 1544131965 264 D0 /
      DATA WGK ( 26) / 0.0615808180 6783293507 8759824240 055 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK51
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(13)*FC
      RESK = WGK(26)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,12
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,13
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(26)*ABS(FC-RESKH)
      DO 20 J=1,25
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQK61
      SUBROUTINE DQK61 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  DQK61
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK61-S, DQK61-D)
C***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rule
C        Standard fortran subroutine
C        Double precision version
C
C
C        PARAMETERS
C         ON ENTRY
C           F      - Double precision
C                    Function subprogram defining the integrand
C                    function F(X). The actual name for F needs to be
C                    declared E X T E R N A L in the calling program.
C
C           A      - Double precision
C                    Lower limit of integration
C
C           B      - Double precision
C                    Upper limit of integration
C
C         ON RETURN
C           RESULT - Double precision
C                    Approximation to the integral I
C                    RESULT is computed by applying the 61-point
C                    Kronrod rule (RESK) obtained by optimal addition of
C                    abscissae to the 30-point Gauss rule (RESG).
C
C           ABSERR - Double precision
C                    Estimate of the modulus of the absolute error,
C                    which should equal or exceed ABS(I-RESULT)
C
C           RESABS - Double precision
C                    Approximation to the integral J
C
C           RESASC - Double precision
C                    Approximation to the integral of ABS(F-I/(B-A))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK61
C
      DOUBLE PRECISION A,DABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(30),FV2(30),XGK(31),WGK(31),WG(15)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
C           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
C           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
C                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
C                   GAUSS RULE
C                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
C                   TO THE 30-POINT GAUSS RULE
C
C           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
C
C           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0079681924 9616660561 5465883474 674 D0 /
      DATA WG  (  2) / 0.0184664683 1109095914 2302131912 047 D0 /
      DATA WG  (  3) / 0.0287847078 8332336934 9719179611 292 D0 /
      DATA WG  (  4) / 0.0387991925 6962704959 6801936446 348 D0 /
      DATA WG  (  5) / 0.0484026728 3059405290 2938140422 808 D0 /
      DATA WG  (  6) / 0.0574931562 1761906648 1721689402 056 D0 /
      DATA WG  (  7) / 0.0659742298 8218049512 8128515115 962 D0 /
      DATA WG  (  8) / 0.0737559747 3770520626 8243850022 191 D0 /
      DATA WG  (  9) / 0.0807558952 2942021535 4694938460 530 D0 /
      DATA WG  ( 10) / 0.0868997872 0108297980 2387530715 126 D0 /
      DATA WG  ( 11) / 0.0921225222 3778612871 7632707087 619 D0 /
      DATA WG  ( 12) / 0.0963687371 7464425963 9468626351 810 D0 /
      DATA WG  ( 13) / 0.0995934205 8679526706 2780282103 569 D0 /
      DATA WG  ( 14) / 0.1017623897 4840550459 6428952168 554 D0 /
      DATA WG  ( 15) / 0.1028526528 9355884034 1285636705 415 D0 /
C
      DATA XGK (  1) / 0.9994844100 5049063757 1325895705 811 D0 /
      DATA XGK (  2) / 0.9968934840 7464954027 1630050918 695 D0 /
      DATA XGK (  3) / 0.9916309968 7040459485 8628366109 486 D0 /
      DATA XGK (  4) / 0.9836681232 7974720997 0032581605 663 D0 /
      DATA XGK (  5) / 0.9731163225 0112626837 4693868423 707 D0 /
      DATA XGK (  6) / 0.9600218649 6830751221 6871025581 798 D0 /
      DATA XGK (  7) / 0.9443744447 4855997941 5831324037 439 D0 /
      DATA XGK (  8) / 0.9262000474 2927432587 9324277080 474 D0 /
      DATA XGK (  9) / 0.9055733076 9990779854 6522558925 958 D0 /
      DATA XGK ( 10) / 0.8825605357 9205268154 3116462530 226 D0 /
      DATA XGK ( 11) / 0.8572052335 4606109895 8658510658 944 D0 /
      DATA XGK ( 12) / 0.8295657623 8276839744 2898119732 502 D0 /
      DATA XGK ( 13) / 0.7997278358 2183908301 3668942322 683 D0 /
      DATA XGK ( 14) / 0.7677774321 0482619491 7977340974 503 D0 /
      DATA XGK ( 15) / 0.7337900624 5322680472 6171131369 528 D0 /
      DATA XGK ( 16) / 0.6978504947 9331579693 2292388026 640 D0 /
      DATA XGK ( 17) / 0.6600610641 2662696137 0053668149 271 D0 /
      DATA XGK ( 18) / 0.6205261829 8924286114 0477556431 189 D0 /
      DATA XGK ( 19) / 0.5793452358 2636169175 6024932172 540 D0 /
      DATA XGK ( 20) / 0.5366241481 4201989926 4169793311 073 D0 /
      DATA XGK ( 21) / 0.4924804678 6177857499 3693061207 709 D0 /
      DATA XGK ( 22) / 0.4470337695 3808917678 0609900322 854 D0 /
      DATA XGK ( 23) / 0.4004012548 3039439253 5476211542 661 D0 /
      DATA XGK ( 24) / 0.3527047255 3087811347 1037207089 374 D0 /
      DATA XGK ( 25) / 0.3040732022 7362507737 2677107199 257 D0 /
      DATA XGK ( 26) / 0.2546369261 6788984643 9805129817 805 D0 /
      DATA XGK ( 27) / 0.2045251166 8230989143 8957671002 025 D0 /
      DATA XGK ( 28) / 0.1538699136 0858354696 3794672743 256 D0 /
      DATA XGK ( 29) / 0.1028069379 6673703014 7096751318 001 D0 /
      DATA XGK ( 30) / 0.0514718425 5531769583 3025213166 723 D0 /
      DATA XGK ( 31) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0013890136 9867700762 4551591226 760 D0 /
      DATA WGK (  2) / 0.0038904611 2709988405 1267201844 516 D0 /
      DATA WGK (  3) / 0.0066307039 1593129217 3319826369 750 D0 /
      DATA WGK (  4) / 0.0092732796 5951776342 8441146892 024 D0 /
      DATA WGK (  5) / 0.0118230152 5349634174 2232898853 251 D0 /
      DATA WGK (  6) / 0.0143697295 0704580481 2451432443 580 D0 /
      DATA WGK (  7) / 0.0169208891 8905327262 7572289420 322 D0 /
      DATA WGK (  8) / 0.0194141411 9394238117 3408951050 128 D0 /
      DATA WGK (  9) / 0.0218280358 2160919229 7167485738 339 D0 /
      DATA WGK ( 10) / 0.0241911620 7808060136 5686370725 232 D0 /
      DATA WGK ( 11) / 0.0265099548 8233310161 0601709335 075 D0 /
      DATA WGK ( 12) / 0.0287540487 6504129284 3978785354 334 D0 /
      DATA WGK ( 13) / 0.0309072575 6238776247 2884252943 092 D0 /
      DATA WGK ( 14) / 0.0329814470 5748372603 1814191016 854 D0 /
      DATA WGK ( 15) / 0.0349793380 2806002413 7499670731 468 D0 /
      DATA WGK ( 16) / 0.0368823646 5182122922 3911065617 136 D0 /
      DATA WGK ( 17) / 0.0386789456 2472759295 0348651532 281 D0 /
      DATA WGK ( 18) / 0.0403745389 5153595911 1995279752 468 D0 /
      DATA WGK ( 19) / 0.0419698102 1516424614 7147541285 970 D0 /
      DATA WGK ( 20) / 0.0434525397 0135606931 6831728117 073 D0 /
      DATA WGK ( 21) / 0.0448148001 3316266319 2355551616 723 D0 /
      DATA WGK ( 22) / 0.0460592382 7100698811 6271735559 374 D0 /
      DATA WGK ( 23) / 0.0471855465 6929915394 5261478181 099 D0 /
      DATA WGK ( 24) / 0.0481858617 5708712914 0779492298 305 D0 /
      DATA WGK ( 25) / 0.0490554345 5502977888 7528165367 238 D0 /
      DATA WGK ( 26) / 0.0497956834 2707420635 7811569379 942 D0 /
      DATA WGK ( 27) / 0.0504059214 0278234684 0893085653 585 D0 /
      DATA WGK ( 28) / 0.0508817958 9874960649 2297473049 805 D0 /
      DATA WGK ( 29) / 0.0512215478 4925877217 0656282604 944 D0 /
      DATA WGK ( 30) / 0.0514261285 3745902593 3862879215 781 D0 /
      DATA WGK ( 31) / 0.0514947294 2945156755 8340433647 099 D0 /
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           DABSC  - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 30-POINT GAUSS RULE
C           RESK   - RESULT OF THE 61-POINT KRONROD RULE
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
C                    OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK61
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(B+A)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
      FC = F(CENTR)
      RESK = WGK(31)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,15
        JTW = J*2
        DABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-DABSC)
        FVAL2 = F(CENTR+DABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J=1,15
        JTWM1 = J*2-1
        DABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-DABSC)
        FVAL2 = F(CENTR+DABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  15    CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(31)*ABS(FC-RESKH)
      DO 20 J=1,30
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQPSRT
      SUBROUTINE DQPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
C***BEGIN PROLOGUE  DQPSRT
C***SUBSIDIARY
C***PURPOSE  This routine maintains the descending ordering in the
C            list of the local error estimated resulting from the
C            interval subdivision process. At each call two error
C            estimates are inserted using the sequential search
C            method, top-down for the largest error estimate and
C            bottom-up for the smallest error estimate.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QPSRT-S, DQPSRT-D)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Ordering routine
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - Integer
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - Integer
C                       Number of error estimates currently in the list
C
C              MAXERR - Integer
C                       MAXERR points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - Double precision
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - Double precision
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - Integer
C                       Vector of dimension LAST, the first K elements
C                       of which contain pointers to the error
C                       estimates, such that
C                       ELIST(IORD(1)),...,  ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - Integer
C                       MAXERR = IORD(NRMAX)
C
C***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQPSRT
C
      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(*),IORD(*)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
C           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
C           SUBDIVISIONS STILL ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
      double precision function d1mach(i)
c
c***begin prologue  d1mach
c***date written   750101    (yymmdd)
c***revision date  801001    (yymmdd)
c***category no.  q
c***keywords  machine constants
c***author  fox, p. a., (bell labs)
c          hall, a. d., (bell labs)
c          schryer, n. l., (bell labs)
c***purpose  returns double precision machine dependent constants
c***description
c     d1mach can be used to obtain machine-dependent parameters
c     for the local machine environment.  it is a function
c     subprogram with one (input) argument, and can be called
c     as forllows, fo example
c
c     d = d1mach(i)
c
c     where i=1,...,5.  the (output) value of d above is
c     determined by the (input) value of i.  the results for
c     various values of i are discussed below.
c
c  double-precision machine constants
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c  d1mach( 5) = log10(b)
c
c***references   fox p.a., hall a.d., schryer n.l., *framework for a
c               portable library*, acm transactions on mathematical
c               software, vol. 4, no. 2, june 1978, pp. 177-188.
c***end prologue  d1mach
c
c     .. scalar arguments ..
      integer i
c     ..
c     .. local scalars ..
      double precision eps,epsneg,xmax,xmin
      integer ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
c     ..
c     .. external subroutines ..
      external machar
c     ..
c     .. intrinsic functions ..
      intrinsic dble,dlog10
c     ..
      call machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,maxexp,
     1            eps,epsneg,xmin,xmax)
c
      if (i.eq.1) d1mach = xmin
      if (i.eq.2) d1mach = xmax
      if (i.eq.3) d1mach = eps
      if (i.eq.4) d1mach = dble(ibeta)*eps
      if (i.eq.5) d1mach = dlog10(dble(ibeta))
c
      return
c
      end
c      algorithm 665, collected algorithms from acm.
c      this work published in transactions on mathematical software,
c      vol. 14, no. 4, p.303.
      subroutine machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,
     1                  maxexp,eps,epsneg,xmin,xmax)
c-----------------------------------------------------------------------
c  this fortran 77 subroutine is intended to determine the parameters
c   of the floating-point arithmetic system specified below.  the
c   determination of the first three uses an extension of an algorithm
c   due to m. malcolm, cacm 15 (1972), pp. 949-951, incorporating some,
c   but not all, of the improvements suggested by m. gentleman and s.
c   marovich, cacm 17 (1974), pp. 276-277.  an earlier version of this
c   program was published in the book software manual for the
c   elementary functions by w. j. cody and w. waite, prentice-hall,
c   englewood cliffs, nj, 1980.
c
c  the program as given here must be modified before compiling.  if
c   a single (double) precision version is desired, change all
c   occurrences of cs (cd) in columns 1 and 2 to blanks.
c
c  parameter values reported are as follows:
c
c       ibeta   - the radix for the floating-point representation
c       it      - the number of base ibeta digits in the floating-point
c                 significand
c       irnd    - 0 if floating-point addition chops
c                 1 if floating-point addition rounds, but not in the
c                   ieee style
c                 2 if floating-point addition rounds in the ieee style
c                 3 if floating-point addition chops, and there is
c                   partial underflow
c                 4 if floating-point addition rounds, but not in the
c                   ieee style, and there is partial underflow
c                 5 if floating-point addition rounds in the ieee style,
c                   and there is partial underflow
c       ngrd    - the number of guard digits for multiplication with
c                 truncating arithmetic.  it is
c                 0 if floating-point arithmetic rounds, or if it
c                   truncates and only  it  base  ibeta digits
c                   participate in the post-normalization shift of the
c                   floating-point significand in multiplication;
c                 1 if floating-point arithmetic truncates and more
c                   than  it  base  ibeta  digits participate in the
c                   post-normalization shift of the floating-point
c                   significand in multiplication.
c       machep  - the largest negative integer such that
c                 1.0+float(ibeta)**machep .ne. 1.0, except that
c                 machep is bounded below by  -(it+3)
c       negeps  - the largest negative integer such that
c                 1.0-float(ibeta)**negeps .ne. 1.0, except that
c                 negeps is bounded below by  -(it+3)
c       iexp    - the number of bits (decimal places if ibeta = 10)
c                 reserved for the representation of the exponent
c                 (including the bias or sign) of a floating-point
c                 number
c       minexp  - the largest in magnitude negative integer such that
c                 float(ibeta)**minexp is positive and normalized
c       maxexp  - the smallest positive power of  beta  that overflows
c       eps     - the smallest positive floating-point number such
c                 that  1.0+eps .ne. 1.0. in particular, if either
c                 ibeta = 2  or  irnd = 0, eps = float(ibeta)**machep.
c                 otherwise,  eps = (float(ibeta)**machep)/2
c       epsneg  - a small positive floating-point number such that
c                 1.0-epsneg .ne. 1.0. in particular, if ibeta = 2
c                 or  irnd = 0, epsneg = float(ibeta)**negeps.
c                 otherwise,  epsneg = (ibeta**negeps)/2.  because
c                 negeps is bounded below by -(it+3), epsneg may not
c                 be the smallest number that can alter 1.0 by
c                 subtraction.
c       xmin    - the smallest non-vanishing normalized floating-point
c                 power of the radix, i.e.,  xmin = float(ibeta)**minexp
c       xmax    - the largest finite floating-point number.  in
c                 particular  xmax = (1.0-epsneg)*float(ibeta)**maxexp
c                 note - on some machines  xmax  will be only the
c                 second, or perhaps third, largest number, being
c                 too small by 1 or 2 units in the last digit of
c                 the significand.
c
c     latest revision - april 20, 1987
c
c     author - w. j. cody
c              argonne national laboratory
c
c-----------------------------------------------------------------------
c     .. scalar arguments ..
      double precision eps,epsneg,xmax,xmin
      integer ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
c     ..
c     .. local scalars ..
      double precision a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,
     1                 y,z,zero
      integer i,itemp,iz,j,k,mx,nxres
c     ..
c     .. intrinsic functions ..
      intrinsic abs,dble,int
c     ..
c     .. statement functions ..
      double precision conv
c     ..
c     .. statement function definitions ..
c-----------------------------------------------------------------------
      conv(i) = dble(i)
c     ..
      one = conv(1)
      two = one + one
      zero = one - one
c-----------------------------------------------------------------------
c  determine ibeta, beta ala malcolm.
c-----------------------------------------------------------------------
      a = one
   10 a = a + a
      temp = a + one
      temp1 = temp - a
      if (temp1-one.eq.zero) go to 10
      b = one
   20 b = b + b
      temp = a + b
      itemp = int(temp-a)
      if (itemp.eq.0) go to 20
      ibeta = itemp
      beta = conv(ibeta)
c-----------------------------------------------------------------------
c  determine it, irnd.
c-----------------------------------------------------------------------
      it = 0
      b = one
   30 it = it + 1
      b = b*beta
      temp = b + one
      temp1 = temp - b
      if (temp1-one.eq.zero) go to 30
      irnd = 0
      betah = beta/two
      temp = a + betah
      if (temp-a.ne.zero) irnd = 1
      tempa = a + beta
      temp = tempa + betah
      if ((irnd.eq.0) .and. (temp-tempa.ne.zero)) irnd = 2
c-----------------------------------------------------------------------
c  determine negep, epsneg.
c-----------------------------------------------------------------------
      negep = it + 3
      betain = one/beta
      a = one
      do 40 i = 1,negep
          a = a*betain
   40 continue
      b = a
   50 temp = one - a
      if (temp-one.ne.zero) go to 60
      a = a*beta
      negep = negep - 1
      go to 50
c
   60 negep = -negep
      epsneg = a
      if ((ibeta.eq.2) .or. (irnd.eq.0)) go to 70
      a = (a* (one+a))/two
      temp = one - a
      if (temp-one.ne.zero) epsneg = a
c-----------------------------------------------------------------------
c  determine machep, eps.
c-----------------------------------------------------------------------
   70 machep = -it - 3
      a = b
   80 temp = one + a
      if (temp-one.ne.zero) go to 90
      a = a*beta
      machep = machep + 1
      go to 80
c
   90 eps = a
      temp = tempa + beta* (one+eps)
      if ((ibeta.eq.2) .or. (irnd.eq.0)) go to 100
      a = (a* (one+a))/two
      temp = one + a
      if (temp-one.ne.zero) eps = a
c-----------------------------------------------------------------------
c  determine ngrd.
c-----------------------------------------------------------------------
  100 ngrd = 0
      temp = one + eps
      if ((irnd.eq.0) .and. (temp*one-one.ne.zero)) ngrd = 1
c-----------------------------------------------------------------------
c  determine iexp, minexp, xmin.
c
c  loop to determine largest i and k = 2**i such that
c         (1/beta) ** (2**(i))
c  does not underflow.
c  exit from loop is signaled by an underflow.
c-----------------------------------------------------------------------
      i = 0
      k = 1
      z = betain
      t = one + eps
      nxres = 0
  110 y = z
      z = y*y
c-----------------------------------------------------------------------
c  check for underflow here.
c-----------------------------------------------------------------------
      a = z*one
      temp = z*t
      if ((a+a.eq.zero) .or. (abs(z).ge.y)) go to 120
      temp1 = temp*betain
      if (temp1*beta.eq.z) go to 120
      i = i + 1
      k = k + k
      go to 110
c
  120 if (ibeta.eq.10) go to 130
      iexp = i + 1
      mx = k + k
      go to 160
c-----------------------------------------------------------------------
c  this segment is for decimal machines only.
c-----------------------------------------------------------------------
  130 iexp = 2
      iz = ibeta
  140 if (k.lt.iz) go to 150
      iz = iz*ibeta
      iexp = iexp + 1
      go to 140
c
  150 mx = iz + iz - 1
c-----------------------------------------------------------------------
c  loop to determine minexp, xmin.
c  exit from loop is signaled by an underflow.
c-----------------------------------------------------------------------
  160 xmin = y
      y = y*betain
c-----------------------------------------------------------------------
c  check for underflow here.
c-----------------------------------------------------------------------
      a = y*one
      temp = y*t
      if (((a+a).eq.zero) .or. (abs(y).ge.xmin)) go to 170
      k = k + 1
      temp1 = temp*betain
      if (temp1*beta.ne.y) go to 160
      nxres = 3
      xmin = y
  170 minexp = -k
c-----------------------------------------------------------------------
c  determine maxexp, xmax.
c-----------------------------------------------------------------------
      if ((mx.gt.k+k-3) .or. (ibeta.eq.10)) go to 180
      mx = mx + mx
      iexp = iexp + 1
  180 maxexp = mx + minexp
c-----------------------------------------------------------------
c  adjust irnd to reflect partial underflow.
c-----------------------------------------------------------------
      irnd = irnd + nxres
c-----------------------------------------------------------------
c  adjust for ieee-style machines.
c-----------------------------------------------------------------
      if ((irnd.eq.2) .or. (irnd.eq.5)) maxexp = maxexp - 2
c-----------------------------------------------------------------
c  adjust for non-ieee machines with partial underflow.
c-----------------------------------------------------------------
      if ((irnd.eq.3) .or. (irnd.eq.4)) maxexp = maxexp - it
c-----------------------------------------------------------------
c  adjust for machines with implicit leading bit in binary
c  significand, and machines with radix point at extreme
c  right of significand.
c-----------------------------------------------------------------
      i = maxexp + minexp
      if ((ibeta.eq.2) .and. (i.eq.0)) maxexp = maxexp - 1
      if (i.gt.20) maxexp = maxexp - 1
      if (a.ne.y) maxexp = maxexp - 2
      xmax = one - epsneg
      if (xmax*one.ne.xmax) xmax = one - beta*epsneg
      xmax = xmax/ (beta*beta*beta*xmin)
      i = maxexp + minexp + 3
      if (i.le.0) go to 200
      do 190 j = 1,i
          if (ibeta.eq.2) xmax = xmax + xmax
          if (ibeta.ne.2) xmax = xmax*beta
  190 continue
  200 return
c---------- last card of machar ----------
      end
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
C***BEGIN PROLOGUE  I1MACH
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   I1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C   where I=1,...,16.  The (output) value of K above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C   Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C   Integers:
C     assume integers are represented in the S-digit, base-A form
C
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C   Floating-Point Numbers:
C     Assume floating-point numbers are represented in the T-digit,
C     base-B form
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base.
C
C   Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C   Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)
C   891012  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
C           (RWC)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added Convex -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
C           options.  (DWL, RWC and WRB).
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        129 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /          7 /
C     DATA IMACH( 2) /          2 /
C     DATA IMACH( 3) /          2 /
C     DATA IMACH( 4) /          2 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -256 /
C     DATA IMACH(13) /        255 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /       -256 /
C     DATA IMACH(16) /        255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /        -50 /
C     DATA IMACH(16) /         76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /     -32754 /
C     DATA IMACH(16) /      32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -4095 /
C     DATA IMACH(13) /       4094 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -4095 /
C     DATA IMACH(16) /       4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /    6LOUTPUT/
C     DATA IMACH( 5) /         60 /
C     DATA IMACH( 6) /         10 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /       -929 /
C     DATA IMACH(13) /       1070 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /       -929 /
C     DATA IMACH(16) /       1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / Z'7FFFFFFF' /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16383 /
C     DATA IMACH(16) /      16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -pd8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         46 /
C     DATA IMACH( 9) / 1777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /         11 /
C     DATA IMACH( 2) /         12 /
C     DATA IMACH( 3) /          8 /
C     DATA IMACH( 4) /         10 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         24 /
C     DATA IMACH( 6) /          3 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         23 /
C     DATA IMACH( 9) /    8388607 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         38 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /         43 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         63 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         39 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         55 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          7 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1015 /
C     DATA IMACH(16) /       1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) /  Z7FFFFFFF /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         54 /
C     DATA IMACH(15) /       -101 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         62 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1021 /
C     DATA IMACH(13) /       1024 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16381 /
C     DATA IMACH(16) /      16384 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          1 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /      -1024 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /          1 /
C     DATA IMACH( 2) /          1 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
C****************************************************************************
C
C		 	EPS09.f
C
C An interface for the scale dependent nuclear modifications
C 		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
C where f_A is the distribution of the parton flavour f for a PROTON in a
C nucleus A, and f_p is the corresponding parton distribution in the 
C free proton.
C  
C When using this interface, please refer to:
C  
C K.J. Eskola, H. Paukkunen and C.A. Salgado,
C "EPS09 - a New Generation of NLO and LO Nuclear Parton Distribution Functions,"
C Published as JHEP04(2009) 065.
C Eprint: arXiv:0902.4154 [hep-ph].
C
C Questions & comments to:
C   hannu.paukkunen@phys.jyu.fi
C   kari.eskola@phys.jyu.fi
C   carlos.salgado@usc.es
C 
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
C
C The error bands for absolute cross-sections and for
C their nuclear ratios should be computed as explained
C in Secs. 2.5 and 4 of arXiv:0902.4154 [hep-ph]. For
C the absolute cross sections, both the errors in the
C free-proton PDFs f_p(x,Q) and the errors in
C the modifications R_f^A(x,Q) should be accounted for.
C For the nuclear ratios, it is sufficient to account only
C for the errors in the modifications R_f^A(x,Q).
C
C *********************************************************
C *********************************************************



      Subroutine EPS09(order, pset, AAA, xxx, QQQ,
     &                   ruv, rdv, ru, rd, rs, rc, rb, rg)

      Implicit none
      Double precision :: ruv, rdv, ru, rd, rs, rc, rb, rg, QQQ, xxx
      Double precision :: LSTEP, x, Q, Q2, allvalues(1:31,1:8,0:50,0:50)
      Double precision :: x_i=0.000001, arg(4), fu(4), res, fg(3)
      Double precision :: result(9), dummy
      Double precision :: realQ, Q2min=1.69, Q2max=1000000.0, Qsteps=50.0
      Double precision :: n_x, zero=0.0

      Character (Len=50) filenimi

      Integer :: xlinsteps=25, xlogsteps=25, startline, lineno
      Integer :: k, p, t, Qpoint, xpoint, pset, iovar
      Integer :: setnumber,j, A, openchannel, order, AAA
      Integer :: psetlast = -10, Alast = -10, orderlast = -10

      save Alast
      save psetlast
      save orderlast
      save allvalues

C *********************************************
C Stop if the set specifications are wrong ones
C *********************************************

      If (order .NE. 1 .and. order .NE. 2) then
      Write(*,*) 'Wrong order!'
      Write(*,*) 'LO : order = 1'
      Write(*,*) 'NLO: order = 2'
      Stop
      End If

      If (pset .LT. 1 .or. pset .GT. 31) then
      Write(*,*) 'Wrong set!'
      Write(*,*) 'Central set: pset = 1'
      Write(*,*) 'Error sets : pset = 2...31'
      Stop
      End If

C ********************************
C Make sure not to change any
C specifications given by the user
C ********************************

      A  = AAA
      x  = xxx
      Q  = QQQ
      Q2 = Q*Q 

C *******************************
C Freeze x if it's < 10E-6 or > 1
C *******************************

      If (x .LT. x_i) Then
      x = x_i
      End If
      If (x .GT. 1) Then
      x = 1.0
      End If

C ************************************
C Freeze Q^2 if it's < 1.69 or > 10E+6
C ************************************

      If (Q2 .LT. Q2min) Then
      Q2 = Q2min
      End If
      If (Q2 .GT. Q2max) Then
      Q2 = Q2max
      End If

C If the set specifications have been changed, read the tables again

      If (A .NE. Alast .or.
     &   order .NE. orderlast) Then

C      Write(*,*) 'Set changed!'

C Read the table

      If (order .EQ. 1) then

        If (A < 10) Then
        Write(filenimi,'("EPS09LOR_", I1)'), A
        Else If (A < 100) Then
        Write(filenimi,'("EPS09LOR_", I2)'), A
        Else If (A < 1000) Then
        Write(filenimi,'("EPS09LOR_", I3)'), A
        End If

      Else

        If (A < 10) Then
        Write(filenimi,'("EPS09NLOR_", I1)'), A
        Else If (A < 100) Then
        Write(filenimi,'("EPS09NLOR_", I2)'), A
        Else If (A < 1000) Then
        Write(filenimi,'("EPS09NLOR_", I3)'), A
        End If

      End If

      Call NextUnit(openchannel)

      OPEN (openchannel, file = filenimi, status='OLD', IOSTAT=iovar)

      If (iovar .NE. 0) Then
      Write(*,*) 'Missing file: ',filenimi
      stop
      End If

      Do setnumber = 1, 31

      Do k = 0,50

      Read(openchannel,*) dummy

      Do t = 0,50

      Read(openchannel,*) (allvalues(setnumber,p,k,t), p=1,8)

      End Do
      End Do

      End Do

      Close(openchannel)

      psetlast  = pset
      Alast     = A
      orderlast = order

      End If

C Find out the position in the loglog Q^2-grid

      realQ  = Qsteps * (log(log(Q2)/log(Q2min)))/
     &                  (log(log(Q2max)/log(Q2min)))
      Qpoint = Aint(realQ)

      If (Qpoint .LE. 0) Then
         Qpoint = 1
      End If
      If (Qpoint .GE. Anint(Qsteps)-1) Then
         Qpoint = Anint(Qsteps)-1
      End If

      LSTEP = (1.0/(xlogsteps)) * LOG(0.1/x_i)

C *********************
C Interpolate the grids 
C *********************

      Do t=1,8

C Find the position in the x-grid

      If (x .LE. 0.1) then
         n_x  = ((1.0/LSTEP) * Log(x/x_i))
       xpoint = Aint(n_x)
      Else
       n_x    = ((x-0.1)*xlinsteps/(1.0-0.1) + xlogsteps)
       xpoint = Aint(n_x)
      End If

      If (xpoint .LE. 0) Then
        xpoint = 1
      End If

      If (t .EQ. 1 .or. t .EQ. 2) Then
         If (xpoint .GE. (xlinsteps+xlogsteps)-4) Then
         xpoint = (xlinsteps+xlogsteps)-4
         End If
      End If

      If (t .EQ. 3 .or. t .EQ. 4 .or. t .EQ. 5 .or.
     &    t .EQ. 6 .or. t .EQ. 7) Then
        If (xpoint .GE. (xlinsteps+xlogsteps)-7) Then
        xpoint = (xlinsteps+xlogsteps)-7
        End If
      End If
	
      If (t .EQ. 8) Then
        If (xpoint .GE. (xlinsteps+xlogsteps)-4) Then
        xpoint = (xlinsteps+xlogsteps)-4
        End If
      End If

      Do k = 1, 4
      If (xpoint-2+k .LT. xlogsteps) Then
      arg(k) = (x_i) * exp(LSTEP * (xpoint-2+k))
      Else
      arg(k) = 0.1 + (xpoint-2+k-xlogsteps) * (1-0.1)/xlinsteps
      End If
      End Do
 
      Do j=1,3

      fu(1) = allvalues(pset,t,Qpoint-2+j,xpoint-1)
      fu(2) = allvalues(pset,t,Qpoint-2+j,xpoint)
      fu(3) = allvalues(pset,t,Qpoint-2+j,xpoint+1)
      fu(4) = allvalues(pset,t,Qpoint-2+j,xpoint+2)
      Call luovi(fu,arg,4,x,res)
      fg(j) = res

C *****************************************
C *****************************************

      End Do

      arg(1) = Qpoint-1
      arg(2) = Qpoint
      arg(3) = Qpoint+1

      Call luovi(fg,arg,3,realQ,res)
  
      result(t) = res

      End Do

      ruv = max(result(1),zero)
      rdv = max(result(2),zero)
      ru  = max(result(3),zero)
      rd  = max(result(4),zero)
      rs  = max(result(5),zero)
      rc  = max(result(6),zero)
      rb  = max(result(7),zero)
      rg  = max(result(8),zero)

200   Continue

      End Subroutine EPS09

! ********************************
! Modified version of Cern Library
! interpolation routine E100
! ********************************

      SUBROUTINE luovi(F,ARG,MMM,Z,SUM)

      Implicit none
      INTEGER :: MMM
      Double precision  :: F(MMM), ARG(MMM), COF(MMM), SUM, Z
      INTEGER :: M, MM, I, J, JNDEX, INDEX

      MM = MIN(MMM, 20)
      M = MM - 1
      DO 1780 I= 1, MM
      COF(I) = F(I)
 1780 Continue
      DO 1800 I= 1, M
      DO 1790 J= I, M
         JNDEX = MM - J
         INDEX = JNDEX + I
         COF(INDEX) = (COF(INDEX)-COF(INDEX-1))/(ARG(INDEX)-ARG(JNDEX))
 1790 CONTINUE
 1800 CONTINUE
      SUM = COF(MM)
      DO 1810 I= 1, M
         INDEX = MM - I
         SUM = (Z-ARG(INDEX))*SUM + COF(INDEX)
 1810 CONTINUE

      End SUBROUTINE luovi

! **********************
! Find the open I/O unit
! **********************

      Subroutine NextUnit(firstopen)

      Logical EX
      Integer firstopen, N

      Do 10 N = 10, 300
C        write(*,*) N
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            firstopen = N
            Goto 20 
        Endif
10    Continue
      Stop ' There is no available I/O unit. '
20    Continue
      End Subroutine NextUnit
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

