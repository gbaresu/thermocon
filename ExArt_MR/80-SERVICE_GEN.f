      subroutine RAYL_GEN(lam,ab_N2,ab_CO2,ab_H2O,SIG)
      !take in input lambda [micron]
      !take in input N2/CO2/H2O abundances [0-1]
      !give out SIGMA in cm2/molecule
       use P_ATM,only: CC_I
       implicit none
       real*8,intent(in) :: lam,ab_N2,ab_CO2,ab_H2O
       real*8,intent(out) :: SIG
       real*8 :: delN2,delCO2,delH2O
       real*8 :: AN2,ACO2,AH2O,BN2,BCO2,BH2O
       real*8 :: sigN2,sigCO2,sigH2O,lam2
       
       AN2 =29.06e-5
       ACO2=43.9e-5
       !AH2O=
       
       BN2 =7.7e-3
       BCO2=6.4e-3
       !BH2O
       
       delN2=1.05
       delCO2=1.15
       lam2=0.532
       
       sigN2=4.577e-21*delN2/lam**4*(AN2*(1+BN2/lam**2))**2
       sigCO2=4.577e-21*delCO2/lam**4*(ACO2*(1+BCO2/lam**2))**2
       !Taken from Sutton 04
       sigH2O=0.75*sigN2

       !sig is in cm2/molecule
       SIG=ab_N2*sigN2+ab_CO2*sigCO2+ab_H2O*sigH2O

       return
      end subroutine
      
      subroutine INIT_PLANET
       !***************************************************************!
       !                                                               ! 
       ! When z is not given, it is calculated from p and t            !
       ! assuming hydrostatic equilibrium                              !
       ! Here we also calculate the gravitational acceleration         !
       ! [GGRAV]=cm/s2                                                 !
       ! and the speed of sound                                        !
       ! [Cs]=cm/s                                                     !
       !                                                               !
       !***************************************************************!
       use P_INP,only: NLYR
       use P_ATM,only: PP_I,TT_I,ZZ_I,CC_I
       use P_CONST,only: KBOL
       use P_PLA,only: GGRAV,CS,MU
       
       implicit none
       real*8 :: HH,muCO2,muN2,abN2
       integer :: i1
       
       call CALC_GGCS(CS,GGRAV)
       HH=CS**2/GGRAV
       
       abN2=1.00-CC_I(1)
       muN2=28
       muCO2=44
       MU=abN2*muN2+CC_I(1)*muCO2

       !UNLOCK ONLY TO CALCULATE Z (km)
        do i1=1,NLYR
         ZZ_I(i1)=HH*log(PP_I(NLYR)/PP_I(i1))
        enddo
       
       return
      end subroutine
      
      subroutine CALC_GGCS(CS,GGRAV)
        use P_INP,only: MPL,RPL,NLYR
        use P_ATM,only: TT_I
        use P_CONST,only: MEARTH,REARTH,MH,KBOL,GRAV
        
        implicit none
        real*8,intent(out) :: CS,GGRAV
        real*8 :: MM,RR
        real*8 :: gam_MARS,mu_MARS
        
        MM=MPL*MEARTH
        RR=RPL*REARTH
        
        gam_MARS=1.31
        mu_MARS=43.2
        
        CS=(gam_MARS*KBOL*TT_I(NLYR)/(mu_MARS*MH))**0.5
        GGRAV=GRAV*MM/RR**2
        
        return
      end subroutine
      
      subroutine spline_linear_val(ndata,tdata,ydata,tval,yval,ypval)
      !
      !*******************************************************************************
      !
      !  Parameters:
      !
      !    Input, integer NDATA, the number of data points defining the spline.
      !
      !    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
      !    and dependent variables at the data points.  The values of TDATA should
      !    be distinct and increasing.
      !
      !    Input, real TVAL, the point at which the spline is to be evaluated.
      !
      !    Output, real YVAL, YPVAL, the value of the spline and its first
      !    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
      !    equal to TDATA(I) for some I.
      !
      implicit none
      !
      integer ndata
      !
      integer :: left,right
      real*8 :: tdata(ndata)
      real*8 :: ydata(ndata)
      real*8 :: tval,ypval,yval
      
      !
      !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
      !  nearest to, TVAL.
      !
      call rvec_bracket ( ndata, tdata, tval, left, right )
      !
      !  Now evaluate the piecewise linear function.
      !
      ypval = (ydata(right)-ydata(left))/(tdata(right)-tdata(left))
      
      yval = ydata(left) +  ( tval - tdata(left) ) * ypval
      
      return
      end

      subroutine rvec_bracket ( n, x, xval, left, right )
      !
      !*******************************************************************************
      !
      !  Parameters:
      !
      !    Input, integer N, length of input array.
      !
      !    Input, real X(N), an array sorted into ascending order.
      !
      !    Input, real XVAL, a value to be bracketed.
      !
      !    Output, integer LEFT, RIGHT, the results of the search.
      !    Either:
      !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
      !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
      !    or
      !      X(LEFT) <= XVAL <= X(RIGHT).
      !
        implicit none
      !
        integer :: n,i,left,right
        real*8 :: x(n)
        real*8 :: xval
      !
        do i = 2, n - 1
      
          if ( xval < x(i) ) then
            left = i - 1
            right = i
            return
          end if
      
         end do
      
        left = n - 1
        right = n
      
        return
      end

      subroutine CHECK_INPUT
       use P_INP
       implicit none
       if(RUN_MODE.lt.0 .or. RUN_MODE.gt.4) then
        call WARNING
        write(*,'(1X,A,I2,A)') "RUN_MODE=",RUN_MODE," NOT ALLOWED, IT "
        write(*,'(1X,A)') "MUST BE CONTAINED WITHIN 1 AND 4"
        stop
       endif
       if(INPUT_SP.lt.0 .or. INPUT_SP.gt.18) then
        call WARNING
        write(*,'(1X,A,I2,A)') "RUN_MODE=",INPUT_SP," NOT ALLOWED, IT "
        write(*,'(1X,A)') "MUST BE CONTAINED WITHIN 1 AND 18"
        stop
       endif
       if(PL_DIST.le.0) then
        call WARNING
        write(*,'(1X,A)') "PLANET DISTANCE MUST BE GREATER THAN 0"
        stop
       endif
       if(INIT_ATM.le.0 .or. INIT_ATM.gt.3) then
        call WARNING
        write(*,'(1X,A)') "INIT_ATM CAN ONLY BE 1, 2 OR 3"
        stop
       endif
       if(TINIT.le.0 .or. TINIT.gt.1000) then
        call WARNING
        write(*,'(1X,A)') "TINIT NEGATIVE OR 0!"
        stop
       endif
       if(SCALE_HEIGHT.le.0) then
        call WARNING
        write(*,'(1X,A)') "SCALE_HEIGHT NEGATIVE OR 0!"
        stop
       endif
       if(LAM_MIN.lt.0.09) then
        call WARNING
        write(*,*) LAM_MIN
        write(*,'(1X,A)') "LAM_MIN MUST BE EQUAL/HIGHER THAN 0.1"
        stop
       endif
       if(LAM_MAX.lt.0 .or. LAM_MAX.le.LAM_MIN) then
        call WARNING
        write(*,'(1X,A)') 
     +   "LAM_MAX MUST BE POSITIVE AND HIGHER THAN LAM_MIN"
        stop
       endif
       if(LAM_RES.lt.1e-5) then
        call WARNING
        write(*,'(1X,A)') "LAM_RES MUST BE HIGHER THAN 1e-5 mic"
        stop
       endif
       if(NLYR.lt.2 .or. NLYR.gt.101) then
        call WARNING
        write(*,'(1X,A)') "NLYR MUST BE CONTAINED WITHIN [2:101]"
        stop
       endif
      end subroutine

      subroutine GET_YV(NLAM,YV)
        use P_INP,only: TST,RST,PL_DIST,ND_VAR,FLUX_SCALE
        use P_CONST,only: CL,HPLANCK,KBOL,PI,AU,RSUN
        use P_ATM,only: NSPMAX
        use P_SPE,only: LAM_S
       
        implicit none
        real*8,dimension(500),intent(out) :: YV
        integer,intent(in) :: NLAM
        integer :: i1,i2,ix,NSP
        integer,parameter :: NMAX=20000
        real*8,dimension(NLAM) :: BAND,nu,band_cm,III
        real*8 :: nu_i,nu_f,dn
        real*8,dimension(NMAX) :: nn,BB
        real*8 :: BASE,HEIGHT,RR,IT

         do i1=1,NLAM-1
          BAND(i1)=LAM_S(i1)
         enddo
         NSP=NLAM-1
         
         nu=0.0
         !Lam in cm
         band_cm=BAND/1e4
         !Get nu [Hz]
         do ix=1,NSP
          nu(ix)=CL/band_cm(ix)
         enddo

         nu_i=CL/(200.0/1e4)
         nu_f=CL/(0.24/1e4)

         BB(1)=0.0
         nn(1)=nu_i
         do i1=2,NMAX
          dn=(nu_f-nu_i)/NMAX
          nn(i1)=nn(i1-1)+dn
          BB(i1)=2.0*HPLANCK*nn(i1)**3/CL**2
          BB(i1)=BB(i1)/(exp(HPLANCK*nn(i1)/(KBOL*TST))-1)
          BB(i1)=BB(i1)*4.0*PI !erg/s/cm2/Hz
          !write(900,*) nn(i1),BB(i1) 
         enddo

         do ix=1,NSP
          III(ix)=0
          do i1=1,NMAX-1
           BASE=nn(i1+1)-nn(i1)
           HEIGHT=0.5*(BB(i1+1)+BB(i1))
           if(nn(i1).gt.nu(ix+1).and.nn(i1+1).lt.nu(ix)) then
            III(ix)=III(ix)+BASE*HEIGHT
           endif
          enddo       
         enddo 
         
         III=PI*III*(RST*RSUN)**2/(4.0*PI*PL_DIST**2*AU**2)/1e7*1e4
         YV(1:NSP)=III(1:NSP)*ND_VAR*FLUX_SCALE
         !do ix=1,NSP
          !write(901,*) BAND(ix),YV(ix)
          !IT=IT+III(ix)
         !enddo
      end subroutine
             
      subroutine GET_WEIGHTS
       use P_SPE,only: WQ
       ! These are the weigths for the k-distribution
       ! 32 quadrature point from 0 to 1
       implicit none
       integer :: i1
       real*8 :: a
       open(unit=10,file="32_mis.inp",status="old")
       do i1=1,32
        read(10,*) WQ(i1)
       enddo
       close(10)
      end subroutine

      subroutine WARNING
       write(*,'(1X,A)') "***** WARNING !!!"
       write(*,*)" CHECK VAR_NAMES.info to get info about each variable"
       write(*,*) 
      end subroutine

      subroutine SPECTRAL_GRID
       use P_ATM
       use P_CON
       implicit none
       integer :: i1
       real :: LCMMIN,LCMMAX,DL
       
       NP=2000
       LCMMAX=1/LLMIN*1e4
       LCMMIN=1/LLMAX*1e4
       
       DL=log10(LCMMAX/LCMMIN)/(NP-1)

       !write(*,*) LCMMIN,LCMMAX,DL
       do i1=NP,1,-1
        LCM(NP-i1+1)=LCMMIN*10**(DL*(NP-i1))
       enddo
      end subroutine 

      subroutine SPECTRAL_GRID2
       use P_ATM
       use P_CON
       implicit none
       integer :: i1
       real :: LCMMIN,LCMMAX,DL
       
       LCMMAX=30.0
       LCMMIN=0.25
       
       DL=0.05
       !log10(LCMMAX/LCMMIN)/(NP-1)
       LCM(1)=LCMMIN
       !write(*,*) LCMMIN,LCMMAX,DL
       do i1=1,9999
        LCM(i1+1)=LCM(i1)+0.05
        if(LCM(i1+1).gt.30) exit
       enddo
       LCM(i1+1)=30.00
       NP=i1+1
      end subroutine 

      subroutine GET_INT_BANDS(LAM,FLUX)
       use P_INP,only: RSTAR,PL_DIST,ND_VAR,FLUX_SCALE
       use P_CONST,only: AU
       use P_SPE,only: LAM_S,NLAM
       use P_ATM,only: NSPMAX
       implicit none
       real*8,dimension(NSPMAX),intent(in) :: LAM,FLUX
       real*8,dimension(NSPMAX) :: LAM_U,FLUX_U
       real*8,dimension(NSPMAX) :: YVARR
       real*8,dimension(8) :: BBCLASS
       integer :: i1,i2,NEFF
       real*8 :: FACT,BASE,HEIGHT,YY,YYP,ITOT,IIIX
       real*8,dimension(7) :: III,III2
       
       FACT=ND_VAR*FLUX_SCALE*(RSTAR**2)/(AU*PL_DIST)**2
       
       BBCLASS=(/0.24,0.4,0.8,1.31,1.86,2.48,3.24,4.5/)

       !open(unit=10,file="OUT_TH/USE_SPEC.dat",status="replace")
       do i1=1,NSPMAX
        if(LAM(i1).gt.0) then
         LAM_U(i1)=LAM(i1)
         FLUX_U(i1)=FLUX(i1)*FACT
         !write(10,*) LAM_U(i1),FLUX_U(i1)
         !write(10,*) LAM(i1),FLUX(i1)
        else
         exit
        endif
       enddo
       !close(10)
       NEFF=i1-1
       
       III=0
       do i1=1,NEFF-1
        if(LAM_U(i1).ge.BBCLASS(1) .and. LAM_U(i1).lt.BBCLASS(2)) then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(1)=III(1)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(2).and.LAM_U(i1).lt.BBCLASS(3))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(2)=III(2)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(3).and.LAM_U(i1).lt.BBCLASS(4))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(3)=III(3)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(4).and.LAM_U(i1).lt.BBCLASS(5))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(4)=III(4)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(5).and.LAM_U(i1).lt.BBCLASS(6))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(5)=III(5)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(6).and.LAM_U(i1).lt.BBCLASS(7))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(6)=III(6)+BASE*HEIGHT
        elseif(LAM_U(i1).ge.BBCLASS(7).and.LAM_U(i1).lt.BBCLASS(8))then
         BASE=LAM_U(i1+1)-LAM_U(i1)
         HEIGHT=0.5*(FLUX_U(i1+1)+FLUX_U(i1))
         III(7)=III(7)+BASE*HEIGHT
        endif
       enddo
       
       !Integrate
       !YVARR=0.0
       !open(unit=10,file="OUTPUT/ADJ_SPEC.dat",status="replace")
       !do i1=1,NLAM-1
       ! call plint(NSPMAX,FLUX_U,LAM_U,NEFF,LAM_S(i1),LAM_S(i1+1),YY)
       ! YVARR(i1)=YY
       ! write(10,*) LAM_S(i1),YVARR(i1)
       !enddo
       !close(10)
       !Integrate2
       YVARR=0.0
       !open(unit=10,file="OUT_TH/ADJ_SPEC.dat",status="replace")
       do i1=1,NLAM-1
        IIIX=0
        do i2=1,NEFF
         if(LAM_U(i2).ge.LAM_S(i1).and.LAM_U(i2).lt.LAM_S(i1+1)) then
          BASE=LAM_U(i2+1)-LAM_U(i2)
          HEIGHT=0.5*(FLUX_U(i2)+FLUX_U(i2+1))
          IIIX=IIIX+BASE*HEIGHT
         elseif(LAM_U(i2).gt.LAM_S(i1+1)) then
          exit
         endif
        enddo
        YVARR(i1)=IIIX
        !write(10,*) LAM_S(i1),YVARR(i1),LAM_S(i1+1)-LAM_S(i1)
       enddo
       !close(10)
       !Check that integrals are similar
       III2=0
       do i1=1,NLAM
        if(LAM_S(i1).ge.BBCLASS(1) .and. LAM_S(i1).lt.BBCLASS(2)) then
         III2(1)=III2(1)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(2).and.LAM_S(i1).lt.BBCLASS(3))then
         III2(2)=III2(2)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(3).and.LAM_S(i1).lt.BBCLASS(4))then
         III2(3)=III2(3)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(4).and.LAM_S(i1).lt.BBCLASS(5))then
         III2(4)=III2(4)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(5).and.LAM_S(i1).lt.BBCLASS(6))then
         III2(5)=III2(5)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(6).and.LAM_S(i1).lt.BBCLASS(7))then
         III2(6)=III2(6)+YVARR(i1)
        elseif(LAM_S(i1).ge.BBCLASS(7).and.LAM_S(i1).lt.BBCLASS(8))then
         III2(7)=III2(7)+YVARR(i1)
        endif
       enddo
       
       do i1=1,7
        ITOT=ITOT+III2(i1)
        !write(*,*) III(i1),III2(i1)
       enddo
       !!MARS SOLAR CONSTANT, uncomment and print
       !write(*,'(5X,A,E13.3,A)') 
     > ! " ---> INTEGRATED FLUX: [W/m2]",ITOT
       return

      end subroutine
      

      subroutine GET_BB(lam,TT,pp)
       use P_CONST,only: KBOL,HPLANCK,CL,PI
       implicit none
       real*8,intent(in) :: lam,TT
       real*8,intent(out) :: pp
       real*8 :: CC,ll,EE,nu
       
       ll=lam/1e4
       CC=2.0*HPLANCK*CL**2/ll**5
       EE=exp((HPLANCK*CL)/(ll*KBOL*TT))
       
       pp=PI*CC/(EE-1.0)
       pp=pp*1e4/1e7/1e4
       ! write(998,*) lam,CC,EE,pp,ll*KBOL*TT
       ! stop
        return
      end subroutine
      
      
      SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )

c        Calculate phase function Legendre expansion coefficients
c        in various special cases


c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c                      4 : Haze L as specified by Garcia/Siewert
c                      5 : Cloud C.1 as specified by Garcia/Siewert

c              GG      Asymmetry factor for Henyey-Greenstein case

c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)

c      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
c                         (be sure to dimension '0:maxval' in calling
c                          program)

c      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results
c                     in Radiative Transfer, Transp. Theory and Stat.
c                     Physics 14, 437-484, Tables 10 And 17
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   IPHAS, NMOM
      REAL      GG
c     ..
c     .. Array Arguments ..

      REAL      PMOM( 0:NMOM )
c     ..
c     .. Local Scalars ..

      INTEGER   K
c     ..
c     .. Local Arrays ..

      REAL      CLDMOM( 299 ), HAZELM( 82 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..

      DATA HAZELM /  2.41260, 3.23047, 3.37296, 3.23150, 2.89350,
     A               2.49594, 2.11361, 1.74812, 1.44692, 1.17714,
     B               0.96643, 0.78237, 0.64114, 0.51966, 0.42563,
     C               0.34688, 0.28351, 0.23317, 0.18963, 0.15788,
     D               0.12739, 0.10762, 0.08597, 0.07381, 0.05828,
     E               0.05089, 0.03971, 0.03524, 0.02720, 0.02451,
     F               0.01874, 0.01711, 0.01298, 0.01198, 0.00904,
     G               0.00841, 0.00634, 0.00592, 0.00446, 0.00418,
     H               0.00316, 0.00296, 0.00225, 0.00210, 0.00160,
     I               0.00150, 0.00115, 0.00107, 0.00082, 0.00077,
     J               0.00059, 0.00055, 0.00043, 0.00040, 0.00031,
     K               0.00029, 0.00023, 0.00021, 0.00017, 0.00015,
     L               0.00012, 0.00011, 0.00009, 0.00008, 0.00006,
     M               0.00006, 0.00005, 0.00004, 0.00004, 0.00003,
     N               0.00003, 3*0.00002, 8*0.00001 /

      DATA  ( CLDMOM(K), K = 1, 159 ) /
     A  2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859,
     B  8.494,  9.286,  9.856, 10.615, 11.229, 11.851, 12.503, 13.058,
     C 13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934,
     D 17.325, 17.673, 17.999, 18.329, 18.588, 18.885, 19.103, 19.345,
     E 19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401,
     F 20.444, 20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310,
     G 20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486, 19.311,
     H 19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659,
     I 17.428, 17.174, 16.931, 16.668, 16.415, 16.144, 15.883, 15.606,
     J 15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378,
     K 13.098, 12.816, 12.536, 12.257, 11.978, 11.703, 11.427, 11.156,
     L 10.884, 10.618, 10.350, 10.090,  9.827,  9.574,  9.318,  9.072,
     M  8.822, 8.584, 8.340, 8.110, 7.874, 7.652, 7.424, 7.211, 6.990,
     N  6.785, 6.573, 6.377, 6.173, 5.986, 5.790, 5.612, 5.424, 5.255,
     O  5.075, 4.915, 4.744, 4.592, 4.429, 4.285, 4.130, 3.994, 3.847,
     P  3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, 2.866, 2.766,
     Q  2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, 1.940,
     R  1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344,
     S  1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909 /
      DATA  ( CLDMOM(K), K = 160, 299 ) /
     T  0.873, 0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610,
     U  0.581, 0.557, 0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401,
     V  0.384, 0.364, 0.349, 0.331, 0.317, 0.301, 0.288, 0.273, 0.262,
     W  0.248, 0.238, 0.225, 0.215, 0.204, 0.195, 0.185, 0.177, 0.167,
     X  0.160, 0.151, 0.145, 0.137, 0.131, 0.124, 0.118, 0.112, 0.107,
     Y  0.101, 0.097, 0.091, 0.087, 0.082, 0.079, 0.074, 0.071, 0.067,
     Z  0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, 0.044, 0.042,
     A  0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026,
     B  0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016,
     C  0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009,
     D  0.009, 3*0.008, 2*0.007, 3*0.006, 4*0.005, 4*0.004, 6*0.003,
     E  9*0.002, 18*0.001 /


      IF ( IPHAS.LT.1 .OR. IPHAS.GT.5 )
     &     CALL ERRMSG( 'GETMOM--bad input variable IPHAS',.TRUE.)

      IF ( IPHAS.EQ.3 .AND. (GG.LE.-1.0 .OR. GG.GE.1.0) )
     &     CALL ERRMSG( 'GETMOM--bad input variable GG',.TRUE.)

      IF ( NMOM.LT.2 )
     &     CALL ERRMSG( 'GETMOM--bad input variable NMOM',.TRUE.)


      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
         PMOM(K) = 0.0
   10 CONTINUE


      IF ( IPHAS.EQ.2 )  THEN
c                                       ** Rayleigh phase function
         PMOM(2) = 0.1

      ELSE IF ( IPHAS.EQ.3 ) THEN
c                                       ** Henyey-Greenstein phase fcn
         DO  20  K = 1, NMOM
            PMOM(K) = GG**K
   20    CONTINUE

      ELSE IF ( IPHAS.EQ.4 ) THEN
c                                        ** Haze-L phase function
         DO  30  K = 1, MIN(82,NMOM)
            PMOM(K) = HAZELM(K) / ( 2*K+1 )
   30    CONTINUE

      ELSE IF ( IPHAS.EQ.5 ) THEN
c                                        ** Cloud C.1 phase function
         DO  40  K = 1, MIN(298,NMOM)
            PMOM(K) = CLDMOM(K) / ( 2*K+1 )
40       CONTINUE

      END IF

      END

    !------------------------------------------------------------------------
      subroutine plint ( ndim, ftab, xtab, ntab, a, b, result )
    !------------------------------------------------------------------------
      !
      !  PLINT approximates the integral of unequally spaced data.
      !
      !
      !  Discussion:
      !
      !    The method uses piecewise linear interpolation.
      !
      !  Reference:
      !
      !    Philip Davis and Philip Rabinowitz,
      !    Methods of Numerical Integration,
      !    Blaisdell Publishing, 1967.
      !
      !  Modified:
      !
      !    30 October 2000
      !
      !  Parameters:
      !
      !    Input, real FTAB(NDIM), the function values, FTAB(I) = F(XTAB(I)).
      !
      !    Input, real XTAB(NDIM), the abscissas at which the
      !    function values are given.  The XTAB's must be distinct
      !    and in ascending order.
      !
      !    Input, integer NTAB, the number of entries in FTAB and
      !    XTAB.  NTAB must be at least 2.
      !
      !    Input, real A, the lower limit of integration.  A should
      !    be, but need not be, near one endpoint of the interval
      !    (X(1), X(NTAB)).
      !
      !    Input, real B, the upper limit of integration.  B should
      !    be, but need not be, near one endpoint of the interval
      !    (X(1), X(NTAB)).
      !
      !    Output, real RESULT, the approximate value of the integral.
      
      implicit none
      integer,intent(in) :: ndim,ntab
      real*8,intent(in),dimension(ndim) :: ftab,xtab 
      
      real*8 :: a,b,fa,fb
      integer i,ihi,ilo,ind
      real*8 :: result
      real*8 :: slope
      real*8 :: syl
      
      if (a==b) then
       result = 0.0E+00
       return
      end if
      
      if ( ntab < 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PLINT - Fatal error!'
       write ( *, '(a,i6)' ) '  NTAB < 2, NTAB = ',ntab
       stop
      end if
      
      do i = 2, ntab
      if ( xtab(i) <= xtab(i-1) ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PLINT - Fatal error!'
       write ( *, '(a)' ) '  Nodes not in strict increasing order.'
       write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
       write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
       write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
       stop
      end if
      end do
      !
      !  If A > B, temporarily switch A and B, and store sign.
      !
      if ( a > b ) then
       syl = b
       b = a
       a = syl
       ind = -1
      else
       syl = a
       ind = 1
      end if
      !
      !  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
      !  with the possible exception that A and B may be in the same
      !  interval, or completely to the right or left of the XTAB's.
      !
      ilo = 1
      ihi = ntab
      do i = 1, ntab
       if ( a <= xtab(i) ) then
        exit
       end if
       ilo = ilo+1
      end do
      
      do i = 1, ntab
       if ( b >= xtab(i) ) then
        exit
       end if
       ihi = ihi-1
      end do
      !
      !  Treat special cases where A, B lie both to left or both to right
      !  of XTAB interval, or inbetween same pair of XTAB's.
      !
      if ( ihi == 0 ) then
       slope = (ftab(2)-ftab(1))/(xtab(2)-xtab(1))
       fa = ftab(1) + slope*(a-xtab(1))
       fb = ftab(1) + slope*(b-xtab(1))
       result = 0.5 * (b-a) * (fa+fb)
       go to 110
      else if ( ilo == ntab+1 ) then
       slope = (ftab(ntab)-ftab(ntab-1))/(xtab(ntab)-xtab(ntab-1))
       fa = ftab(ntab-1)+slope*(a-xtab(ntab-1))
       fb = ftab(ntab-1)+slope*(b-xtab(ntab-1))
       result = 0.5 * (b-a) * (fa+fb)
       go to 110
      else if ( ihi+1 == ilo ) then
       slope = (ftab(ilo)-ftab(ihi))/(xtab(ilo)-xtab(ihi))
       fa = ftab(ihi)+slope*(a-xtab(ihi))
       fb = ftab(ihi)+slope*(b-xtab(ihi))
       result = 0.5 * (b-a) * (fa+fb)
       go to 110
      end if
      !
      !  Carry out approximate integration.  We know that ILO is no greater
      !  than IHI-1, but equality is possible; A and B may be on either side
      !  of a single XTAB(I).  That's OK, then the loop below won't be executed
      !  at all.
      !
      result = 0.0E+00
      do i = ilo, ihi-1
       result = result + 0.5 * (xtab(i+1)-xtab(i))*(ftab(i)+ftab(i+1))
      end do
      !
      !  Add contribution from A-ILO and IHI-B.
      !  Still have to watch out if ILO = 1 or IHI=NTAB...
      !
      if ( ilo == 1 ) then
       slope = (ftab(2)-ftab(1)) / (xtab(2)-xtab(1))
       fa = ftab(1) + slope*(a-xtab(1))
       result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
      else
       slope = (ftab(ilo)-ftab(ilo-1)) / (xtab(ilo)-xtab(ilo-1))
       fa = ftab(ilo-1) + slope*(a-xtab(ilo-1))
       result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
      end if
      
      if ( ihi == ntab ) then
       slope = (ftab(ntab)-ftab(ntab-1)) / (xtab(ntab)-xtab(ntab-1))
       fb = ftab(ntab-1) + slope*(b-xtab(ntab-1))
       result = result + 0.5*(b-xtab(ntab))*(fb+ftab(ntab))
      else
       slope = (ftab(ihi+1)-ftab(ihi)) / (xtab(ihi+1)-xtab(ihi))
       fb = ftab(ihi) + slope*(b-xtab(ihi))
       result = result + 0.5*(b-xtab(ihi))*(fb+ftab(ihi))
      end if
      !
      !  Restore original values of A and B, reverse sign of integral
      !  because of earlier switch.
      !
110   continue
      
      if ( ind /= 1 ) then
       ind = 1
       syl = b
       b = a
       a = syl
       result = -result
      end if
      
      return
      end
      
      
      
