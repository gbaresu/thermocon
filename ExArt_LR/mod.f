
      module DISORT_PAR
       INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI,
     &          MXSQT
       PARAMETER ( MXCLY = 104, MXULV = 105, MXCMU = 40, MXUMU = 40,
     &          MXPHI = 40, MI = MXCMU / 2, MI9M2 = 9*MI - 2,
     &          NNLYRI = MXCMU*MXCLY, MXSQT = 1000 )
      end module 

      module P_DIS
       !DISORT properties
       use DISORT_PAR
       integer :: NTAU,NUMU,NPHI,IBCND,NSTR
       character(len=127) :: HEADER
       character(len=3) :: BLANKS
       integer :: MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU
       PARAMETER ( MAXCLY = MXCLY, MAXCMU = MXCMU, MAXPHI = MXPHI, 
     +  MAXULV = MXULV,MAXUMU = MXUMU )
       logical :: USRTAU,USRANG,LAMBER,PLANK,ONLYFL,DELTAM,PRNT(7)
       logical :: AZMAVG
       real :: ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     +         HL( 0:MAXCMU ), PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ),
     +         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     +         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXCLY )

       real :: FLDN(MAXULV),FLDR(MAXULV),RFLDIR(MAXULV),RFLDN(MAXULV),
     +         FLUPD( MAXULV ),FLUP( MAXULV ),FLDIR( MAXULV),
     +         DFDT( MAXULV ), UAVG( MAXULV ), U0U( MAXUMU, MAXULV ),
     +         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     +         TRNMED( MAXUMU )

      end module

      module P_SPE
       !Final set of P and T to build the spectra
       integer,parameter :: NSPMAX=10000
       integer,parameter :: NLMAX=200
       integer,parameter :: PQ=32
       integer :: NLYR,WM
       integer,dimension(15) :: BAND
       real*8,dimension(15) :: BANDV
       real*8,dimension(15,PQ) :: yyPT
       integer,dimension(NLMAX) :: LABP,LABT,INDXP,INDXT
       integer,dimension(NLMAX) :: LABC,LABH,INDXC,INDXH
       real*8,dimension(NLMAX,NSPMAX,PQ) :: TTQ,SSQ,RRQ
       real*8,dimension(NLMAX,PQ) :: FLDNQ,FLDRQ,FLUPQ
       real*8,dimension(NLMAX) :: PPP,TTT,ZZZ,ZZX,CCC,HHH,NNN
       real*8,dimension(NLMAX) :: tsca,psca,csca,hsca      
       real*8,dimension(NLMAX) :: tint,pint,cint,hint      
       !INTERPOLATION
       real*8,dimension(NLMAX) :: PPX,TTX,CCX,HHX
       integer :: NBB,NWW
      end module
       
      module P_INP
       !Stellar parameters
       integer :: MODE                                         !Running Mode
       integer :: WRITE_OUT                                    !Write Output
       integer :: INPUT_SP                                     !Input Spectrum
       real*8 :: ZANG                                     !Zenith Angle
       real*8 :: TST                                     !Stellar T
       real*8 :: RST                                     !Radius T
      end module
       
      module P_ST
       !Stellar parameters
       real*8 :: YV(14)                                         !Fluxes [W/m2] in the 14 bands
      end module

      module P_PL
       !Planet parameters
       real*8 :: gg                                             !gravity acceleration
       real*8 :: alb                                            !bond albedo (between 0.24 and 4.5 micron)
       real*8 :: MU                                             !molecular mean weight
       real*8 :: PL_DIST                                        !planet distance in AU
       real*8 :: FL_SCALE                                       !scale flux
       real*8 :: ND_VAR                                         !day/night variation
       real*8 :: PL_RAD                                         !planet mass in earth masses
       real*8 :: PL_MAS                                         !planet radius in earth radii
       real*8 :: ALBI                                           !planet bond albedo in the visible
       integer :: LAMBI                                          !lambertian reflection
      end module
      
      module P_CONST
       real*8 :: PI=3.14159
       real*8 :: RSUN=6.955e10
       real*8 :: AU=1.496e13
       real*8 :: KBOL=1.38e-16
       real*8 :: MH=1.660538e-24
       real*8 :: MEARTH=5.972e27
       real*8 :: REARTH=6.371e8
       real*8 :: GRAV=6.67259e-8
       real*8 :: HPLANCK=6.62e-27
       real*8 :: CL=2.99e10
      end module

