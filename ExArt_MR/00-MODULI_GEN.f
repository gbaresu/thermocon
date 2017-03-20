      module P_INP
       integer :: RUN_MODE,INPUT_SP,INIT_ATM,NLYR,INTERP_KT
       real*8 :: TINIT,SCALE_HEIGHT,LAM_MIN,LAM_MAX,LAM_RES,ZANG
       real*8 :: RPL,MPL,PL_DIST,ND_VAR,FLUX_SCALE,RSTAR,ALBI
       logical :: HYDRO
       character(len=1) :: INTERP_SP,TAU_OUT
       character(len=20) :: MNAME
       real*8 :: SCALE_PRESS,SCALE_CO2,SHIFT_TEMP,H2O_MODE,TST,RST
      end module
      
      module P_ATM
       !Atmospheric and planet/star properties
       integer,parameter :: NLMAX=120
       integer,parameter :: NSPMAX_SM=1000
       integer,parameter :: NSPMAX=100000
       real*8,dimension(NLMAX) :: TT,ZZ,PP
       real*8,dimension(NLMAX) :: TT_I,ZZ_I,PP_I,HH_I,CC_I,O3_I
       real*8 :: LLMIN,LLMAX,LLRES
       real*8,dimension(NLMAX) :: H2O,O3
       real*8,dimension(NSPMAX) :: LAM_EXT,FLUX_EXT
       real*8,dimension(NSPMAX) :: LAM_INTRP,FLUX_INTRP
       real*8 :: CO2,N2,CH4
      end module

      module P_MOD
       integer,parameter :: MMAX=20000
       character(len=3) :: IV 
      end module 
      
      module P_PLA
       real*8 :: CS,GGRAV,MU
      end module
      
      module P_TRA
       use P_ATM,only: NLMAX,NSPMAX,NSPMAX_SM
       integer,parameter :: NTT=9
       real*8,dimension(NLMAX) :: TT_T,ZZ_T,PP_T,CC_T,HH_T
       real*8,dimension(NSPMAX,NTT) :: TAU_LOS 
      
      end module
      
      module P_SPE
       !Final set of P and T to build the spectra
       use P_ATM,only: NLMAX,NSPMAX,NSPMAX_SM
       use P_INP,only: NLYR
       integer :: NLAM,NP
       integer,parameter :: PQ=32
       real*8,dimension(NLMAX,NSPMAX) :: TAU_S,WREAL,FDO,FUP
       real*8,dimension(NLMAX) :: TT_S,ZZ_S,PP_S,CC_S,HH_S,FDOT,FUPT
       real*8,dimension(NSPMAX) :: LAM_S
       real*8,dimension(PQ) :: WQ
       real*8,dimension(NLMAX,PQ) :: FLDNQ,FLDRQ,FLUPQ
       real*8,dimension(NLMAX,PQ,NSPMAX_SM) :: TTQ,SSQ
       real*8,dimension(NLMAX,PQ,NSPMAX_SM) :: TTR,TTC
       real*8,dimension(NLMAX,PQ) :: TTQ2,SSQ2
       real*8,dimension(NLMAX) :: tsca,psca,csca,hsca
       
       !INTERPOLATION
       real*8,dimension(NLMAX) :: PPX,TTX,CCX,HHX
       integer,dimension(NLMAX) :: LABP,LABT,LABC,LABH
       integer,dimension(NLMAX) :: INDXP,INDXT,INDXC,INDXH
      end module
       
      module DISORT_PAR
       INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI,
     &          MXSQT
       PARAMETER ( MXCLY = 102, MXULV = 102, MXCMU = 40, MXUMU = 40,
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
     +         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )
       real :: PR( 0:MAXCMU, MAXCLY ),PC( 0:MAXCMU, MAXCLY )
       real :: FLDN(MAXULV),FLDR(MAXULV),RFLDIR(MAXULV),RFLDN(MAXULV),
     +         FLUPD( MAXULV ),FLUP( MAXULV ),FLDIR( MAXULV),
     +         DFDT( MAXULV ), UAVG( MAXULV ), U0U( MAXUMU, MAXULV ),
     +         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     +         TRNMED( MAXUMU )

      end module
      
      module P_CON
       !VARIOUS PARAMETERS
       integer,parameter :: NLMAX=200,NSMAX=9999
       integer :: H2Osw,O3sw,NP
       real,dimension(NSMAX) :: FLAM,LCM
       real,dimension(NSMAX,NLMAX) :: TAU,ALB
       real :: INTSP,LMIN,LMAX
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

