      module PAR_INP
       integer,parameter :: NLMAX=200
       
       real*8 :: MPL_INP                    ![MEarth] Planet mass
       real*8 :: RPL_INP                    ![REarth] Planet radius
       real*8 :: DPL_INP                    ![AU] Planet distance from star
       real*8 :: TST_INP                    ![K] Stellar temperature input
       real*8 :: RST_INP                    ![Rsun] Stellar radius
       real*8 :: ALB_INP                    ![] Bond Albedo
       real*8 :: SF_INP                     ![] Stellar Flux scale factor
       real*8 :: SZA_INP                    ![DEG] Solar zenith angle
       real*8 :: CLRMOD_INP                 ![] CLR switch 
       real*8 :: CLR_INP                    ![K Km-1] Critical lapse rate
	   
       real*8 :: I_TIME_INP	                ![hours] Integration time
       real*8 :: C_CRI_INP	                ![] Convergence criterion [(Tnew-Told) per day all layers]
       integer :: M_INP 	                ![] Use moderate or low res model 
       integer :: I_MAX_INP 	            ![] Maximum number of iteration
       integer :: CH_CON_INP		        ![] (1) use convective correction (0) radiative only

       real*8 :: PFACT_INP                  ![] Scale PRE at each level
       real*8 :: CO2FACT_INP                ![] Scale CO2 at each level
       real*8 :: TSHIFT_INP                 ![K] Shifts TEMP at each level
       
       real*8,dimension(NLMAX) :: PRE_INP   ![mbar] Pressure values  
       real*8,dimension(NLMAX) :: TEMP_INP  ![K] Temperature values  
       real*8,dimension(NLMAX) :: CO2_INP   ![ratio] CO2 concentration 
       real*8,dimension(NLMAX) :: H2O_INP   ![ratio] Water concentration
       
       integer :: WAT_MODE_INP              ![] Water switch
       real*8 :: WAT_HUM_INP                ![perc] Water humidity at the ground
       
      end module PAR_INP

      module PAR_PLA
      
       real*8 :: MPL                        ![g] Planet mass
       real*8 :: RPL                        ![cm] Planet radius
       real*8 :: DPL                        ![cm] Planet distance from star
       real*8 :: GACC                       ![cm s-2] Gravitational acceleration
       real*8 :: Cp                         ![erg g-1 K-1] Specific heat constant pressure
       real*8 :: Cv                         ![erg g-1 K-1] Specific heat constant volume
       real*8 :: GAMMA                      ![] Gamma
       real*8 :: MU                         ![] Mean molecular mass
       real*8 :: CS                         ![cm s-1] Sound speed
       real*8 :: ALB                        ![] Bond Albedo
       real*8 :: SZA                        ![rad] Solar zenit angle
       real*8 :: CSZA                       ![] Cosine of the solar zenit angle
       real*8 :: NDVAR                      ![] Flux factor, night day var
       real*8 :: CLR                        ![K km-1] Critical Lapse Rate
       real*8 :: CLR_DRY                    ![K km-1] Critical Lapse Rate Wet
       real*8 :: CLR_WET                    ![K km-1] Critical Lapse Rate Dry
       
      end module PAR_PLA
      
      module PAR_STA
      
       real*8 :: RST                        ![cm] Star radius
       real*8 :: TST                        ![K] Star temperature
       
      end module PAR_STA
      
      module PAR_ATM
       use PAR_INP,only: NLMAX
       integer :: NL                         ![] Number of atmospheric layers
       real*8,dimension(NLMAX) :: CO2         ![] CO2 concentration
       real*8,dimension(NLMAX) :: N2          ![] N2 concentration
       real*8,dimension(NLMAX) :: H2O         ![] H2O concentration
       real*8,dimension(NLMAX) :: ZZZ         ![cm] Height
       real*8,dimension(NLMAX) :: PRE_USE1    ![barye] Pressure
       real*8,dimension(NLMAX) :: Tinit       ![K] Temperature
       real*8,dimension(NLMAX) :: Ttau        ![K] Temperature
       real*8,dimension(NLMAX) :: T0          ![K] Temperature 
       real*8,dimension(NLMAX) :: T1          ![K] Temperature
       real*8,dimension(NLMAX) :: T2          ![K] Temperature
       real*8 :: HSCA                       ![cm] Scale Height
       real*8 :: WAT_HUM                    ![] Water humidity at the ground
       
       real*8,dimension(NLMAX) :: FL_EX_TOT_DWN              ![W/m2] Flux from exart tot down
       real*8,dimension(NLMAX) :: FL_EX_DIR_DWN              ![W/m2] Flux from exart direct down
       real*8,dimension(NLMAX) :: FL_EX_DIF_UPP              ![W/m2] Flux from exart diff up
       real*8,dimension(NLMAX) :: FL_EX_DWN_MAN              ![W/m2] Net VIS Flux Plus IR DWN Flux
       real*8,dimension(NLMAX) :: dTdt      ![erg s-1] Heating function
       real*8 :: TGR_EQ                     ![] Ground Temperature imposing rad eq
       real*8 :: TGR_MAN                    ![] Ground Temperature imposing manabe initial condition
       
      end module PAR_ATM

      module PAR_CODE
       use PAR_INP,only: NLMAX
       
       integer :: N_IT                       ![sec] Iteration number
       real*8 :: I_TIME                     ![sec] Integration time
       real*8 :: C_CRI                      ![] Convergence criterion
       integer :: I_MAX                     ![] Maximum number of iterations
       integer :: CH_CON                    ![] Convective correction switch
       integer :: WAT_MODE                  ![] Water switch
       integer :: COUNT_CL                  ![] Count Convection Loops
       integer :: CHECK_CC                  ![] CHECK Convection Loops
       character,dimension(NLMAX) :: IS_CRIT![] CHECK Convection Loops
       
      end module PAR_CODE
      
      module CONST
       !Constants
       real*8 :: KBOL=1.380658e-16          ![erg K-1]
       real*8 :: SIGBOL=5.67051e-5          ![erg cm-2 s-1 K-4]
       real*8 :: PI=3.14159265              ![] Pi
       real*8 :: AU=1.49597870e13           ![cm] Astronomical Unit
       real*8 :: MEARTH=5.972e27            ![g] Earth mass
       real*8 :: REARTH=6.371e8             ![cm] Earth radius
       real*8 :: GGRAV=6.67259e-8           ![cm3 g-1 s-2] Grav constant
       real*8 :: CvN2=743e4                 ![erg g-1 K-1] N2 specific heat
       real*8 :: CvCO2=655e4                ![erg g-1 K-1] CO2 specific heat 
       real*8 :: mH=1.6733e-24              ![g] hydrogen mass
       real*8 :: mN2=28                     ![uma] molecular nitrogen mass
       real*8 :: mCO2=44                    ![uma] CO2 mass
       real*8 :: Rw=461.5e4                 ![erg K-1 g-1]
       real*8 :: Rd=287.0e4                 ![erg K-1 g-1]
       real*8 :: Lvap=2.5e10                ![erg g-1]
       end module CONST
