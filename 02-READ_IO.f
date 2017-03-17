      subroutine READ_INPUT
       use PAR_INP
       use PAR_ATM,only: WAT_HUM
       use PAR_PLA,only: MPL,RPL,DPL,ALB,SZA,CSZA,NDVAR
       use PAR_STA,only: TST,RST
       use PAR_CODE,only: I_MAX,CH_CON,I_TIME,C_CRI,WAT_MODE
       use CONST,only: REARTH,MEARTH,PI,AU
       implicit none
       
       open(unit=10,file="INP_TH/input_th.inp",status="old")
        read(10,*) 
        read(10,*) 
        read(10,*) MPL_INP
        read(10,*) RPL_INP
        read(10,*) DPL_INP
        read(10,*) TST_INP
        read(10,*) RST_INP
        read(10,*) ALB_INP
        read(10,*) SF_INP
        read(10,*) SZA_INP
        read(10,*) CLRMOD_INP
        read(10,*) CLR_INP
        read(10,*) WAT_MODE_INP
        read(10,*) WAT_HUM_INP
        read(10,*)
        read(10,*) PFACT_INP
        read(10,*) CO2FACT_INP
        read(10,*) TSHIFT_INP
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) M_INP
        read(10,*) I_TIME_INP
        read(10,*) C_CRI_INP
        read(10,*) I_MAX_INP
        read(10,*) CH_CON_INP
       close(10)
       
       MPL=MPL_INP*MEARTH
       RPL=RPL_INP*REARTH
       DPL=DPL_INP*AU
       TST=TST_INP
       RST=RST_INP
       ALB=ALB_INP
       NDVAR=SF_INP
       SZA=SZA_INP!/180.0*PI
       CSZA=cos(SZA)
       I_MAX=I_MAX_INP
       C_CRI=C_CRI_INP
       I_TIME=I_TIME_INP*3600.0
       CH_CON=CH_CON_INP
       WAT_MODE=WAT_MODE_INP
       WAT_HUM=WAT_HUM_INP/100.0
       
      end subroutine
  
      subroutine READ_INIT_ATM
       use PAR_INP,only: PRE_INP,TEMP_INP,CO2_INP,H2O_INP,PFACT_INP,
     >  CO2FACT_INP,TSHIFT_INP
       use PAR_ATM,only: NL,N2,CO2,H2O,PRE_USE1,Tinit
       implicit none
       integer :: i1
       
       i1=1
       open(unit=10,file="INP_TH/atms_th.inp",status="old")
        do 
         read(10,*,END=1000) PRE_INP(i1),TEMP_INP(i1),
     >    CO2_INP(i1),H2O_INP(i1) 
         i1=i1+1
        enddo
1000   close(10)
       NL=i1-1
       !From mbar to cgs
       PRE_USE1=PRE_INP*PFACT_INP
       Tinit=TEMP_INP+TSHIFT_INP
       CO2=CO2_INP*CO2FACT_INP
       H2O=H2O_INP
       N2=1.0-CO2
       ! TODO add shift and scale
      end subroutine

      subroutine WRITE_ATM_EXART
       use PAR_ATM,only: NL,PRE_USE1,Tinit,Ttau,CO2,H2O
       implicit none
       integer :: i1
       
       open(unit=10,file="INP_TH/atms_ex.inp",status="replace")
       do i1=1,NL
        write(10,1000) PRE_USE1(i1),Tinit(i1),CO2(i1),H2O(i1)
       enddo
       close(10)
1000   format(4(E13.6))
       Ttau=Tinit
      end subroutine

      subroutine WRITE_ATM_EXART_LOOP
       use PAR_ATM,only: NL,PRE_USE1,Ttau,T2,CO2,H2O
       use PAR_CODE,only: N_IT
       implicit none
       character(len=100) :: fname
       integer :: i1
       
       write(fname,'(A,I4.4,A)') "OUT_EX/ATM",N_IT,".out"
       open(unit=11,file=trim(fname),status="replace")
       open(unit=10,file="INP_TH/atms_ex.inp",status="replace")
       do i1=1,NL
        write(10,1000) PRE_USE1(i1),T2(i1),CO2(i1),H2O(i1)
        write(11,1000) PRE_USE1(i1),T2(i1),CO2(i1),H2O(i1)
       enddo
       close(10)
       close(11)
1000   format(4(E13.6))
       Ttau=T2
      end subroutine

      subroutine WRITE_INP_EXART
       use PAR_ATM,only: NL
       use PAR_STA,only: TST,RST
       use PAR_PLA,only: DPL,RPL,MPL,NDVAR,CSZA,ALB,SZA
       use CONST,only: AU,MEARTH,REARTH
       implicit none
       integer :: i1
       
       open(unit=10,file="INP_TH/input_ex.inp",status="replace")
        write(10,*) "17            : INPUT SP" 
        write(10,*) TST,"         : STAR TEMP" 
        write(10,*) RST,"         : STAR RAD" 
        write(10,*) DPL/AU,"         : PL DIST" 
        write(10,*) NDVAR,"         : ND VAR" 
        write(10,*) "1.00         : SOL FACT" 
        write(10,*) SZA,"         : ZANG" 
        write(10,*) 
        write(10,*) RPL/REARTH,"     : RPL"
        write(10,*) MPL/MEARTH,"     : MPL"
        write(10,*) 
        write(10,*) 0.24,"     LMIN"
        write(10,*) 100.5,"     LMAX"
        write(10,*) 
        write(10,*) "N       : INTERP"
        write(10,*) "2       : INTERP K"
        write(10,*) "N       : TAU OUT"
        write(10,*) 
        write(10,*) NL,"    : NLAYER"
        write(10,*) ALB,"    : ALBEDO BOND VISIBLE"
       close(10)
       
      end subroutine
      
      subroutine READ_EXART
       use PAR_ATM,only: NL,FL_EX_DIF_UPP,FL_EX_DIR_DWN,FL_EX_TOT_DWN,
     +  FL_EX_DWN_MAN
       use PAR_INP,only: M_INP
       implicit none
       integer :: i1
       real*8 :: PX
       
       if(M_INP.eq.1) then
        open(unit=10,file="OUT_TH/FL_EXART_MR.dat",status="old")
       elseif(M_INP.eq.2) then
        open(unit=10,file="OUT_TH/FL_EXART_LR.dat",status="old")       
       endif
       do i1=1,NL
        read(10,*) PX,FL_EX_TOT_DWN(i1),FL_EX_DIR_DWN(i1),
     >   FL_EX_DIF_UPP(i1),FL_EX_DWN_MAN(i1)
       enddo
       close(10)
       !Fromj J/s/m2 to erg/s/cm2
       FL_EX_TOT_DWN=FL_EX_TOT_DWN*1e7*1e-4
       FL_EX_DIR_DWN=FL_EX_DIR_DWN*1e7*1e-4
       FL_EX_DIF_UPP=FL_EX_DIF_UPP*1e7*1e-4
       FL_EX_DWN_MAN=FL_EX_DWN_MAN*1e7*1e-4
      end subroutine
