      program THERMOCON
       use PAR_CODE,only: I_MAX,N_IT,CH_CON,COUNT_CL,CHECK_CC,IS_CRIT
       implicit none
       integer :: i1
       real*8 :: TX,cp
       
       call system("cat /dev/null > tground.txt")

       !Read input file and input atmosphere
       call READ_INPUT
       call READ_INIT_ATM
       !Calculate various atmospheres properties
       call PREP_ATM
       !Write out first atmosphere and input file to ExArt
       call WRITE_INP_EXART
       call WRITE_ATM_EXART
       
       do i1=1,I_MAX
        write(*,*) i1
        N_IT=i1
        call RUN_EXART
        call READ_EXART

        call CALC_HEAT
        call TIME_STEP
        
        call INIT_MANABE
        
        IS_CRIT=""
        if(CH_CON.eq.1) then
         !call CHECK_CONV(1)
         call CHECK_CONV_OLD(1)
         COUNT_CL=1
          write(*,*) "ITERATION OUT",COUNT_CL," ",CHECK_CC
         do while(COUNT_CL.lt.99999 .and. CHECK_CC.eq.0)
          COUNT_CL=COUNT_CL+1
          !call CHECK_CONV(2)
          call CHECK_CONV_OLD(2)
         enddo
         write(*,*) "ITERATION IN",COUNT_CL
        else
         call SKIP_CONV
        endif
        !RESET TEMPERATURE and WATER values for next iteration
        call PRINT_OUT
        call SET_LOOP
        call PREP_ATM_LOOP
        call WRITE_ATM_EXART_LOOP
       enddo
      end
      
      subroutine PRINT_OUT
       use PAR_ATM,only: T2,PRE_USE1,NL,ZZZ,TGR_MAN,TGR_EQ,Ttau
       use PAR_PLA,only: CLR
       use PAR_CODE,only: IS_CRIT,COUNT_CL,N_IT
       implicit none
       integer :: i1
       real*8 :: DZ,DTDZ,CC
       real*8,dimension(NL) :: DTCHECK
       do i1=1,NL
        DTCHECK(i1)=abs(T2(i1)-Ttau(i1))
        DZ=(ZZZ(i1-1)-ZZZ(i1))/1e5
        DTDZ=(T2(i1)-T2(i1-1))/DZ
        write(*,*) ZZZ(i1)/1e5,PRE_USE1(i1)/1e3,
     >   T2(i1),DTDZ,DTCHECK(i1),IS_CRIT(i1)
       enddo
       CC=0
       do i1=1,NL
        if(DTCHECK(i1).lt.1e-3) CC=CC+1
       enddo
       !if((CC-NL)/NL.lt.0.05) stop
       write(*,*) 
       write(*,'(A,F10.4,A,I5.5)') "CLR: ",CLR," COUNT: ",COUNT_CL
       write(*,'(2(A,F10.4))') "TGR_MAN: ",TGR_MAN," TGR_EQ: ",TGR_EQ
       write(*,*) TGR_EQ-TGR_MAN
       open(unit=10,file="tground.txt",access="append",status="old")
        write(10,*) N_IT,T2(NL),DTCHECK(NL) 
       close(10)
       !write(*,'(A,I3.3)') "CONVERGED LAYERS: ",CC
      end subroutine
      
      subroutine SKIP_CONV
       use PAR_ATM,only: T0,T1,T2      
       T2=T1
      end subroutine
      
      subroutine CHECK_CONV(XX)
       use PAR_ATM,only: T0,T1,T2,PRE_USE1,NL,ZZZ
       use PAR_PLA,only: Cp,GACC,CLR
       use PAR_CODE,only: I_TIME,CHECK_CC,IS_CRIT
       use PAR_INP,only: NLMAX
       use CONST,only: SIGBOL
       implicit none
       integer,intent(in) :: XX
       integer :: i1,kk
       real*8 :: LRG,AA,BB,CC,DD,EE,DPN,DPNM,C1,C2,C3,TIN,TT,DZ,DPK,DPKM
       real*8,dimension(NLMAX) :: LR_LAY,TW2,TW1,TW0
       !!
       ! T0 temperature after adding heating
       ! T1 temperature after computing manabe initial condition
       ! TW working temperature in the convective loop
       ! T2 temperature set corrected for convection
       !!
       C1=SIGBOL*I_TIME
       C2=Cp*DPN/GACC
       C3=Cp*DPNM/GACC
       
       TW0=T0
       TW1=T1
       TW2=0.0
       TW2(1)=TW1(1)
       do kk=1,9999
       
       !do i1=1,NL
       ! write(*,*) "a",TW0(i1),TW1(i1),TW2(i1)
       !enddo
       
       CHECK_CC=0
       !Check ground layer first
       DZ=(ZZZ(NL-1)-ZZZ(NL))/1e5 !In Km
       LR_LAY(NL-1)=(TW1(NL)-TW0(NL-1))/DZ
       if(LR_LAY(NL-1).gt.CLR) then
        DPN=PRE_USE1(NL)-PRE_USE1(NL-1)
        DPNM=PRE_USE1(NL-1)-PRE_USE1(NL-2)
        
        AA=SIGBOL*I_TIME
        BB=Cp/GACC*(DPN+DPNM)
        CC=0.0
        DD=0.0
        EE=C1*TW1(NL)**4+C2*TW1(NL)+
     >   C3*CLR*DZ+C3*TW0(NL-1)
     
        TIN=TW0(NL)
        
        call solve_quartic(AA,BB,CC,DD,EE,TIN,TT)
        TW2(NL)=TT
        TW1(NL-1)=TW2(NL)-CLR*DZ
        IS_CRIT(NL)="*"
       else 
        TW2(NL)=TW1(NL)
        TW1(NL-1)=TW0(NL-1)
        CHECK_CC=1
       endif
       !Restore first layer lapse rate
       LR_LAY(NL-1)=(TW2(NL)-TW1(NL-1))/((ZZZ(NL-1)-ZZZ(NL))/1e5)
       
       !Recalculate other lapse rates (only NL-2 is "new"
       do i1=1,NL-2
        DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5
        LR_LAY(i1)=(TW1(i1+1)-TW0(i1))/DZ
       enddo

       do i1=NL-2,1,-1
        if((LR_LAY(i1)-CLR).gt.1e-3) then
         CHECK_CC=0
         DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5
         DPK=PRE_USE1(i1+1)-PRE_USE1(i1)
         DPKM=PRE_USE1(i1)-PRE_USE1(i1-1)
         
         TW2(i1+1)=1/(DPK+DPKM)*(DPK*TW1(i1+1)+DPKM*TW0(i1)+
     >     DPKM*CLR*DZ)
         TW1(i1)=TW2(i1+1)-CLR*DZ

         LR_LAY(i1)=(TW2(i1+1)-TW1(i1))/DZ
         IS_CRIT(i1+1)="*"
        else
         TW2(i1+1)=TW1(i1+1)
         TW1(i1)=TW0(i1)
         CHECK_CC=1
        endif
       enddo
       !do i1=1,NL
       ! write(*,*) kk,LR_LAY(i1),TW2(i1),TW1(i1),TW0(i1)
       !enddo

       do i1=1,NL-1
        DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5
        LR_LAY(i1)=(TW2(i1+1)-TW1(i1))/DZ
       enddo
       do i1=1,NL
        DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5        
        !write(*,*) kk,DZ,TW2(i1),TW1(i1),TW0(i1),LR_LAY(i1)
       enddo
       TW0=TW2
       enddo
       !write(*,*)
       T2=TW2

      end subroutine
      
      subroutine CHECK_CONV_OLD(XX)
       use PAR_ATM,only: T0,T1,T2,PRE_USE1,NL,ZZZ
       use PAR_PLA,only: Cp,GACC,CLR
       use PAR_CODE,only: I_TIME,CHECK_CC,IS_CRIT
       use PAR_INP,only: NLMAX
       use CONST,only: SIGBOL
       implicit none
       integer,intent(in) :: XX
       integer :: i1,kk,COUNT_CC
       real*8 :: LRG,AA,BB,CC,DD,EE,DPN,DPNM,C1,C2,C3,TIN,TT,DZ,DPK,DPKM
       real*8,dimension(NLMAX) :: LR_LAY,TW2,TW1,TW0
       !!
       ! T0 temperature after adding heating
       ! T1 temperature after computing manabe initial condition
       ! TW working temperature in the convective loop
       ! T2 temperature set corrected for convection
       !!
       
       TW0=T0
       TW1=T1
       TW2=0.0
       TW2(1)=TW1(1)
       !write(*,*)
       !do i1=NL-10,NL
       ! write(*,*) T2(i1),T1(i1)
       !enddo       
       !do i1=1,NL
       ! write(*,*) "a",TW0(i1),TW1(i1),TW2(i1)
       !enddo
       
       !Check ground layer first
       DZ=(ZZZ(NL-1)-ZZZ(NL))/1e5 !In Km
       LR_LAY(NL-1)=(TW1(NL)-TW0(NL-1))/DZ

       if(LR_LAY(NL-1).gt.CLR) then
        !write(*,*) "SUPER",LR_LAY(NL-1),CLR
        CHECK_CC=0
        DPN=PRE_USE1(NL)-PRE_USE1(NL-1)
        DPNM=PRE_USE1(NL-1)-PRE_USE1(NL-2)

        C1=SIGBOL*I_TIME
        C2=Cp*DPN/GACC
        C3=Cp*DPNM/GACC
        
        AA=SIGBOL*I_TIME
        BB=0.0
        CC=0.0
        DD=Cp/GACC*(DPN+DPNM)
        EE=C1*TW1(NL)**4+C2*TW1(NL)+C3*CLR*DZ+C3*TW0(NL-1)
        TIN=TW0(NL)
        
        call solve_quartic(AA,BB,CC,DD,EE,TIN,TT)
        TW2(NL)=TT
        TW1(NL-1)=TW2(NL)-CLR*DZ
        CHECK_CC=1
        IS_CRIT(NL)="*"
        !Restore
        LR_LAY(NL-1)=(TW2(NL)-TW1(NL-1))/((ZZZ(NL-1)-ZZZ(NL))/1e5)
       else 
        TW2(NL)=TW1(NL)
        CHECK_CC=1
       endif
       !Restore first layer lapse rate

       !do i1=NL-10,NL
       ! write(*,*) T2(i1),T1(i1),LR_LAY(i1)
       !enddo
       
       do i1=NL-2,1,-1
        DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5
        LR_LAY(i1)=(TW1(i1+1)-TW0(i1))/DZ
        if((LR_LAY(i1)-CLR).gt.1e-3) then
         CHECK_CC=0
         DPK=PRE_USE1(i1+1)-PRE_USE1(i1)
         DPKM=PRE_USE1(i1)-PRE_USE1(i1-1)
         
         TW2(i1+1)=1/(DPK+DPKM)*(DPK*TW1(i1+1)+DPKM*TW0(i1)+
     >     DPKM*CLR*DZ)
         TW1(i1)=TW2(i1+1)-CLR*DZ
         !Restore
         LR_LAY(i1)=(TW2(i1+1)-TW1(i1))/DZ
         IS_CRIT(i1+1)="*"
        else
         TW2(i1+1)=TW1(i1+1)
         CHECK_CC=1
        endif
       enddo
       !Store results
       T0=TW1
       T1=TW2
       T2=TW2
       !Recheck
       COUNT_CC=0
       do i1=1,NL-1
        DZ=(ZZZ(i1)-ZZZ(i1+1))/1e5
        LR_LAY(i1)=(T2(i1+1)-T2(i1))/DZ
        if((LR_LAY(i1)-CLR).gt.1e-3) then
         COUNT_CC=COUNT_CC+1
        endif
       enddo
       if(COUNT_CC.eq.0) then 
        CHECK_CC=1
       else
        CHECK_CC=0
       endif
       !do i1=1,NL
       ! write(*,*) T2(i1),T1(i1),LR_LAY(i1),LR_LAY(i1)-CLR,COUNT_CC
       !enddo

      end subroutine
      
      subroutine SET_LOOP
       use PAR_ATM,only: T2,Ttau
       implicit none
       
       Ttau=T2
      end subroutine
      
      subroutine TIME_STEP
       use PAR_ATM,only: Ttau,T0,Tinit,PRE_USE1,NL,dTdt
       implicit none
       real*8 :: dTave
       integer :: i1

       do i1=1,NL
        !dTave=(dTdt(i1+1)+dTdt(i1))*0.5
        T0(i1)=Ttau(i1)+dTdt(i1)!ave
       enddo
       return
      end subroutine
      
      subroutine CALC_HEAT
       use PAR_ATM,only: FL_EX_DIF_UPP,FL_EX_TOT_DWN,NL,PRE_USE1,dTdt
       use PAR_PLA,only: GACC,CP
       use PAR_CODE,only: N_IT,I_TIME
       implicit none
       integer :: i1
       character(len=100) :: fname
       real*8 :: FINN,FOUT,DF,DP
       
       dTdt=0.0
       !write(fname,'(A,I4.4,A)') "OUT_EX/HEAT",N_IT,".out"
       !open(unit=100,file=trim(fname),status="replace")
       do i1=1,NL
        FINN=FL_EX_DIF_UPP(i1+1)+FL_EX_TOT_DWN(i1)
        FOUT=FL_EX_DIF_UPP(i1)+FL_EX_TOT_DWN(i1+1)
        DF=FINN-FOUT
        !From mbar to cgs
        DP=(PRE_USE1(i1+1)-PRE_USE1(i1))
        dTdt(i1)=GACC/CP*DF/DP
        dTdt(i1)=dTdt(i1)*I_TIME
        !write(100,*) PRE_USE1(i1),dTdt(i1)
       enddo
       !close(100)
      end subroutine

      subroutine RUN_EXART
       use PAR_INP,only: M_INP
       implicit none
       character(len=200) :: cmd
       
       if(M_INP.eq.1) then
        write(cmd,'(A)') "./exart_gen > out_exart.txt"
       elseif(M_INP.eq.2) then
        write(cmd,'(A)') "./exart_lr > out_exart.txt"       
       endif
       call system(trim(cmd))
       
      end subroutine
      
