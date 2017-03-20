      program CONRAD
       use P_INP,only: RUN_MODE
       real*8 :: t0

       !c- Read general input file: INPUT/INPUT01
       call READ_INPUT
       !c- Read input spectrum specified in INPUT01 file
       !call READ_SPECTRA
       !c- Read atmospheric layers
       call READ_ATM
       !c- Get Speed of sound and gravitational acceleration
       call INIT_PLANET
       !c- Prepare gas optical depth at each layer
       !stop
       !c- Run main routine
       call RUNNER(RUN_MODE)
       call CPU_TIME(t0)
       !write(*,*) t0
      end
      
      subroutine RUNNER(RM)
       integer,intent(in) :: RM
       !write(*,'(1X,A)') "RUNNING MODE INFOS"
       if(RM.eq.1) then
        !write(*,'(5X,A)') " ---> SPECTRAL MODE"
        call READ_TAU(RM,0)
        call ASSIGN_VARIABLES
        call RUN_SPECTRA
        call RUN_DISORT
       elseif(RM.eq.2) then
        write(*,'(5X,A)') " ---> THERMOCONVECTIVE MODEL"
       elseif(RM.eq.3) then
        write(*,'(5X,A)') " ---> TRANSMISSION MODE"
        call TR_MODE
       elseif(RM.eq.4) then
        write(*,'(5X,A)') " ---> THERMOCON+SPECTRA"
       endif
      end subroutine

      subroutine RUN_SPECTRA
       use P_ATM
       use P_SPE
       use P_PLA,only: GGRAV,MU
       use P_INP,only: LAM_MIN,LAM_MAX,LAM_RES,TAU_OUT,INTERP_KT,MNAME
       use P_CONST,only: MH,KBOL
       implicit none
       integer :: i1,i2,i3,ix
       character(len=100) :: fname
       real*8 :: lam,tra(32),K_GAS(32)
       real*8 :: TII,K_RAY,K_CL,b_gas,b_ray,ttt,gg,mm,kbb
       real*8,dimension(NLYR,NSPMAX) :: TXX,TRR,TCC
       real*8,dimension(NSPMAX) :: LAM_T
       
       call GET_WEIGHTS
       
       !c- Get GGRAV in m/s2 and mh in kms
       gg=GGRAV!/1e2
       mm=MH!/1e3
       kbb=KBOL!/1e7
       
       !c- Fill arrays for the calculation of the spectra
       PP_S=PP
       ZZ_S=ZZ
       TT_S=TT
       LAM_S=0.0
       SSQ=0.0
       TTQ=0.0
       
       !write(*,'(5X,A,I3,A,I3,A)') 
     > ! " ---> ",NLYR," BOUNDARY LAYERS, ",NLYR-1," LAYERS"
     
       do i1=1,NLYR
        !Write out tau
        
        !write(*,'(A1,5X,A,f8.2,A)',advance="no")
     >  !char(13)," ---> READING K-DIST: ",real(i1-1)/(NLYR-1)*1e2,"%"

        b_gas=(PP_S(i1)-PP_S(i1-1))/gg/MU
        b_ray=(PP_S(i1)-PP_S(i1-1))/gg/MU

        if(i1.eq.1) b_gas=(PP_S(2)-PP_S(1))/gg/MU
        if(i1.eq.1) b_ray=(PP_S(2)-PP_S(1))/gg/MU
        
        write(fname,'(A,I3.3,A)') "OUT_EX/",i1,"_k_GEN.dat"
        open(unit=10,file=trim(fname),status='old')
        lam=LAM_MIN
        i3=0
        do while(lam.ge.0 .and. lam.le.LAM_MAX)
         K_RAY=0
         if(lam.ge.LAM_MIN) i3=i3+1
         if(lam.ge.0.39) then
          read(10,*,END=1010) lam,K_GAS(1:32)
          K_GAS(1:PQ)=K_GAS(1:PQ)*b_gas/mm
          do ix=1,PQ
           if(K_GAS(ix).gt.100.0) K_GAS(ix)=100.0
          enddo
         endif
         TII=0
         !! RAYLEIGH scattering
         call RAYL_GEN(lam,1-CC_I(i1),CC_I(i1),HH_I(i1),K_RAY)
         !From cm-1 to tau
         K_RAY=K_RAY*b_ray/mm
         !! OTHER
         if(i1.eq.70) then
          K_CL=0.0
         else
          K_CL=0.0
         endif
         !no line absorption below 0.39
         if(lam.ge.0.40) then
          do i2=1,PQ
           if(K_GAS(i2).lt.0.00) K_GAS(i2)=0.00
           ttt=exp(-K_GAS(i2))
           !Avoid tau=infty and then NaN
           !if(ttt.lt.tiny(one)) ttt=tiny(one)
           TII=TII+WQ(i2)*ttt
           TTQ(i1,i2,i3)=K_GAS(i2)+K_RAY+K_CL
           if(ttt.eq.0) TTQ(i1,i2,i3)=0.0
           SSQ(i1,i2,i3)=(K_RAY+K_CL)/TTQ(i1,i2,i3)
          enddo
         else
          TTQ(i1,1:PQ,i3)=K_RAY+K_CL
          SSQ(i1,1:PQ,i3)=1.0
          TII=1.00
          lam=lam+0.01
         endif
         LAM_T(i3)=lam
         TXX(i1,i3)=TII
         TRR(i1,i3)=K_RAY
         TCC(i1,i3)=K_CL
         TTR(i1,1:PQ,i3)=K_RAY
         TTC(i1,1:PQ,i3)=K_CL
        enddo
1010    close(10)
       enddo
       NLAM=i3
       do i3=1,NLAM
        LAM_S(i3)=LAM_T(i3)
       enddo
       if(TAU_OUT.eq."Y") then
        do i1=1,NLYR
         write(fname,'(A,I3.3,A)') "OUT_TH/L",i1,"_GEN.out"
         open(unit=100+i1,file=trim(fname),status="replace")
         do i3=1,NLAM
          write(100+i1,*) LAM_S(i3),-log(TXX(i1,i3)),
     +     TRR(i1,i3),TCC(i1,i3)
         enddo
         close(100+i1)
        enddo
       endif

      end subroutine
      
      subroutine RUN_DISORT
       use P_INP,only: NLYR,PL_DIST,ND_VAR,FLUX_SCALE,LAM_MIN,LAM_MAX,
     +  RSTAR,MNAME,ALBI,ZANG
       use P_ATM,only: NLMAX,FLUX_EXT,LAM_EXT,NSPMAX,TT_I,ZZ_I,PP_I
       use P_SPE,only: TAU_S,LAM_S,WREAL,NLAM,NP,FDO,FUP,FDOT,FUPT,PQ,
     +  TTQ,SSQ,WQ,FLDNQ,FLDRQ,FLUPQ,TTC,TTR
       use P_DIS
       use P_CONST
       use P_MOD,only: IV
       implicit none
       integer :: i1,i2,i3,ix,iy
       character(len=100) :: fname
       real*8 :: III1,III2,III3,BASE,HEIGHT,HEIGHT2,fx,DT,TT,ITX
       real :: DL,PX(0:7)
       real,dimension(NSPMAX,NLMAX) :: F_LAM1,F_LAM2,F_LAM3,F_LAM4
       real*8,dimension(200) :: YVARR
       real*8,dimension(NLMAX) :: FFF1,FFF2,FFF3,FFF4
       real*8,dimension(NSPMAX) :: LMID
              
       DTAUC=0.0
       !open(unit=101,file="OUTPUT/F_TOT_"//IV//".txt",status="replace")
       !open(unit=102,file="OUTPUT/F_DIR_"//IV//".txt",status="replace")
       !open(unit=103,file="OUTPUT/F_DIF_"//IV//".txt",status="replace")
       !open(unit=104,file="OUTPUT/F_UPD_"//IV//".txt",status="replace")
       
       
       !call GET_INT_BANDS(LAM_EXT,FLUX_EXT)
       call GET_YV(NLAM,YVARR)
       ITX=0.0
       do i1=1,NLAM-1
        !write(*,*) LAM_S(i1),YVARR(i1)
        ITX=ITX+YVARR(i1)
       enddo
       !write(*,*) ITX
       !stop


       do i1=1,NLAM-1
        !write(*,'(A1,A,f8.2,A)',advance="no")
     >  !char(13),"RUNNING DISORT ",real(i1)/(NLAM-1)*100,"%"
        DL=(LAM_S(i1+1)-LAM_S(i1))
        FBEAM = YVARR(i1)
        if(LAM_S(i1).gt.4.5) FBEAM=0.0

        do i2=1,PQ
         PMOM = 0.0
         PC=0.0
         PR=0.0
         PX=0.0
         NSTR = 8
         CALL  GETMOM( 3, -0.99, NSTR, PX )
         !write(*,*) PX
         !stop
         do ix=1,NLYR
          PC(0:NSTR-1,ix)=PX(0:NSTR-1)
         enddo
         PR(0,1:NLYR)=1.0
         PR(2,1:NLYR)=0.1
         !Average 
         do iy=0,NSTR-1
         do ix=1,NLYR
          PMOM(iy,ix)=
     +    (PR(iy,ix)*TTR(ix,i2,i1)+PC(iy,ix)*TTC(ix,i2,i1))
     +     /(TTC(ix,i2,i1)+TTR(ix,i2,i1))
         enddo
         enddo
         
         USRTAU    = .FALSE.
         NTAU      = 0
         UTAU      = 0.0
         USRANG    = .FALSE.
         NUMU      = 0
         UMU       = 0.0
         NPHI      = 0
         PHI       = 0.0
         IBCND     = 0
         UMU0      = cos(ZANG*PI/180.0)
         PHI0      = 0.0
         LAMBER    = .TRUE.
         HL        = 0.0
         PLANK     = .TRUE.
         TEMIS     = 0.0
         
         if(LAM_S(i1).lt.4.5) then
          ALBEDO    = ALBI
         elseif(LAM_S(i1).ge.4.5) then
          ALBEDO = 0.0
         endif
         DELTAM    = .TRUE.
         ONLYFL    = .TRUE.
         FISOT     = 0.0!TT_I(1)
         ACCUR     = 0.001
         TEMPER(1:NLYR)    = TT_I(1:NLYR)
         TTEMP     = TT_I(1)
         BTEMP     = TT_I(NLYR)
         
        PRNT=(/.false.,.false.,.false.,.false.,.false.,.false.,.false./)
  
         WVNMLO    = 1/LAM_S(i1+1)*1e4
         WVNMHI    = 1/LAM_S(i1)*1e4
        do i3=1,NLYR
         DTAUC(i3)=TTQ(i3,i2,i1)
         SSALB(i3)=SSQ(i3,i2,i1)
         if(isnan(SSALB(i1))) SSALB(i1)=0.0
         if(SSALB(i1).lt.0.0) SSALB(i1)=0.0
         if(SSALB(i1).gt.1.0) SSALB(i1)=1.0
         if(DTAUC(i3).lt.1e-5) DTAUC(i3)=0.0
         if(DTAUC(i3).gt.1e3) DTAUC(i3)=1e3
        enddo
        do i3=1,NLYR
         if(SSALB(i3).gt.1) then
          SSALB(i3)=1.00
         endif
         if(SSALB(i3).lt.0) SSALB(i3)=SSALB(i3-1)
         if(DTAUC(i3).lt.0) DTAUC(i3)=DTAUC(i3-1)
        enddo

         if(0.eq.0) then
         call DISORT(NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     +                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     +                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     +                   FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     +                   DELTAM, PLANK, ONLYFL, ACCUR, PRNT, HEADER,
     +                   MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, RFLDIR,
     +                   RFLDN, FLUP, DFDT, UAVG, UU, U0U, ALBMED,
     +                   TRNMED)
         endif
         !COLLECT FLUXES [Given by disort in W/m2]
         do i3=1,NLYR+1
          FLDNQ(i3,i2)=RFLDIR(i3)
          FLDRQ(i3,i2)=RFLDN(i3)
          FLUPQ(i3,i2)=FLUP(i3)
         enddo
        enddo
        FFF1(1:NLYR+1)=0
        FFF2(1:NLYR+1)=0
        FFF3(1:NLYR+1)=0
        do i3=1,NLYR+1
         do i2=1,PQ
          FFF1(i3)=FFF1(i3)+FLDNQ(i3,i2)*WQ(i2)
          FFF2(i3)=FFF2(i3)+FLDRQ(i3,i2)*WQ(i2)
          FFF3(i3)=FFF3(i3)+FLUPQ(i3,i2)*WQ(i2)
         enddo
        enddo
        F_LAM1(i1,1:NLYR+1)=FFF1(1:NLYR+1)+FFF2(1:NLYR+1)
        F_LAM2(i1,1:NLYR+1)=FFF1(1:NLYR+1)
        F_LAM3(i1,1:NLYR+1)=FFF3(1:NLYR+1)
        !write(101,*) LAM_S(i1),(FFF1(1:NLYR+1)+FFF2(1:NLYR+1))/DL
        !write(102,*) LAM_S(i1),FFF1(1:NLYR+1)/DL
        !write(103,*) LAM_S(i1),FFF2(1:NLYR+1)/DL
        !write(104,*) LAM_S(i1),FFF3(1:NLYR+1)/DL  
      enddo
      !close(101)
      !close(102)
      !close(103)        
      !close(104)        
      
      call PRINT_FLUX_INT(F_LAM1,F_LAM2,F_LAM3)

      write(*,*)
      end subroutine
      
      subroutine PRINT_FLUX_INT(F_LAM1,F_LAM2,F_LAM3)
       use P_INP,only: NLYR,MNAME
       use P_SPE,only: NLAM,LAM_S
       use P_ATM,only: NSPMAX,NLMAX,ZZ_I,PP_I
       use P_MOD,only: IV
       implicit none
       integer :: i1,i2,i3
       real,intent(in),dimension(NSPMAX,NLMAX) :: F_LAM1,F_LAM2,F_LAM3
       real*8,dimension(NSPMAX):: F_LAM_LAY1,F_LAM_LAY2,F_LAM_LAY3,LMID
       real*8,dimension(NLMAX) :: PPPR
       real*8 :: III1,III2,III3,III4,HEIGHT,BASE,DL,IIUPV
       PPPR(2:NLYR+1)=PP_I(1:NLYR)
       
       open(unit=100,file="OUT_TH/FL_EXART_MR.dat",status="replace")

       do i1=2,NLYR+1
        III1=0.0
        III2=0.0
        III3=0.0
        IIUPV=0.0
        do i2=1,NLAM-1
         DL=LAM_S(i1+1)-LAM_S(i1)
         III1=III1+F_LAM1(i2,i1)
         III2=III2+F_LAM2(i2,i1)
         III3=III3+F_LAM3(i2,i1)
         if(LAM_S(i2).lt.4.5) then
          IIUPV=IIUPV+F_LAM3(i2,i1)
         endif
         F_LAM_LAY1(1:NLAM)=F_LAM1(1:NLAM,i1)/DL
         F_LAM_LAY2(1:NLAM)=F_LAM2(1:NLAM,i1)/DL
         F_LAM_LAY3(1:NLAM)=F_LAM3(1:NLAM,i1)/DL
        enddo  
        write(100,*) PPPR(i1),III1,III2,III3,III1-IIUPV
        !write(*,*) PPPR(i1)/100,III1,III2,III3,III1-IIUPV
       enddo 
       close(100)       
      end subroutine
      

      subroutine ASSIGN_VARIABLES
       use P_ATM
       PP=PP_I       
       ZZ=ZZ_I       
       TT=TT_I       
      end subroutine


