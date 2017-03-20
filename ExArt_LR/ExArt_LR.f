      program ExArt
       ! Radiative transfer in a CO2/H2O/N2 atmosphere 
       ! All quantities must be expressed in cgs 
       ! [apart from the fluxes which are in W/m2]
       use P_INP,only: MODE,WRITE_OUT,ZANG
       use P_DIS
       use P_ST,only: YV
       use P_PL,only: MU,gg,ALBI,LAMBI
       use P_SPE,only: FLDNQ,FLDRQ,FLUPQ,PPP,TTT,CCC,HHH,TTQ,SSQ,NLYR,
     >  NLMAX,psca,tsca,PQ,BAND,BANDV,NNN,WM,NBB,NWW
       use P_CONST,only: PI
       
       implicit none
       integer :: i1,i2,i3
       real*8,dimension(NLMAX,15,PQ) :: tauKK,tauRAYL
       real*8 :: III1,III2,III3,BASE,HEIGHT,wl1,wl2,wlm,WQ(PQ)
       real*8 :: lam,DL,dummy,SIG,SIG2,DP
       character(len=100) :: fname
       real*8,dimension(NLMAX) :: FFF1,FFF2,FFF3,PPPR
       real*8,dimension(15,NLMAX) :: F_LAM1,F_LAM2,F_LAM3
       real*8,dimension(NLMAX) :: III1vis, III1ir, III3vis, III3ir
       real*8,dimension(NLMAX) :: III2vis, III2ir
       real*8 :: t0
      
       call READ_INPUT
       MODE=1
       WRITE_OUT=0
       NBB=14
       NWW=32
       !D: Collect weights for the K_DISTRIBUTION
       open(unit=10,file="32_mis.dat",status="old")
       do i1=1,NWW
        read(10,*) WQ(i1),dummy
       enddo
       close(10)
              
       !D: Read in atmosphere
       ! Pressure in Barye, Temp in K
      open(unit=10,file="INP_TH/atms_ex.inp",status="old")
       i1=0
       do
        i1=i1+1
        read(10,*,END=1000) PPP(i1),TTT(i1), CCC(i1), HHH(i1)
        !write(*,*) PPP(i1),TTT(i1)+50, CCC(i1), HHH(i1)
        NNN(i1)=1.d0-CCC(i1)-HHH(i1)
       enddo
1000   close(10)
       NLYR=i1-1
       
       
       call GET_MU
       call GET_BANDS
       call GET_YV
       call READ_TAU

             
       do i1=1,NLYR
        write(fname,'(A,I3.3,A)') "OUT_EX_LR/TAU/",i1,".out"
     
        open(unit=1000,file=trim(fname),status="old")
        do i2=1,NBB
         read(1000,*) dummy,tauKK(i1,i2,1:NWW)
        enddo
        close(1000)
       enddo

       do i1=1,NLYR
        DP=(PPP(i1)-PPP(i1-1))*psca(i1)*tsca(i1)
        do i2=1,NBB
         wl1=1.0/BANDV(i2)*1e4
         wl2=1.0/BANDV(i2+1)*1e4
         wlm=0.5*(wl1+wl2)
         lam=1/wlm*1e4
         call RAYL_GEN(lam,NNN(i1),CCC(i1),HHH(i1),SIG)
         tauRAYL(i1,i2,1:NWW)=SIG*DP
        enddo
        do i2=1,NBB
         do i3=1,NWW
          TTQ(i1,i2,i3)=tauKK(i1,i2,i3)*DP+tauRAYL(i1,i2,i3)
          SSQ(i1,i2,i3)=tauRAYL(i1,i2,i3)/TTQ(i1,i2,i3)
         enddo
         if(i1.eq.98) write(444,*) BANDV(i2),tauKK(i1,i2,1:NWW)*DP
        enddo
       enddo

       PPPR(1)=1e-6
       PPPR(2:NLYR+1)=PPP(1:NLYR)

        do i3=1,NBB
         if(WRITE_OUT.eq.1) write(*,*) "BAND",i3
         DL=BANDV(i3+1)-BANDV(i3)
         do i2=1,PQ
          FBEAM=YV(i3)
          HEADER=""
          PMOM = 0.0
          NSTR = 8
          !D: legendre coefficients for Rayleigh scattering
          PMOM(0,1:NLYR)=1.0
          PMOM(1,1:NLYR)=0.0
          PMOM(2,1:NLYR)=0.1

          USRTAU    = .FALSE.
          NTAU      = 0
          UTAU      = 0.0
          USRANG    = .FALSE.
          NUMU      = 0
          UMU       = 0.0
          NPHI      = 8
          PHI       = 0.0
          IBCND     = 0
          UMU0      = cos(ZANG*PI/180.0)
          PHI0      = 0.0
          HL        = 0.0
          if(BANDV(i3).lt.4.5) then
           PLANK     = .TRUE.
          else
           PLANK     = .TRUE.
          endif
          TEMIS     = 0.0
!        
          if(i3.le.7) ALBEDO    = ALBI
          if(i3.ge.8) ALBEDO    = 0.0
          if(LAMBI.eq.1) then 
           LAMBER=.TRUE.
          elseif(LAMBI.eq.0) then
           LAMBER=.FALSE.
          endif
!       
          DELTAM    = .TRUE.
          ONLYFL    = .TRUE.
          FISOT     = 0.0
          ACCUR     = 0.001
          TEMPER(1:NLYR)    = TTT(1:NLYR)
          TTEMP     = TTT(1)
          BTEMP     = TTT(NLYR)
        PRNT=(/.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./)
!   
          WVNMLO    = 1/BANDV(i3+1)*1e4
          WVNMHI    = 1/BANDV(i3)*1e4

          DTAUC=0.0
          SSALB=0.0
          do i1=1,NLYR
            DTAUC(i1)=TTQ(i1,i3,i2)
           if(isnan(DTAUC(i1))) DTAUC(i1)=0.0
           if(DTAUC(i1).lt.1e-5) DTAUC(i1)=0.0
           if(DTAUC(i1).gt.1e3) DTAUC(i1)=1e3
           SSALB(i1)=SSQ(i1,i3,i2)
           !SSALB(i1)=(SSQ(i1,i3,i2)+SSQ(i1-1,i3,i2))
           if(isnan(SSALB(i1))) SSALB(i1)=0.0
           if(SSALB(i1).lt.0.0) SSALB(i1)=0.0
           if(SSALB(i1).gt.1.0) SSALB(i1)=1
          enddo
          
          call DISORT(NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     +                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     +                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     +                   FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     +                   DELTAM, PLANK, ONLYFL, ACCUR, PRNT, HEADER,
     +                   MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, RFLDIR,
     +                   RFLDN, FLUP, DFDT, UAVG, UU, U0U, ALBMED,
     +                   TRNMED)
          do i1=1,NLYR+1
           if(RFLDIR(i1).lt.0) RFLDIR(i1)=0.0
           if(RFLDN(i1).lt.0) RFLDN(i1)=0.0
           if(FLUP(i1).lt.0) FLUP(i1)=0.0
           FLDNQ(i1,i2)=RFLDIR(i1)/DL
           FLDRQ(i1,i2)=RFLDN(i1)/DL
           FLUPQ(i1,i2)=FLUP(i1)/DL
          enddo
         enddo

         FFF1(1:NLYR+1)=0
         FFF2(1:NLYR+1)=0
         FFF3(1:NLYR+1)=0
         do i1=1,NLYR+1
          do i2=1,PQ
           FFF1(i1)=FFF1(i1)+FLDNQ(i1,i2)*WQ(i2)
           FFF2(i1)=FFF2(i1)+FLDRQ(i1,i2)*WQ(i2)
           FFF3(i1)=FFF3(i1)+FLUPQ(i1,i2)*WQ(i2)
          enddo
         enddo
         F_LAM1(i3,1:NLYR+1)=(FFF1(1:NLYR+1)+FFF2(1:NLYR+1))
         F_LAM2(i3,1:NLYR+1)=FFF1(1:NLYR+1)
         F_LAM3(i3,1:NLYR+1)=FFF3(1:NLYR+1)
        enddo
        
         !open(unit=200, file="OUTPUT/fl_vis.dat", status="replace")
         !open(unit=201, file="OUTPUT/fl_ir.dat", status="replace")

        open(unit=800, file="OUT_TH/FL_EXART_LR.dat", status="replace")
        do i1=2,NLYR+1
         III1=0
         III2=0
         III3=0
         do i2=1,7
          BASE=BANDV(i2+1)-BANDV(i2)
          !BASE=1
          !HEIGHT=(F_LAM1(i2,i1)+F_LAM1(i2+1,i1))
          HEIGHT=F_LAM1(i2,i1)
          III1=III1+BASE*HEIGHT
          !HEIGHT=(F_LAM2(i2,i1)+F_LAM2(i2+1,i1))
          HEIGHT=F_LAM2(i2,i1)
          III2=III2+BASE*HEIGHT
          !HEIGHT=(F_LAM3(i2,i1)+F_LAM3(i2+1,i1))
          HEIGHT=F_LAM3(i2,i1)
          III3=III3+BASE*HEIGHT
         enddo  
         
         III1vis(i1)=III1
         III2vis(i1)=III2
         III3vis(i1)=III3
         
         write(*,*) PPPR(i1),III1,III2,III3
         !write(555,*) PPPR(i1)*1e2,III1,III2,III3
        enddo

        do i1=2,NLYR+1
         III1=0
         III2=0
         III3=0
         do i2=8,NBB
          BASE=BANDV(i2+1)-BANDV(i2)
          HEIGHT=F_LAM1(i2,i1)
          III1=III1+BASE*HEIGHT
          HEIGHT=F_LAM2(i2,i1)
          III2=III2+BASE*HEIGHT
          HEIGHT=F_LAM3(i2,i1)
          III3=III3+BASE*HEIGHT
         enddo  
         III1ir(i1)=III1
         III2ir(i1)=III2
         III3ir(i1)=III3
         write(556,*) PPPR(i1),III1,III2,III3
         write(800,*) PPPR(i1),
     >    (III1vis(i1)+III1ir(i1)),
     >    (III2vis(i1)+III2ir(i1)),
     >    (III3vis(i1)+III3ir(i1)),
     >    III1vis(i1)-III3vis(i1)+III1ir(i1)
        enddo 
        close(800)
        !close(200)
        !close(201)
        call CPU_TIME(t0)
       write(*,*) t0

      end program


       subroutine READ_TAU
         !***************************************************************
         !                                                              *
         ! Collect all the available tau_k                              * 
         !                                                              *
         !***************************************************************
         use P_INP,only: MODE,WRITE_OUT
         use P_SPE
         implicit none
         character(len=200) :: fname
         integer :: i1,i2,TTV,PPV,CCV,HHV,NP,NQ
         real*8,dimension(NSPMAX,PQ) :: yaaaa,yaaab,yaaba,yaabb
         real*8,dimension(NSPMAX,PQ) :: yabaa,yabab,yabba,yabbb
         real*8,dimension(NSPMAX,PQ) :: ybaaa,ybaab,ybaba,ybabb
         real*8,dimension(NSPMAX,PQ) :: ybbaa,ybbab,ybbba,ybbbb,yy
         real*8,dimension(NSPMAX) :: lam
         real*8,dimension(NLMAX) :: ppuse,ttuse,ccuse,hhuse

         call GET_INDEXES(PPV,TTV,CCV,HHV)

         ppuse=PPP
         ttuse=TTT
         ccuse=CCC
         hhuse=HHH
         NQ=PQ
         !D: Interpolate log values
         PPX=log10(PPX)
         ppuse=log10(ppuse)
         if(WRITE_OUT.eq.1) then
          open(unit=20,file="OUT_EX_LR/interp.res",status="replace")
           !write(*,*) "*** INTERPOLATING ABSORPTION COEFFICIENTS"
         endif
         
         do i1=1,NLYR
          if(WRITE_OUT.eq.1) then
           write(20,*) 
           write(20,'(A,I3.3,4(F15.4))') "LAYER ",i1,
     >      ppuse(i1),ttuse(i1),ccuse(i1),hhuse(i1)
          endif
          do i2=1,PPV-1
           if(ppuse(i1).ge.PPX(i2).and.ppuse(i1).le.PPX(i2+1)) then
            !If within interval, done, exit
            INDXP(i1)=i2
            pint(i1)=(ppuse(i1)-PPX(i2))/(PPX(i2+1)-PPX(i2))
            psca(i1)=1
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4),1X,I2)') 'PRE CASE A: ',
     >       LABP(INDXP(i1)),PPX(i2),ppuse(i1),PPX(i2+1),pint(i1)
            endif
            exit
           elseif(ppuse(i1).lt.PPX(1)) then
            !If lower than minimum, done, exit
            INDXP(i1)=1
            pint(i1)=0.00
            !psca(i1)=10**ppuse(i1)/10**PPX(1)
            psca(i1)=ppuse(i1)/PPX(1)
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4),1X,I2)') 'PRE CASE B: ',
     >       LABP(INDXP(i1)),PPX(i2),ppuse(i1),PPX(i2+1),psca(i1)
            endif
            exit
           elseif(ppuse(i1).gt.PPX(PPV)) then
            !If higher than maximum, done, exit
            INDXP(i1)=PPV-1
            pint(i1)=1.00
            psca(i1)=ppuse(i1)/PPX(PPV)
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4),1X,I2)') 'PRE CASE C: ',
     >       LABP(INDXP(i1)),PPX(i2),ppuse(i1),PPX(i2+1),psca(i1)
            endif
            exit
           endif
          enddo
          do i2=1,TTV-1
           if(ttuse(i1).ge.TTX(i2).and.ttuse(i1).le.TTX(i2+1)) then
            INDXT(i1)=i2
            tint(i1)=(ttuse(i1)-TTX(i2))/(TTX(i2+1)-TTX(i2))
            tsca(i1)=1
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE A: ',
     >        INDXT(i1),TTX(i2),ttuse(i1),TTX(i2+1),tint(i1)
            endif
            exit
           elseif(ttuse(i1).lt.TTX(1)) then
            INDXT(i1)=1
            tint(i1)=0.00
            tsca(i1)=ttuse(i1)/TTX(1)
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE B: ',
     >        INDXT(i1),TTX(i2),ttuse(i1),TTX(i2+1),tsca(i1)
            endif
            exit
           elseif(ttuse(i1).gt.TTX(TTV)) then
            INDXT(i1)=TTV-1
            tint(i1)=1.00
            tsca(i1)=ttuse(i1)/TTX(TTV)
            if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE C: ',
     >        INDXT(i1),TTX(i2),ttuse(i1),TTX(i2+1),tsca(i1)
            endif
            exit
          endif
          enddo
         do i2=1,CCV-1
          if(ccuse(i1).ge.CCX(i2).and.ccuse(i1).le.CCX(i2+1)) then
           INDXC(i1)=i2
           cint(i1)=(ccuse(i1)-CCX(i2))/(CCX(i2+1)-CCX(i2))
           if(CCX(i2+1).eq.CCX(i2)) cint(i1)=1.0
           csca(i1)=1
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'CO2 CASE C: ',
     >        INDXC(i1),CCX(i2),ccuse(i1),CCX(i2+1),cint(i1)
           endif
           exit
          elseif(ccuse(i1).lt.CCX(1)) then
           INDXC(i1)=1
           cint(i1)=0.00
           csca(i1)=ccuse(i1)/CCX(1)
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'CO2 CASE B: ',
     >        INDXC(i1),CCX(i2),ccuse(i1),CCX(i2+1),csca(i1)
           endif
           exit
          elseif(ccuse(i1).gt.CCX(CCV)) then
           INDXC(i1)=CCV-1
           cint(i1)=1.00
           csca(i1)=ccuse(i1)/CCX(CCV)
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'CO2 CASE C: ',
     >        INDXC(i1),CCX(i2),ccuse(i1),CCX(i2+1),csca(i1)
           endif
           exit
         endif
         enddo
         do i2=1,HHV-1
          if(hhuse(i1).ge.HHX(i2).and.hhuse(i1).le.HHX(i2+1)) then
           INDXH(i1)=i2
           hint(i1)=(hhuse(i1)-HHX(i2))/(HHX(i2+1)-HHX(i2))
           if(HHX(i2+1).eq.HHX(i2)) hint(i1)=1.00
           hsca(i1)=1
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE A: ',
     >        INDXH(i1),HHX(i2),hhuse(i1),HHX(i2+1),hint(i1)
           endif 
           exit
          elseif(hhuse(i1).lt.HHX(1)) then
           INDXH(i1)=1
           hint(i1)=0.00
           hsca(i1)=hhuse(i1)/HHX(1)
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE A: ',
     >        INDXH(i1),HHX(i2),hhuse(i1),HHX(i2+1),hsca(i1)
           endif 
           exit
          elseif(hhuse(i1).gt.HHX(HHV)) then
           INDXH(i1)=HHV-1
           hint(i1)=1.00
           hsca(i1)=hhuse(i1)/HHX(HHV)
           if(WRITE_OUT.eq.1) then
             write(20,'(A,I2.2,4(E13.4))') 'H2O CASE A: ',
     >        INDXH(i1),HHX(i2),hhuse(i1),HHX(i2+1),hsca(i1)
           endif
           exit
         endif
         enddo
        
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)),
     >     LABH(INDXH(i1)),yaaaa,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)),
     >     LABH(INDXH(i1)+1),yaaab,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)+1),
     >     LABH(INDXH(i1)),yaaba,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)+1),
     >     LABH(INDXH(i1)+1),yaabb,lam,i1,NQ)

         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),yabaa,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),yabab,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),yabba,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),yabbb,lam,i1,NQ)


         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),ybaaa,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),ybaab,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),ybaba,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),ybabb,lam,i1,NQ)

         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),ybbaa,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),ybbab,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),ybbba,lam,i1,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),ybbbb,lam,i1,NQ)
     
         yy(1:NSPMAX,1:NQ)=
     >     (1-pint(i1))*(1-tint(i1))*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    yaaaa(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*                    yaaab(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    yaaba(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*                    yaabb(1:NSPMAX,1:NQ)+
     >
     >     (1-pint(i1))*tint(i1)*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    yabaa(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*                    yabab(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    yabba(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*                    yabbb(1:NSPMAX,1:NQ)+
     >
     >     pint(i1)*(1-tint(i1))*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    ybaaa(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*                    ybaab(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    ybaba(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*                    ybabb(1:NSPMAX,1:NQ)+
     >
     >     pint(i1)*tint(i1)*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    ybbaa(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*                    ybbab(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*                    ybbba(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*                    ybbbb(1:NSPMAX,1:NQ)
              
          !Write out
          write(fname,'(A,I3.3,A)') "OUT_EX_LR/TAU/",i1,".out"
          open(unit=12,file=trim(fname),status="replace")
          do i2=1,NBB
           if(isnan(yy(i2,1))) yy(i2,1:NQ)=yy(i2-1,1:NQ)
!!         k is in cm2/mol
           write(12,*) lam(i2),yy(i2,1:NQ)
          enddo
!          close(11)
          close(12)
         enddo
         if(WRITE_OUT.eq.1) then
          close(20)
         endif
         !D: Restore linear P
         ppuse=10**ppuse
        end subroutine
!        
        subroutine GET_TAU(j,k,l,m,yy,lam,NL,NW)
         use P_INP,only: MODE
         use P_PL,only: mu,gg
         use P_SPE,only: NSPMAX,NWW,WM
         use P_CONST,only: MH
         implicit none
         integer,intent(in) :: j,k,l,m,NW
         real*8 :: tt(NWW),KBOL,PM,TM
         integer :: i1,i2
         real*8,intent(out) :: yy(NSPMAX,NWW)
         real*8,intent(out) :: lam(NSPMAX)
         integer,intent(in) :: NL
         character(len=100) :: fname
         
         !D: need to scale absorption coefficients
         write(fname,'(A,4(I2.2),A)') "K_ABS_LR/MODEL",j,k,l,m,"_LR.out"
         i1=1
         !k are given in cm2
         open(unit=10,file=trim(fname),status="old")
         do 
          i1=i1+1
          read(10,*,END=1001) lam(i1),tt(1:NWW)
          tt(1:NWW)=tt(1:NWW)/MU/MH/gg!/1e4
          do i2=1,NWW
           if(tt(i2).le.0.0) tt(i2)=1e-30
          enddo
          yy(i1,1:NWW)=tt(1:NWW)
         enddo
 1001    close(10)
         yy(1,1:NWW)=0.0
         return
        end subroutine
!        
       subroutine GET_INDEXES(PPV,TTV,CCV,HHV)
        use P_INP,only: MODE
        use P_SPE,only: PPX,TTX,HHX,CCX,LABP,LABT,LABH,LABC,WM
        implicit none
        integer,intent(out) :: PPV,TTV,CCV,HHV
        integer :: i1,i2
        character(len=50) :: fname
        
        if(MODE.eq.1) write(fname,'(A)') "CALC_K_LR.dat"
        
        open(unit=10,file=trim(fname),status="old")
         read(10,*) 
         read(10,*) PPV
         do i1=1,PPV
          read(10,*) PPX(i1),LABP(i1)
         enddo
         read(10,*) TTV
         do i1=1,TTV
          read(10,*) TTX(i1),LABT(i1)
         enddo
         read(10,*) CCV
         do i1=1,CCV
          read(10,*) CCX(i1),LABC(i1)
         enddo
         read(10,*) HHV
         do i1=1,HHV
          read(10,*) HHX(i1),LABH(i1)
         enddo
        close(10)
        !From atm to Barye
        PPX=PPX*1013250
        return
       end subroutine
        

      subroutine RAYL_GEN(lam,ab_N2,ab_CO2,ab_H2O,SIG)
      !take in input lambda [micron]
      !take in input N2/CO2/H2O abundances [0-1]
      !take in input press and temp [atm and K]
      !give out SIGMA in cm2/molecule
       use P_PL,only: gg,mu
       use P_CONST,only: MH,KBOL,PI
       implicit none
       real*8,intent(in) :: lam,ab_N2,ab_CO2,ab_H2O
       real*8,intent(out) :: SIG
       real*8 :: delN2,delCO2,delH2O
       real*8 :: AN2,ACO2,AH2O,BN2,BCO2,BH2O
       real*8 :: sigN2,sigCO2
       real*8 :: Dh2o,rh2o,sigH2O,rair,NN
       
       AN2 =29.06e-5
       ACO2=43.9e-5
       AH2O=516e-5
       BN2 =7.7e-3
       BCO2=6.4e-3
       BH2O=0.0e-3
       
       delN2=1.05
       delCO2=1.15
       delH2O=1.35 !D=0.17
       
       sigN2=4.577e-21*delN2/lam**4*(AN2*(1+BN2/lam**2))**2
       sigCO2=4.577e-21*delCO2/lam**4*(ACO2*(1+BCO2/lam**2))**2
       sigH2O=0.75*sigN2
       

       !sig is in cm2/molecule
       SIG=ab_N2*sigN2+ab_CO2*sigCO2+ab_H2O*sigH2O
       SIG=SIG/MU/MH/gg
       return
      end subroutine

      subroutine READ_INPUT
       !READ FILE input_ex.inp
       use P_PL, only: PL_DIST,ND_VAR,FL_SCALE,gg,PL_RAD,PL_MAS,
     +   LAMBI,ALBI
       use P_ST, only: YV
       use P_INP,only: MODE,INPUT_SP,WRITE_OUT,TST,RST,ZANG
       use P_CONST,only: GRAV,MEARTH,REARTH
       implicit none
       real*8 :: MPL,RPL
       real*8 :: LMIN,LMAX
       
       open(unit=10,file="INP_TH/input_ex.inp",status="old")
        read(10,*) INPUT_SP
        read(10,*) TST
        read(10,*) RST
        read(10,*) PL_DIST
        read(10,*) ND_VAR
        read(10,*) FL_SCALE
        read(10,*) ZANG
        read(10,*) 
        read(10,*) PL_RAD
        read(10,*) PL_MAS
        read(10,*)
        read(10,*) LMIN
        read(10,*) LMAX
        read(10,*) 
        read(10,*) 
        read(10,*) 
        read(10,*) 
        read(10,*) 
        read(10,*)
        read(10,*) ALBI
       close(10)
       LAMBI=1
       MPL=MEARTH*PL_MAS
       RPL=REARTH*PL_RAD
       gg=GRAV*MPL/RPL**2
       
       !write(*,*) MPL,RPL,gg
       !stop
      end subroutine
      
      subroutine GET_YV
       use P_INP,only: TST,RST
       use P_PL,only: PL_DIST,PL_RAD,FL_SCALE,ND_VAR
       use P_ST, only: YV
       use P_SPE,only: BANDV
       use P_CONST,only: CL,HPLANCK,KBOL,PI,RSUN,AU,REARTH
       
       implicit none
       integer :: i1,i2,ix,NSP
       real*8,dimension(15) :: nu,band_cm,III
       real*8,dimension(10000) :: BB,nn
       real*8 :: nu_i,nu_f,dn,BASE,HEIGHT,RR,IT
       
       NSP=10000
       !Lam in cm
       band_cm=BANDV/1e4
       !Get nu [Hz]
       nu=CL/band_cm
       nu_i=CL/(1000.0/1e4)
       nu_f=CL/(0.24/1e4)
       
       BB(1)=0.0
       nn(1)=nu_i
       do i1=2,NSP
        dn=(nu_f-nu_i)/NSP
        nn(i1)=nn(i1-1)+dn
        BB(i1)=2.0*HPLANCK*nn(i1)**3/CL**2
        BB(i1)=BB(i1)/(exp(HPLANCK*nn(i1)/(KBOL*TST))-1)
        BB(i1)=BB(i1)*4.0*PI !erg/s/cm2/Hz
        write(900,*) nn(i1),BB(i1) 
       enddo
       
       do ix=1,15
        III(ix)=0
        do i1=1,NSP-1
         BASE=nn(i1+1)-nn(i1)
         HEIGHT=0.5*(BB(i1+1)+BB(i1))
         if(nn(i1).gt.nu(ix+1).and.nn(i1+1).lt.nu(ix)) then
          III(ix)=III(ix)+BASE*HEIGHT
         endif
        enddo       
       enddo       
       RR=RST*RSUN
       III=PI*III*RR**2/(4.0*PI*PL_DIST**2*AU**2)/1e7*1e4
       do i1=1,14
        IT=IT+III(i1)
       enddo
       !This fluxes are given in W/m2
       YV(1:14)=III(1:14)*ND_VAR*FL_SCALE
      end subroutine
      
      subroutine GET_YV2
       use P_INP, only: WRITE_OUT
       use P_ST, only: YV
       use P_PL, only: PL_DIST,FL_SCALE,ND_VAR
       implicit none
       integer :: i1
       real*8 :: YVTOT
       
       YV(1)=25.61
       YV(2)=102.38
       YV(3)=56.69
       YV(4)=20.05
       YV(5)=8.25
       YV(6)=3.93
       YV(7)=2.27
       YV(8)=0.0
       YV(9)=0.0
       YV(10)=0.0
       YV(11)=0.0
       YV(12)=0.0
       YV(13)=0.0
       YV(14)=0.0
      
       YV(1)=53.70
       YV(2)=332.15
       YV(3)=179.28
       YV(4)=73.71
       YV(5)=27.58
       YV(6)=12.23
       YV(7)=6.41
       YV(8)=0.0
       YV(9)=0.0
       YV(10)=0.0
       YV(11)=0.0
       YV(12)=0.0
       YV(13)=0.0
       YV(14)=0.0
      
       YV=YV/PL_DIST**2*(2.0*ND_VAR)*FL_SCALE
       YVTOT=0.0
       do i1=1,7
        YVTOT=YVTOT+YV(i1)
       enddo
       if(WRITE_OUT.eq.1) write(*,'(A,F13.4,A)') 
     +  "*** TOTAL FLUX: ",YVTOT," W/m2" 
      end subroutine
      
      subroutine GET_MU
       use P_SPE,only: CCC,NNN,HHH
       use P_CONST,only: MH
       use P_PL,only: MU
       implicit none
       real*8 :: mCO2,mN2,mH2O
       
       mH2O=18.00
       mCO2=44.00
       mN2=28.00
       
       MU=CCC(1)*mCO2+NNN(1)*mN2+HHH(1)*mH2O

      end subroutine
      
      subroutine GET_BANDS
      use P_SPE,only: BANDV
      use P_INP,only: MODE
      implicit none

      BANDV(1)= 0.24
      BANDV(2)= 0.40
      BANDV(3)= 0.80
      BANDV(4)= 1.31
      BANDV(5)= 1.86
      BANDV(6)= 2.48
      BANDV(7)= 3.24
      BANDV(8)= 4.50
      BANDV(9)= 8.00
      BANDV(10)=12.00
      BANDV(11)=14.00
      BANDV(12)=16.00
      BANDV(13)=24.00
      BANDV(14)=60.0
      BANDV(15)=1000.0
      
      if(MODE.eq.3) then
       BANDV(8)= 4.50
       BANDV(9)= 8.00
       BANDV(10)=16.00
       BANDV(11)=32.00
       BANDV(12)=60.00       
       BANDV(13)=1000.00       
      endif

      end subroutine

      
