       subroutine READ_INPUT
        !***************************************************************
        !                                                              *
        ! Read input file INPUT01 where several parameters are         *
        ! specified by the user.                                       *
        ! See INPUT/GUIDE.info to know more                            *            
        !                                                              *    
        !***************************************************************
        use P_INP
        use P_MOD,only: IV
        implicit none
        
        open(unit=10,file="INP_TH/input_ex.inp",status="old")
         read(10,*) INPUT_SP
         read(10,*) TST
         read(10,*) RST
         read(10,*) PL_DIST
         read(10,*) ND_VAR
         read(10,*) FLUX_SCALE
         read(10,*) ZANG
         read(10,*) 
         read(10,*) RPL
         read(10,*) MPL
         read(10,*) 
         read(10,*) LAM_MIN
         read(10,*) LAM_MAX
         read(10,*) 
         read(10,*) INTERP_SP
         read(10,*) INTERP_KT
         read(10,*) TAU_OUT
         read(10,*) 
         read(10,*) NLYR
         read(10,*) ALBI
        close(10)


        RUN_MODE=1
        SCALE_HEIGHT=10
        LAM_RES=1e-2
        SCALE_PRESS=1.0
        SCALE_CO2=1.0
        H2O_MODE=1
        SHIFT_TEMP=0.0
        INIT_ATM=3
        TINIT=300 !not used
        HYDRO=.true. !not used
        call CHECK_INPUT
        if(LAM_MIN.lt.4.5) write(IV,'(A3)') "VIS"
        if(LAM_MAX.gt.4.5) write(IV,'(A3)') "IRR"
       end subroutine
       
       subroutine READ_TAU(MODE,NX)
        !***************************************************************
        !                                                              *
        ! Collect all the available tau_k                              * 
        !                                                              *
        !***************************************************************
        use P_SPE,only: PPX,TTX,CCX,HHX,LABP,LABT,LABC,LABH,
     >   INDXT,INDXP,INDXC,INDXH,psca,tsca,csca,hsca,PQ
        use P_ATM,only: PP_I,TT_I,CC_I,HH_I,NSPMAX,NLMAX
        use P_INP,only: NLYR,MNAME
        use P_TRA
        implicit none
        integer,intent(in) :: MODE,NX
        character(len=200) :: cmd,fname
        integer :: i1,i2,TTV,PPV,CCV,HHV,NP
        integer :: NQ,NNN
        real*8 :: NPAR,KBOL
        real*8,dimension(NLYR) :: pint,tint,cint,hint
        real*8,dimension(NSPMAX,PQ) :: yaaaa,yaaab,yaaba,yaabb
        real*8,dimension(NSPMAX,PQ) :: yabaa,yabab,yabba,yabbb
        real*8,dimension(NSPMAX,PQ) :: ybaaa,ybaab,ybaba,ybabb
        real*8,dimension(NSPMAX,PQ) :: ybbaa,ybbab,ybbba,ybbbb,yy
        real*8,dimension(NSPMAX) :: lam
        real*8,dimension(NLMAX) :: ppuse,ttuse,ccuse,hhuse
        
        NQ=PQ
        
        if(MODE.eq.1) then
          ppuse=PP_I
          ttuse=TT_I
          ccuse=CC_I
          hhuse=HH_I
          NNN=NLYR
        elseif(MODE.eq.3) then
          ppuse=PP_T
          ttuse=TT_T
          ccuse=CC_T
          hhuse=HH_T
          NNN=21
        endif
        
        call GET_INDEXES(PPV,TTV,CCV,HHV)
        
        PPX=log10(PPX)
        ppuse=log10(ppuse)
        
        do i1=1,NNN
        ! write(*,'(A,I3.3,A,I3.3)')
c     >  !"INTERPOLATING LAYER ",i1," OF ",NNN
         do i2=1,PPV-1
          if(ppuse(i1).ge.PPX(i2).and.ppuse(i1).le.PPX(i2+1)) then
           !If within interval, done, exit
           INDXP(i1)=i2
           pint(i1)=(ppuse(i1)-PPX(i2))/(PPX(i2+1)-PPX(i2))
           psca(i1)=1
           !write(*,*) 'PP A',INDXP(i1),ppuse(i1),psca(i1),pint(i1)
           exit
          elseif(ppuse(i1).lt.PPX(1)) then
           !If lower than minimum, done, exit
           INDXP(i1)=1
           pint(i1)=0.00
           psca(i1)=ppuse(i1)/PPX(1)
           !write(*,*) 'PP B',INDXP(i1),ppuse(i1),psca(i1),pint(i1)
           exit
          elseif(ppuse(i1).gt.PPX(PPV)) then
           !If higher than maximum, done, exit
           INDXP(i1)=PPV-1
           pint(i1)=1.00
           psca(i1)=10**ppuse(i1)/10**PPX(PPV)
           !write(*,*) 'PP C',INDXP(i1),ppuse(i1),psca(i1),pint(i1)
           exit
          endif
         enddo
         do i2=1,TTV-1
          if(ttuse(i1).ge.TTX(i2).and.ttuse(i1).le.TTX(i2+1)) then
           INDXT(i1)=i2
           tint(i1)=(ttuse(i1)-TTX(i2))/(TTX(i2+1)-TTX(i2))
           tsca(i1)=1
         !write(*,*) 'TT A',INDXT(i1),ttuse(i1),tsca(i1),tint(i1)
           exit
          elseif(ttuse(i1).lt.TTX(1)) then
           INDXT(i1)=1
           tint(i1)=0.00
           tsca(i1)=ttuse(i1)/TTX(1)
          ! write(*,*) 'TT B',INDXT(i1),ttuse(i1),tsca(i1),tint(i1)
           exit
          elseif(ttuse(i1).gt.TTX(TTV)) then
           INDXT(i1)=TTV-1
           tint(i1)=1.00
           tsca(i1)=ttuse(i1)/TTX(TTV)
           !write(*,*) 'TT C',INDXT(i1),ttuse(i1),tsca(i1),tint(i1)
           exit
         endif
         enddo
         do i2=1,CCV-1
          if(ccuse(i1).ge.CCX(i2).and.ccuse(i1).le.CCX(i2+1)) then
           INDXC(i1)=i2
           cint(i1)=(ccuse(i1)-CCX(i2))/(CCX(i2+1)-CCX(i2))
           if(CCX(i2+1).eq.CCX(i2)) cint(i1)=1.0
           csca(i1)=1
           !write(*,*) 'CO2 A',INDXC(i1),ccuse(i1),csca(i1),cint(i1)
           exit
          elseif(ccuse(i1).lt.CCX(1)) then
           INDXC(i1)=1
           cint(i1)=1.00
           csca(i1)=ccuse(i1)/CCX(1)
           !write(*,*) 'CO2 B',INDXC(i1),ccuse(i1),csca(i1),cint(i1)
           exit
          elseif(ccuse(i1).gt.CCX(CCV)) then
           INDXC(i1)=CCV-1
           cint(i1)=1.00
           csca(i1)=ccuse(i1)/CCX(CCV)
           !write(*,*) 'CO2 C',INDXC(i1),ccuse(i1),csca(i1),cint(i1)
           exit
         endif
         enddo
         do i2=1,HHV-1
          if(hhuse(i1).ge.HHX(i2).and.hhuse(i1).le.HHX(i2+1)) then
           INDXH(i1)=i2
           hint(i1)=(hhuse(i1)-HHX(i2))/(HHX(i2+1)-HHX(i2))
           if(HHX(i2+1).eq.HHX(i2)) hint(i1)=1.00
           hsca(i1)=1
           !write(*,*) 'H2O A',INDXH(i1),hhuse(i1),hsca(i1),hint(i1)
           exit
          elseif(hhuse(i1).lt.HHX(1)) then
           INDXH(i1)=1
           hint(i1)=0.00
           hsca(i1)=hhuse(i1)/HHX(1)
           !write(*,*) 'H2O B',INDXH(i1),hhuse(i1),hsca(i1),hint(i1)
           exit
          elseif(hhuse(i1).gt.HHX(HHV)) then
           INDXH(i1)=HHV-1
           hint(i1)=1.00
           hsca(i1)=hhuse(i1)/HHX(HHV)
           !write(*,*) 'H2O C',INDXH(i1),hhuse(i1),hsca(i1),hint(i1)
           exit
         endif
         enddo

         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)),
     >     LABH(INDXH(i1)),yaaaa,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)),
     >     LABH(INDXH(i1)+1),yaaab,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)+1),
     >     LABH(INDXH(i1)),yaaba,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)),LABC(INDXC(i1)+1),
     >     LABH(INDXH(i1)+1),yaabb,lam,NP,NQ)

         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),yabaa,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),yabab,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),yabba,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),yabbb,lam,NP,NQ)


         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),ybaaa,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),ybaab,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),ybaba,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),ybabb,lam,NP,NQ)

         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)),ybbaa,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)),LABH(INDXH(i1)+1),ybbab,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)),ybbba,lam,NP,NQ)
         call GET_TAU(LABP(INDXP(i1)+1),LABT(INDXT(i1)+1),
     >     LABC(INDXC(i1)+1),LABH(INDXH(i1)+1),ybbbb,lam,NP,NQ)

         yy(1:NSPMAX,1:NQ)=
     >     (1-pint(i1))*(1-tint(i1))*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yaaaa(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yaaab(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yaaba(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*(1-tint(i1))*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yaabb(1:NSPMAX,1:NQ)+
     >
     >     (1-pint(i1))*tint(i1)*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yabaa(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yabab(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yabba(1:NSPMAX,1:NQ)+
     >     (1-pint(i1))*tint(i1)*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  yabbb(1:NSPMAX,1:NQ)+
     >
     >     pint(i1)*(1-tint(i1))*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybaaa(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybaab(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybaba(1:NSPMAX,1:NQ)+
     >     pint(i1)*(1-tint(i1))*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybabb(1:NSPMAX,1:NQ)+
     >
     >     pint(i1)*tint(i1)*(1-cint(i1))*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybbaa(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*(1-cint(i1))*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybbab(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*cint(i1)*(1-hint(i1))*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybbba(1:NSPMAX,1:NQ)+
     >     pint(i1)*tint(i1)*cint(i1)*hint(i1)*
     >       csca(i1)*hsca(i1)*psca(i1)*tsca(i1)*  ybbbb(1:NSPMAX,1:NQ)
              
         if(MODE.eq.1) then
          write(fname,'(A,I3.3,A)') "OUT_EX/",i1,"_k_GEN.dat"
         elseif(MODE.eq.3) then
          write(fname,'(A,2(I3.3),A)') "OUT_EX/",NX,i1,"_k_GEN.dat"
         endif
         
         open(unit=12,file=trim(fname),status="replace")
         do i2=1,NP-1
          if(isnan(yy(i2,1))) yy(i2,1:NQ)=yy(i2-1,1:NQ)
          write(12,*) lam(i2),yy(i2,1:NQ)
         enddo
         close(12)
        enddo
        
        PPX=10**PPX
        ppuse=10**ppuse
        !Write out psca and tsca
        open(unit=10,file="OUT_EX/scale.dat",status="replace")
         do i1=1,NLYR
          write(10,*) psca(i1),tsca(i1)
         enddo
        close(10)

       end subroutine
       
       subroutine GET_TAU(j,k,l,m,yy,lam,NP,NQ)
        use P_ATM,only: NSPMAX
        use P_SPE,only: PPX,TTX
        implicit none
        integer,intent(in) :: j,k,l,m,NQ
        real*8 :: tt(NQ),NPAR
        integer :: i1,CHECK,IND
        real*8,intent(out) :: yy(NSPMAX,NQ)
        real*8,intent(out) :: lam(NSPMAX)
        integer,intent(out) :: NP
        character(len=100) :: fname2,fname,str
        
        call CHECK_0(j,k,l,m,CHECK,str)
c        write(fname,'(A19,4(I2.2),A)') "INPUT/TAU/OLD_GRID/",
c     >    j,k,l,m,".dat"
        if(CHECK.eq.1) then
         write(*,*) trim(str)
         yy=0
        else
         write(fname,'(A,4(I2.2),A)') 
     +     "ExArt_MR/K_ABS_MARS/MODEL",j,k,l,m,"_MR.out"
         IND=INDEX(trim(fname),"1308_MR.out")
         if(IND.gt.0) then
          write(fname,'(A,2(I2.2),A)') 
     +     "ExArt_MR/K_ABS_MARS/MODEL",j,k,"1307_MR.out"
         endif
         open(unit=10,file=trim(fname),status="old")
         !write(*,*) trim(fname)
         i1=0
         do
          i1=i1+1
          read(10,*,END=1001) lam(i1),tt(1:NQ)
          !absorption coefficient
          yy(i1,1:NQ)=tt(1:NQ)
         enddo
1001     close(10)
         NP=i1
        endif
        return
       end subroutine
       
       subroutine GET_INDEXES(PPV,TTV,CCV,HHV)
        use P_SPE,only: PPX,TTX,HHX,CCX,LABP,LABT,LABH,LABC
        implicit none
        integer,intent(out) :: PPV,TTV,CCV,HHV
        integer :: i1,i2
        character(len=50) :: fname
        
        write(fname,'(A)') "ExArt_MR/CALC_K.dat"
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
       
       subroutine READ_SPECTRA
        !***************************************************************
        !                                                              *
        ! Read the stellar spectrum specified by the user, then:        * 
        !  1-Calculate Stellar Luminosity and Flux                     *
        !  2-Interpolate with the user specified resolution            *
        !  3-[if INTERP_SP is Y] Write spectrum in OUTPUT/INT_SPEC     *
        !                                                              *
        !***************************************************************
        use P_INP,only: INPUT_SP,PL_DIST,LAM_MIN,LAM_MAX,LAM_RES,
     +   INTERP_SP,ND_VAR,FLUX_SCALE,RSTAR
        use P_SPE,only: NP
        use P_ATM,only: LAM_EXT,FLUX_EXT,LAM_INTRP,FLUX_INTRP
        use P_CONST,only: RSUN,AU,PI
        implicit none
        integer,parameter :: LMAX=100000
        character(len=10),dimension(18) :: STSP_NAME
        character(len=100) :: FNAME
        integer :: i1,NP2
        real*8,dimension(LMAX) :: LAM,FLUX,BBFLUX
        real*8,dimension(200) :: YV
        real*8 :: HEIGHT,BASE,IX,LSTAR,F_PL,LLMIN,LLMAX,LLRES
        real*8 :: YV2,DLAM,LAML,SFACL,SFACF,LL
        
        !write(*,'(A)') " SPECTRAL INFOS"
        
        LLMIN=LAM_MIN
        LLMAX=LAM_MAX
        LLRES=LAM_RES
        
        call GET_YV(YV)
        
       end subroutine
       
       subroutine READ_ATM
        !***************************************************************
        !                                                              *
        ! Read the atmosphere specified by the user [INPUT03]          * 
        !  Quota (km), Pressure (pa), Temp (K), H2O (VMR), O3 (VMR)    *
        !                                                              *
        !***************************************************************
        use P_INP,only: NLYR,INIT_ATM,MNAME,SCALE_PRESS,
     >   SCALE_CO2,SHIFT_TEMP,H2O_MODE
        use P_ATM,only: TT_I,PP_I,CC_I,ZZ_I,HH_I,O3_I
        implicit none
        integer :: i1
        real*8 :: sig
        real*8,dimension(NLYR) :: PPT,TTT
        
        if(INIT_ATM.eq.1) then
         write(*,*) "NOT READY YET, CHOOSE INIT_ATM 2"
        elseif(INIT_ATM.eq.2) then
         open(unit=10,file="INPUT/ATM_IN_T1",status="old")
         do i1=1,NLYR
          read(10,*) ZZ_I(i1),PP_I(i1),TT_I(i1),HH_I(i1),O3_I(i1)
         enddo
         close(10)
        elseif(INIT_ATM.eq.3) then
         open(unit=10,file="INP_TH/atms_ex.inp",status="old")
         HH_I=0.0
         O3_I=0.0
         !Pressures, Temperatures and optical depths are defined at layer interfaces
         do i1=1,NLYR
          read(10,*) PP_I(i1),TT_I(i1),CC_I(i1),HH_I(i1)
         enddo
         !From atm to pas
         PP_I=PP_I*SCALE_PRESS
         CC_I=CC_I*SCALE_CO2
         TT_I=TT_I+SHIFT_TEMP
         if(H2O_MODE.eq.1) then
          !call HUM_PROF(PP_I,TT_I,HH_I,CC_I)
         else
          HH_I=0.0
         endif
         
         
         !open(unit=20,file="OUT_TH/atms.dat",status="replace")
         !do i1=1,NLYR
         ! write(20,*) PP_I(i1),TT_I(i1),CC_I(i1),HH_I(i1)
         !enddo
         !close(20)
         close(10)
        endif

       end subroutine
       
       subroutine HUM_PROF(PP,TT,HH,CCC)
        use P_INP,only: NLYR
        implicit none
        integer :: i1
        real*8 :: QQ,AA,BB,CC,psat,tot
        real*8,dimension(NLYR),intent(in) :: PP,TT
        real*8,dimension(NLYR),intent(inout) :: CCC
        real*8,dimension(NLYR),intent(out) :: HH  
        
        AA=8.07131
        BB=1730.63
        CC=233.426
        
        do i1=1,NLYR
          QQ=PP(i1)/PP(NLYR)
	      psat=133.32*10**(AA-BB/(CC+TT(i1)-273.15))
	      HH(i1)=0.6*((QQ-0.02)/(1-0.02))
          HH(i1)=HH(i1)*psat/PP(i1)
	      if(HH(i1).lt.1e-10) HH(i1)=0.0
          if((CCC(i1)+HH(i1)).gt.1.00) CCC(i1)=1.0-HH(i1)
        enddo
        return
       end subroutine
       
       subroutine CHECK_0(j,k,l,m,CHECK,str)
        implicit none
        integer,intent(in) :: j,k,l,m
        integer,intent(out) :: CHECK
        character(len=100),intent(out) :: str
        
        CHECK=0
        if(j.eq.0) then
         write(str,'(A)') "WARNING! FOUND 0 in PRESSURE INDEX"
        CHECK=1
        elseif(k.eq.0) then
         write(str,'(A)') "WARNING! FOUND 0 in TEMPERATURE INDEX"
        CHECK=1
        elseif(l.eq.0) then
         write(str,'(A)') "WARNING! FOUND 0 in CO2 INDEX"
        CHECK=1
        elseif(m.eq.0) then
         write(str,'(A)') "WARNING! FOUND 0 in H2O INDEX"
        CHECK=1
        endif
        return
       end subroutine
       
       
