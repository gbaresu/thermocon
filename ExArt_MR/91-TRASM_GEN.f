      subroutine TR_MODE
       use P_ATM,only: NSPMAX,ZZ_I
       use P_INP,only: NLYR,RPL,MPL
       use P_TRA
       use P_CONST,only: MEARTH,REARTH,RSUN
       implicit none
       integer :: i1,i2,NLAM
       real*8,dimension(NTT) :: alt
       real*8 :: BASE,HEIGHT,RPU,MPU,III,yy
       real*8,dimension(NSPMAX) :: lam
       real*8,dimension(NTT) :: zz,FUNC,XTAB
       
       RPU=RPL*REARTH/1e5
       MPU=MPL*MEARTH
       zz(1)=0.0
       zz(2)=1.0
       zz(3)=2.0
       zz(4)=3.0
       zz(5)=5.0
       zz(6)=10.0
       zz(7)=20.0
       zz(8)=50.0
       zz(9)=ZZ_I(2)/1e5

       alt=(/zz(1),zz(2),zz(3),zz(4),zz(5),zz(6),zz(7),zz(8),zz(9)/)

       call GET_WEIGHTS

       do i1=1,NTT
        write(*,*) alt(i1)
        call BUILD_LOS(alt(i1))
        
        call READ_TAU(3,i1)
        write(*,*)
        call INTEGRAL_1(i1,lam,NLAM)
       enddo
       
       do i2=1,NLAM
        III=0
        do i1=1,NTT
         FUNC(i1)=(RPU+alt(i1))*(1.-exp(-TAU_LOS(i2,i1)))
         XTAB(i1)=alt(i1)
        enddo
        call plint(NTT,FUNC,XTAB,NTT,alt(1),alt(NTT),yy)
        write(1000,*) lam(i2),
     >   (RPU**2)/(RSUN/1e5)**2,
     >   (RPU**2+2*yy)/(RSUN/1e5)**2,
     >   (RPU**2+yy)**0.5/RPU,
     >   (RPU**2+2*zz(9)*RPU+zz(9)**2)/(RSUN/1e5)**2
       enddo
      end
      
      subroutine INTEGRAL_1(IX,lam,NLAM)
       use P_TRA
       use P_CONST,only: KBOL
       use P_ATM,only: NSPMAX,CC_I,HH_I
       use P_SPE,only: WQ,PQ
       implicit none
       integer,intent(in) :: IX
       integer,intent(out) :: NLAM
       integer :: i1,i2,i3
       real*8 :: BASE,HEIGHT,TTT,III,TAU_RAY,N2,H2O
       real*8,dimension(21) :: NNTOT
       character(len=50) :: fname
       real*8,dimension(NSPMAX,32) :: K_GAS
       real*8,dimension(NSPMAX,21) :: TAU_GAS
       real*8,intent(out),dimension(NSPMAX) :: lam
       
       N2=1-CC_I(1)
       H2O=HH_I(1)

       
       do i1=1,21
        NNTOT(i1)=PP_T(i1)/(KBOL/1e7)/TT_T(i1)/1e6
        write(fname,'(A,2(I3.3),A)') 
     >   "OUTPUT/",IX,i1,"_k_GEN.dat"
        open(unit=10,file=trim(fname),status="old")
        i2=0
        do
         TTT=0
         i2=i2+1
         read(10,*,END=1000) lam(i2),K_GAS(i2,1:32)
         call RAYL_GEN(lam(i2),N2,CC_I(1),H2O,TAU_RAY)
         !K GAS in cm2/mol as well as K_RAYL
         K_GAS(i2,1:32)=K_GAS(i2,1:32)+TAU_RAY
         !write(*,*) NNTOT(i1),PP_T(i1)/TT_T(i1)/KBOL
         !stop
         do i3=1,PQ
          TTT=TTT+WQ(i3)*exp(-1.0*K_GAS(i2,i3)*NNTOT(i1))
         enddo
         TAU_GAS(i2,i1)=-log(TTT)
        enddo
1000    close(10)
       enddo
       NLAM=i2

       do i2=1,NLAM
       III=0.0
        do i1=1,20
         BASE=(ZZ_T(i1+1)-ZZ_T(i1))*1e3
         HEIGHT=0.5*(TAU_GAS(i2,i1+1)+TAU_GAS(i2,i1))
         if(BASE.lt.0) stop
         III=III+BASE*HEIGHT
        enddo 
        TAU_LOS(i2,IX)=III
        write(300+IX,*) lam(i2),III
       enddo
       
      
      end subroutine
      
      subroutine BUILD_LOS(zz)
       use P_INP,only: RPL,MPL,NLYR
       use P_ATM,only: ZZ_I,PP_I,TT_I,CC_I,HH_I
       use P_TRA       
       use P_CONST,only: MEARTH,REARTH
       implicit none
       integer :: i1,ind
       integer,parameter :: NPX=20
       real*8,intent(in) :: zz
       real*8 :: zx,ll,RPU,MPU,xxmax,llmax,Hmax
       real*8 :: pval,tval,cval,hval,yval2,xnmax,xnmin
       real*8,dimension(NLYR) :: zzuse,ppuse,ttuse,ccuse,hhuse
       real*8,dimension(NPX) :: lll,ppl,ttl,ccl,hhl
       real*8,dimension(NPX) :: xx 
       
       Hmax=ZZ_I(1)/1e5
       
       RPU=RPL*REARTH/1e5
       MPU=MPL*MEARTH
       
       xx=0.0

       xxmax=acos((RPU+zz)/(RPU+Hmax))
       llmax=(RPU+Hmax)*sin(xxmax)
       
       do i1=NLYR,1,-1
        zzuse(i1)=ZZ_I(NLYR-i1+1)/1e5
        ppuse(i1)=PP_I(NLYR-i1+1)
        ttuse(i1)=TT_I(NLYR-i1+1)
        ccuse(i1)=CC_I(NLYR-i1+1)
        hhuse(i1)=HH_I(NLYR-i1+1)
       enddo
       lll=0.0
       i1=0
       xnmax=0.0
       xnmin=1e10
       do i1=1,NPX
        xx(i1)=10**(real(i1-1)/20.)-1.0
        if(xx(i1).gt.xnmax) xnmax=xx(i1)
        if(xx(i1).lt.xnmin) xnmin=xx(i1)
       enddo
       xx=(xx)/xnmax*0.99*xxmax

       do i1=1,NPX
        zx=(RPU+zz)/cos(xx(i1))-RPU
        ll=(RPU+zx)*sin(xx(i1))
        ind=NPX-i1+1
        lll(ind)=llmax-ll
        !Pressure
        call spline_linear_val(NLYR,zzuse,ppuse,zx,pval,yval2)
        ppl(ind)=pval
        !Temperature
        call spline_linear_val(NLYR,zzuse,ttuse,zx,tval,yval2)
        ttl(ind)=tval
        !CO2
        call spline_linear_val(NLYR,zzuse,ccuse,zx,cval,yval2)
        ccl(ind)=cval
        !H2O
        call spline_linear_val(NLYR,zzuse,hhuse,zx,hval,yval2)
        hhl(ind)=hval
       enddo
       
       do i1=1,NPX
        ZZ_T(i1)=lll(i1)
        PP_T(i1)=ppl(i1)
        TT_T(i1)=ttl(i1)
        CC_T(i1)=ccl(i1)
        HH_T(i1)=hhl(i1)
       enddo
       
       do i1=NPX+1,2*NPX-1
        ZZ_T(i1)=ZZ_T(i1-1)+lll(2*NPX-i1)
        PP_T(i1)=ppl(2*NPX-i1)
        TT_T(i1)=ttl(2*NPX-i1)
        CC_T(i1)=ccl(2*NPX-i1)
        HH_T(i1)=hhl(2*NPX-i1)
       enddo
       
       do i1=1,2*NPX-1
        write(*,'(I2,5(f13.4))') i1,ZZ_T(i1),PP_T(i1),TT_T(i1),
     >   CC_T(i1),HH_T(i1)
       enddo
       end subroutine
