      subroutine INIT_COND_MACKAY
       use PAR_ATM,only: temp,tav,fdown,nz,dt,pre,fup
       use PAR_INP,only: int_time
       use COST,only: water,Cpw,rhow,dh,sig,Cpa,ga
       implicit none
       real*8 :: tt,a,b,c,d,e,dp
       water=Cpw*rhow*dh
       dt=1/water*(fdown(nz)-sig*temp(nz)**4)*int_time
       tav(nz)=tav(nz)+dt
      end subroutine

      subroutine INIT_COND_MAN
       use PAR_ATM,only: temp,tav,fdown,nz,dt,pre,fup,told,pav,dTdt
       use PAR_INP,only: int_time,fact
       use COST,only: sig,Cpa,ga
       use PAR_GEN,only: check,foldg,foldt
       implicit none
       real*8 :: tt,a,b,c,d,e,dp,tgr,ttoa,ttn,hh
       
       tgr=((fdown(nz))/sig)**0.25
       ttoa=((fdown(1))/sig)**0.25
       dp=pre(nz)-pre(nz-1)
       
       a=sig*int_time
       b=0
       c=0
       d=Cpa/ga*dp
       e=sig*tgr**4*int_time+Cpa/ga*dp*temp(nz)
       call solve_quartic(a,b,c,d,e,temp(nz),ttn)
       
       dp=pre(2)-pre(1)
       a=sig*int_time
       b=0
       c=0
       d=Cpa/ga*dp
       e=sig*ttoa**4*int_time+Cpa/ga*dp*temp(1)
       call solve_quartic(a,b,c,d,e,temp(1),tt)
       
       temp(1)=temp(2)!tt
       temp(nz)=ttn
       
       write(*,*) "INITIAL CONDITIONS: teq [GR] :",tgr
       write(*,*) "INITIAL CONDITIONS: teq [TOA]:",ttoa
      end subroutine

      subroutine DEF_INIT
       use PAR_ATM,only: temp,pav,pre,nz,tav,quo,told
       implicit none
       integer :: i1
       real :: a1,a2,psum
       told=temp
       do i1=1,nz-1
        pav(i1)=(pre(i1)+pre(i1+1))/2        
       enddo
      end subroutine

      subroutine READ_SBDART(it)
       use PAR_ATM,only: nz,fdown,fup,dTdt,heat,pre,temp,told,Cpt,pav,
     >   quo
       use PAR_INP,only: fact,chh
       use COST,only: ga,Cpa
       use PAR_SYS,only: DPL
       implicit none
       integer,intent(in) :: it
       integer :: i1
       real*8 :: fdir,rad,dp,df,dfin,dfout,pp,quos,hh,hhp,hhh,dTnz,xx
       real*8,dimension(100) :: dTdtt,dTave
       character(len=50) :: cmd
       open(unit=10,file="F_INT_TOT.txt",status="old")
       do i1=1,nz
        read(10,*) pp,fdown(i1),xx,fup(i1)
        !dTdtt(i1-1)=hhh
       enddo

       open(unit=600,file="fort.600",status="replace")
       fdown=fdown!/fact/Dpl**2
       fup=fup!/fact/Dpl**2
       if(chh.eq.1) dTdtt=dTdtt/Dpl**2/fact
       do i1=1,nz-1
        dp=pre(i1+1)-pre(i1)
        dfin =fdown(i1)-fup(i1)
        dfout=fdown(i1+1)-fup(i1+1)
        df=dfin-dfout
        heat(i1)=ga/Cpt(i1)*df/dp
        write(600,*) pre(i1),heat(i1),ga,Cpt(i1)
       enddo
       close(600)

       do i1=1,nz-1
        if(chh.eq.1) then
         hh=dTdtt(i1)/24/3600
         hhp=dTdtt(i1)/24/3600
        elseif(chh.eq.2) then 
         hh=heat(i1)
         hhp=heat(i1)
        endif 
        dTave(i1)=hh*pre(i1)+hhp*pre(i1+1)
        dTave(i1)=dTave(i1)/(pre(i1)+pre(i1+1))
       enddo
       close(10)
       do i1=1,nz
        if(chh.eq.1) dTdt(i1)=dTave(i1)
        if(chh.eq.2) dTdt(i1)=heat(i1)
       enddo
       
      write(cmd,'(A,I0,A4)')"cp F_INT_TOT.txt OUT_CONRAD/o",it,".dat"
      call SYSTEM(cmd)
      write(cmd,'(A,I0,A4)')"cp INPUT/atms.dat OUT_ATMS/a",it,".dat"
       call SYSTEM(cmd)
      end subroutine
      
      subroutine ERASE_OLD
       character(len=50) :: cmd
       cmd="rm OUT_CONRAD/*"
       call SYSTEM(cmd)
       cmd="rm OUT_ATMS/*"
       call SYSTEM(cmd)
      end subroutine


      subroutine GEN_ATM
       use PAR_ATM,only: quo,temp,pre,nz,n2,co2,h2o
       use PAR_SYS, only: RPL,MPL
       use COST,only: ga, Ggrav
       real*8 :: HSCA
       
       mu=44.*co2(nz)+28.*n2(nz)+18.*h2o(nz)
       HSCA=1.38e-23*temp(nz)/1.67e-27/44.0/3.72/1e3
       ga=0.01*Ggrav*Mpl/Rpl**2
       write(*,*) ga, Mpl, Rpl
       stop
       do i1=1,nz
        quo(i1)=HSCA*log(pre(nz)/pre(i1))
        !write(*,*) quo(i1),pre(i1),temp(i1),co2(i1),h2o(i1)
       enddo

      end subroutine

      subroutine GET_CP
       use PAR_ATM,only: quo,temp,cpt,hum,pre,nz,h2o
       use COST,only: Rgas,Cpa
       implicit none
       integer i1
       real*8,dimension(100) :: h2oloc
       real*8 :: lheat,dair,q1,q2,dq,dT
       h2oloc=h2o
       h2oloc=h2oloc/1e3
       do i1=1,nz
        lheat=2510-2.38*(temp(i1)-273.15)
        dair=pre(i1)/Rgas/temp(i1)
        q1=h2oloc(i1)/dair
        q2=h2oloc(i1+1)/dair
        dq=q2-q1
        dT=temp(i1+1)-temp(i1)
        if(dT.eq.0) then 
         Cpt(i1)=Cpa
        else
         Cpt(i1)=Cpa*(1+lheat/Cpa*(dq/dT))
        endif
        !write(*,*) Cpt(i1)
       enddo
      end subroutine

      subroutine GET_PINT
      !Prepares pressure profile MACKAY style
      use PAR_ATM,only: nz,premc,quo,pre
      implicit none
      integer :: i1
      real*8,dimension(100) :: sig
      real*8 :: ds,dp,p
      
      ds=1./nz
      sig(1)=1.0/nz
      do i1=2,nz-1
       sig(i1)=sig(i1-1)+ds
      enddo
      do i1=1,nz-1
       dp=6*(sig(i1)-sig(i1)**2)*ds
       p=sig(i1)**2*(3-2*sig(i1))
       !write(*,*) quo(i1),p*1010
      premc(i1)=p*pre(nz)
      enddo
      premc(nz)=pre(nz)
      !do i1=1,nz
      ! write(*,*) premc(i1)
      !enddo
      end subroutine 

      subroutine MAKE_BB
      use PAR_ATM, only: INTEGRAL
      use PAR_SYS, only: TSTAR
      use COST, only: Ssun
      use PAR_INP, only: sza, SOLFAC
      implicit none
      integer,parameter :: NMAX=1000000
      integer,parameter :: NB=15
      real*8 :: iii,hpl,cl,kbol,nu,kk,sig,ITOT,ITOT2,IX,ITOTCH
      integer :: i1,i2,NN,ixx
      real*8 :: BASE,HEIGHT,yval
      real*8,dimension(NMAX) :: kkk,bbb,kkx,bbx
      real*8,dimension(NB) :: band,II
      character(len=100) :: cmd
      
      band(15)=1/10.*1e8
      band(14)=1/166.*1e8
      band(13)=1/416.*1e8
      band(12)=1/625.*1e8
      band(11)=1/710.*1e8
      band(10)=1/833.*1e8
      band(9)=1/1250.*1e8
      band(8)=1/2222.22*1e8
      band(7)=1/3087.37*1e8
      band(6)=1/4030.63*1e8
      band(5)=1/5370.57*1e8
      band(4)=1/7651.11*1e8
      band(3)=1/12500.0*1e8
      band(2)=1/25000.0*1e8
      band(1)=1/41666.67*1e8

      hpl=6.62e-27
      kbol=1.38e-16
      sig=5.67e-5
      cl=3e10
      open(456, file="temp", status="replace")
      !WRITE OUT angstrom  erg/cm2/s/ang
      do i1=10,50000
       kk=real(i1)
       nu=kk*cl
       iii=(2*hpl*nu**3/cl**2)/(exp((hpl*nu)/(kbol*Tstar))-1)
       write(456,*) (cl/nu)*1e8,iii*nu**2/cl/1e8
      enddo
      close(456)
      
      write(cmd,'(A)') "tac temp > bb.dat"
      call system(trim(cmd))
      
      ITOT=sig*Tstar**4/3.14
      open(unit=10,file="bb.dat",status="old")
      i1=0
      do
       i1=i1+1
       read(10,*,END=199) kkk(i1),bbb(i1)
       kkk(i1)=kkk(i1)
      enddo
199   close(10)
      NN=i1
      II=0
      do i2=1,NB-1
       do i1=1,NN-1
        if(kkk(i1).ge.band(i2).and.kkk(i1).le.band(i2+1)) then
         BASE=kkk(i1+1)-kkk(i1)
         HEIGHT=0.5*(bbb(i1+1)+bbb(i1))
         II(i2)=II(i2)+BASE*HEIGHT
        endif
       enddo
!       write(*,*) II(i2)/ITOT
      enddo
      
      call plint(NN-1,bbb,kkk,NN-1,band(1),band(NB),yval)
      ITOT2=yval
!      write(*,*) ITOT,ITOT2,ITOT/ITOT2
      ITOTCH=0.d0
      do i2=1,NB-1
       ixx=0
       kkx=0
       bbx=0
       do i1=1,NN-1
        if(kkk(i1).ge.band(i2).and.kkk(i1).le.band(i2+1)) then
         ixx=ixx+1
         kkx(ixx)=kkk(i1)
         bbx(ixx)=bbb(i1)
        endif
       enddo
       call plint(NN,bbx,kkx,ixx,band(i2),band(i2+1),yval)
       INTEGRAL(i2)=yval/ITOT*sza*SOLFAC*Ssun
       ITOTCH=ITOTCH+INTEGRAL(i2)
      enddo
      IX=0
      do i1=1,NB-1
       IX=IX+II(i1)
      enddo
      end subroutine 

	
