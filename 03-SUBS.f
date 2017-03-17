      subroutine INIT_MANABE
       use PAR_ATM,only: T0,T1,PRE_USE1,NL,TGR_EQ,FL_EX_DWN_MAN,TGR_MAN
       use PAR_CODE,only: I_TIME
       use PAR_PLA,only: Cp,GACC
       use CONST,only: SIGBOL
       implicit none
       integer :: i1
       real*8 :: COEFF_A,COEFF_B,DPN
       real*8 :: a,b,c,d,e,init,tt
       
       TGR_EQ=(FL_EX_DWN_MAN(NL)/SIGBOL)**0.25
       
       DPN=PRE_USE1(NL)-PRE_USE1(NL-1)
       COEFF_A=SIGBOL*I_TIME
       COEFF_B=Cp/GACC*DPN
       
       a=COEFF_A
       b=0.0
       c=0.0
       d=COEFF_B
       e=COEFF_A*TGR_EQ**4+COEFF_B*T0(NL)
       init=T0(NL)
       
       call solve_quartic(a,b,c,d,e,init,tt)
       
       T1=T0
       T1(NL)=tt
       TGR_MAN=tt
       
      end subroutine
      
      subroutine PREP_ATM
       use PAR_ATM,only: H2O,CO2,N2,Tinit,PRE_USE1,NL,HSCA,WAT_HUM,
     >  ZZZ
       use PAR_INP,only: CLRMOD_INP,CLR_INP
       use PAR_PLA,only: CS,GACC,MPL,RPL,Cp,Cv,GAMMA,MU,CLR,
     >  CLR_WET,CLR_DRY
       use PAR_CODE,only: WAT_MODE
       use CONST,only: GGRAV,CvN2,CvCO2,mH,mCO2,MN2,KBOL,Rd,Rw,Lvap
       implicit none
       integer :: i1
       real*8 :: QQ,HUM,PSAT,AA,BB,CC,cpN2,cpCO2,cpH2O,WS
       
       !Calculate gravitational acceleration
       GACC=GGRAV*MPL/RPL**2
       !Calculate average Cp and Cv
       call fit_cpN2(Tinit(NL),cpN2)
       call fit_cpCO2(Tinit(NL),cpCO2)
       call fit_cpH2O(Tinit(NL),cpH2O)
       
       Cp=CpN2*N2(NL)+CpCO2*CO2(NL)
       Cv=CvN2*N2(1)+CvCO2*CO2(1)
       GAMMA=Cp/Cv
       ! TODO include temperature dependence on these values?
       !Calculate mean molecular mass
       MU=N2(1)*mN2+CO2(1)*mCO2
       !Calculate sound speed
       CS=sqrt(GAMMA*KBOL*Tinit(NL)/MU/mH)
       !Calculate Scale Height
       HSCA=CS**2/GACC
       !Calculate humidity [if switch is 1]
       AA=8.07131
       BB=1730.63
       CC=233.426
       if(WAT_MODE.eq.2) then
        write(*,*) H2O
        stop
       elseif(WAT_MODE.eq.1) then
        do i1=NL,1,-1
         QQ=PRE_USE1(i1)/PRE_USE1(NL)
         HUM=WAT_HUM*(QQ-0.02)/(1.00-0.02)
         PSAT=133.32*10**(AA-BB/(CC+Tinit(i1)-273.15))
         H2O(i1)=HUM*PSAT/(PRE_USE1(i1)/1e1)
         if(QQ.lt.0.02) H2O(i1)=H2O(i1+1)
         if(H2O(i1).gt.1.00) H2O(i1)=1.00
        enddo
       elseif(WAT_MODE.eq.0) then
        H2O=0.0
       endif
       !Calculate Critical Lapse Rate
       if(CLRMOD_INP.eq.0) then
        WS=18.00/MU*H2O(NL)
        CLR_DRY=GACC/Cp*1e5
        CLR_WET=CLR_DRY*(1+(Lvap*WS)/(Rd*Tinit(NL)))/
     +   (1+(WS*Lvap**2)/(Cp*Rw*Tinit(NL)**2))
        CLR=WAT_HUM*CLR_WET+(1-WAT_HUM)*CLR_DRY
       elseif(CLRMOD_INP.eq.1) then
        CLR=CLR_INP
       endif
       !Calculate heights
       do i1=1,NL
        ZZZ(i1)=HSCA*log(PRE_USE1(NL)/PRE_USE1(i1))
       enddo

       return
      end subroutine
      
      
      subroutine PREP_ATM_LOOP
       use PAR_ATM,only: H2O,CO2,N2,T2,PRE_USE1,NL,HSCA,WAT_HUM,
     >   ZZZ
       use PAR_PLA,only: CS,GACC,MPL,RPL,Cp,Cv,GAMMA,MU,CLR,
     >  CLR_WET,CLR_DRY
       use PAR_INP,only: CLRMOD_INP,CLR_INP
       use PAR_CODE,only: WAT_MODE
       use CONST,only: GGRAV,CvN2,CvCO2,mH,mCO2,MN2,KBOL,
     >  Rd,Rw,Lvap
       implicit none
       integer :: i1
       real*8 :: QQ,HUM,PSAT,AA,BB,CC,cpN2,cpCO2,cpH2O,WS
       
       !Calculate average Cp and Cv
       call fit_cpN2(T2(NL),cpN2)
       call fit_cpCO2(T2(NL),cpCO2)
       call fit_cpH2O(T2(NL),cpH2O)
       
       Cp=CpN2*N2(NL)+CpCO2*CO2(NL)
       Cv=CvN2*N2(1)+CvCO2*CO2(1)
       GAMMA=Cp/Cv
       !Calculate sound speed
       CS=sqrt(GAMMA*KBOL*T2(NL)/MU/mH)
       !Calculate Scale Height
       HSCA=CS**2/GACC
       !Calculate humidity [if switch is 1]
       AA=8.07131
       BB=1730.63
       CC=233.426
       if(WAT_MODE.eq.2) then
        write(*,*) H2O
        stop
       elseif(WAT_MODE.eq.1) then
        do i1=NL,1,-1
         QQ=PRE_USE1(i1)/PRE_USE1(NL)
         HUM=WAT_HUM*(QQ-0.02)/(1-0.02)
         PSAT=133.32*10**(AA-BB/(CC+T2(i1)-273.15))
         H2O(i1)=HUM*PSAT/(PRE_USE1(i1)/10.0)
         !write(*,*) i1,T2(i1),PSAT/(PRE_USE1(i1)/10.0)
         if(QQ.lt.0.02) H2O(i1)=H2O(i1+1)
         if(H2O(i1).gt.1.00) H2O(i1)=1.00
        enddo
       elseif(WAT_MODE.eq.0) then
        H2O=0.0
       endif
       !Calculate Critical Lapse Rate
       if(CLRMOD_INP.eq.0) then
        WS=18.00/MU*H2O(NL)
        CLR_DRY=GACC/Cp*1e5
        CLR_WET=CLR_DRY*(1+(Lvap*WS)/(Rd*T2(NL)))/
     +   (1+(WS*Lvap**2)/(Cp*Rw*T2(NL)**2))
        CLR=WAT_HUM*CLR_WET+(1-WAT_HUM)*CLR_DRY
       elseif(CLRMOD_INP.eq.1) then
        CLR=CLR_INP
       endif
       !Calculate heights
       do i1=1,NL
        ZZZ(i1)=HSCA*log(PRE_USE1(NL)/PRE_USE1(i1))
       enddo
       
       return
      end subroutine
      

      subroutine solve_quartic(a,b,c,d,e,init,tt)
       !Solves 4 degree equation
       implicit none 
       real*8,intent(in) :: a,b,c,d,e,init
       real*8,intent(out) :: tt
       real*8 :: q1,tb,tx,sm,pm
       integer :: i1
       
       q1=0.5
       tb=init
       do i1=1,999
        tx=tb-q1
        pm=a*tx**4+d*tx
        sm=e
        if(abs(pm-sm).lt.0.001) then
         tt=tb
         exit
        else if(pm.gt.sm) then
         tb=tx
         if(q1.lt.0) then
          q1=-q1/2
         else
          q1=q1
         endif
        else if(pm.lt.sm) then
         tb=tx
         if(q1.gt.0) then 
          q1=-q1/2
         else
          q1=q1
         endif
        endif
       enddo
      end subroutine
      
      subroutine fit_cpN2(TT,yv)
       implicit none
       integer,parameter :: NX=12
       real*8,intent(in) :: TT
       real*8,intent(out) :: yv
       real*8,dimension(NX) :: TX,cp
       real*8 :: yvp
       
       TX=(/175,200,225,250,275,300,325,350,375,400,450,500/)
       cp=(/1.039,1.039,1.039,1.039,1.039,1.040,1.040,1.041,
     +   1.042,1.044,1.049,1.056/) 
     
       call spline_linear_val(NX,TX,cp,TT,yv,yvp)
       !From mks to cgs
       yv=yv*1e4*1e3
       return 
      end
      
      subroutine fit_cpCO2(TT,yv)
       implicit none
       integer,parameter :: NX=12
       real*8,intent(in) :: TT
       real*8,intent(out) :: yv
       real*8,dimension(NX) :: TX,cp
       real*8 :: yvp
       
       TX=(/175,200,225,250,275,300,325,350,375,400,450,500/)
       cp=(/0.709,0.735,0.763,0.791,0.819,0.846,0.871,0.895,0.918,
     +  0.939,0.978,1.014/)

       call spline_linear_val(NX,TX,cp,TT,yv,yvp)
       !From mks to cgs
       yv=yv*1e4*1e3
       return 
      end
      
      subroutine fit_cpH2O(TT,yv)
       implicit none
       integer,parameter :: NX=16
       real*8,intent(in) :: TT
       real*8,intent(out) :: yv
       real*8,dimension(NX) :: TX,cp
       real*8 :: yvp
       
       TX=(/273,277,278,283,288,293,298,303,308,313,318,323,328,333,338,
     +  343/)
       cp=(/4.217,4.205,4.202,4.192,4.185,4.182,4.180,4.178,4.178,4.179,
     +   4.181,4.182,4.183,4.185,4.188,4.191/) 
     
       call spline_linear_val(NX,TX,cp,TT,yv,yvp)
       !From mks to cgs
       yv=yv*1e4*1e3
       return 
      end
      
      subroutine spline_linear_val(ndata,tdata,ydata,tval,yval,ypval)
      !
      !*******************************************************************************
      !
      !  Parameters:
      !
      !    Input, integer NDATA, the number of data points defining the spline.
      !
      !    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
      !    and dependent variables at the data points.  The values of TDATA should
      !    be distinct and increasing.
      !
      !    Input, real TVAL, the point at which the spline is to be evaluated.
      !
      !    Output, real YVAL, YPVAL, the value of the spline and its first
      !    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
      !    equal to TDATA(I) for some I.
      !
      implicit none
      !
      integer ndata
      !
      integer :: left,right
      real*8 :: tdata(ndata)
      real*8 :: ydata(ndata)
      real*8 :: tval,ypval,yval
      
      !
      !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
      !  nearest to, TVAL.
      !
      call rvec_bracket ( ndata, tdata, tval, left, right )
      !
      !  Now evaluate the piecewise linear function.
      !
      ypval = (ydata(right)-ydata(left))/(tdata(right)-tdata(left))
      
      yval = ydata(left) +  ( tval - tdata(left) ) * ypval
      
      return
      end

      subroutine rvec_bracket ( n, x, xval, left, right )
      !
      !*******************************************************************************
      !
      !  Parameters:
      !
      !    Input, integer N, length of input array.
      !
      !    Input, real X(N), an array sorted into ascending order.
      !
      !    Input, real XVAL, a value to be bracketed.
      !
      !    Output, integer LEFT, RIGHT, the results of the search.
      !    Either:
      !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
      !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
      !    or
      !      X(LEFT) <= XVAL <= X(RIGHT).
      !
        implicit none
      !
        integer :: n,i,left,right
        real*8 :: x(n)
        real*8 :: xval
      !
        do i = 2, n - 1
      
          if ( xval < x(i) ) then
            left = i - 1
            right = i
            return
          end if
      
         end do
      
        left = n - 1
        right = n
      
        return
      end

      
      
      
      
      
