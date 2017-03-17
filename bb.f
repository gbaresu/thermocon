      program bb
      
      implicit none
      integer,parameter :: NMAX=1000000
      integer,parameter :: NB=15
      real*8 :: iii,hpl,cl,kbol,TT,nu,kk,sig,ITOT,ITOT2,IX
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
      
      write(*,*) band/1e4
      
      TT=5777
      
      hpl=6.62e-27
      kbol=1.38e-16
      sig=5.67e-5
      cl=3e10
    
      !WRITE OUT angstrom  erg/cm2/s/ang
      do i1=10,50000
       kk=real(i1)
       nu=kk*cl
       iii=(2*hpl*nu**3/cl**2)/(exp((hpl*nu)/(kbol*TT))-1)
       write(10,*) (cl/nu)*1e8,iii*nu**2/cl/1e8
      enddo

      write(cmd,'(A)') "tac fort.10 > bb.dat"
      call system(trim(cmd))
      
      ITOT=sig*TT**4/3.14
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
       write(*,*) II(i2)/ITOT
      enddo
      
      call plint(NN-1,bbb,kkk,NN-1,band(1),band(NB),yval)
      ITOT2=yval
      write(*,*)
      write(*,*) ITOT,ITOT2,ITOT/ITOT2
      write(*,*)
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
       write(*,*) yval/ITOT,yval/ITOT2
      enddo

      do i2=1,NB-1
       do i1=1,NN-1
        if(kkk(i1).ge.band(i2).and.kkk(i1).le.band(i2+1)) then
         write(20,*) kkk(i1),II(i2)/(band(i2+1)-band(i2))*3.14
        endif
       enddo
      enddo

      IX=0
      do i1=1,NB-1
       IX=IX+II(i1)
       !write(*,*) band(i1)
      enddo


      end program
