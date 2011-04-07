      PROGRAM compensator
      !
      ! a dispersion compensator for pulsed electron ?????
      !  
      !
      use msimsl


      IMPLICIT NONE
      common/par/ pi,me,e,length,vini
      real*8 pi,me,e,length,vini
      common/tole/ tol
      real*8 tol
      common/time/ endoftim
      real*8 endoftim 
      common/worksp/rwksp
      real*8 rwksp(43592)
      call iwkin(43592)
      OPEN(UNIT=11,FILE='disp.dat',type='replace')
      OPEN(UNIT=12,FILE='velocity.dat',type='replace')        
      OPEN(UNIT=13,FILE='times.dat',type='replace')        
      OPEN(UNIT=14,FILE='width.dat',type='replace')        

      !write(6,*) ' stop integration at (micro s):'
      !read(*,*) endoftim

      !endoftim=endoftim*1.E-6 
      tol=.0000001D0

      write(6,*) "version 1.5"
      call parameters
      call result

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      STOP
      END



      subroutine parameters
      implicit none
      common/par/ pi,me,e,length,vini
      real*8 pi,me,e,length,vini
      common/time/ endoftim
      real*8 endoftim 

      PI=dACOS(-1d0)
      me=9.10938188D-31
      e=1.60217646D-19
      length=.32815D0
      !length=2D0
      vini=3.23D7
      endoftim=length/vini

      return
      end



      subroutine fcn(neq,t,y,yprime)
      implicit none
      integer neq
      real*8 t,y(neq),yprime(neq)
      common/par/ pi,me,e,length,vini
      real*8 pi,me,e,length,vini
      real*8 Bx,By,Ex,Ez
      common/trajectory/Bz,Ey
      real*8 Bz,Ey


      Ex=0D0					 
      Ez=0D0
      Bx=0D0
      By=0D0


      yprime(1)=y(4)
      yprime(2)=y(5)
      yprime(3)=y(6)
      yprime(4)=(e*(y(5)*Bz-y(6)*By)+e*Ex)/me  
      yprime(5)=(e*(y(6)*Bx-y(4)*Bz)+e*Ey)/me  
      yprime(6)=(e*(y(4)*By-y(5)*Bx)+e*Ez)/me	 


      return
      end

      !p varies angle, q varies velocity
      subroutine result
      implicit none
      common/trajectory/Bz,Ey
      real*8 Bz,Ey
      common/par/ pi,me,e,length,vini
      real*8 pi,me,e,length,vini
      common/tole/ tol
      real*8 tol
      integer counter,test
      real*8 y1old
      common/time/ endoftim
      real*8 divider,endoftim
      integer mxparm,neq,i,ido,p,q,k,l,z,m,time(6),lowstepcheck,iold
      parameter (mxparm=120,neq=6,p=3,q=3)
      integer r,error, ti, scanl,scanprec
      parameter (scanl=20,scanprec=2)
      parameter (r=5)
      real*8 fcn,param(mxparm),t,tend,tft,y(neq),B,y0,Ed,E0,Bdw,x1,x2,dw
      parameter (B=.031D0,dw=.8D-1)
      real*8 vout(2,1000,p*q),mem(2,1000,p*q),weintime(7,p*q),dataout(7,2,p)
      real*8 pend,pos1,pos2,endtime,step,histep,lowstep,sls,comp(p,q)
      real*8 dw1,dw2,dummy,ddw,wfls,bls,disp(1:p*q),debug(10),yold,theta
      integer step2,step1
      parameter (pend=1D-1,ddw=1D-4)
      real*8 tftold,hiTstep,lowTstep,wfTstep,bTstep,Vg(p*q),scandata(p,q,0:scanl*2),scandx(0:scanl*2)
      real*8 timecheck
      external fcn,divprk,sset
      ! !run test electron to get information about path
      write(6,*) "Initial Electron"
      length=me*vini/(e*B)*4D0*pi+pend
      endoftim=length/vini
      !lowTstep=endoftim/lowstep
      !sls=lowTstep*vini

      !lowstep=length/(((ddw/1D1)/4D0)/10D0) !set how small the steps are
      lowTstep=ddw/(4D1*vini)        !based how wide the shortest boundary is
      if (lowTstep.gt.7.75D-15) then
              lowTstep=7.75D-15
      endif
      bTstep=lowTstep*1D1
      bls=bTstep*vini
      wfTstep=lowTstep*.1D2
      wfls=wfTstep*vini
      sls=lowTstep*vini        !distance infront of border to go to hi-res
      hiTstep=sls/(1D3*vini)   !hi-res step for approaching borders 
      !!!do m=1,r                !decreases the filter width each iteration
      !dw=pend-m*ddw
      !E0=-B*vini/4D0          !sets the value of the Wien filter so c=-1
      E0=-(pend*B*vini)/(dw*4D0)
      pos1=(pend-dw)/4D0 !first B-field border
      pos2=pend-pos1     !second B-field border
      dw1=pend/2D0-dw/2D0 !start of Wien filter
      dw2=pend/2D0+dw/2D0 !end of Wien filter
      !lowstep=length/(pos1/10)
      !midstep=endoftim*vini*1D6/dw
      !histep=1000*endoftim/(sls/vini) 
      !hiTstep=endoftim/histep !time changed for each step
      write(6,51) pos1,pos2 !write info to console
      write(6,52) dw1,dw2
51    format(x,'pos1=',E12.4,x,"pos2=",E12.7)
52    format(x,'dw1=',E12.4,x,'dw2=',E12.7)
      write(6,53) lowTstep,hiTstep
      write(6,54) 'sls=',sls,'lowDstep=',lowTstep*vini,'hiDstep=',hiTstep*vini
53    format(x,'lowTstep=',E12.4,x,'hiTstep=',E12.4)
54    format(x,3(A,E12.4,x))
      write(6,'(x,A,F12.4,x,A,E12.4)') 'length=',length,'end time=',endoftim
      ido=1   !set initial conditions.
      step=lowTstep
      t=0D0
      Bz=0D0
      Ey=0D0
      counter=0
      do i=1,neq
              y(i)=0D0
      enddo
      y(4)=vini	  !cos(.01*(k-(p+1)/2))
      !	write(6,*) "v0=", y(4)
      !y(5)=vini*sin(.01*(k-(p+1)/2))	  !(1+0.01*(j-3))
      !	y(2)=j*0.001   !(start at 11 position separated by 1 mm each)
      !	write(6,*) y(2)
      call sset(mxparm,0.0,param,1)
      param(4)=1000000000
      param(10)=1.0   
      tft=0D0
      i=0
      do while(counter.ne.8.or.y(1).lt.pend)
              !write(6,*) y(1),y(2)
              !READ (*,*) dummy
              i=i+1
              tft=tft+step
              tend=tft					    
              call divprk(ido,neq,fcn,t,tend,tol,param,y)

60    format(x,'counter=',i1,x,'flag=',i1,x,'steps=',i10)
              if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls).and.(counter.eq.0)) then
                      step=hiTstep
                      !write(6,60) counter,0,i
                      lowstepcheck=lowstepcheck+1
              endif
              if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls).and.(counter.eq.2)) then
                      step=hiTstep
                      !write(6,60) counter,2,i
                      lowstepcheck=lowstepcheck+1
              endif

              if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls).and.(counter.eq.2)) then
                      step=lowTstep
                      !write(6,60) counter,2,i
                      lowstepcheck=lowstepcheck+1
              endif
              if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls).and.(counter.eq.1)) then
                      !write(6,60) counter, 1, i
                      step=lowTstep
                      lowstepcheck=lowstepcheck+1
              endif
              if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls).and.(counter.eq.1)) then
                      !write(6,60) counter, 1, i
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
              endif
              if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls).and.counter.eq.3) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter,3, i
              endif
              if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls).and.counter.eq.4) then
                      step=lowTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter,4, i
              endif
              if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls).and.counter.eq.4) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter,4, i
              endif
              if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls).and.(counter.eq.6)) then
                      step=lowTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 5, i
              endif
              if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls).and.counter.eq.5) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 6, i
              endif

              if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls).and.(counter.eq.6)) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 5, i
              endif
              if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls).and.counter.eq.7) then
                      step=lowTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 6, i
              endif
              if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls).and.counter.eq.7) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 6, i
              endif
              !if ((y1old.le.pend-sls).and.(y(1).gt.pend-sls)
              !1.and.(counter.eq.8)) then
              if ((tft.gt.endoftim-lowTstep).and.(tftold.le.endoftim-lowTstep)) then
                      step=hiTstep
                      lowstepcheck=lowstepcheck+1
                      !write(6,60) counter, 7, i
              endif
61    format(x,'counter=',i1,8x,'steps=',i10)
              if ((y1old.le.pos1).and.(y(1).gt.pos1).and.(counter.eq.0)) then
                      counter=1
                      step=bTstep
                      time(1)=i
                      !write(6,61) counter,i
              endif
              if ((y1old.gt.pos1).and.(y(1).le.pos1).and.(counter.eq.1)) then
                      counter=2
                      step=bTstep
                      !write(6,61) counter,i

              endif
              if ((y1old.le.pos1).and.(y(1).gt.pos1).and.(counter.eq.2)) then
                      step=lowTstep
                      counter=3
                      !write(6,61) counter,i
                      time(2)=i
              endif
              if ((y1old.le.dw1).and.(y(1).gt.dw1).and.(counter.eq.3)) then
                      !write(6,*) pend/2-dw/2,pend/2+dw/2
                      counter=4
                      step=wfTstep
                      !write(6,61) counter,i
                      time(3)=i
                      y0=y(2)
              endif
              if ((y1old.le.dw2).and.(y(1).gt.dw2).and.(counter.eq.4)) then
                      step=lowTstep
                      counter=5
                      !write(6,61) counter,i
                      time(4)=i
              endif
              if ((y1old.le.pos2).and.(y(1).gt.pos2).and.(counter.eq.5)) then
                      step=bTstep
                      counter=6
                      !write(6,61) counter,i
                      time(5)=i
              endif
              if ((y1old.gt.pos2).and.(y(1).le.pos2).and.(counter.eq.6)) then
                      step=bTstep
                      counter=7
                      !write(6,61) counter,i
              endif
              if ((y1old.le.pos2).and.(y(1).gt.pos2).and.(counter.eq.7)) then
                      step=lowTstep
                      counter=8
                      !write(6,61) counter,i
                      time(6)=i
                      x1=y(1)
              endif

              if (counter.eq.1) then 
                      Bz=-B
              endif 
              if (counter.eq.2) then 
                      Bz=B
              endif 
              if (counter.eq.3) then 
                      Bz=0D0
              endif 
              if (counter.eq.4) then 
                      !sign of wienfilter should be negative
                      Ey=E0
                      Bz=Ey/vini
              endif
              if (counter.eq.5) then 
                      Bz=0D0
                      Ey=0D0
              endif 
              if (counter.eq.6) then 
                      Bz=B
                      Ey=0D0
              endif 
              if (counter.eq.7) then 
                      Bz=-B
              endif 
              if (counter.eq.8) then 
                      Bz=0D0
              endif
              y1old=y(1)
              tftold=tft
      enddo
                      timecheck=tft
      x2=y(1)
      ido=3
      endtime=tft
      call divprk(ido,neq,fcn,t,tend,tol,param,y)
      write(6,*) "Lowstepcheck=",lowstepcheck
      lowstepcheck=0
      z=i	! this sets the time steps to stop 
      !the v0 electron at correct distance
      !first electron done, do the rest
      Ed=E0 
      divider=endtime/1000D0
      do l=1,q
              do k=1,p
                      test=0
                      ido=1
                      t=0d0
                      ti=0
                      Bz=0.
                      Ey=0.
                      counter=0
                      lowstepcheck=0
                      do i=1,neq
                              y(i)=0d0
                      enddo

                      y(4)=(1.00D0+2D0**1*2.5D-5*(l-(q+1D0)/2D0))*vini*cos(.001D0*(k-(p+1D0)/2D0))
                      y(5)=(1.00D0+2D0**1*2.5D-5*(l-(q+1D0)/2D0))*vini*sin(.001D0*(k-(p+1D0)/2D0))
                      !	y(2)=j*0.001   !(start at 11 position separated by 1 mm each)
                      !	write(6,*) y(2)
                      write(6,50) k,p,l,q,m,r,Ed,y(4)
50    format(x/,x/,x,i1,'/',i1,x,i1,'/',i1,x,i2,'/',x,i2/,' Ed=',E12.2,x,'v0=',E12.2//)
                      call sset(mxparm,0.0,param,1)
                      param(4)=1000000000
                      param(10)=1D0 
                      tft=0D0
                      i=0
                      step=lowTstep
                      do while(tft.le.endtime+scanl*lowTstep*scanprec)
                              if ((test*divider).le.tft.and.test.lt.999) then
                                      test=test+1

                                      mem(1,test,k+p*(l-1))=y(1)
                                      mem(2,test,k+p*(l-1))=y(2)
                                      vout(1,test,k+p*(l-1))=tft
                                      vout(2,test,k+p*(l-1))=sqrt(y(4)*y(4)+y(5)*y(5))
                              endif
                              i=i+1
                              tft=tft+step
                              tend=tft  
                              call divprk(ido,neq,fcn,t,tend,tol,param,y)
                              if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls).and.counter.eq.2) then
                                      step=lowTstep
                                      !write(6,60) counter, 1, i
                                      lowstepcheck=lowstepcheck+1
                              endif
                              if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls).and.(counter.eq.1)) then
                                      step=lowTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 2, i
                              endif
                              if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls).and.(counter.eq.0.or.counter.eq.2)) then
                                      step=hiTstep
                                      !write(6,60) counter, 1, i
                                      lowstepcheck=lowstepcheck+1
                              endif
                              if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls).and.(counter.eq.1)) then
                                      step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 2, i
                              endif
                              if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls).and.counter.eq.3) then
                                      step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter,3, i
                              endif
                              if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls).and.counter.eq.4) then
                                      step=lowTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter,4, i
                              endif
                              if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls).and.counter.eq.4) then
                                      step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter,4, i
                              endif
                              if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls).and.(counter.eq.6)) then
                                      step=lowTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 5, i
                              endif
                              if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls).and.(counter.eq.6)) then
                                      step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 5, i
                              endif
                              if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls).and.counter.eq.7) then
                                      step=lowTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 6, i
                              endif
                              if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls).and.(counter.eq.7.or.counter.eq.5)) then
                                      step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 6, i
                              endif
                              !if ((y1old.le.pend-sls).and.(y(1).gt.pend-sls)
                              !1.and.(counter.eq.8)) then
                              if ((tft.gt.endoftim-lowTstep*2D0).and.(tftold.le.endoftim-lowTstep*2D0)) then
                                      !step=hiTstep
                                      lowstepcheck=lowstepcheck+1
                                      !write(6,60) counter, 7, i
                              endif
                              if ((tft.gt.timecheck+lowTstep*scanprec*(ti-scanl)).and.
     1(tftold.le.timecheck+lowTstep*scanprec*(ti-scanl))) then
                                      scandata(k,l,ti)=y(1)
                                      !write(6,*) ti,counter, tft, scandata(k,l,ti)
                                      if (ti.ne.scanl*2) then
                                              ti=ti+1
                                      endif
                              endif
62    format(18x,A,i10)
                              if ((y1old.le.pos1).and.(y(1).gt.pos1).and.(counter.eq.0)) then
                                      counter=1
                                      !write(6,62) "steps=", i
                                      step=bTstep
                                      !			dataout(1,1,k)=y(1)
                                      !			dataout(2,1,k)=y(2)
                                      !			dataout(3,1,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !			dataout(4,1,k)=tft
                                      weintime(1,k+p*(l-1))=tft
                              endif
                              if ((y1old.gt.pos1).and.(y(1).le.pos1).and.(counter.eq.1)) then
                                      counter=2
                                      !write(6,62) "steps=", i
                                      step=bTstep 
                              endif
                              if ((y1old.le.pos1).and.(y(1).gt.pos1).and.(counter.eq.2)) then
                                      counter=3	
                                      !write(6,62) "steps=", i
                                      step=lowTstep
                                      !			dataout(1,2,k)=y(1)
                                      !			dataout(2,2,k)=y(2)
                                      !			dataout(3,2,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !			dataout(4,2,k)=tft
                                      weintime(2,k+p*(l-1))=tft
                              endif
                              if ((y1old.le.dw1).and.(y(1).gt.dw1).and.(counter.eq.3)) then
                                      !debug(1)=y(4)
                                      !write(6,*) "dv=", e*B*(y(2)-y0)/(4*me)
                                      !write(6,*) "dvg=", Ed*e*(y(2)-y0)/(2*me*vini)
                                      yold=sqrt(y(4)**2+y(5)**2)
                                      theta=ATAN(y(5)/y(4))
                                      y(4)=sqrt(yold**2+2D0*e*Ed*(y(2)-y0)/me)*cos(theta)
                                      y(5)=sqrt(yold**2+2D0*e*Ed*(y(2)-y0)/me)*sin(theta)
                                      debug(2)=y(4)-vini
                                      Vg(k+p*(l-1))=y(4)
                                      !dummy=e*(y(2)-y0)/me*(B/4+Ed/vini)
                                      debug(3)=e*(y(2)-y0)/me*(B/4D0+Ed/((vini)))
                                      write(6,*) debug(2), debug(3), test
                                      counter=4
                                      !write(6,62) "steps=", i
                                      step=wfTstep
                                      !write(6,*) "enter field",e*Ed*(y(2)-y0)/(me*vini)

                                      !		dataout(1,3,k)=y(1)
                                      !		dataout(2,3,k)=y(2)
                                      !		dataout(3,3,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !		dataout(4,3,k)=tft
                                      weintime(3,k+p*(l-1))=tft
                              endif

                              if ((y1old.le.dw2).and.(y(1).gt.dw2).and.(counter.eq.4)) then
                                      yold=sqrt(y(4)**2+y(5)**2)
                                      theta=ATAN(y(5)/y(4))
                                      y(4)=sqrt(yold**2-2D0*e*Ed*(y(2)-y0)/me)*cos(theta)
                                      y(5)=sqrt(yold**2-2D0*e*Ed*(y(2)-y0)/me)*sin(theta)
                                      !write(6,*) "leave field",-e*Ed*(y(2)-y0)/(me*vini)
                                      counter=5
                                      !write(6,62) "steps=", i
                                      step=lowTstep
                                      !			dataout(1,4,k)=y(1)
                                      !			dataout(2,4,k)=y(2)
                                      !			dataout(3,4,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !			dataout(4,4,k)=tft
                                      weintime(4,k+p*(l-1))=tft
                              endif
                              if ((y1old.le.pos2).and.(y(1).gt.pos2).and.(counter.eq.5)) then
                                      counter=6
                                      !write(6,62) "steps=", i
                                      step=bTstep
                                      !			dataout(1,5,k)=y(1)
                                      !			dataout(2,5,k)=y(2)
                                      !			dataout(3,5,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !			dataout(4,5,k)=tft
                                      weintime(5,k+p*(l-1))=tft
                              endif
                              if ((y1old.gt.pos2).and.(y(1).le.pos2).and.(counter.eq.6)) then
                                      counter=7
                                      !write(6,62) "steps=", i
                                      step=bTstep
                              endif
                              if ((y1old.le.pos2).and.(y(1).gt.pos2).and.(counter.eq.7)) then
                                      counter=8
                                      !write(6,62) "steps=", i
                                      step=lowTstep
                                      !dataout(1,6,k)=y(1)
                                      !dataout(2,6,k)=y(2)
                                      !dataout(3,6,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !dataout(4,6,k)=tft
                                      step1=i
                                      weintime(6,k+p*(l-1))=tft
                              endif
                              if((tftold.le.endtime).and.(tft.gt.endtime)) then
                                      comp(k,l)=y(1)
                              endif

                              if (counter.eq.1) then 
                                      Bz=-B
                              endif 
                              if (counter.eq.2) then 
                                      Bz=B
                              endif 
                              if (counter.eq.3) then 
                                      Bz=0D0
                              endif 
                              if (counter.eq.4) then 
                                      !sign of wienfilter should be negative
                                      Ey=Ed
                                      !Bz=Ey/(sqrt((vini+(y(2)-y0)*e*B/(4*me))**2+e*Ed*(y(2)-y0)/me))
                                      Bz=Ed*(1/vini-(y(2)-y0)*(e*(B*vini+4D0*Ed))/(4D0*me*vini**3)
     1+((y(2)-y0)*(e*(B*vini+4D0*Ed))/(4D0*me*vini**3))**2)

                                      !Bz=E0/y(4)
                                      !Bz=Ey/vini
                              endif
                              if (counter.eq.5) then 
                                      Bz=0D0
                                      Ey=0D0
                              endif 
                              if (counter.eq.6) then 
                                      Bz=B
                                      Ey=0D0
                              endif 
                              if (counter.eq.7) then 
                                      Bz=-B
                              endif 
                              if (counter.eq.8) then 
                                      Bz=0D0
                              endif
                              y1old=y(1)
                              tftold=tft

                              if (time(6).eq.i) then
                                      !write(6,*) m,i,time(m)
                                      !dataout(1,1,k)=y(1)
                                      !dataout(2,1,k)=y(2)
                                      !dataout(3,1,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                                      !dataout(4,1,k)=dATAN(y(5)/y(4))
                                      !dataout(5,1,k)=tft
                                      !dataout(6,1,k)=y(4)-vini
                                      !dataout(7,1,k)=y(1)-x1
                              endif                
                      enddo
                      step2=i
                      write(6,*) "steps=", step2-step1
                      if(lowstepcheck.ne.14) then
                              error=1
                      endif
                      lowstepcheck=0
                      write(6,*) "end step=", i
                      !mem(1,test,k+3*(l-1))=y(1)
                      !mem(2,test,k+3*(l-1))=y(2)
                      !dataout(1,2,k)=y(1)
                      !dataout(2,2,k)=y(2)
                      !dataout(3,2,k)=sqrt(y(4)*y(4)+y(5)*y(5))
                      !dataout(4,2,k)=dATAN(y(5)/y(4))
                      !dataout(5,2,k)=tft
                      !dataout(6,2,k)=y(4)-vini
                      !dataout(7,2,k)=y(1)-x2
                      weintime(7,k+p*(l-1))=tft
                      ido=3
                      mem(1,1000,k+p*(l-1))=y(1)
                      mem(2,1000,k+p*(l-1))=y(2)
                      vout(1,1000,k+p*(l-1))=tft
                      vout(2,1000,k+p*(l-1))=sqrt(y(4)*y(4)+y(5)*y(5))

                      !	write(6,*) "test2" 
                      call divprk(ido,neq,fcn,t,tend,tol,param,y)

                      tft=0d0

              enddo
      enddo
      write(6,'(A,E12.4)') ' dx=',maxval(comp(1:p,1:q))-minval(comp(1:p,1:q))
      !if (error.eq.1) then
      !write(6,*) "Invalid Test", lowstepcheck
      !write(12,101) percent, 
      !maxval(comp(1:p,1:q))-minval(comp(1:p,1:q))
      !else
      !write(11,101) percent, 
      !maxval(comp(1:p,1:q))-minval(comp(1:p,1:q))
      !endif
      lowstepcheck=0
      !!!enddo
      do test=0,scanl*2
              scandx(test)=maxval(scandata(1:p,1:q,test))-minval(scandata(1:p,1:q,test))
              !write(6,*) test, scandx(test)
              write(14,104) scandata(1:p,1:q,test)-scandata(2,2,test), scandx(test)
      enddo
      if(minval(scandx(0:scanl*2)).eq.scandx(0).or.minval(scandx(0:scanl*2)).eq.scandx(scanl*2)) then
              write(6,*) "WARNING!! False Minimum. increase scanl or scanprec"
      endif


      write(6,*) "dx is also=", minval(scandx(0:scanl*2))

      !do test=1,1000
      !        write(11,102) mem(1:2,test,1:p*q)
      !        write(12,102) mem(1,test,(p*q+1)/2), vout(2,test,1:p*q)
      !        !disp(1:p*q)=mem(1,test,1:p*q)-mem(1,test,(p*q+1)/2)
      !        write(14,102) mem(1,test,(p*q+1)/2), mem(1,test,3)-mem(1,test,2)!maxval(mem(1,test,1:p*q))-minval(mem(1,test,1:p*q))
      !enddo
      do test=1,7
              write(13,102) weintime(test,1:p*q)
      enddo
      do test=1,p*q
              write(6,*) Vg(test)
      enddo
      !write(6,*) "counter=",counter

101   format(2E15.7)
102  	format(18E15.7)
103   format(7(E12.7,','))
104   format(10E15.7)
      return
      end
