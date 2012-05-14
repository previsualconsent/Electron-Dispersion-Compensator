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
      OPEN(UNIT=12,FILE='width.dat',type='replace')        
      OPEN(UNIT=13,FILE='endpoint.dat',type='replace')        
      OPEN(UNIT=15,FILE='scan.dat',type='replace')        
      OPEN(UNIT=16,FILE='histogram.dat',type='replace')        

      !write(6,*) ' stop integration at (micro s):'
      !read(*,*) endoftim

      !endoftim=endoftim*1.E-6 
      tol=1D-15

      write(6,*) char(7)
      write(6,*) "version 1.5"
      call parameters
      call result

      write(6,*) char(7)
      close(11)
      close(12)
      close(13)
      close(15)
      close(16)
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
      vini=3.23D1
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
      real*8 B,y0
      integer counter
      common/func/B,y0,counter


      Ex=0D0					 
      Ez=0D0
      Bx=0D0
      By=0D0

      !Bz is constant unless counter=4
      if(counter.eq.4) then
         Bz= Ey/(vini+(B/4+Ey/vini)*e*(y(2)-y0)/me)
      endif

      yprime(1)=y(4)
      yprime(2)=y(5)
      yprime(3)=y(6)
      yprime(4)=((y(5)*Bz-y(6)*By)+Ex)*e/me
      yprime(5)=((y(6)*Bx-y(4)*Bz)+Ey)*e/me
      yprime(6)=((y(4)*By-y(5)*Bx)+Ez)*e/me


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
      integer mxparm,neq,i,ido,p,q,k,l,boundarycheck,scani
      parameter (mxparm=120,neq=6,p=31,q=31)
      integer error, ti, scanl,midk,midl,n_bins
      parameter (scanl=40,n_bins=40)
      integer hist(2,0:scanl,1:n_bins)
      real*8 fcn,param(mxparm),t,tend,tft,y(neq),E0,dw,dE,dtheta
      real*8 mem(2,1000,p*q)
      real*8 B,y0,bin_size,minimum,maximum,sums,sumsq,pulsecenter
      common/func/B,y0,counter

      parameter (dw=.5D5,dE=1D-3,dtheta=2D-3)
      real*8 pend,pos1,pos2,endtime,step,endTstep,sls
      real*8 dw1,dw2,ddw,wfls,bls,yold,theta,theta2,scanprec
      parameter (pend=1D5,ddw=1D-4,scanprec=2.5D-1)
      real*8 tftold,hiTstep,lowTstep,wfTstep,bTstep,scandata(p,q,0:scanl),scandx(0:scanl*2)
      real*8 timecheck,xf(p,q),yf(p,q),ET0,ET,endpoints(2,p,q),sd
      logical scanning
      
      external fcn,divprk,sset

      scanning=.false.
      B=3.1D-15 !B-Field value for turning fields 

      !initial calculations
      midk=(p+1)/2
      midl=(q+1)/2
      length=me*vini/(e*B)*4D0*pi+pend   !total length based on path
      endoftim=length/vini          !total time for travel

      lowTstep=endoftim/1D5      !low precision time step
      sls=lowTstep*vini             !average distance for 1 lowTstep


      bTstep=lowTstep*5D-1           !time step for B-fields
      bls=bTstep*vini               !average distance of 1 bTstep

      wfTstep=lowTstep*1D1         !time step for WF
      wfls=wfTstep*vini             !average distance for 1 wfTstep

      hiTstep=lowTstep*1D-5          !hi-res step for approaching borders 
      endTstep=lowTstep*1D-5    !set step size for final section
      E0=-(pend*B*vini)/(dw*4D0)    !set E field for WF based on theory
      write(6,*) me*vini/(e*B)
      

      !some distance calculations
      pos1=(pend-dw)/4D0 !first B-field border
      pos2=pend-pos1     !second B-field border
      dw1=pend/2D0-dw/2D0 !start of Wien filter
      dw2=pend/2D0+dw/2D0 !end of Wien filter

51    format(x,'pos1=',E12.4,x,"pos2=",E12.7)
52    format(x,'dw1=',E12.4,x,'dw2=',E12.7)
53    format(x,'lowTstep=',E12.4,x,'hiTstep=',E12.4)
54    format(x,3(A,E12.4,x))

      ido=1   !set initial conditions.

      step=lowTstep
      Bz=0D0
      Ey=0D0
      t=0D0
      counter=0
      do i=1,neq
         y(i)=0D0
      enddo
      y(4)=vini

      call sset(mxparm,0.0,param,1)

      param(4)=1000000000
      param(10)=1.0   
      tft=0D0
      i=0

      do while(counter.ne.8.or.y(1).lt.pend) !loop until the electron reaches the target

         i=i+1             !increment the step
         tft=tft+step      !increment the time
         tend=tft          !set new tend  

         call divprk(ido,neq,fcn,t,tend,tol,param,y) !call integrater

60    format(x,'counter=',i1,x,'flag=',i2,x,'steps=',i10,x,'Tstep=',E12.4)
61    format(x,'counter=',i1,9x,'steps=',i10,x,'Tstep=',E12.4)

         select case (counter) !switch on counter- this reduced long if statements
                                    !each case checks to see if the particle is getting
                                    !close to the next barrier. If so, reduce the step
                                    !size. When it crosses the barrier, hand off to 
                                    !next case. boundary check is used for error checking
            case (0)  !Free space before device
               Bz=0D0
               Ey=0D0
               if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                  counter=1
                  step=bTstep
                  Bz=-B
               endif
            case (1) !First Turning Field
               if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.gt.pos1).and.(y(1).le.pos1)) then
                  counter=2
                  Bz=B
                  step=bTstep
               endif
            case (2) !Second Turning Field
               if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                  step=lowTstep
                  counter=3
                  Bz=0D0
               endif

            case (3) !empty space between turning fields and 
               if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.dw1).and.(y(1).gt.dw1)) then ! enter thw WF
                  counter=4
                  Ey=E0
                  Bz=Ey/vini
                  step=wfTstep
                  y0=y(2)                   ! y0 has the center of the WF
                  timecheck=tft
               endif

            case (4)

               if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls)) then
                  step=lowTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.dw2).and.(y(1).gt.dw2)) then
                  step=lowTstep
                  counter=5
                  Bz=0D0
                  Ey=0D0
               endif

            case (5)
               if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                  step=bTstep
                  counter=6
                  Bz=B
                  Ey=0D0
               endif

            case (6)
               if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.gt.pos2).and.(y(1).le.pos2)) then
                  step=bTstep
                  counter=7
                  Bz=-B
               endif
            case (7)
               if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
               endif
               if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                  step=lowTstep
                  counter=8
                  Bz=0D0
               endif

            end select
            !lower step size while approaching final region.
            if ((tft.gt.endoftim-lowTstep).and.(tftold.le.endoftim-lowTstep)) then
               step=hiTstep
               boundarycheck=boundarycheck+1
            endif

            !store info from last run
            y1old=y(1)
            tftold=tft
         enddo

          ido=3
          endtime=tft
          write(6,*) endtime-endoftim
          write(6,*) "debug", pend,tft, y(1)
          call divprk(ido,neq,fcn,t,tend,tol,param,y) !close out ODE workspace



          !first electron done, do the rest

          divider=endtime/1000D0
          do l=1,q
              do k=1,p
                  !reset variables for next run
                  test=0
                  ido=1
                  t=0D0
                  ti=0
                  Bz=0D0
                  Ey=0D0
                  counter=0
                  scanning=.false.
                  scani=0
                  boundarycheck=0
                  do i=1,neq
                      y(i)=0D0
                  enddo

                  !sets energy spread to dE and angle spread to dtheta
                  y(4)=(1.00D0+dE/(4D0*(q-1)/2D0)*(l-midl))*vini*dcos(dtheta/2D0*(k-midk)/((p-1)/2D0))
                  y(5)=(1.00D0+dE/(4D0*(q-1)/2D0)*(l-midl))*vini*dsin(dtheta/2D0*(k-midk)/((p-1)/2D0))

                  ET0=.5D0*me*(y(4)**2+y(5)**2)

                  write(6,50) k,p,l,q,dE
50    format(x/,x,i2,'/',i2,x,i2,'/',i2,x," dEm=",E15.7)
                  call sset(mxparm,0.0,param,1)
                  param(4)=1000000000
                  param(10)=1D0 
                  tft=0D0
                  i=0
                  step=lowTstep

62    format(18x,A,i10)

                  ! loop until the scan has finished
                  do while(scani.le.scanl)
                     if ((test*divider).le.tft.and.test.lt.999) then
                        test=test+1
                        mem(1,test,k+p*(l-1))=y(1)
                        mem(2,test,k+p*(l-1))=y(2)
                     endif

                     ET=.5D0*me*(y(4)**2+y(5)**2)+Ey*q*y(2)
                      i=i+1
                      tft=tft+step
                      tend=tft  
                      if ((ET-ET0)/ET0.gt.1D-10) then
                         write(6,*) "SOMETHING IS WRONG! ET: ",ET," ET0: ", ET0
                         write(6,*) counter
                         ET0=ET
                      endif

                      if(MOD(i,100000).eq.0) then
66    format(i3,'%')
                         write(6,66) int(tft/endtime*100)
                      endif
                      call divprk(ido,neq,fcn,t,tend,tol,param,y)


                      !same switch as above
                      select case (counter)
                      case (0)
                          Bz=0D0
                          Ey=0D0
                          if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                              counter=1
                              Bz=-B
                              step=bTstep
                          endif
                      case (1)
                          if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)) then
                              !step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.gt.pos1).and.(y(1).le.pos1)) then
                              counter=2
                              step=bTstep
                              Bz=B
                          endif
                      case (2)
                          if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                              !step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                              step=lowTstep
                              counter=3
                              Bz=0D0
                          endif

                      case (3)
                          if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.dw1).and.(y(1).gt.dw1)) then
                              !ENTERING WIEN FILTER
                              !calculate old velocity and angle
                              yold=y(4)
                              theta=ATAN(y(5)/y(4))

                              !set new velocity based on calculated energy change
                              y(4)=dsqrt(yold**2+2D0*E0*e*(y(2)-y0)/me)
                              !write(6,*) "yold=",yold,"vx_0 in filter=", y(4)
                              !write(6,*) "y=", y(2)-y0                  
                              !write(6,*) "y_a=", 

                              counter=4
                              Ey=E0
                              !Matched B field. For particle traveling strait through, the field should be balanced
                              !Bz=E0*(1/vini-((y(2)-y0)*(e*(B*vini+4D0*E0))/(4D0*me*vini**3)+((y(2)-y0)*(e*(B*vini+4D0*E0))/(4D0*me))**2*(1/vini**6))*1D0)
                              Bz= E0/(vini+(B/4+E0/vini)*e*(y(2)-y0)/me)
                              !Bz= E0/(dsqrt((vini+(y(2)-y0)*e*B/(4D0*me))**2+2D0*E0*e*(y(2)-y0)/me))
                              step=wfTstep
                          endif

                      case (4)
                         Bz= E0/(vini+(B/4+E0/vini)*e*(y(2)-y0)/me)
                         if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls)) then
                            step=lowTstep
                            boundarycheck=boundarycheck+1
                         endif
                         if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls)) then
                            step=hiTstep
                            boundarycheck=boundarycheck+1
                         endif
                          if ((y1old.le.dw2).and.(y(1).gt.dw2)) then
                              !calculate velocity change for leaving WF
                              yold=y(4)
                              theta=ATAN(y(5)/y(4))
                              theta2=theta


                              y(4)=dsqrt(yold**2-2D0*E0*e*(y(2)-y0)/me)
                              counter=5
                              Bz=0D0
                              Ey=0D0
                              step=lowTstep
                          endif

                      case (5)
                          if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                              step=bTstep
                              counter=6
                              Bz=B
                              Ey=0D0
                          endif

                      case (6)
                          if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)) then
                              !step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.gt.pos2).and.(y(1).le.pos2)) then
                              step=bTstep
                              counter=7
                              Bz=-B
                          endif
                      case (7)
                          if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls)) then
                              step=lowTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                              step=lowTstep
                              counter=8
                              Bz=0D0
                          endif

                       case(8)
                          if(scanning) then
                             scandata(k,l,scani)=y(1)
                             if(k.eq.midk-1.and.l.eq.midl-1) write(6,*) y(1)
                             scani=scani+1
                             step=lowTstep*scanprec
                             !write(6,*) scani,scanl,tft,endtime
                          endif

                          if(.not.scanning) then
                             step=(endtime-lowTstep*scanprec*scanl/2D0)-tft
                             scanning=.true.
                             scani=0
                          endif
                      end select
                      if ((tft.gt.endoftim-lowTstep*2D0).and.(tftold.le.endoftim-lowTstep*2D0)) then
                          !step=hiTstep
                          boundarycheck=boundarycheck+1
                      endif

                      !write out predicted locatoin of focus
                      if((tftold.lt.endtime).and.(tft.gt.endtime)) then
                          xf(k,l)=y(1)
                          yf(k,l)=y(2)
                          endpoints(1,k,l)=xf(k,l)-pend
                          endpoints(2,k,l)=yf(k,l)
                      endif

                      !store info from last run
                      y1old=y(1)
                      tftold=tft

                  enddo

                  !if not all the boundaries were hit, throw an error
                  if(boundarycheck.ne.14) then
                     error=1
                  endif
                  boundarycheck=0

                  mem(1,1000,k+p*(l-1))=y(1)
                  mem(2,1000,k+p*(l-1))=y(2)

                  !close out ODE workspace
                  ido=3
                  call divprk(ido,neq,fcn,t,tend,tol,param,y)

                  tft=0d0

          enddo

          !calculate minimum width


          if (error.eq.1) write(6,*) "Invalid Test"
      enddo
      write(6,*) "dx=", maxval(xf(1:p,1:q))-minval(xf(1:p,1:q))
      write(6,*) "dE=",dE
      write(6,*) "dv=", dE/4D0*vini
      write(6,*) "dtheta=",dtheta
      if (error.eq.1) write(6,*) "Invalid Test", boundarycheck

      do test=0,scanl
         sums=0D0
         sumsq=0D0
         pulsecenter=scandata(midk,midl,test)
         do l=1,q
            do k=1,p
               scandata(k,l,test)=scandata(k,l,test)-pulsecenter
               sumsq=sumsq+scandata(k,l,test)**2D0
               sums=sums+scandata(k,l,test)
            enddo
         enddo

         call stddev(scandata(1:p,1:q,test),p*q,sd)
         scandx(test)=sd

         write(15,104) lowTstep*scanprec*(test-scanl/2D0)*vini,(scandata(1:p,1:q,test)-scandata((p+1)/2,(q+1)/2,test))/vini
         minimum=minval(scandata(1:p,1:q,test))
         maximum=maxval(scandata(1:p,1:q,test))
         bin_size=(maximum-minimum)/n_bins
         do i=1,n_bins
            hist(1,test,i)=count(scandata(1:p,1:q,test).le.bin_size*i+minimum.and.scandata(1:p,1:q,test).gt.bin_size*(i-1)+minimum)
            hist(2,test,i)=minimum+bin_size*i
         enddo

         write(12,101) lowTstep*scanprec*(test-scanl/2D0)*vini,scandx(test)/vini
      enddo
      do test=1,n_bins
         write(16,106) hist(1:2,:,test)
      enddo

      if(minval(scandx(0:scanl)).eq.scandx(0).or.minval(scandx(0:scanl)).eq.scandx(scanl)) then
         write(6,*) "WARNING!! False Minimum. increase scanl or scanprec"
      endif


      !write(6,*) "dx is also=", minval(scandx(0:scanl*2))

      do test=1,1000
          write(11,102) mem(1:2,test,1:p*q)
      enddo

      do l=1,q
         do k=1,p
            write(13,104) endpoints(1:2,k,l)
            !write(13,101) xf(k,l)-pend,yf(k,l)
         enddo
      enddo

      do test=0,scanl

      enddo


101   format(2E20.12)
102   format(18E20.12)
103   format(7(E12.7,','))
104   format(962E15.7)
106   format(201i5)
      return
      end

      subroutine stddev(x,npts,sd)
      integer npts,i
      real*8 x(npts),avg,sd,sums,sumsq
      sums = 0D0
      sumsq = 0D0

      sums=sum(x)
      avg=sums/npts
      do i=1,npts
         sumsq = sumsq + (x(i)-avg)**2
      enddo
      sd =  dsqrt(sumsq/npts)
      end
