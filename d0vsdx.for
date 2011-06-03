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
      OPEN(UNIT=11,FILE='tuning.dat',type='replace')
      OPEN(UNIT=14,FILE='width.dat',type='replace')        

      !write(6,*) ' stop integration at (micro s):'
      !read(*,*) endoftim

      !endoftim=endoftim*1.E-6 
      tol=.0000001D0

      write(6,*) "version 1.5"
      write(6,*) char(7)
      call parameters
      call result
      write(6,*) char(7)

      close(11)
      close(14)
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


      Ex=0D0					 
      Ez=0D0
      Bx=0D0
      By=0D0


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
      real*8 endoftim
      integer mxparm,neq,i,ido,p,q,r,k,l,m,boundarycheck
      parameter (mxparm=120,neq=6,p=3,q=3,r=11)
      integer error, ti, scanl,scanprec,midk,midl
      parameter (scanl=400,scanprec=15)
      real*8 fcn,param(mxparm),t,tend,tft,y(neq),B,y0,E0,x1,dw,dE,dEm,dtheta,dthetam
      parameter (B=3.1D-14,dw=.5D5,dE=2.5D-5,dtheta=2D-4)
      real*8 pend,pos1,pos2,endtime,step,sls
      real*8 dw1,dw2,ddw,wfls,bls,yold,theta
      parameter (pend=1D5,ddw=1D-4)
      real*8 tftold,hiTstep,lowTstep,wfTstep,bTstep,endTstep,scandata(p,q,0:scanl*2),scandx(1:r,0:scanl*2)
      real*8 timecheck,percent,outdata(3,r)
      logical debugout
      parameter (debugout=.false.)
      
      external fcn,divprk,sset
      !initial calculations
      midk=(p+1)/2
      midl=(q+1)/2
      length=me*vini/(e*B)*4D0*pi+pend   !total length based on path
      endoftim=length/vini          !total time for travel

      lowTstep=endoftim/1D5      !low precision time step
      sls=lowTstep*vini             !average distance for 1 lowTstep

      write(6,*) "Accuracy:",sls*scanprec
      bTstep=lowTstep*1D1           !time step for B-fields
      bls=bTstep*vini               !average distance of 1 bTstep

      wfTstep=lowTstep*1D2         !time step for WF
      wfls=wfTstep*vini             !average distance for 1 wfTstep

      hiTstep=lowTstep*1D-3          !hi-res step for approaching borders 
      endTstep=lowTstep*1D-2    !set step size for final section

      do m=1,r
          dthetam=dtheta**(1+(m-(r+1D0)/2D0)*1D-2)
          E0=-(pend*B*vini)/(dw*4D0)    !set E field for WF based on theory

          !some distance calculations
          pos1=(pend-dw)/4D0 !first B-field border
          pos2=pend-pos1     !second B-field border
          dw1=pend/2D0-dw/2D0 !start of Wien filter
          dw2=pend/2D0+dw/2D0 !end of Wien filter

51    format(x,'pos1=',E12.4,x,"pos2=",E12.7)
52    format(x,'dw1=',E12.4,x,'dw2=',E12.7)
53    format(x,'lowTstep=',E12.4,x,'hiTstep=',E12.4)
54    format(x,3(A,E12.4,x))
          if (debugout) then

              write(6,51) pos1,pos2 !write info to console
              write(6,52) dw1,dw2
              write(6,53) lowTstep,hiTstep
              write(6,54) 'sls=',sls,'lowDstep=',lowTstep*vini,'hiDstep=',hiTstep*vini
              write(6,'(x,A,F12.4,x,A,E12.4)') 'length=',length,'end time=',endoftim

          endif   

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

60    format(x,'counter=',i1,x,'flag=',i1,x,'steps=',i10)
61    format(x,'counter=',i1,8x,'steps=',i10)

              select case (counter) !switch on counter- this reduced long if statements
                                    !each case checks to see if the particle is getting
                                    !close to the next barrier. If so, reduce the step
                                    !size. When it crosses the barrier, hand off to 
                                    !next case. boundary check is used for error checking
              case (0) 
                  Bz=0D0
                  Ey=0D0
                  if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                      step=hiTstep
                      if (debugout) write(6,60) counter,0,i
                      boundarycheck=boundarycheck+1
                  endif
                  if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                      if (debugout) write(6,61) counter,i
                      counter=1
                      step=bTstep
                      Bz=-B
                  endif
              case (1)
                  if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls)) then
                      if (debugout) write(6,60) counter, 1, i
                      step=lowTstep
                      boundarycheck=boundarycheck+1
                  endif
                  if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)) then
                      if (debugout) write(6,60) counter, 1, i
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                  endif
                  if ((y1old.gt.pos1).and.(y(1).le.pos1)) then
                      counter=2
                      Bz=B
                      step=bTstep
                      if (debugout) write(6,61) counter,i
                  endif
              case (2)
                  if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                      step=hiTstep
                      if (debugout) write(6,60) counter,2,i
                      boundarycheck=boundarycheck+1
                  endif
                  if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls)) then
                      step=lowTstep
                      if (debugout) write(6,60) counter,2,i
                      boundarycheck=boundarycheck+1
                  endif
                  if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                      step=lowTstep
                      counter=3
                      Bz=0D0
                      if (debugout) write(6,61) counter,i
                  endif

              case (3)
                  if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls)) then
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter,3, i
                  endif
                  if ((y1old.le.dw1).and.(y(1).gt.dw1)) then
                      if (debugout) write(6,*) pend/2-dw/2,pend/2+dw/2
                      counter=4
                      Ey=E0
                      Bz=Ey/vini
                      if (debugout) write(6,*) Ey,Bz*y(4)
                      step=wfTstep
                      if (debugout) write(6,61) counter,i
                      y0=y(2)                   ! y0 has the center of the WF
                      if (debugout) write(6,*) "y0=",y0
                  endif

              case (4)

                  if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls)) then
                      step=lowTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter,4, i
                  endif
                  if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls)) then
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter,4, i
                  endif
                  if ((y1old.le.dw2).and.(y(1).gt.dw2)) then
                      step=lowTstep
                      counter=5
                      Bz=0D0
                      Ey=0D0
                      if (debugout) write(6,61) counter,i
                  endif

              case (5)
                  if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter, 6, i
                  endif
                  if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                      step=bTstep
                      counter=6
                      Bz=B
                      Ey=0D0
                      if (debugout) write(6,61) counter,i
                  endif

              case (6)
                  if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls)) then
                      step=lowTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter, 5, i
                  endif
                  if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)) then
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter, 5, i
                  endif
                  if ((y1old.gt.pos2).and.(y(1).le.pos2)) then
                      step=bTstep
                      counter=7
                      Bz=-B
                      if (debugout) write(6,61) counter,i
                  endif
              case (7)
                  if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls)) then
                      step=lowTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter, 6, i
                  endif
                  if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                      step=hiTstep
                      boundarycheck=boundarycheck+1
                      if (debugout) write(6,60) counter, 6, i
                  endif
                  if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                      step=lowTstep
                      counter=8
                      Bz=0D0
                      timecheck=tft+lowTstep*scanprec*3 !grabs a timestamp from the beginning of region 8
                      if (debugout) write(6,61) counter,i
                  endif
              end select
              !lower step size while approaching final region.
              if ((tft.gt.endoftim-lowTstep).and.(tftold.le.endoftim-lowTstep)) then
                  step=hiTstep
                  boundarycheck=boundarycheck+1
                  if (debugout) write(6,60) counter, 7, i
              endif

              !store info from last run
              y1old=y(1)
              tftold=tft
          enddo

          ido=3
          endtime=tft
          call divprk(ido,neq,fcn,t,tend,tol,param,y) !close out ODE workspace

          if(debugout) write(6,*) "boundarycheck=",boundarycheck


          !first electron done, do the rest
          do l=1,q
              do k=1,p
                  !reset variables for next run
                  test=0
                  ido=1
                  t=0D0
                  ti=0
                  Bz=0.
                  Ey=0.
                  counter=0
                  boundarycheck=0
                  do i=1,neq
                      y(i)=0d0
                  enddo

                  !sets energy spread to dE and angle spread to dtheta
                  y(4)=(1.00D0+dE/4D0*(l-midl))*vini*cos(dthetam/2D0*(k-midk))
                  y(5)=(1.00D0+dE/4D0*(l-midl))*vini*sin(dthetam/2D0*(k-midk))

                  write(6,50) k,p,l,q,m,r
50    format(x/,x,i1,'/',i1,x,i1,'/',i1,x,i3,'/',i3)
                  call sset(mxparm,0.0,param,1)
                  param(4)=1000000000
                  param(10)=1D0 
                  tft=0D0
                  i=0
                  step=lowTstep

62    format(18x,A,i10)

                  ! loop until the scan has finished
                  do while(tft.le.timecheck+2D0*scanl*lowTstep*scanprec)
                      i=i+1
                      tft=tft+step
                      tend=tft  
                      call divprk(ido,neq,fcn,t,tend,tol,param,y)

                      !same switch as above
                      select case (counter)
                      case (0)
                          Bz=0D0
                          Ey=0D0
                          if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                              step=hiTstep
                              if (debugout) write(6,60) counter,0,i
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                              if (debugout) write(6,61) counter,i
                              counter=1
                              Bz=-B
                              step=bTstep
                          endif
                      case (1)
                          if ((y1old.ge.pos1+bls).and.(y(1).lt.pos1+bls)) then
                              if (debugout) write(6,60) counter, 1, i
                              step=lowTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)) then
                              if (debugout) write(6,60) counter, 1, i
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.gt.pos1).and.(y(1).le.pos1)) then
                              counter=2
                              step=bTstep
                              Bz=B
                              if (debugout) write(6,61) counter,i
                          endif
                      case (2)
                          if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)) then
                              step=hiTstep
                              if (debugout) write(6,60) counter,2,i
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1-bls).and.(y(1).gt.pos1-bls)) then
                              step=lowTstep
                              if (debugout) write(6,60) counter,2,i
                              boundarycheck=boundarycheck+1
                          endif
                          if ((y1old.le.pos1).and.(y(1).gt.pos1)) then
                              step=lowTstep
                              counter=3
                              if (debugout) write(6,61) counter,i
                              Bz=0D0
                          endif

                      case (3)
                          if ((y1old.le.dw1-sls).and.(y(1).gt.dw1-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter,3, i
                          endif
                          if ((y1old.le.dw1).and.(y(1).gt.dw1)) then
                              !if (debugout) write(6,*) "dv=", e*B*(y(2)-y0)/(4*me)
                              !if (debugout) write(6,*) "dvg=", E0*e*(y(2)-y0)/(2*me*vini)

                              !ENTERING WIEN FILTER
                              !calculate old velocity and angle
                              yold=sqrt(y(4)**2+y(5)**2)
                              theta=ATAN(y(5)/y(4))

                              !set new velocity based on calculated energy change
                              y(4)=sqrt(yold**2+2D0*E0*e*(y(2)-y0)/me)*cos(theta)
                              y(5)=sqrt(yold**2+2D0*E0*e*(y(2)-y0)/me)*sin(theta)

                              counter=4
                              Ey=E0
                              !Matched B field. For particle traveling strait through, the field should be balanced
                              Bz=E0*(1/vini-((y(2)-y0)*(e*(B*vini+4D0*E0))/(4D0*me*vini**3)+((y(2)-y0)*(e*(B*vini+4D0*E0))/(4D0*me))**2*(1/vini**6))*1D0)
                              if (debugout) write(6,62) "steps=", i
                              step=wfTstep
                              if (debugout) write(6,*) "enter field",e*E0*(y(2)-y0)/(me*vini)
                          endif

                      case (4)
                          if ((y1old.le.dw2-wfls).and.(y(1).gt.dw2-wfls)) then
                              step=lowTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter,4, i
                          endif
                          if ((y1old.le.dw2-sls).and.(y(1).gt.dw2-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter,4, i
                          endif
                          if ((y1old.le.dw2).and.(y(1).gt.dw2)) then
                              !calculate velocity change for leaving WF
                              yold=sqrt(y(4)**2+y(5)**2)
                              theta=ATAN(y(5)/y(4))
                              y(4)=sqrt(yold**2-2D0*E0*e*(y(2)-y0)/me)*cos(theta)
                              y(5)=sqrt(yold**2-2D0*E0*e*(y(2)-y0)/me)*sin(theta)
                              if (debugout) write(6,*) "leave field",-e*E0*(y(2)-y0)/(me*vini)
                              counter=5
                              Bz=0D0
                              Ey=0D0
                              if (debugout) write(6,62) "steps=", i
                              step=lowTstep
                          endif

                      case (5)
                          if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter, 6, i
                          endif
                          if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                              step=bTstep
                              counter=6
                              Bz=B
                              Ey=0D0
                              if (debugout) write(6,61) counter,i
                          endif

                      case (6)
                          if ((y1old.ge.pos2+bls).and.(y(1).lt.pos2+bls)) then
                              step=lowTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter, 5, i
                          endif
                          if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter, 5, i
                          endif
                          if ((y1old.gt.pos2).and.(y(1).le.pos2)) then
                              step=bTstep
                              counter=7
                              Bz=-B
                              if (debugout) write(6,61) counter,i
                          endif
                      case (7)
                          if ((y1old.le.pos2-bls).and.(y(1).gt.pos2-bls)) then
                              step=lowTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter, 6, i
                          endif
                          if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)) then
                              step=hiTstep
                              boundarycheck=boundarycheck+1
                              if (debugout) write(6,60) counter, 6, i
                          endif
                          if ((y1old.le.pos2).and.(y(1).gt.pos2)) then
                              step=endTstep
                              counter=8
                              Bz=0D0
                              if (debugout) write(6,61) counter,i
                              x1=y(1)
                          endif

                      end select
                      if ((tft.gt.endoftim-lowTstep*2D0).and.(tftold.le.endoftim-lowTstep*2D0)) then
                          !step=hiTstep
                          boundarycheck=boundarycheck+1
                          if (debugout) write(6,60) counter, 7, i
                      endif

                      !taking time snapshots to calculate widths
                      if ((tft.gt.timecheck+lowTstep*scanprec*ti).and.(tftold.le.timecheck+lowTstep*scanprec*ti)) then

                          scandata(k,l,ti)=y(1)

                          !if (debugout) write(6,*) ti,counter, tft, scandata(k,l,ti)
                          if (ti.ne.scanl*2) then
                              ti=ti+1
                          endif
                      endif

                      !write out predicted locatoin of focus
                      if((tftold.le.endtime).and.(tft.gt.endtime)) then
                          if (debugout) write(6,*) y(1)
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
                  if (debugout) write(6,*) "end step=", i

                  !close out ODE workspace
                  ido=3
                  call divprk(ido,neq,fcn,t,tend,tol,param,y)

                  tft=0d0

              enddo
          enddo

          !calculate width of beam at diff timestamps
          do test=0,scanl*2
              scandx(m,test)=maxval(scandata(1:p,1:q,test))-minval(scandata(1:p,1:q,test))
              !if (debugout) write(6,*) test, scandx(m,test)
          enddo
          !if minimum occurs at end of scan, print error
          if(minval(scandx(m,0:scanl*2)).eq.scandx(m,scanl*2)) then
              write(6,*) "WARNING!! False Minimum. increase scanl or scanprec"
          endif

          !if minimum occurs at the begninning of scan, print error
          if(minval(scandx(m,0:scanl*2)).eq.scandx(m,0)) then
              write(6,*) "Does not converge"
          endif

          !calculate minimum width
          write(6,*) "dx=", minval(scandx(m,0:scanl*2))
          outdata(1,m)=dthetam
          outdata(2,m)=minval(scandx(m,0:scanl*2))

          !find the location of the minimum width
          do test=0,scanl*2
              if (minval(scandx(m,0:scanl*2)).eq.scandx(m,test)) then
                  outdata(3,m)=scandata(midk,midl,test)
              endif
          enddo
          

          if (error.eq.1) write(6,*) "Invalid Test", boundarycheck
      enddo

      !write out the data
      do test=1,r
          write(11,101) outdata(1:3,test) !percent, dx, pos
      enddo
      do test=0,scanl*2
          write(14,104) scandx(1:r,test) !shows how the pulse width of each run
      enddo                              !changed over time


101   format(3E15.7)
102   format(18E15.7)
103   format(7(E12.7,','))
104   format(101E15.7)
      return
      end
