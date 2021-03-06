	PROGRAM compensator
c
c       a dispersion compensator for pulsed electron ?????
c        
c
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
	OPEN(UNIT=11,FILE='fast-dv-dtheta.dat',type='replace')       

c	write(6,*) ' stop integration at (micro s):'
c	read(*,*) endoftim

c	endoftim=endoftim*1.E-6 
	tol=.0000001D0
        
        write(6,*) "version 1.5"
	call parameters
	call result

	close(11)
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
	
	subroutine result
	implicit none
	common/trajectory/Bz,Ey
	real*8 Bz,Ey
	common/par/ pi,me,e,length,vini
	real*8 pi,me,e,length,vini
	common/tole/ tol
	real*8 tol
	integer count,test
	real*8 y1old
	common/time/ endoftim
	real*8 endoftim
	integer mxparm,neq,i,ido,p,q,k,l,z,m,time(6),lowstepcheck
	parameter (mxparm=120,neq=6,p=3,q=3)
	integer divider,r
        parameter (divider=9000,r=50)
        real*8 fcn,param(mxparm),t,tend,tft,y(neq),B,y0,Ed,E0,x1,x2,dw
	parameter (B=.031D0)
	!real*8 mem(2,1000,p*q),dataout(7,2,p),vout(1000,p,2),weintime(7,p)
        real*8 mem(2,r)
        real*8 pend,pos1,pos2,endtime,step,histep,lowstep,sls,comp(p,q)
        real*8 dw1,dw2
        parameter (pend=1D0)
	real*8 tftold,hiTstep,lowTstep,midstep
        external fcn,divprk,sset
c       !run test electron to get information about path
	write(6,*) "Initial Electron"
        length=me*vini/(e*B)*4D0*pi+pend
        endoftim=length/vini
        !lowTstep=endoftim/lowstep
        !sls=lowTstep*vini
	
        !lowstep= length/(((pend-dw)/4D0)/10)
        do m=1,r
        dw=pend*(1-m*1D-5)
        E0=-B*vini/4D0
                pos1=(pend-dw)/4D0
                pos2=pend-pos1
                dw1=pend/2D0-dw/2D0
                dw2=pend/2D0+dw/2D0
                lowstep=length/(pos1/10)
                midstep=endoftim*vini*1D6/dw
                lowTstep=endoftim/lowstep
                sls=lowTstep*vini
                histep=1000*endoftim/(sls/vini)
                hiTstep=endoftim/histep
	write(6,51) pos1,pos2
	write(6,52) dw1,dw2
51    format(x,'pos1=',E12.4,x,"pos2=",E12.7)
52    format(x,'dw1=',E12.4,x,'dw2=',E12.7)
        write(6,53) lowstep,histep
	write(6,54) 'sls=',sls,'lowDstep=',
     1lowTstep*vini,'hiDstep=',hiTstep*vini
53    format(x,'lowstep=',E12.4,x,'histep=',E12.4)
54    format(x,3(A,E12.4,x))
        write(6,'(x,A,F12.4,x,A,E12.4)') 
     1'length=',length,'end time=',endoftim
        ido=1
		step=lowstep
                t=0D0
	        Bz=0D0
	        Ey=0D0
		count=0
		do i=1,neq
		   y(i)=0D0
		enddo
		y(4)=vini	  !cos(.01*(k-(p+1)/2))
c		write(6,*) "v0=", y(4)
		!y(5)=vini*sin(.01*(k-(p+1)/2))	  !(1+0.01*(j-3))
c		y(2)=j*0.001   !(start at 11 position separated by 1 mm each)
c		write(6,*) y(2)
		call sset(mxparm,0.0,param,1)
		param(4)=1000000000
		param(10)=1.0   
		tft=0D0
		i=0
		do while(count.ne.8.or.y(1).lt.pend)
                !write(6,*) y(1),y(2)
			i=i+1
			tft=tft+endoftim/step
			tend=tft					    
			call divprk(ido,neq,fcn,t,tend,tol,param,y)
			
                        if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)
     1.and.(count.eq.0.or.count.eq.2)) then
                                step=histep
                                write(6,*) count, 1, i
                                lowstepcheck=lowstepcheck+1
                        endif
               if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)
     1.and.(count.eq.1)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 2, i
                        endif
               if ((y1old.le.dw1-sls)
     1.and.(y(1).gt.dw1-sls).and.count.eq.3) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count,3, i
                        endif
	       if ((y1old.le.dw2-sls)
     1.and.(y(1).gt.dw2-sls).and.count.eq.4) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count,4, i
                        endif
	         if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)
     1.and.(count.eq.6)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 5, i
                        endif
       		if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)
     1.and.(count.eq.7.or.count.eq.5)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 6, i
                        endif
c               	if ((y1old.le.pend-sls).and.(y(1).gt.pend-sls)
c     1.and.(count.eq.8)) then
        if ((tft.gt.endoftim-lowTstep).and.
     1(tftold.le.endoftim-lowTstep)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                               write(6,*) count, 7, i
                        endif

                        if ((y1old.le.pos1).and.(y(1).gt.pos1)
     1.and.(count.eq.0)) then
				count=1
                                step=lowstep
				time(1)=i
                                write(6,*) count,i
			endif
			if ((y1old.gt.pos1).and.(y(1).le.pos1)
     1.and.(count.eq.1)) then
				count=2
                                step=lowstep
				write(6,*) count

			endif
			if ((y1old.le.pos1).and.(y(1).gt.pos1)
     1.and.(count.eq.2)) then
				step=lowstep
				count=3
                                write(6,*) count
				time(2)=i
			endif
		if ((y1old.le.dw1).and.(y(1).gt.dw1)
     1.and.(count.eq.3)) then
                                write(6,*) pend/2-dw/2,pend/2+dw/2
					count=4
                                        step=lowstep
				write(6,*) count
					time(3)=i
					y0=y(2)
			endif
		if ((y1old.le.dw2).and.(y(1).gt.dw2)
     1.and.(count.eq.4)) then
					step=lowstep
				count=5
                                        write(6,*) count
					time(4)=i
			endif
			if ((y1old.le.pos2).and.(y(1).gt.pos2)
     1.and.(count.eq.5)) then
					step=lowstep
				count=6
                                        write(6,*) count
					time(5)=i
			endif
			if ((y1old.gt.pos2).and.(y(1).le.pos2)
     1.and.(count.eq.6)) then
					step=lowstep
				count=7
                                        write(6,*) count
			endif
			if ((y1old.le.pos2).and.(y(1).gt.pos2)
     1.and.(count.eq.7)) then
					step=lowstep
				count=8
                                        write(6,*) count
					time(6)=i
                                        x1=y(1)
			endif
							 
			if (count.eq.1) then 
				Bz=-B
			endif 
			if (count.eq.2) then 
				Bz=B
			endif 
			if (count.eq.3) then 
				Bz=0D0
			endif 
			if (count.eq.4) then 
				   !sign of wienfilter should be negative
				Ey=E0
				Bz=Ey/vini
			endif
			if (count.eq.5) then 
				Bz=0D0
				Ey=0D0
			endif 
			if (count.eq.6) then 
				Bz=B
				Ey=0D0
			endif 
			if (count.eq.7) then 
				Bz=-B
			endif 
			if (count.eq.8) then 
				Bz=0D0
			endif
			y1old=y(1)
                        tftold=tft
		enddo
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
        do l=1,q
	do k=1,p
		test=0
		ido=1
		t=0d0
	    Bz=0.
	    Ey=0.
		count=0
                lowstepcheck=0
		do i=1,neq
		   y(i)=0d0
		enddo
	y(4)=(1.00D0+2.5D-5*(l-(q+1)/2))*vini*cos(.001D0*(k-(p+1)/2)) !+vini*(k-(p+1)/2)*.005D0
	y(5)=(1.00D0+2.5D-5*(l-(q+1)/2))*vini*sin(.001D0*(k-(p+1)/2))	  !(1+0.01*(j-3))
c		y(2)=j*0.001   !(start at 11 position separated by 1 mm each)
c		write(6,*) y(2)
		write(6,50) k,p,l,q,m,r,Ed
50    format(x,i1,'/',i1,x,i1,'/',i1,x,i2,'/',i2,x,'Ed=',E12.2)
		call sset(mxparm,0.0,param,1)
		param(4)=1000000000
		param(10)=1D0 
		tft=0D0
		write(6,*) "Ed=",Ed,"v0=", y(4)
                i=0
                step=lowstep
                do while(tft.le.endtime)
                i=i+1
			tft=tft+endoftim/step
			tend=tft  
			call divprk(ido,neq,fcn,t,tend,tol,param,y)
               if ((y1old.le.pos1-sls).and.(y(1).gt.pos1-sls)
     1.and.(count.eq.0.or.count.eq.2)) then
                                step=histep
                                write(6,*) count, 1, i
                                lowstepcheck=lowstepcheck+1
                        endif
               if ((y1old.ge.pos1+sls).and.(y(1).lt.pos1+sls)
     1.and.(count.eq.1)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 2, i
                        endif
               if ((y1old.le.dw1-sls)
     1.and.(y(1).gt.dw1-sls).and.count.eq.3) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count,3, i
                        endif
	       if ((y1old.le.dw2-sls)
     1.and.(y(1).gt.dw2-sls).and.count.eq.4) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count,4, i
                        endif
	         if ((y1old.ge.pos2+sls).and.(y(1).lt.pos2+sls)
     1.and.(count.eq.6)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 5, i
                        endif
       		if ((y1old.le.pos2-sls).and.(y(1).gt.pos2-sls)
     1.and.(count.eq.7.or.count.eq.5)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                                write(6,*) count, 6, i
                        endif
c               	if ((y1old.le.pend-sls).and.(y(1).gt.pend-sls)
c     1.and.(count.eq.8)) then
        if ((tft.gt.endoftim-lowTstep).and.
     1(tftold.le.endoftim-lowTstep)) then
                                step=histep
                                lowstepcheck=lowstepcheck+1
                               write(6,*) count, 7, i
                        endif
               	if ((y1old.le.pos1).and.(y(1).gt.pos1)
     1.and.(count.eq.0)) then
				count=1
                                write(6,*) "step=", i
                                step=lowstep
c				dataout(1,1,k)=y(1)
c				dataout(2,1,k)=y(2)
c				dataout(3,1,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c				dataout(4,1,k)=tft
                                !weintime(1,k)=tft
			endif
			if ((y1old.gt.pos1).and.(y(1).le.pos1).
     1and.(count.eq.1)) then
				count=2
                                write(6,*) "step=", i
                                step=lowstep 
			endif
			if ((y1old.le.pos1).and.(y(1).gt.pos1).
     1and.(count.eq.2)) then
				count=3	
                                write(6,*) "step=", i
                                step=lowstep
c				dataout(1,2,k)=y(1)
c				dataout(2,2,k)=y(2)
c				dataout(3,2,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c				dataout(4,2,k)=tft
                                !weintime(2,k)=tft
			endif
			if ((y1old.le.dw1).and.(y(1).gt.
     1dw1).and.(count.eq.3)) then
			y(4)=y(4)+e*Ed*(y(2)-y0)/(me*vini)
					count=4
                                 write(6,*) "step=", i
                                        step=lowstep
               !write(6,*) "enter field",e*Ed*(y(2)-y0)/(me*vini)

c			dataout(1,3,k)=y(1)
c			dataout(2,3,k)=y(2)
c			dataout(3,3,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c			dataout(4,3,k)=tft
                        !weintime(3,k)=tft
			endif
			if ((y1old.le.dw2).and.(y(1).gt.
     1dw2).and.(count.eq.4)) then
			y(4)=y(4)-e*Ed*(y(2)-y0)/(me*vini)	
		!write(6,*) "leave field",-e*Ed*(y(2)-y0)/(me*vini)
					count=5
                                  write(6,*) "step=", i
                                       step=lowstep
c				dataout(1,4,k)=y(1)
c				dataout(2,4,k)=y(2)
c				dataout(3,4,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c				dataout(4,4,k)=tft
                                !weintime(4,k)=tft
			endif
			if ((y1old.le.pos2).and.(y(1).gt.pos2).
     1and.(count.eq.5)) then
					count=6
                                 write(6,*) "step=", i
                                        step=lowstep
c				dataout(1,5,k)=y(1)
c				dataout(2,5,k)=y(2)
c				dataout(3,5,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c				dataout(4,5,k)=tft
                                !weintime(5,k)=tft
			endif
			if ((y1old.gt.pos2).and.(y(1).le.pos2).
     1and.(count.eq.6)) then
					count=7
                                 write(6,*) count,"step=", i
                                        step=lowstep
			endif
			if ((y1old.le.pos2).and.(y(1).gt.pos2)
     1.and.(count.eq.7)) then
					count=8
                                 write(6,*) "step=", i
                                        step=lowstep
c				dataout(1,6,k)=y(1)
c				dataout(2,6,k)=y(2)
c				dataout(3,6,k)=sqrt(y(4)*y(4)+y(5)*y(5))
c				dataout(4,6,k)=tft
                                !weintime(6,k)=tft
			endif
		
			if (count.eq.1) then 
				Bz=-B
			endif 
			if (count.eq.2) then 
				Bz=B
			endif 
			if (count.eq.3) then 
				Bz=0D0
			endif 
			if (count.eq.4) then 
                                !sign of wienfilter should be negative
				Ey=Ed
				Bz=Ey/vini
			endif
			if (count.eq.5) then 
				Bz=0D0
				Ey=0D0
			endif 
			if (count.eq.6) then 
				Bz=B
				Ey=0D0
			endif 
			if (count.eq.7) then 
				Bz=-B
			endif 
			if (count.eq.8) then 
				Bz=0D0
			endif
			y1old=y(1)
                        tftold=tft
			if (i/divider.ge.test.and.test.lt.1000) then
				test=test+1
	
c				write (6,*) y(1),y(2),test,k,count
				!mem(1,test,k+3*(l-1))=y(1)
				!mem(2,test,k+3*(l-1))=y(2)
                                !vout(test,k,1)=y(4)
                                !vout(test,k,2)=y(5)
			endif
	
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
                if(lowstepcheck.ne.9) then
                        write(6,*) "Invalid Test", lowstepcheck
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
                !weintime(7,k)=tft
                ido=3
                comp(k,l)=y(1)
c		write(6,*) "test2" 
                call divprk(ido,neq,fcn,t,tend,tol,param,y)
	
	    tft=0d0

    	enddo
        enddo
        write(6,*) maxval(comp(1:p,1:q))-minval(comp(1:p,1:q))
        write(11,101) pend-dw, 
     1maxval(comp(1:p,1:q))-minval(comp(1:p,1:q))
        enddo
c	do test=1,p
c	!	write(12,102) dataout(1:7,1,test)
c	
c                write(13,102) dataout(1:7,2,test)
c                write(15,102) weintime(1:7,test)
c	enddo

c      write(12,103) mem(1,1000,1:p)
c	write(6,*) "count=",count

101   format(2E15.7)
102  	format(14E15.7)
103   format(7(E12.7,','))
	return
	end
