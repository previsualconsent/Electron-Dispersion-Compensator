      PROGRAM test
      IMPLICIT NONE
      real*8 me,e,B0,dw,dtot,vx,vini,E0,B,Bz,y
      integer i

      OPEN(UNIT=11,FILE='acc.dat',type='replace')
      me=9.10938188D-31
      e=1.60217646D-19

      vini=3.23D1
      vx=32.280885206547300D0

      dtot=1D5
      dw=.5D5

      B=1.5D-14
      E0=-(dtot*B*vini)/(dw*4D0)
      B0=E0/vini


      write(6,*) e*(-vx*B0+E0)/me

      write(6,*) e,vx,B0,E0,me

      do i=-20,20
         y=i*1D0
         Bz=E0/(vini+(B/4+E0/vini)*e*(y)/me)
         write(6,101) y,((-vx*Bz)+E0)*e/me
         write(11,101) y,((-vx*Bz)+E0)*e/me
      enddo

      close(11)
101   format(2E15.7)
      STOP
      END
