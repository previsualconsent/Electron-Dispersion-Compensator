      PROGRAM test
      !
      ! a dispersion compensator for pulsed electron ?????
      !  
      !
      use msimsl


      IMPLICIT NONE
      integer i
      write(6,*) "version 1.5"
      i=1
      select case (i)
        case (1)
            write(*,*) 1
            i=2
        case (2)
            write(*,*) 2
      end select
      write(*,*) i

      STOP
      END


