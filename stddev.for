      Program test
      implicit none
      integer npts, i
      parameter (npts=100)
      real*8 x(100),sd
      
      do i =1,npts
         x(i)=sin(i/2D0)
      enddo

      call stddev(x,npts,sd)
      write(*,*) sd
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
      sd =  sqrt(sumsq/npts)
      end


