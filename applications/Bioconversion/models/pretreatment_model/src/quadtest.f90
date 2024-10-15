program quadtest

      use quadrature_module

      implicit none

      integer :: npoints
      real*8,allocatable :: x(:)
      real*8,allocatable :: w(:)
      integer :: i

      npoints=1

      call gll_gen(x,w,npoints)

      do i=1,npoints+1
      	
      		print *,x(i),w(i)

	enddo

end program quadtest
