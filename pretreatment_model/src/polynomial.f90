module polynomial_module

use globalvars_module
implicit none


contains

!========================================================
subroutine polyinitialize(pl,n,coeffarray)

	      integer :: n
	      real*8 :: coeffarray(:)
	      integer :: i

	      type(polynomial) :: pl

	      pl%n=n
	      allocate(pl%coeffs(n+1))

	      do i=1,n+1
	      	pl%coeffs(i)=coeffarray(i)
	      enddo

end subroutine
!========================================================
subroutine zeropolyinitialize(pl,n)

		integer,intent(in) :: n
		type(polynomial) :: pl

		integer :: i
		real*8,allocatable :: coeffarray(:)

		allocate(coeffarray(n+1))

		do  i=1,n+1
			coeffarray(i)=0.d0
		enddo

		call polyinitialize(pl,n,coeffarray)	

end subroutine zeropolyinitialize
!========================================================
subroutine polymultiply(p1,p2,pout)

		type(polynomial), intent(in) :: p1,p2
		type(polynomial), intent(out) :: pout
		integer :: i,j

		call zeropolyinitialize(pout,p1%n+p2%n)

		do i=1,p1%n+1
		do j=1,p2%n+1
			pout%coeffs(i+j-1)=pout%coeffs(i+j-1)&
					+p1%coeffs(i)*p2%coeffs(j)
		enddo
		enddo

end subroutine polymultiply
!========================================================
subroutine polyderivative(pin,pout)

	type(polynomial),intent(in) :: pin
	type(polynomial),intent(out) :: pout

	integer :: i

	pout%n=pin%n-1
	allocate(pout%coeffs(pout%n+1))

	do i=1,pout%n+1
		pout%coeffs(i)=pin%coeffs(i+1)*i
		!print *,"coeff of ",i,":",pout%coeffs(i)
	enddo	

end subroutine polyderivative
!========================================================
subroutine printpoly(pl)

		type(polynomial) :: pl
		integer :: i

		print *,"order of pl:",pl%n

		do i=1,pl%n
			write(*,'(F4.2,A,I2,A)',advance='no')pl%coeffs(i),"x^",(i-1),"+"
		enddo
		write(*,'(F8.2,A,I2)')pl%coeffs(i),"x^",(i-1)

end subroutine printpoly
!========================================================
subroutine createlagrangepoly(x,np,polyarray)
		
	integer :: i,j,np
	real*8 :: x(np)
	real*8 :: dnmtr
	type(polynomial) :: polyarray(np)

	real*8 :: initcoeff(1)
	real*8 :: termcoeff(2)
	real*8,allocatable :: lagrcoeffs(:)
	type(polynomial) :: term
	type(polynomial) :: temp
	type(polynomial) :: current

	call zeropolyinitialize(term,1)

	do i=1,np
		
		dnmtr=1.d0
		do j=1,np
			if(j .ne. i) then
			dnmtr=dnmtr*(x(i)-x(j))
			endif
		enddo

		call zeropolyinitialize(polyarray(i),0)
		polyarray(i)%coeffs(1)=1.d0/dnmtr
		do j=1,np

			if(j .ne. i) then
				
				term%coeffs(1)=-x(j)
				term%coeffs(2)=1.d0
				call polymultiply(polyarray(i),term,temp)

				polyarray(i)=temp
			endif
		enddo

		print *,"order of polyarray",i,"is",polyarray(i)%n

	enddo
			

end subroutine createlagrangepoly
!========================================================
function polyfindval(pl,x) result(total)

	      real*8 :: x,total,pdt
	      type(polynomial) :: pl
	      integer :: i

	      total=0.0
	      pdt=1.0

	      do i=1,pl%n+1

	      	total=total+pl%coeffs(i)*pdt
		pdt=pdt*x

	      enddo

end function
!========================================================
function polyfindderval(pl,x) result(total)

	      real*8 :: x,total,pdt
	      type(polynomial) :: pl
	      type(polynomial) :: plder
	      integer :: i

	      call polyderivative(pl,plder)
	      total=0.0
	      pdt=1.0

	      do i=1,plder%n+1

	      	total=total+plder%coeffs(i)*pdt
		pdt=pdt*x

	      enddo

end function
!========================================================

end module polynomial_module
