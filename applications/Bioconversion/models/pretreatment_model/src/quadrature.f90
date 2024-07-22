module quadrature_module

      use polynomial_module
      implicit none

      contains

!===================================================================
subroutine trapzpointsandweights(x,w,np)
	
	real*8,allocatable,intent(out) :: w(:)
	real*8,allocatable,intent(out) :: x(:)

	integer,intent(out) :: np
	real*8 :: dx
	integer :: i

	np=10 !for now

	dx=2.d0/(np-1)

	allocate(w(np))
	allocate(x(np))

	w(1)=0.5*dx
	x(1)=-1.0

	w(np)=0.5*dx
	x(np)=1.0

	do i=2,np-1
		w(i)=dx
		x(i)=-1.d0+(i-1)*dx
	enddo
	
end subroutine trapzpointsandweights
!===================================================================
subroutine gausspointsandweights(x,w,npoints)

		integer,intent(in) :: npoints
		real*8, allocatable, intent(out) :: w(:)
		real*8, allocatable, intent(out) :: x(:)

		allocate(w(npoints))
		allocate(x(npoints))

		if(npoints .eq. 1) then
			w(1)=2.0
			x(1)=0.0
		endif

		if(npoints .eq. 2) then
			w(1)=1.0
			w(2)=1.0

			x(1)=-1.0/sqrt(3.0)
			x(2)=1.0/sqrt(3.0)
		endif

		if(npoints .eq. 3) then
			w(1)=5.0/9.0
			w(2)=8.0/9.0
			w(3)=5.0/9.0

			x(1)=-sqrt(3.0/5.0)
			x(2)=0.0
			x(3)=sqrt(3.0/5.0)
		endif

		if(npoints .eq. 4) then
			w(1) = (18.0-sqrt(30.0))/36.0
			w(2) = (18.0+sqrt(30.0))/36.0
			w(3) = (18.0+sqrt(30.0))/36.0
			w(4) = (18.0-sqrt(30.0))/36.0

			x(1) = -sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)
			x(2) = -sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)
			x(3) =  sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)
			x(4) =  sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)
		endif

end subroutine gausspointsandweights
!===================================================================
subroutine gll_gen(x,w,N)

      integer, intent(in) :: N
      real*8,allocatable, intent(out) :: x(:),w(:)

      real*8 :: dleg(N+1)
      real*8 :: pi,tol,x_it,xold
      integer :: maxit,N1,i,j,k

      allocate(x(N+1))
      allocate(w(N+1))

      pi = acos(-1.)

      tol = 1d-15

      N1 = N+1

      maxit = 1000   ! max iterations for newton-raphson
  
      x(1) = -1.d0
      x(N1) = 1.d0
 
      do i = 1, N+1

        x_it = -cos(pi * float(i-1) / N) ! initial guess - chebyshev points

        do j = 1, maxit 

    	  xold = x_it
          dleg(1) = 1.d0
          dleg(2) = x_it

	  do k = 2,N  

            dleg(k+1) = (  (2.d0*float(k) - 1.d0) * dleg(k) * x_it &
                    - (float(k)-1.d0)*dleg(k-1) ) / float(k)
          
          enddo

          x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                (float(N1) * dleg(N1) ) 

          if (abs(x_it - xold) .lt. tol) goto 10

        enddo
        write(*,*) 'max itertions reached'
10      continue

        x(i) = x_it
        w(i) = 2.d0 / (float(N * N1) * dleg(N1)**2 )

      enddo

      return

end subroutine gll_gen
!===================================================================
 subroutine gl_gen(x,w,n)
      integer n
      double precision x1,x2
      double precision eps
      parameter (eps=3.d-14)
      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1
      real*8,allocatable :: x(:),w(:)

      allocate(x(n))
      allocate(w(n))

      m=(n+1)/2

      x1 = -1.d0
      x2 = +1.d0

      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592653589793d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end
!===================================================================
subroutine globattopointsandweights(x,w,npoints)

		real*8,allocatable,intent(out) :: w(:),x(:)
		integer, intent(in) :: npoints

		allocate(w(npoints))
		allocate(x(npoints))

		if(npoints .eq. 2) then
			w(1)=1.0
			w(2)=1.0

			x(1)=-1.0
			x(2)=1.0
		endif

		if(npoints .eq. 3) then

			w(1)=1.0/3.0
			w(2)=4.0/3.0
			w(3)=1.0/3.0

			x(1)=-1.0
			x(2)=0.0
			x(3)=1.0
		endif

		if(npoints .eq. 4) then

			w(1)=1.0/6.0
			w(2)=5.0/6.0
			w(3)=5.0/6.0
			w(4)=1.0/6.0

			x(1)=-1.0
			x(2)=-sqrt(1.0/5.0)
			x(3)=sqrt(1.0/5.0)
			x(4)=1.0
		endif

		if(npoints .eq. 5) then

			w(1)=1.0/10.0
			w(2)=49.0/90.0
			w(3)=32.0/45.0
			w(4)=49.0/90.0
			w(5)=1.0/10.0

			x(1)=-1.0
			x(2)=-sqrt(3.0/7.0)
			x(3)=0.0
			x(4)=sqrt(3.0/7.0)
			x(5)=1.0
		endif


end subroutine globattopointsandweights
!===================================================================
subroutine getuniformlagrangepoints(n,x)

		integer :: n
		real*8 :: x(:)

		integer :: i
		real*8 :: delx
		real*8 :: xmin,xmax

		xmin = -1.d0
		xmax =  1.d0
		delx= (xmax-xmin)/n

		do i=1,n+1
			x(i)=xmin+(i-1)*delx
		enddo

end subroutine getuniformlagrangepoints
!===================================================================
subroutine convolvetrapz(p1,p2,a,b,npoints,answer)

	      type(polynomial), intent(in) :: p1,p2
	      real*8, intent(in) :: a,b
	      integer, intent(in) :: npoints
	      real*8 , intent(out) :: answer

	      integer :: i
	      real*8 :: delx,integval
	      real*8 :: h1,h2

	      delx=(b-a)/(npoints-1)
	      integval=0.0

	      do i=1,npoints-1

	      	h1 = polyfindval(p1,a+(i-1)*delx)*polyfindval(p2,a+(i-1)*delx)
		h2 = polyfindval(p1,a+i*delx)*polyfindval(p2,a+i*delx)

		integval = integval + 0.5*delx*(h1+h2)

	      enddo

	      answer=integval

end subroutine convolvetrapz
!===================================================================

end module quadrature_module
