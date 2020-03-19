program testpoly

      use polynomial_module
!      use quadrature_module
      implicit none

      type(polynomial) :: p1
      type(polynomial) :: p2
      type(polynomial) :: p1der1,p1der2,p1der3
      real*8 :: coeffarray(3),coeffs(2)
      real*8 :: cmul1(2),cmul2(2)
      real*8 :: val
      type(polynomial) :: pmult
      integer :: i
      type(polynomial) :: pmul1,pmul2
      type(polynomial) :: temp
      type(polynomial) :: parray(4)
      real*8 :: endp(4)

      coeffarray(1)=1.0
      coeffarray(2)=2.0
      coeffarray(3)=3.0

      call polyinitialize(p1,2,coeffarray)
      
      coeffs(1)=-5.0
      coeffs(2)=2.0
      
      call polyinitialize(p2,1,coeffs)


      cmul1(1)=1.0
      cmul1(2)=-1.0

      cmul2(1)=1.0
      cmul2(2)=1.0


      call polyinitialize(pmul1,1,cmul1)
      call polyinitialize(pmul2,1,cmul2)
      !call zeropolyinitialize(pmult,2)

      print *,"init pmult"
      call polymultiply(pmul1,pmul2,pmult)
      print *,pmult%coeffs(1:pmult%n+1)
      
      !call zeropolyinitialize(temp,3)
      print *,"init temp"
      call polymultiply(pmult,p2,temp)
      print *,"multiplied"

      pmult=temp

      print *,pmult%coeffs(1:pmult%n+1)
      

      endp(1)=-1.d0
      endp(2)=-0.33
      endp(3)=0.33
      endp(4)=1.d0

      print *,"lagrange polynomials"
      call createlagrangepoly(endp,4,parray)

      do i=1,4
        print *,"order =",parray(i)%n
      	print *,parray(i)%coeffs(1:parray(i)%n+1)
      enddo

    !  call convolvetrapz(p1,p2,0.d0,1.d0,10,val)
    !  print *,"integral = ",val
    !  call convolvetrapz(p1,p2,0.d0,1.d0,100,val)
    !  print *,"integral = ",val

      call polyderivative(p1,p1der1)
      call polyderivative(p1der1,p1der2)
      call polyderivative(p1der2,p1der3)

      print *,"p1der order :",p1der3%n
      print *,"val at 5=",polyfindval(p1der3,5.d0)

end program testpoly
