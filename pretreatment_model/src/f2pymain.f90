!/* -*- Mode: F90; tab-width: 3 -*- */ 
! I suggest the use of an editor that follows default fortran indentation with
! spaces rather than tab characters, JJS 1/21/14

! this is indeed a useful way to define and access variables (scalar and
! arrays) from both python and fortran; however, it is not preferred practice
! for obtained an output array
! module ptmain

! 	real*8,allocatable :: pointsoln(:,:)
! 	real*8,allocatable :: interpsoln(:,:)

! contains

subroutine main(interpsolnout,m,n, inputfilename)

    use polynomial_module
    use inputs_module
    use exprparser_module 
    use element_module
    use outputs_module
    use timestepper_module
	
    implicit none

   real*8, intent(out) :: interpsolnout(m,n)
   integer, intent(in) :: m,n
   character(*), intent(in) :: inputfilename

    integer :: i
    real*8 :: delr,endp(2)

   real*8,allocatable :: pointsoln(:,:)
   real*8,allocatable :: interpsoln(:,:)

   call readinpfile(inpobj, inputfilename)

   allocate(elementarray(inpobj%nelements))

   delr=(inpobj%maxr-inpobj%minr)/(inpobj%nelements)

   do i=1,inpobj%nelements

   endp(1)=inpobj%minr+(i-1)*delr
   endp(2)=inpobj%minr+i*delr

   if(i .eq. 1) then
       call setelementattribs(elementarray(i),endp,inpobj,1)
   else if(i .eq. inpobj%nelements) then
       call setelementattribs(elementarray(i),endp,inpobj,2)
   else
       call setelementattribs(elementarray(i),endp,inpobj)
   endif
   enddo

   if(inpobj%timeorder .eq. 0) then 
       call timestepping_Implicit(inpobj,elementarray,pointsoln,interpsoln)
   else
       call timestepping_explicit(inpobj,elementarray,pointsoln,interpsoln)
   endif

   write(*,'(A)',advance='yes')
   print *,"======================================="

   ! put values in output array
   interpsolnout = interpsoln

   ! deallocate `elementarray` to prevent python crash?  JJS 1/13/21
   deallocate(elementarray)
   
   return
end subroutine main

! end module ptmain
