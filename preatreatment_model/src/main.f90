program pretreatmentmodeling

	use polynomial_module
	use inputs_module
	use exprparser_module 
	use element_module
	use outputs_module
	use timestepper_module
	
	implicit none

	integer :: i
	real*8 :: delr,endp(2)
	real*8,allocatable :: pointsoln(:,:)
	real*8,allocatable :: interpsoln(:,:)

	call readinpfile(inpobj,'pretreat_defs.inp')

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

end program pretreatmentmodeling
