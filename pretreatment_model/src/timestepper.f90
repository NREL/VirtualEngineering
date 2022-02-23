module timestepper_module

  use inputs_module
  use outputs_module
  use element_module
  use solvergmres_module

  implicit none

contains
  !========================================================================================================
  subroutine timestepping_Explicit(inpobj,elementarray,pointsoln,interpsoln)

    type(inputdata), intent(inout)  :: inpobj
    type(element), intent(inout) :: elementarray(:)
    real*8                       :: l2norm(globalnumeq),totaltime,printtime
    integer                      :: iterations,df,i,neq,printit
    logical                      :: isitRK
    integer                      :: numstages
    integer                      :: rkstage
    real*8                       :: multiplier
    real*8,allocatable           :: pointsoln(:,:)
    real*8,allocatable           :: interpsoln(:,:)


    if(inpobj%timeorder > 1) then
	    isitRK = .true.
	    numstages = 5
    else
    	    isitRK = .false.
            numstages = 1
    endif

    call printsolution(elementarray,inpobj%nelements,0,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,pointsoln)
    call printpolysolution(elementarray,inpobj%nelements,0,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,interpsoln)
    totaltime  = 0.d0
    printtime  = 0.d0
    iterations = 0
    printit    = 1

    do while(totaltime < inpobj%finaltime)

       totaltime  = totaltime  + inpobj%timestep
       printtime  = printtime  + inpobj%timestep       
       iterations = iterations + 1
       
       call applydirichletbc(elementarray(inpobj%nelements),inpobj)

       !save prvs soln if RK
       do i=1,inpobj%nelements
	      call saveprvssoln(elementarray(i))
       enddo

	if(isitRK) then

		do rkstage=1,numstages
			multiplier=rkcoeffs(rkstage)
			call compute_explicitrhs(elementarray,inpobj,totaltime)
			call dotimestep_explicit(elementarray,inpobj,multiplier)
			print *,"rkstage:",rkstage
		enddo

	else
		multiplier=1.d0
		call compute_explicitrhs(elementarray,inpobj,totaltime)
		call dotimestep_explicit(elementarray,inpobj,multiplier)
	endif

       call calc_l2norm_of_residual(elementarray,inpobj%nelements,l2norm)
   
       do i=1,globalnumeq
          write(*,'(F18.11)',advance='no')l2norm(i)
       enddo
       write(*,'(F18.11)',advance='yes')

       if(printtime*inpobj%tscale .ge. inpobj%printtime) then
    	
	call printsolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
			inpobj%xscale,inpobj%scalefactors,pointsoln)

    	call printpolysolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
			inpobj%xscale,inpobj%scalefactors,interpsoln)
	
	print *
	print *,"saving solution at time:",totaltime
	printtime=0.d0
	printit=printit+1

	endif

       print *,"totaltime:",totaltime*inpobj%tscale,"iterations:",iterations
       print *,"======================================="

    enddo

    print *,"final time:",totaltime*inpobj%tscale
    call printsolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,pointsoln)
    call printpolysolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,interpsoln)


  end subroutine timestepping_Explicit
  !========================================================================================================
  subroutine compute_explicitrhs(elementarray,inpobj,totaltime)
    
    type(inputdata), intent(inout)  :: inpobj
    type(element), intent(inout) :: elementarray(:)
    real*8, intent(in) :: totaltime
    integer :: i,neq

      !boundary contributions
       do i=1,inpobj%nelements

          call computemassmatX(elementarray(i),inpobj,totaltime)
          call computediffmatX(elementarray(i),elementarray(i)%solnvec,inpobj,totaltime)

       enddo
       !      found residuals and mass matrix

       ! update residuals at internal boundaries
       do i=1,inpobj%nelements-1

          do neq=1,globalnumeq

             elementarray(i)%residual(elementarray(i)%p+1,neq) = elementarray(i)%residual(elementarray(i)%p+1,neq) + &
                  elementarray(i+1)%residual(1,neq)
             elementarray(i+1)%residual(1,neq) = elementarray(i)%residual(elementarray(i)%p+1,neq)
             elementarray(i)%massmatdiag(elementarray(i)%p+1,neq) = elementarray(i)%massmatdiag(elementarray(i)%p+1,neq) + &
                  elementarray(i+1)%massmatdiag(1,neq)
             elementarray(i+1)%massmatdiag(1,neq) = elementarray(i)%massmatdiag(elementarray(i)%p+1,neq)

          enddo

       enddo
       !     assembled matrices

       !boundary contributions
       call boundarycontribution(elementarray(1),elementarray(1)%solnvec,inpobj)
       call boundarycontribution(elementarray(inpobj%nelements),elementarray(inpobj%nelements)%solnvec,inpobj)

  end subroutine compute_explicitrhs
  !========================================================================================================
  subroutine dotimestep_explicit(elementarray,inpobj,multiplier)
    
    type(inputdata), intent(inout)  :: inpobj
    type(element), intent(inout) :: elementarray(:)
    real*8 :: multiplier
    integer :: i,neq,df

		  do i=1,inpobj%nelements

          		do neq=1,globalnumeq

             			do df=1,elementarray(i)%p+1

               				 elementarray(i)%solnvec(df,neq) = elementarray(i)%solnvec_prev(df,neq) +&
                     				multiplier*inpobj%timestep*elementarray(i)%residual(df,neq)/elementarray(i)%massmatdiag(df,neq)
             			enddo
          		enddo
       		enddo
    
  
  end subroutine dotimestep_explicit
  !========================================================================================================
  subroutine timestepping_Implicit(inpobj,elementarray,pointsoln,interpsoln)
    
    type(inputdata), intent(inout)  :: inpobj
    type(element), intent(inout) :: elementarray(:)
    real*8                       :: l2norm(globalnumeq),totaltime,printtime
    integer                      :: iterations,df,i,neq,printit
    integer                      :: numstages
    integer                      :: rkstage
    real*8                       :: multiplier
    real*8,allocatable           :: pointsoln(:,:)
    real*8,allocatable           :: interpsoln(:,:)

    real*8, allocatable           :: X(:)
    real*8, allocatable           :: X0(:)
    real*8, allocatable           :: b(:)
    real*8, allocatable           :: testvec(:)

    integer :: ndf,nelem,p,N,el

    nelem = inpobj%nelements
    p     = elementarray(1)%p

    ndf   = p*nelem + 1
    N     = ndf*globalnumeq

    allocate( X(n))
    allocate(X0(n))
    allocate( b(n))
    allocate(testvec(n))

    testvec = 1.d0


    call printsolution    (elementarray,inpobj%nelements,0,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors, pointsoln)

    call printpolysolution(elementarray,inpobj%nelements,0,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,interpsoln)
    totaltime  = 0.d0
    printtime  = 0.d0
    printit    = 1
    iterations = 0
    
    do while(totaltime < inpobj%finaltime)

       totaltime  = totaltime  + inpobj%timestep	
       printtime  = printtime  + inpobj%timestep       
       iterations = iterations + 1
      
       do i=1,inpobj%nelements
       	elementarray(i)%currenttime = totaltime
       enddo

       call applydirichletbc(elementarray(inpobj%nelements),inpobj)
       
       do i=1,inpobj%nelements

          call computemassmatX(elementarray(i),inpobj,totaltime)
          call computediffmatX(elementarray(i),elementarray(i)%solnvec,inpobj,totaltime)

       enddo
       
       !boundary contributions
       call boundarycontribution(elementarray(1),elementarray(1)%solnvec,inpobj)
       call boundarycontribution(elementarray(inpobj%nelements),elementarray(inpobj%nelements)%solnvec,inpobj)
       
       call findBvec(b,N)
       
       call findX0(X0,elementarray,inpobj,n)
     
       call performgmres(b,X0,X,10,N,600,findAX,noprecond)
       
       call updatelocalsoln(X,elementarray,inpobj,n)
       
       write(*,'(F18.11)',advance='yes')

       if(printtime*inpobj%tscale .ge. inpobj%printtime) then
    	
	call printsolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
			inpobj%xscale,inpobj%scalefactors,pointsoln)

    	call printpolysolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
			inpobj%xscale,inpobj%scalefactors,interpsoln)

	print *
	print *,"saving solution at time:",totaltime
	printtime=0.d0
	printit=printit+1

	endif
       
       print *,"totaltime:",totaltime*inpobj%tscale,"iterations:",iterations
       print *,"======================================="

    enddo
    
    print *,"final time:",totaltime*inpobj%tscale
    call printsolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,pointsoln)
    call printpolysolution(elementarray,inpobj%nelements,printit,inpobj%printpointsperelement,&
		    inpobj%xscale,inpobj%scalefactors,interpsoln)
  

  end subroutine timestepping_Implicit
  !========================================================================================================
  subroutine findX0(X0,elementarray,inpobj,n)

	integer, intent(in)          :: n
	real*8, intent(out)          :: X0(n)
	type(element),intent(inout)  :: elementarray(:)
	type(inputdata),intent(in)   :: inpobj

	integer :: i,nelem,p,ndf,neq,df,counter

	nelem = inpobj%nelements
	p     = elementarray(1)%p
	ndf   = p*nelem+1

	counter=0
	do neq=1,globalnumeq
		
		do i=1,nelem
			
			do df=1,p
				counter=counter+1
				X0(counter)=elementarray(i)%solnvec(df,neq)
			enddo

		enddo
		counter=counter+1
		X0(counter)=elementarray(nelem)%solnvec(p+1,neq)
	enddo
				
				
  end subroutine findX0
  !========================================================================================================
  subroutine updatelocalsoln(X,elementarray,inpobj,n)
	
	integer, intent(in) :: n
	real*8, intent(out) :: X(n)
	type(element),intent(inout) :: elementarray(:)
	type(inputdata),intent(in) :: inpobj
	
	integer :: i,nelem,p,ndf,neq,df,counter

	nelem = inpobj%nelements
	p     = elementarray(1)%p
	ndf   = p*nelem+1

	counter=0
	do neq=1,globalnumeq

		do i=1,nelem

			do df=1,p
				counter=counter+1
				elementarray(i)%solnvec(df,neq)=X(counter)
			enddo
			elementarray(i)%solnvec(p+1,neq)=X(counter+1)
		enddo
		counter=counter+1
	enddo
	

  end subroutine updatelocalsoln
  !========================================================================================================
  subroutine findAX(AX,X,n)
  
	integer, intent(in)   :: n
	real*8, intent(inout) :: X(n)
	real*8, intent(out)   :: AX(n)

	integer :: p
	real*8 :: localsolnvec(elementarray(1)%p+1,globalnumeq)
	integer :: counter,ind,i
	integer :: df,neq,ndf
	real*8  :: delt
      
	p = elementarray(1)%p
	delt = inpobj%timestep

	do i=1,inpobj%nelements

	  call getlocalsolnvecfromglobal(X,localsolnvec,inpobj%nelements,i,p,globalnumeq)
          call computemassmatX(elementarray(i),inpobj,elementarray(i)%currenttime)
          call computediffmatX(elementarray(i),localsolnvec,inpobj,elementarray(i)%currenttime)

        enddo
       
	!boundary contributions
       
       call getlocalsolnvecfromglobal(X,localsolnvec,inpobj%nelements,1,p,globalnumeq)
       call boundarycontribution(elementarray(1),localsolnvec,inpobj)
       call getlocalsolnvecfromglobal(X,localsolnvec,inpobj%nelements,inpobj%nelements,p,globalnumeq)
       call boundarycontribution(elementarray(inpobj%nelements),localsolnvec,inpobj)
       
	!fill from forward
	counter=0
        do neq=1,globalnumeq
	
	   do i=1,inpobj%nelements
		
	   	do df=1,p
	   		counter=counter+1
	   		AX(counter) = elementarray(i)%KX(df,neq)
		enddo

          enddo
	 counter=counter+1
         AX(counter) = elementarray(inpobj%nelements)%KX(p+1,neq)

       enddo

       ndf = (inpobj%nelements*p+1)

       !connection terms
       do neq=1,globalnumeq
       
       		do i=1,inpobj%nelements-1
			
			ind = (neq-1)*ndf+i*p+1
			AX(ind) = AX(ind)   + &
					elementarray(i)%KX(p+1,neq)

		enddo
	enddo

	!add mass matrix contribution
	counter=0
	do neq=1,globalnumeq
		
		do i=1,inpobj%nelements

			do df=1,p
				counter=counter+1
				AX(counter) = AX(counter) - elementarray(i)%massmatdiag(df,neq)*X(counter)/delt
			enddo
		enddo
		counter=counter+1
		AX(counter)=AX(counter) - elementarray(inpobj%nelements)%massmatdiag(p+1,neq)*X(counter)/delt
	enddo

	!connection mass matrix terms
        do neq=1,globalnumeq
       
       		do i=1,inpobj%nelements-1
			
			ind = (neq-1)*ndf + i*p+1
			AX(ind) = AX(ind)  - &
				elementarray(i)%massmatdiag(p+1,neq)*X(ind)/delt

		enddo
	enddo
			
			

  end subroutine findAX
  !========================================================================================================
  subroutine findBvec(b,n)
  	
	integer, intent(in) :: n
	real*8, intent(out)  :: b(n)

	integer :: nelem,p,ind,df,i
	integer :: ndf,neq,counter
	real*8 :: delt

	nelem = inpobj%nelements
	p     = elementarray(1)%p

	ndf = nelem*p + 1
	delt = inpobj%timestep

        !without connection terms
	counter=0
	do neq=1,globalnumeq
		
		do i=1,nelem

			do df=1,p
				counter=counter+1
				b(counter)=elementarray(i)%F(df,neq)
			enddo
		enddo
		counter=counter+1
		b(counter) = elementarray(nelem)%F(p+1,neq)
	enddo

	!connection terms
	do neq=1,globalnumeq
		
		do i=1,nelem-1
			ind=(neq-1)*ndf + i*p + 1
			b(ind) = b(ind) + elementarray(i)%F(p+1,neq)
		enddo

	enddo

	!mass matrix terms
	counter=0
	do neq=1,globalnumeq
		
		do i=1,nelem

			do df=1,p
				counter=counter+1
				b(counter) = b(counter) + elementarray(i)%massmatdiag(df,neq)&
				*elementarray(i)%solnvec(df,neq)/delt
			enddo
		enddo
		counter=counter+1
		b(counter)=b(counter) + elementarray(inpobj%nelements)%massmatdiag(p+1,neq)&
					*elementarray(nelem)%solnvec(p+1,neq)/delt
	enddo

	!connection mass matrix terms
       do neq=1,globalnumeq
       
       		do i=1,nelem-1
			
			ind = (neq-1)*ndf + i*p+1
			b(ind) = b(ind)  + &
				elementarray(i)%massmatdiag(p+1,neq)*&
				elementarray(i)%solnvec(p+1,neq)/delt

		enddo
	enddo

	!b is -F - M/delt u^n
	b(:)=-1.d0*b(:)


  end subroutine findBvec
  !========================================================================================================
  subroutine noprecond(MinvX,X,n)
	
        integer, intent(in) :: n
	real*8, intent(inout) :: MinvX(n)
	real*8, intent(inout) :: X(n)
	
	MinvX = X
  
  end subroutine noprecond
  !========================================================================================================
  subroutine jacobiprecond(MinvX,X,n)
  
	integer, intent(in) :: n
	real*8, intent(inout) :: MinvX(n)
	real*8, intent(inout) :: X(n)
	
	real*8 :: diag(n)
	real*8 :: M,delt
	integer :: counter,neq,i,df
	integer :: ind,p,nelem,ndf

	diag = 0.d0
	delt = inpobj%timestep
	p     = elementarray(1)%p
	nelem = inpobj%nelements
	ndf = nelem*p + 1
	
	!find Diagonal term
	counter=0
	do neq=1,globalnumeq
		
		do i=1,nelem

			do df=1,p
				counter=counter+1
				M = elementarray(i)%massmatdiag(df,neq)
				diag(counter)=diag(counter)+elementarray(i)%diag(df,neq)-M/delt
			enddo
		enddo
		counter=counter+1
		M = elementarray(nelem)%massmatdiag(p+1,neq)
		diag(counter)=diag(counter)+elementarray(nelem)%diag(p+1,neq)-M/delt
	enddo

	!connection mass matrix terms
        do neq=1,globalnumeq
       
       		do i=1,nelem-1
			
			ind       = (neq-1)*ndf + i*p+1
			M         = elementarray(i)%massmatdiag(p+1,neq)
			diag(ind) = diag(ind) + &
				elementarray(i)%diag(p+1,neq) - M/delt

		enddo
	enddo

	MinvX = X/diag

  end subroutine jacobiprecond
  !========================================================================================================
  subroutine getlocalsolnvecfromglobal(globalvec,localvec,nelem,elnum,p,neq)
  
	integer, intent(in) :: p,neq,elnum,nelem
	real*8, intent(in) :: globalvec(:)
	real*8, intent(out) :: localvec(p+1,neq)

	integer :: i

	integer :: ndf,eq

	ndf = nelem*p+1

	do eq=1,neq
		localvec(:,eq)=globalvec(  (eq-1)*ndf + (elnum-1)*p + 1 : (eq-1)*ndf + &
				elnum*p + 1  )
	enddo

  end subroutine getlocalsolnvecfromglobal
  !========================================================================================================

end module timestepper_module
