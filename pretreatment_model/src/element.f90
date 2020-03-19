module element_module

  use polynomial_module
  use inputs_module
  use quadrature_module
  use transport_module

  implicit none



contains
  !=============================================================
  subroutine setelementattribs(el,endp,inpobj,belementside)

    type(element)    :: el
    type(inputdata)  :: inpobj
    integer,optional :: belementside
    integer          :: i,neq,p
    real*8           :: x
    real*8           :: arrayvar(pretreatvars)
    real*8           :: initval
    real*8           :: endp(2)

    el%boundaryside=0 
    !whether the element is at the boundary
    !1 indicates left and 2 indicates right

    if(present(belementside)) then
       el%boundaryside = belementside
    endif

    !scaling x dimensions
    el%x(1) = endp(1)/inpobj%xscale
    el%x(2) = endp(2)/inpobj%xscale

    !scaled size
    el%dx = el%x(2)-el%x(1)

    !polynomial order
    p = inpobj%porder
    el%p = p

    allocate(   el%basisfunc(p+1))
    allocate(el%basisfuncder(p+1))
    allocate(el%elementnodes(p+1))

    allocate( el%massmatdiag((p+1),globalnumeq))
    allocate(    el%residual((p+1),globalnumeq))
    allocate(     el%solnvec((p+1),globalnumeq))
    allocate(el%solnvec_prev((p+1),globalnumeq))
    
    allocate(el%KX( (p+1), globalnumeq) )
    allocate(el%F((p+1), globalnumeq) )
    allocate(el%diag((p+1), globalnumeq) )

    !call getuniformlagrangepoints(p,el%elementnodes)
    !call globattopointsandweights(el%elementnodes,el%globattoweights,p+1)
    call gll_gen(el%elementnodes,el%globattoweights,p)

    !get gauss lagrange quadrature points for 1st element to avoid r=0
    call gl_gen(el%gaussquadx,el%gaussquadw,p+1)

    !create lagrange interpolation polynomials at nodes
    call createlagrangepoly(el%elementnodes,p+1,el%basisfunc)

    !store derivatives before hand
    do i=1,p+1
          call polyderivative(el%basisfunc(i),el%basisfuncder(i))
    enddo
   
     !Initial conditions   
     do i=1,p+1
	  el%solnvec(i,Steam_eq)       = 0.d0
	  el%solnvec(i,Liquid_eq)      = inpobj%el0/inpobj%scalefactors(Liquid_eq)
	  el%solnvec(i,Temperature_eq) = inpobj%T0/inpobj%scalefactors(Temperature_eq)
	  el%solnvec(i,Xylan_eq)       = inpobj%fX0/inpobj%scalefactors(Xylan_eq)&
			  		 *(1-inpobj%ep0)
	  el%solnvec(i,Xylog_eq)       = 0.d0
	  el%solnvec(i,Xylose_eq)      = 0.d0
	  el%solnvec(i,Furfural_eq)    = 0.d0
     enddo

  end subroutine setelementattribs
  !=============================================================
  subroutine saveprvssoln(el)

		  type(element) :: el
		  integer :: i,neq

		  do neq=1,globalnumeq
		  	do i=1,el%p+1
			    el%solnvec_prev(i,neq) = el%solnvec(i,neq)
			enddo
		  enddo

  end subroutine saveprvssoln
  !=============================================================
  subroutine computediffmatX(el,X,inpobj,time) 

    type(inputdata), intent(inout)    :: inpobj
    type(element), intent(inout)   :: el
    real*8, intent(in)             :: time
    real*8, intent(in)             :: X(:,:)

    real*8             :: dcoeff,rcoeff
    real*8             :: sourcemult,convcoeff
    type(polynomial)   :: testp,trialp
    type(polynomial)   :: testpder,trialpder
    real*8,allocatable :: quadw(:),quadx(:)
    real*8             :: integvaldiff,integvalconv
    real*8             :: integvalreact,integvalsource
    real*8             :: dxdzeta,dzetadx
    real*8             :: signdiff,signreact,signsource,signconv
    real*8             :: xscale,tscale
    real*8,allocatable :: convcoeffs(:,:),diffcoeffs(:,:)
    real*8,allocatable :: reactcoeffs(:,:),sourceterms(:,:)
    real*8             :: diffscaling,convscaling,reactscaling
    real*8             :: solnvec(globalnumeq)
    real*8             :: solnvecder(globalnumeq)
    integer            :: i,j,neq,quadnp,q

    dxdzeta    =  0.5*el%dx
    dzetadx    =  2.d0/el%dx
    
    !diffusion term shows up on RHS and -ve from weakening
    signdiff   =  -1.d0 
    signreact  =   1.d0
    signsource =   1.d0
    signconv   =   1.d0 !shows up on LHS
    xscale     =  inpobj%xscale
    tscale     =  inpobj%tscale
    
    diffscaling  = tscale/xscale/xscale
    convscaling  = tscale/xscale
    reactscaling = tscale

    quadnp=el%p+1
    allocate( quadx(quadnp) )
    allocate( quadw(quadnp) )
    
    if(el%boundaryside .eq. 1) then
	    quadx = el%gaussquadx
	    quadw = el%gaussquadw
    else
	quadx = el%elementnodes
	quadw = el%globattoweights
    endif

    !print *,size(quadw),size(quadx) 

    allocate(  convcoeffs(quadnp,globalnumeq) )
    allocate(  diffcoeffs(quadnp,globalnumeq) )
    allocate( reactcoeffs(quadnp,globalnumeq) )
    allocate( sourceterms(quadnp,globalnumeq) )

    do q=1,quadnp
    	
        call getsolnvecatlocation(el,el%solnvec,getphyscoord(quadx(q),el%x(1),el%x(2)),solnvec)
        call getsolnvecderatlocation(el,el%solnvec,getphyscoord(quadx(q),el%x(1),el%x(2)),solnvecder)
    	
	call calcdiffandconvcoeffs(solnvec,solnvecder,getphyscoord(quadx(q),el%x(1),el%x(2)),&
				inpobj,diffcoeffs(q,:),convcoeffs(q,:))
    	
	call calcreactandsourceterms(solnvec,getphyscoord(quadx(q),el%x(1),el%x(2)),&
			     time,inpobj,reactcoeffs(q,:),sourceterms(q,:))
   enddo

    !equation loop
    do neq=1,globalnumeq

    	!print *,"dcoeff of neq:",neq,"is ",diffcoeffs(1,neq)
        !test function loop
       do i=1,el%p+1

          el%residual(i,neq) = 0.d0
	  el%KX(i,neq)       = 0.d0
	  el%F(i,neq)        = 0.d0
    	  el%diag(i,neq)     = 0.d0
          testp    = el%basisfunc(i)
          testpder = el%basisfuncder(i)

          !trial function loop
          do j=1,el%p+1

             trialp    = el%basisfunc(j)
             trialpder = el%basisfuncder(j)

             integvaldiff  =  0.d0
             integvalreact =  0.d0
	     integvalconv  =  0.d0
             
           do q=1,quadnp

	   	dcoeff    =  diffcoeffs(q,neq)
		rcoeff    =  reactcoeffs(q,neq)
		convcoeff =  convcoeffs(q,neq)

		integvaldiff = integvaldiff + quadw(q) &
                        *dcoeff&
                        *polyfindval(testpder,quadx(q))&
                        *polyfindval(trialpder,quadx(q))
               
                integvalreact=integvalreact + quadw(q) &
                           *rcoeff&
                           *polyfindval(testp,quadx(q))&
                           *polyfindval(trialp,quadx(q))

		integvalconv=integvalconv + quadw(q)  &
			  *convcoeff&
			  *polyfindval(testpder,quadx(q))&
		        *polyfindval(trialp,quadx(q))

             enddo
             
             integvaldiff  =  diffscaling*integvaldiff*(dzetadx)*(dzetadx)*(dxdzeta)
             integvalreact =  reactscaling*integvalreact*(dxdzeta)
	     integvalconv  =  convscaling*integvalconv*dzetadx*dxdzeta

	     if(i .eq. j) then
		     el%diag(i,neq) = el%diag(i,neq) + signdiff  * integvaldiff
		     el%diag(i,neq) = el%diag(i,neq) + signreact * integvalreact
		     el%diag(i,neq) = el%diag(i,neq) + signconv  * integvalconv
	     endif

	     el%KX(i,neq) = el%KX(i,neq) + signdiff  *  integvaldiff  *   X(j,neq)
             el%KX(i,neq) = el%KX(i,neq) + signreact *  integvalreact *   X(j,neq)
             el%KX(i,neq) = el%KX(i,neq) + signconv  *  integvalconv  *   X(j,neq)
             
	     el%residual(i,neq) = el%KX(i,neq)

          enddo
          
          integvalsource=0.d0   
          do q=1,quadnp
	     
         	sourcemult     =  sourceterms(q,neq)      

                integvalsource =  integvalsource + quadw(q) &
                                  *sourcemult&
                                  *polyfindval(testp,quadx(q))
        enddo
                
                el%F(i,neq)    =     (tscale/inpobj%scalefactors(neq))*&
				     signsource*integvalsource*dxdzeta

		el%residual(i,neq) = el%residual(i,neq) + el%F(i,neq)
       enddo
 

    enddo


  end subroutine computediffmatX
  !=============================================================
  subroutine computemassmatX(el,inpobj,time) 

    !Making use of globatto weights

    type(inputdata) :: inpobj
    type(element) :: el
    integer :: i,j,neq,quadnp,q
    type(polynomial) :: testp,trialp
    real*8 :: arrayvar(pretreatvars)
    real*8 :: time
    real*8 :: integval
    real*8 :: dxdzeta,dzetadx
    real*8 :: xscale,tscale
    real*8 :: delt

    dxdzeta=0.5*el%dx
    dzetadx=2.d0/el%dx

    xscale = inpobj%xscale
    tscale = inpobj%tscale
    delt   = inpobj%timestep

    !call globattopointsandweights(quadw,quadx,quadnp) 

    do neq=1,globalnumeq

       do i=1,el%p+1
          el%massmatdiag(i,neq) = el%globattoweights(i)*(dxdzeta)
       enddo

    enddo

  end subroutine computemassmatX
  !=============================================================
  subroutine applydirichletbc(el,inpobj)
    
    type(element) :: el
    type(inputdata) :: inpobj
    
    if(el%boundaryside .eq. 2) then
	      el%solnvec(el%p+1,Steam_eq) = inpobj%c_steam_bulk/inpobj%scalefactors(Steam_eq)
    endif


  end subroutine applydirichletbc
  !=============================================================
  subroutine boundarycontribution(el,solnX,inpobj)

    type(element) :: el
    type(inputdata) :: inpobj
    real*8,intent(in) :: solnX(:,:)

    integer :: neq
    real*8 :: kval
    real*8 :: signvar
    real*8 :: xscale
    real*8 :: tscale

    if(el%boundaryside .eq. 2) then

             el%residual(el%p+1,Steam_eq)  = 0.d0

	     el%KX(el%p+1,Steam_eq)        = 0.d0
	     el%F (el%p+1,Steam_eq)        = 0.d0
	  
	     el%diag(el%p+1,Steam_eq)      = 0.d0
             
	     el%solnvec(el%p+1,Steam_eq)  = inpobj%c_steam_bulk/inpobj%scalefactors(Steam_eq)

             signvar=1.0

             el%residual(el%p+1,Temperature_eq) = el%residual(el%p+1,Temperature_eq) + signvar*inpobj%h*inpobj%T_steam*&
			     					(inpobj%tscale/inpobj%scalefactors(Temperature_eq)/inpobj%xscale)
             
	     el%F(el%p+1,Temperature_eq)        = el%F(el%p+1,Temperature_eq) + signvar*inpobj%h*inpobj%T_steam*&
			     					(inpobj%tscale/inpobj%scalefactors(Temperature_eq)/inpobj%xscale)
             
	     el%residual(el%p+1,Temperature_eq) = el%residual(el%p+1,Temperature_eq) - &
	     					  signvar*el%solnvec(el%p+1,Temperature_eq)*inpobj%h*(inpobj%tscale/inpobj%xscale)
	     
	     el%KX(el%p+1,Temperature_eq)       = el%KX(el%p+1,Temperature_eq) - &
	     					  signvar*solnX(el%p+1,Temperature_eq)*inpobj%h*(inpobj%tscale/inpobj%xscale)
	     
	     el%diag(el%p+1,Temperature_eq)     = signvar*inpobj%h*(inpobj%tscale/inpobj%xscale)

    endif

  end subroutine boundarycontribution
  !=============================================================
  subroutine getinterpolatedsolution(el,x,values,npoints)

		  type(element), intent(in) :: el
		  integer, intent(in) :: npoints
		  real*8,allocatable,intent(out) :: x(:)
		  real*8,allocatable,intent(out) :: values(:,:)

		  real*8 :: delzeta
		  real*8 :: minzeta
		  real*8 :: maxzeta
		  real*8 :: zeta
		  integer :: i,p,neq

		  minzeta=-1.d0
		  maxzeta=1.d0

		  allocate(x(npoints))
		  allocate(values(npoints,globalnumeq))

		  delzeta=(maxzeta-minzeta)/(npoints-1.d0)

		  do i=1,npoints
			
		  	zeta=minzeta+(i-1)*delzeta
			x(i)=getphyscoord(zeta,el%x(1),el%x(2))

			!print *,"zeta=",zeta,el%x(1),el%x(2)
			
			do neq=1,globalnumeq

			values(i,neq)=0.0

			do p=1,el%p+1
				values(i,neq)=values(i,neq)+el%solnvec(p,neq)*polyfindval(el%basisfunc(p),zeta)
			enddo

			enddo
		enddo

  end subroutine getinterpolatedsolution
  !=============================================================
  subroutine getsolnvecatlocation(el,nodalsoln,x,values)

	type(element), intent(in)  :: el
	real*8, intent(in)         :: x
	real*8, intent(out)        :: values(globalnumeq)
	real*8, intent(in)         :: nodalsoln(el%p+1,globalnumeq)
	
	integer :: p,neq
	real*8 :: zeta
	real*8 :: val
	zeta=getmastercoord(x,el%x(1),el%x(2))

	do neq=1,globalnumeq
		
		values(neq)=0.d0

		do p=1,el%p+1

			values(neq) = values(neq) + &
			nodalsoln(p,neq)*polyfindval(el%basisfunc(p),zeta)

		enddo
	enddo
  
  end subroutine getsolnvecatlocation
  !=============================================================
  subroutine getsolnvecderatlocation(el,nodalsoln,x,values)

	type(element), intent(in)  :: el
	real*8, intent(in)         :: x
	real*8, intent(out)        :: values(globalnumeq)
	real*8, intent(in)         :: nodalsoln(el%p+1,globalnumeq)
	
	integer :: p,neq
	real*8 :: zeta
	real*8 :: val
	real*8 :: dzetadx

	zeta=getmastercoord(x,el%x(1),el%x(2))
	dzetadx=2.d0/el%dx

	do neq=1,globalnumeq
		
		values(neq)=0.d0

		do p=1,el%p+1

			values(neq) = values(neq) + &
			nodalsoln(p,neq)*polyfindderval(el%basisfunc(p),zeta)*dzetadx

		enddo
	enddo
  
  end subroutine getsolnvecderatlocation
  !=============================================================
  function getmastercoord(physx,x1,x2) result(zeta)

    real*8, intent(in) :: x1,x2,physx
    real*8 :: zeta

    zeta=(2.d0*physx-(x1+x2))/(x2-x1)

  end function getmastercoord
  !=============================================================
  function getphyscoord(zeta,x1,x2) result(physx)

    real*8, intent(in) :: zeta,x1,x2
    real*8 :: physx

    physx=0.5*(1-zeta)*x1+0.5*(1+zeta)*x2

  end function getphyscoord
  !=============================================================
end module element_module
