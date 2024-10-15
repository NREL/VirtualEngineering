module transport_module

   use globalvars_module
   use inputs_module

   implicit none

  interface setvaluesforvars
	  module procedure setvaluesforvars_default
          module procedure setvaluesforvars_all
  end interface setvaluesforvars

      contains

!=====================================================================
  subroutine setvaluesforvars_default(arrayvar,rcoord,time)

    real*8 :: arrayvar(pretreatvars)
    real*8 :: rcoord
    real*8 :: time

    arrayvar(1)=rcoord
    arrayvar(2)=time

    !for time being (oct 24)
    !later on add the variables from solnvec
    arrayvar(3)=0.0
    arrayvar(4)=0.0
    arrayvar(5)=0.0
    arrayvar(6)=0.0
    arrayvar(7)=0.0
    arrayvar(8)=0.0
    arrayvar(9)=0.0

  end subroutine setvaluesforvars_default
  !=============================================================
  subroutine setvaluesforvars_all(arrayvar,rcoord,time,solnvec,inpobj)

    type(inputdata),intent(in) :: inpobj
    real*8,intent(inout) :: arrayvar(pretreatvars)
    real*8,intent(in) :: rcoord
    real*8,intent(in) :: time
    real*8, intent(in) :: solnvec(globalnumeq)
    integer :: i

    arrayvar(1)=rcoord*inpobj%xscale
    arrayvar(2)=time*inpobj%tscale

    do i=1,globalnumeq
    	arrayvar(2+i)=solnvec(i)*inpobj%scalefactors(i)
    enddo

  end subroutine setvaluesforvars_all
!============================================================================
subroutine calcdcoeffs(T,inpobj)

		real*8, intent(in) :: T
		type(inputdata), intent(inout) :: inpobj

		real*8 :: nmtom,gmtokg,sqmtosqcm
		real*8 :: dpore,rxo,rxy,rf
		real*8 :: Tsteam

		nmtom = 1.0e-9
		gmtokg = 1.0e-3
		sqmtosqcm = 10000.0
		Tsteam = inpobj%T_steam

		rxo   = 0.5*inpobj%d_xylog *nmtom
		rxy   = 0.5*inpobj%d_xylose*nmtom
		rf    = 0.5*inpobj%d_furf  *nmtom
		dpore = inpobj%d_pore      *nmtom
		
		inpobj%D_S0=0.3333*dpore &
		*sqrt(8.0*Runiv*Tsteam/PI/(inpobj%M_w*gmtokg))*sqmtosqcm

		inpobj%D_xylog0  = kB*T/(6.0*PI*eta*rxo)*sqmtosqcm
		inpobj%D_xylose0 = kB*T/(6.0*PI*eta*rxy)*sqmtosqcm
		inpobj%D_F0      = kB*T/(6.0*PI*eta*rf)*sqmtosqcm

end subroutine calcdcoeffs
!============================================================================
subroutine calcdiffandconvcoeffs(solnvec,solnvecder,rcoord,inpobj,diffcoeffs,convcoeffs)

	      real*8,intent(in) :: rcoord
	      type(inputdata),intent(inout) :: inpobj
	      real*8 :: fX,ep,el
	      real*8 :: tau_l,tau_g,T,rhoeff,keff,Ceff,rhoCeff
	      real*8 :: solnvec(globalnumeq)
	      real*8 :: solnvecder(globalnumeq)

	      real*8,intent(out) :: diffcoeffs(globalnumeq)
	      real*8,intent(out) :: convcoeffs(globalnumeq)

	      real*8 :: E_dXylog,E_dXylose,E_dF,E_dS
	      integer :: i

	      fX = solnvec(Xylan_eq)*inpobj%scalefactors(Xylan_eq) 
	      ep = calc_ep(inpobj%fX0,inpobj%ep0,fX)

	      tau_l = inpobj%tau_l
	      tau_g = inpobj%tau_g
	      el    = solnvec(Liquid_eq)*inpobj%scalefactors(Liquid_eq)
	      T     = solnvec(Temperature_eq)*inpobj%scalefactors(Temperature_eq)

	      E_dS  = inpobj%E_dS
	      E_dXylog = inpobj%E_dXylog
	      E_dXylose = inpobj%E_dXylose
	      E_dF = inpobj%E_dF

	      call calcdcoeffs(T,inpobj)

	      diffcoeffs(Steam_eq)=(ep-el)/tau_g * &
			inpobj%D_S0 * exp(-E_dS/Runiv/T)
	
	      diffcoeffs(Liquid_eq) = 0.d0

	      rhoeff=calc_rhoeff(inpobj%rho_s,inpobj%rho_l,&
		      inpobj%rho_g,ep,el)
	      keff = calc_keff(inpobj%k_s,inpobj%k_l, &
			inpobj%k_g,ep,el)
	      Ceff = calc_Ceff(inpobj%C_s,inpobj%C_l, &
			inpobj%C_g,ep,el)

	      !rhoCeff = rhoeff*Ceff
	      rhoCeff = calc_rhoCeff(inpobj%rho_s,inpobj%rho_l,inpobj%rho_g,&
	      		      inpobj%C_s,inpobj%C_l,inpobj%C_g,ep,el)

	      diffcoeffs(Temperature_eq)=keff/rhoCeff

	      diffcoeffs(Xylan_eq) = 0.d0

	      diffcoeffs(Xylog_eq)=(el/tau_l)*inpobj%D_xylog0*&
		      		exp(-E_DXylog/Runiv/T)/el
	      
	      diffcoeffs(Xylose_eq)=(el/tau_l)*inpobj%D_xylose0*&
		      		exp(-E_DXylose/Runiv/T)/el

	      diffcoeffs(Furfural_eq)=(el/tau_l)*inpobj%D_F0*&
		      		exp(-E_DF/Runiv/T)/el


	      if(inpobj%coordsystem .eq. 1) then
	      
		!spherical
		do i=1,globalnumeq
	      		convcoeffs(i)=-2.d0*diffcoeffs(i)/(rcoord*inpobj%xscale)
	      	enddo

      	      else
	     
	       !cylindrical or linear
	      	do i=1,globalnumeq
	      		convcoeffs(i)=0.d0
	      	enddo

      	      endif

	      !print *,"solnvecder of liquid eq:",solnvecder(Liquid_eq)

	      convcoeffs(Xylog_eq)    = diffcoeffs(Xylog_eq)   /el  *  solnvecder(Liquid_eq) 
	      convcoeffs(Xylose_eq)   = diffcoeffs(Xylose_eq)  /el  *  solnvecder(Liquid_eq)
	      convcoeffs(Furfural_eq) = diffcoeffs(Furfural_eq)/el  *  solnvecder(Liquid_eq) 

	      !print *,"diffcoeff at ",rcoord,":",diffcoeffs
	      !print *,"convcoeff at ",rcoord,":",convcoeffs

end subroutine calcdiffandconvcoeffs
!=============================================================================
subroutine calcreactandsourceterms(solnvec,rcoord,time,inpobj,reactcoeffs,sourceterms)
	      
	real*8,intent(in) ::rcoord,time
	type(inputdata),intent(in) :: inpobj
	
	real*8 :: solnvec(globalnumeq)
	real*8,intent(out) :: sourceterms(globalnumeq)
        real*8,intent(out) :: reactcoeffs(globalnumeq)

	real*8 :: arrayvar(pretreatvars)
	real*8 :: ep,el,kcond,kxylog,kxyl1,kxyl2,kfur,fX,cacid
	real*8 :: csteam,cxylog,cxylose,rhoeff,Ceff,rhoCeff
	real*8 :: tempscale
	real*8 :: kevap
	real*8 :: el0
	real*8 :: diameter

	real*8 :: M_xylan,M_xylog,M_furf,M_xylose !molecular weights

	M_xylan  = inpobj%M_X
	M_xylog  = inpobj%M_XO
	M_xylose = inpobj%M_XY
	M_furf   = inpobj%M_F


	diameter = (inpobj%maxr-inpobj%minr)/10; 
	!aspect ratio is typically 20 
	!(Roche et al., biotech and bioeng., 104,2,2009) 

        call setvaluesforvars(arrayvar,rcoord,time,solnvec,inpobj)
	tempscale=inpobj%scalefactors(Temperature_eq)

	kcond   = getcondensationrate(inpobj%k_cond,  tempscale*solnvec(Temperature_eq))
	kxylog  = getarrhreactionrate(inpobj%k_xylog, tempscale*solnvec(Temperature_eq))
	kxyl1   = getarrhreactionrate(inpobj%k_xyl1,  tempscale*solnvec(Temperature_eq))
	kxyl2   = getarrhreactionrate(inpobj%k_xyl2,  tempscale*solnvec(Temperature_eq))
	kfur    = getarrhreactionrate(inpobj%k_fur,   tempscale*solnvec(Temperature_eq))

	el = solnvec(Liquid_eq)*inpobj%scalefactors(Liquid_eq)
	fX = solnvec(Xylan_eq)*inpobj%scalefactors(Xylan_eq)
	ep = calc_ep(inpobj%fX0,inpobj%ep0,fX)
	cacid  = inpobj%c_acid0*inpobj%el0/el
	csteam = solnvec(Steam_eq)*inpobj%scalefactors(Steam_eq)
	cxylog = solnvec(Xylog_eq)*inpobj%scalefactors(Xylog_eq)
	cxylose = solnvec(Xylose_eq)*inpobj%scalefactors(Xylose_eq)
	
	rhoeff=calc_rhoeff(inpobj%rho_s,inpobj%rho_l,&
		      inpobj%rho_g,ep,el)
	Ceff = calc_Ceff(inpobj%C_s,inpobj%C_l, &
			inpobj%C_g,ep,el)

	!rhoCeff = rhoeff*Ceff
	rhoCeff = calc_rhoCeff(inpobj%rho_s,inpobj%rho_l,inpobj%rho_g,&
			      inpobj%C_s,inpobj%C_l,inpobj%C_g,ep,el)
	kevap   = getevaporationrate(inpobj%k_evap, tempscale*solnvec(Temperature_eq), el)
	el0     = inpobj%k_evap(3)
	
	reactcoeffs(Steam_eq)  = -kcond*(ep-el)
	reactcoeffs(Liquid_eq) = -kcond*csteam*inpobj%M_w/inpobj%rho_l - kevap
	reactcoeffs(Temperature_eq) = -4.d0*inpobj%h/diameter
	reactcoeffs(Xylan_eq)  =  -(kxylog+kxyl1)*cacid*ep
	reactcoeffs(Xylog_eq)  =  -kxyl2*cacid
	reactcoeffs(Xylose_eq) = -kfur*cacid
	reactcoeffs(Furfural_eq) = 0.d0


	sourceterms(Steam_eq)   = 0.d0 + kevap*(el-el0)*inpobj%rho_l/inpobj%M_w

	sourceterms(Liquid_eq)  = kcond*csteam*ep*inpobj%M_w/inpobj%rho_l + kevap*el0
	
	sourceterms(Temperature_eq) = inpobj%L_cond*kcond*(ep-el)*csteam*inpobj%M_w/&
			(rhoCeff) - inpobj%L_cond*kevap*(el-el0)*inpobj%rho_l/(rhoCeff)&
			+ 4.d0*inpobj%h/diameter*inpobj%T_steam;

	sourceterms(Xylan_eq) = 0.d0
	
	sourceterms(Xylog_eq) = kxylog*inpobj%rho_s*fX*cacid/inpobj%M_X*ep  * (M_xylan/M_xylog)
	
	sourceterms(Xylose_eq) = kxyl1*inpobj%rho_s*fX*cacid/inpobj%M_X*ep  * (M_xylan/M_xylose) +&
				kxyl2*cxylog*cacid * (M_xylog/M_xylose)
	
	sourceterms(Furfural_eq) = kfur*cxylose*cacid * (M_xylose/M_furf)


	!print *,"sourceterms at ",rcoord,":",sourceterms
	!print *,"react terms at ",rcoord,":",reactcoeffs

end subroutine calcreactandsourceterms
!=====================================================================
function getarrhreactionrate(rateparams,T) result(rateval)
	
	real*8 :: rateval
	real*8,intent(in) :: rateparams(3)
	real*8,intent(in) :: T

	real*8 :: A,E,alpha;

	A     = rateparams(1)
	alpha = rateparams(2)
	E     = rateparams(3)

	rateval=A*(T**(alpha))*exp(-E/Runiv/T)

  end function getarrhreactionrate
!=============================================================================
function getcondensationrate(rateparams,T) result(rateval)

	!the rate is calculated from kinetic theory as (3c/4r)*alpha
	!c is the average thermal speed of steam
	!r is a typical length scale which denotes the ratio of volume to 
	!surface area of region filled by steam
	!alpha is the condensation coefficient
	!The rate constant is 0 for temperature greater than Tsat

	real*8, intent(in) :: rateparams(3)
	real*8, intent(in) :: T

	real*8 :: Tsat,alpha,radius
	real*8 :: rateval,c
	real*8 :: Mw

	real*8 :: tol
	real*8 :: delT

	tol = 0.001

	alpha  = rateparams(1)
	Tsat   = rateparams(2)
	radius = rateparams(3)*0.01 !in metres

	Mw = 18.0e-3 !kg/mol

	rateval=0.0

	delT=tol*Tsat

	if(T .le. Tsat-delT) then
		rateval=alpha
	elseif((T > Tsat-delT) .and. (T .le. (Tsat+delT))) then
		rateval=0.5*alpha*(1-(T-Tsat)/delT)
	else
		rateval=0.d0
	endif

end function getcondensationrate
!=============================================================================
function getevaporationrate(rateparams,T,el) result(rateval)

	real*8, intent(in) :: rateparams(3)
	real*8, intent(in) :: T
	real*8, intent(in) :: el

	real*8 :: rateval
	real*8 :: Tsat,alpha,el0
	real*8 :: tol,delT

	tol = 0.001
	delT = tol*Tsat

	alpha  = rateparams(1)
	Tsat   = rateparams(2)
	el0    = rateparams(3)

	rateval=0.0

	if(el .ge. el0) then
		
		if(T .le. Tsat-delT) then
			rateval=0.d0
		elseif((T > Tsat-delT) .and. (T .le. (Tsat+delT))) then
			rateval=0.5*alpha*(1+(T-Tsat)/delT)
		else
			rateval=alpha
		endif
		
	endif

end function getevaporationrate
!=============================================================================
function calc_ep(fx0,ep0,fX) result(ep)

      real*8,intent(in) :: fx0,ep0,fx
      real*8 :: ep
      real*8 :: a,b,c,D
      real*8 :: r1,r2,dratio,k

      !ep = 1.d0-(1.d0-fX0)*(1.d0-ep0)/(1.d0-fX)
      !ep = fX0*(1-ep0)+ep0-fX

      dratio=2.0

      k = (fX0/(1-fX0) + dratio)/(1-ep0)

      a = k;
      b = -(dratio+k*fX);
      c = dratio*fX - fX;

      D = b*b - 4.d0*a*c

      r1 = 0.5*(-b + sqrt(D))/a
      r2 = 0.5*(-b - sqrt(D))/a

      !print *,"r1,r2:",r1,r2

      if(r1 > 0) then
	      ep=(1-r1)
      else if(r2 > 0) then
	      ep=(1-r2)
      else
	      print *,"error, r1 and r2:",r1,r2
	      ep=ep0
      endif


end function calc_ep
!=============================================================================
function calc_keff(k_s,k_l,k_g,ep,el) result(keff)

      real*8,intent(in):: k_s,k_l,k_g,ep,el
      real*8 :: keff

      keff=(1.d0-ep)*k_s+el*k_l+(ep-el)*k_g

end function calc_keff
!=============================================================================
function calc_rhoeff(rho_s,rho_l,rho_g,ep,el) result(rhoeff)

      real*8, intent(in) :: rho_s,rho_l,rho_g,ep,el
      real*8 :: rhoeff

      rhoeff = (1.d0-ep)*rho_s+el*rho_l+(ep-el)*rho_g

end function calc_rhoeff
!=============================================================================
function calc_Ceff(C_s,C_l,C_g,ep,el) result(Ceff)

      real*8, intent(in) :: C_s,C_l,C_g,ep,el
      real*8 :: Ceff

      Ceff=(1.d0-ep)*C_s+el*C_l+(ep-el)*C_g

end function calc_Ceff
!=============================================================================
function calc_rhoCeff(rho_s,rho_l,rho_g,C_s,C_l,C_g,ep,el) result(rhoCeff)
      
      real*8, intent(in) :: C_s,C_l,C_g,ep,el
      real*8, intent(in) :: rho_s,rho_l,rho_g
      real*8 :: rhoCeff

      rhoCeff=(1.d0-ep)*rho_s*C_s + el*rho_l*C_l + (ep-el)*rho_g*C_g

end function calc_rhoCeff
!=============================================================================

end module transport_module
