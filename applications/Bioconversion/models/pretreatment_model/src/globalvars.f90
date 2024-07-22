module globalvars_module

implicit none
integer, parameter :: maxpostfixsize=1000
integer, parameter :: maxstacksize=1000
integer, parameter :: maxstringsize=1000
integer,parameter :: pretreatvars=9
!for time being (oct 24)
integer, parameter:: globalnumeq=7
character,parameter :: globalvars(pretreatvars)=(/'R',&
	't','S','L','T',&
	'X','Y','Z','F'/)

integer,parameter :: Steam_eq=1
integer,parameter :: Liquid_eq=2
integer,parameter :: Temperature_eq=3
integer,parameter :: Xylan_eq = 4
integer,parameter :: Xylog_eq = 5
integer,parameter :: Xylose_eq = 6
integer,parameter :: Furfural_eq = 7

real*8,parameter :: Runiv=8.314 !J/K/mol
real*8,parameter :: kB=1.381e-23 !J/K/#
real*8,parameter :: PI=3.14159265 !J/K/mol
! Pa s (water viscosity at 100 C, flattens out
!after 100 
!see http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html)
real*8,parameter :: eta=2.8e-4 ! Pa s
real*8,parameter :: rkcoeffs(5)=(/0.0533,0.1263,0.2375,0.4414,1.0/)


type polynomial
	integer :: n
	real*8,allocatable :: coeffs(:)
end type polynomial

type pfixexpr
	character,allocatable :: expr(:)
end type pfixexpr

type inputdata

	real*8 :: minr,maxr
	integer :: nelements
	integer :: porder,timeorder
	real*8 :: timestep,finaltime,printtime
	real*8 :: xscale,tscale
	integer :: printpointsperelement
	integer :: coordsystem

	!Initial and boundary conditions
	real*8 :: fx0,ep0,c_steam_bulk,&
		c_acid0,el0,T0,T_steam

	real*8 :: d_pore, d_xylose, d_xylog, d_furf

	!Thermodynamics and energetics 
	real*8 :: E_dS,E_dXylog,E_dXylose,E_dF,L_cond,&
		k_s,k_l,k_g,C_s,C_l,C_g,h

	!densities,tortuosities,diff coeffs,mol. wts
	real*8 :: rho_s,M_X,M_XO,M_XY,M_F,tau_g,D_S0,M_w,rho_l,&
		tau_l,D_xylog0,D_Xylose0,D_F0,rho_g

	!Rate coefficients
	real*8 :: k_cond(3)
	real*8 :: k_xylog(3) 
	real*8 :: k_xyl1(3)
	real*8 :: k_xyl2(3) 
	real*8 :: k_fur(3) 
	real*8 :: k_evap(3) 

	real*8 :: scalefactors(globalnumeq)

end type inputdata
type element

	real*8 :: x(2)
	real*8 :: dx
	integer :: p
	type(polynomial),allocatable:: basisfunc(:)
	type(polynomial),allocatable:: basisfuncder(:)
	real*8,allocatable :: elementnodes(:)
	real*8,allocatable :: globattoweights(:)
	real*8,allocatable :: massmatdiag(:,:)
	real*8,allocatable :: residual(:,:)
	real*8,allocatable :: solnvec(:,:)
	real*8,allocatable :: solnvec_prev(:,:)
	
	real*8,allocatable :: gaussquadx(:)
	real*8,allocatable :: gaussquadw(:)
	integer :: boundaryside

	real*8, allocatable :: KX(:,:)
	real*8, allocatable ::  F(:,:)
	real*8, allocatable :: diag(:,:)
	real*8 :: currenttime

end type element

!global variable
type(element),allocatable :: elementarray(:)
type(inputdata) :: inpobj

end module globalvars_module
