module inputs_module
use exprparser_module
use globalvars_module
      implicit none


contains

!=========================================================
subroutine readinpfile(inputs, inputfilename)

	character(LEN=maxstringsize) :: temp
	type(inputdata) :: inputs
	character(*), intent(in) :: inputfilename
	integer :: strlen
	character, allocatable :: reactionrateinfix(:)
	character(LEN=maxstringsize) :: reactionrate
	integer :: numeq

	!open(unit=4,file="pretreat_defs.inp")

	! works OK if the file exists, but results in a runtime error if it
	! doesn't; is there a way to check the filename and catch errors? JJS
	! 1/23/14
	open(unit=4,file=inputfilename)
	
	read(4,*)temp
	read(4,*)temp,inputs%coordsystem
	read(4,*)temp,inputs%minr
	read(4,*)temp,inputs%maxr
	read(4,*)temp,inputs%nelements
	read(4,*)temp,inputs%porder
	read(4,*)temp,inputs%timeorder
	read(4,*)temp,inputs%timestep
	read(4,*)temp,inputs%finaltime
	read(4,*)temp,inputs%printtime
	read(4,*)temp,inputs%printpointsperelement
	read(4,*)temp,inputs%xscale
	read(4,*)temp,inputs%tscale
	read(4,*)temp,inputs%scalefactors(Steam_eq)
	read(4,*)temp,inputs%scalefactors(Liquid_eq)
	read(4,*)temp,inputs%scalefactors(Temperature_eq)
	read(4,*)temp,inputs%scalefactors(Xylan_eq)
	read(4,*)temp,inputs%scalefactors(Xylog_eq)
	read(4,*)temp,inputs%scalefactors(Xylose_eq)
	read(4,*)temp,inputs%scalefactors(Furfural_eq)

	inputs%finaltime=inputs%finaltime/inputs%tscale
	inputs%timestep=inputs%timestep/inputs%tscale

	!read(4,*)temp !blank line
	print *,"temp:",temp
	read(4,*)temp !Heading
	read(4,*)temp,inputs%fx0
	read(4,*)temp,inputs%ep0
	read(4,*)temp,inputs%c_steam_bulk
	read(4,*)temp,inputs%c_acid0
	read(4,*)temp,inputs%el0
	read(4,*)temp,inputs%T0
	read(4,*)temp,inputs%T_steam
	print *,"T_steam:",inputs%T_steam


	!read(4,*)temp !blank line
	read(4,*)temp !Heading
	read(4,*)temp,inputs%E_dS
	read(4,*)temp,inputs%E_dXylog
	read(4,*)temp,inputs%E_dXylose
	read(4,*)temp,inputs%E_dF
	read(4,*)temp,inputs%L_cond
	read(4,*)temp,inputs%k_s
	read(4,*)temp,inputs%k_l
	read(4,*)temp,inputs%k_g
	read(4,*)temp,inputs%C_s
	read(4,*)temp,inputs%C_l
	read(4,*)temp,inputs%C_g
	read(4,*)temp,inputs%h

	!read(4,*)temp
	read(4,*)temp
	read(4,*)temp,inputs%rho_s
	read(4,*)temp,inputs%rho_l
	read(4,*)temp,inputs%rho_g
	read(4,*)temp,inputs%M_X
	read(4,*)temp,inputs%M_XO
	read(4,*)temp,inputs%M_XY
	read(4,*)temp,inputs%M_F
	read(4,*)temp,inputs%M_w
	read(4,*)temp,inputs%tau_g
	read(4,*)temp,inputs%tau_l
	read(4,*)temp,inputs%d_pore
	read(4,*)temp,inputs%d_xylog
	read(4,*)temp,inputs%d_xylose
	read(4,*)temp,inputs%d_furf

 	!read(4,*)temp
	read(4,*)temp
	read(4,*)temp,inputs%k_cond(1),inputs%k_cond(2),inputs%k_cond(3)
	read(4,*)temp,inputs%k_evap(1),inputs%k_evap(2),inputs%k_evap(3)
	read(4,*)temp,inputs%k_xylog(1),inputs%k_xylog(2),inputs%k_xylog(3)
	read(4,*)temp,inputs%k_xyl1(1),inputs%k_xyl1(2),inputs%k_xyl1(3)
	read(4,*)temp,inputs%k_xyl2(1),inputs%k_xyl2(2),inputs%k_xyl2(3)
	read(4,*)temp,inputs%k_fur(1),inputs%k_fur(2),inputs%k_fur(3)

	close(4)
end subroutine readinpfile
!=========================================================
subroutine printchararray(chararray)

		character, intent(in) :: chararray(:)
		integer :: length

		length=size(chararray)

		print *,chararray(1:length)

end subroutine printchararray
!=========================================================
end module inputs_module
