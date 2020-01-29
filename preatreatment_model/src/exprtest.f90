program pretreatmentmodeling

	use polynomial_module
	use exprparser_module 
	
	implicit none

	type(polynomial) :: p
	real*8 :: arr(3)
	character,allocatable :: pfix(:)
	character(LEN=maxstacksize) :: expr
	character,allocatable :: infixexpr(:)
	integer :: i,lengexpr
	character :: vars(4)
	real*8 :: varvals(4)
	real*8 :: fnval
	character :: testch

	!functions
	!real*8 :: polyfindval

	!expr='(A+B*(C+D)/(A-B))/(D+5)^2'
	!expr='(A+B*C*D+A+B)/(A+B+C+D+(A+B)*(C+D))'
	!expr='32*U*E(m1.6*U*e(m5)/(3.2*U*e(m5)))'

	!expr='(A-B)/(3.0*U*e(m4))*5.0*U*E(m4.0*U*e(2)/(8314.0*U*300.0))'
	!expr='(A-B)/(3.0*U*e(m4))*5.0'
	!expr='3.0*U*e(m5.0)*A^20.0+4.0*U*E(3.0*U*E(1))'
	!expr='3.0'
	expr='E(m1.0/A)'
	
	lengexpr=stringlen(expr)
	allocate(infixexpr(lengexpr))

	do i=1,lengexpr
		infixexpr(i)=expr(i:i)
	enddo
		

	call convertinfixtopfix(infixexpr,pfix)

	print *,pfix(1:size(pfix))

	vars(1)='A'
	vars(2)='B'
	vars(3)='C'
	vars(4)='D'

	varvals(1)=1.d0
	varvals(2)=1.d0
	varvals(3)=1.d0
	varvals(4)=1.d0

	call evaluateexpr(pfix,vars,varvals,fnval)

	print *,"finalval=",fnval
	
end program pretreatmentmodeling
