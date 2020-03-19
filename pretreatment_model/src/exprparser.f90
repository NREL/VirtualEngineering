module exprparser_module

  use globalvars_module
  implicit none

contains
  !===================================================================
  subroutine strtochararray(str,chararray)

         character(LEN=*) :: str
	 character,allocatable :: chararray(:)

	 integer :: strlength,i

	 strlength=stringlen(str)
	 allocate(chararray(strlength))

	 do i=1,strlength
	 	chararray(i)=str(i:i)
	enddo
  
  end subroutine strtochararray
  !===================================================================
  subroutine push(stack,ch,top)

    character :: stack(:)
    character :: ch
    integer :: top

    ! need to check for max stack size
    top=top+1
    stack(top)=ch

  end subroutine push
  !===================================================================
  subroutine realpush(stack,realval,top)

    real*8 :: stack(:)
    real*8 :: realval
    integer :: top

    ! need to check for max stack size
    top=top+1
    stack(top)=realval

  end subroutine realpush
  !===================================================================
  subroutine pop(stack,ch,top)

    character :: stack(:)
    character :: ch
    integer :: top
    logical :: isempty

    isempty=.false.
    ch=stack(top)
    top=top-1

  end subroutine pop
  !===================================================================
  subroutine realpop(stack,realval,top)

    real*8 :: stack(:)
    real*8 :: realval
    integer :: top

    realval=stack(top)
    top=top-1

  end subroutine realpop
  !===================================================================
  subroutine evaluatenumstring(numberstring,ncount,val1)

    integer :: ncount
    character :: numberstring(:)
    real*8 :: val1,digitval,power
    integer :: numlength,i
    integer :: decpointloc

    val1=0.0
    numlength=ncount

    !print *,"number string :",numberstring(1:numlength)

    do i=1,numlength
       if(numberstring(i) .eq. '.') then
          exit
       endif
    enddo

    decpointloc=i

    !number without decimal points
    if(decpointloc .eq. numlength) then
       decpointloc=numlength+1
    endif

    !left hand side
    do i=1,decpointloc-1
       read(numberstring(i),*)digitval
       power = decpointloc-1-i
       val1=val1+digitval*(10**power)
    enddo

    !right hand side
    do i=decpointloc+1,numlength
       read(numberstring(i),*)digitval
       power = decpointloc-i
       val1=val1+digitval*(10**power)
    enddo

  end subroutine evaluatenumstring
  !===================================================================
  subroutine convertinfixtopfix(infixexpr,postfixexpr)

    character :: infixexpr(:)
    character, allocatable :: postfixexpr(:)
    integer :: infixexprsize
    integer :: prnum,prnumstack,i
    integer :: topvar,topops
    logical :: diffnumflag

    !functions
    !integer :: stringlen,priority
    !logical :: isop

    character :: opstack(maxstacksize)
    character :: varstack(maxstacksize)
    character :: parse,popchar

    do i=1,maxstacksize
       opstack(i)='j'
       varstack(i)='j'
    enddo
    !infixexprsize = stringlen(infixexpr)
    infixexprsize = size(infixexpr)

    topvar=1
    topops=1

    do i=1,infixexprsize

       parse=infixexpr(i)
       !print *,"i=",i,"char=",parse
       !print *,"opstack =",opstack(2:topops)
       !print *,"varstack =",varstack(2:topvar)

       if(isop(parse)) then

          if(parse .ne. ')') then

             prnum=priority(parse)
             prnumstack=priority(opstack(topops))

             do while((prnumstack .gt. prnum) .and. (topops .gt. 1) )

                if( opstack(topops) .ne. '(' ) then
                   call pop(opstack,popchar,topops)
                   call push(varstack,popchar,topvar)
                else
                   exit
                endif

             enddo

             call push(opstack,parse,topops)
          else
             popchar='j'	
             do while(popchar .ne. '(')

                call pop(opstack,popchar,topops)

                if(popchar .ne. '(') then
                   call push(varstack,popchar,topvar)
                endif

             enddo
          endif
       else

       	if(is_numeric(parse)) then
		       diffnumflag=.true.
         else
	       diffnumflag=.false.
         endif

          call push(varstack,parse,topvar)
       endif

       !print *,"i=",i,"varstack=",varstack(2:topvar)

    enddo

    do i=2,topops
       call pop(opstack,popchar,topops)
       if(popchar .ne. '(') then
          call push(varstack,popchar,topvar)
       endif
    enddo

    allocate(postfixexpr(topvar-1))

    do i=2,topvar
       postfixexpr(i-1)=varstack(i)
    enddo

  end subroutine convertinfixtopfix
  !===================================================================
  subroutine evaluateexpr(pfixexpr,vars,varvals,finalval)

    character :: pfixexpr(:)
    character :: vars(:)
    real*8 :: varvals(:)
    real*8 :: finalval
    integer :: topstack
    character :: popchar1,popchar2
    real*8 :: val1,val2,res

    real*8,allocatable :: stack(:)
    character :: parse
    character :: numberstring(10)

    integer :: exprsize
    integer :: numvars,i
    integer :: numberstringcount

    exprsize=size(pfixexpr)
    numvars = size(vars)

    allocate(stack(exprsize+1))
    
    do i=1,exprsize+1
    	stack(i)=0.d0
    enddo
    
    topstack=1
    i=1

    do while(i<=exprsize)

       numberstringcount=0
       parse=pfixexpr(i)

       if(isop(parse)) then

          call realpop(stack,val1,topstack)
          val2=0.0

          if(isopunary(parse)) then

             res=domath(val1,val2,parse)
          else
             call realpop(stack,val2,topstack)
             res = domath(val1,val2,parse)
          endif

          call realpush(stack,res,topstack)
	  i=i+1

       else
          if(is_numeric(parse)) then

	     !print *,"i:",i

	     do while((is_numeric(parse) .eqv. .true.) .or. (parse .eq. '.'))
                
	     	numberstringcount=numberstringcount+1
                numberstring(numberstringcount)=parse
		i=i+1
		
	     	if(i > exprsize) then
			exit
		else
                	parse=pfixexpr(i)
		endif
                
             end do
     
	     !print *,"numberstring :",numberstring
       	     call evaluatenumstring(numberstring,numberstringcount,val1)
       
          else
                  if(parse .eq. 'U') then
			  val1=1.0
	          else
		  	val1=varvals(findlocation(vars,numvars,parse))
		  endif

		  i=i+1

          endif
       
       	  call realpush(stack,val1,topstack)
       
  	endif

    enddo

    finalval=stack(topstack)

  end subroutine evaluateexpr
  !===================================================================
  function domath(val1,val2,op) result(res)

    real*8 :: val1,val2
    character :: op
    real*8 :: res
    real*8 :: ten

    ten=10.d0

    select case(op)

    case('+') 
       res=val2+val1
    case('-')
       res=val2-val1
    case('*')
       res=val2*val1
    case('/')
       res=val2/val1
    case('^')
       res=val2**val1
    case('m')
       res=-val1
    case('e')
       res=ten**(val1)
    case('E')
       res=exp(val1)

    end select

  end function domath
  !===================================================================
  function findlocation(charray,length,ch) result(loc)

    character :: charray(:)
    integer :: length
    integer :: loc,i
    character :: ch

    do i=1,length
       if(charray(i) .eq. ch) then
          exit
       endif
    enddo

    loc=i

  end function findlocation
  !===================================================================
  function stringlen(string) result(leng)

    character(LEN=*) :: string
    integer :: i,leng

    i=1
    do while(string(i:i) .ne. ' ')
       i=i+1
    enddo

    leng=i-1

  end function stringlen
  !===================================================================
  function priority(op) result(priornum)

    character :: op
    integer :: priornum

    select case(op)

    case('(')
       priornum=5
    case('E')
       priornum=4
    case('m')
       priornum=4
    case('e')
       priornum=4
    case('^')
       priornum=4
    case('/')
       priornum=3
    case('*')
       priornum=2
    case('+')
       priornum=1
    case('-')
       priornum=1

    case('j')
       priornum=0

    end select

  end function priority
  !===================================================================
  function isop(ch) result(isit)

    character :: ch
    logical :: isit

    isit=.false.

    select case(ch)

    case('E')
       isit=.true.
    case('m')
       isit=.true.
    case('e')
       isit=.true.
    case('^')
       isit=.true.
    case('*')
       isit=.true.
    case('/')
       isit=.true.
    case('+')
       isit=.true.
    case('-')
       isit=.true.
    case('(')
       isit=.true.
    case(')')
       isit=.true.

    end select

  end function isop
  !===================================================================
  function isopunary(ch) result(isit)

    character :: ch
    logical :: isit

    isit=.false.

    select case(ch)

    case('E')
       isit=.true.
    case('m')
       isit=.true.
    case('e')
       isit=.true.
    case('^')
       isit=.false.
    case('*')
       isit=.false.
    case('/')
       isit=.false.
    case('+')
       isit=.false.
    case('-')
       isit=.false.
    case('(')
       isit=.false.
    case(')')
       isit=.false.

    end select

  end function isopunary
  !===================================================================
  FUNCTION is_numeric(ch) result(isit)

    character :: ch
    logical :: isit
    real :: x
    integer :: e

    read(ch,*,IOSTAT=e) x
    isit = e == 0
    
    if(ch .eq. '/') then
	    isit=.false.
    endif

  END FUNCTION is_numeric
  !===================================================================
end module exprparser_module
