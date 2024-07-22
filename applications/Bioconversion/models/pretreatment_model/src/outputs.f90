module outputs_module

    use element_module
    use globalvars_module

    implicit none
contains

    !===================================================================================	    
    subroutine printsolution(elarray,nelements,timenum,pointsperelement,xscale,scalefactors,pointsoln)

        type(element), intent(in) :: elarray(:)
        integer, intent(in) :: nelements
        integer, intent(in) :: timenum
        integer, intent(in) :: pointsperelement

        integer :: i,neq,j
        character(LEN=maxstringsize) :: filename
        real*8,allocatable :: x,soln
        real*8,allocatable :: elsoln(:,:)
        real*8 :: xscale
        real*8 :: scalefactors(globalnumeq)
        real*8,allocatable :: pointsoln(:,:)
        logical :: writelastoneflag
        integer :: cnt

        if(allocated(pointsoln) .neqv. .true.) then
            allocate(pointsoln(nelements*elarray(1)%p+1,globalnumeq+1))
        endif

        write(filename,'(A,I4.4,A)')"out_",timenum,".dat"

        open(unit=5,file=filename)

        cnt=1
        do i=1,nelements

        do j=1,elarray(i)%p

        x=getphyscoord(elarray(i)%elementnodes(j),&
            elarray(i)%x(1),elarray(i)%x(2))*xscale

        write(5,'(F15.8)',advance='no') x
        pointsoln(cnt,1)=x

        do neq=1,globalnumeq

        soln=elarray(i)%solnvec(j,neq)*scalefactors(neq)
        write(5,'(F15.8)',advance='no')soln
        pointsoln(cnt,neq+1)=soln
        enddo

        cnt=cnt+1

        write(5,'(A)',advance='yes')

        enddo
        enddo

        i=nelements
        j=elarray(i)%p+1

        x=getphyscoord(elarray(i)%elementnodes(j),&
            elarray(i)%x(1),elarray(i)%x(2))*xscale;
        pointsoln(cnt,1)=x
        write(5,'(F15.8)',advance='no')x

        do neq=1,globalnumeq

        soln=elarray(i)%solnvec(j,neq)*scalefactors(neq)
        write(5,'(F15.8)',advance='no')soln
        pointsoln(cnt,neq+1)=soln

        enddo
        write(5,'(A)',advance='yes')

        close(5)

    end subroutine printsolution
    !===================================================================================	   
    subroutine printpolysolution(elarray,nelements,timenum,pointsperelement,xscale,scalefactors,interpsoln)

        type(element), intent(in) :: elarray(:)
        integer, intent(in) :: nelements
        integer, intent(in) :: timenum
        integer, intent(in) :: pointsperelement

        integer :: i,neq,j
        character(LEN=maxstringsize) :: filename
        real*8,allocatable :: x(:)
        real*8,allocatable :: elsoln(:,:)
        real*8 :: xscale
        real*8 :: scalefactors(globalnumeq)
        character :: tabbing(5)
        real*8,allocatable :: interpsoln(:,:)
        integer :: cnt

        if(allocated(interpsoln) .neqv. .true.) then
            allocate(interpsoln(nelements*(pointsperelement-1)+1,globalnumeq+1))
        endif

        write(filename,'(A,I4.4,A)')"outpoly_",timenum,".dat"

        open(unit=5,file=filename)

        write(5,'(A A A A A A A A)')"x(cm)   ","Steam(molcm3)   ","Liquidfrac   ","Temperature(K)   ",&
            "Xylanfrac   ","Xylog(molcm3)   ","Xylose(molcm3)   ","Furfural(molcm3)"

        cnt=1
        do i=1,nelements

        call getinterpolatedsolution(elarray(i),x,elsoln,pointsperelement)

        do j=1,pointsperelement-1

        write(5,'(F15.8)',advance='no')x(j)*xscale
        interpsoln(cnt,1)=x(j)*xscale

        do neq=1,globalnumeq

        write(5,'(E20.10)',advance='no')elsoln(j,neq)*scalefactors(neq)
        interpsoln(cnt,neq+1)=elsoln(j,neq)*scalefactors(neq)

        enddo

        cnt=cnt+1
        write(5,'(A)',advance='yes')

        enddo
        enddo

        i=nelements
        j=pointsperelement

        call getinterpolatedsolution(elarray(i),x,elsoln,pointsperelement)

        write(5,'(F15.8)',advance='no')x(j)*xscale
        interpsoln(cnt,1)=x(j)*xscale

        do neq=1,globalnumeq

        write(5,'(E20.10)',advance='no')elsoln(j,neq)*scalefactors(neq)
        interpsoln(cnt,neq+1)=elsoln(j,neq)*scalefactors(neq)

        enddo
        write(5,'(A)',advance='yes')

        close(5)

    end subroutine printpolysolution
    !===================================================================================	   
    subroutine calc_l2norm_of_residual(elarray,nelements,l2normres)

        integer,intent(in) :: nelements
        type(element),intent(in) :: elarray(nelements)

        real*8,intent(out) :: l2normres(globalnumeq)
        integer :: neq,i,p,pmax

        do  neq=1,globalnumeq

        l2normres(neq)=0.0

        do i=1,nelements

        if(i < nelements) then
            pmax=elarray(i)%p
        else
            pmax=elarray(i)%p+1
        endif

        do p=1,pmax
        l2normres(neq)=l2normres(neq)+elarray(i)%residual(p,neq)*elarray(i)%residual(p,neq)
        enddo

        enddo

        l2normres(neq)=sqrt(l2normres(neq))
        enddo



    end subroutine calc_l2norm_of_residual
    !===================================================================================	   
end module outputs_module
