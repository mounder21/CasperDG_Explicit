Module initializeSolution_module
     use my_kinddefs
     use inputs_module, only: totalModes
     use initializeBasis_module, only: numGaussPts
     use integrate_module

     implicit none 
     real(dp),allocatable :: solCoeffs(:,:,:),spRes(:,:,:)
     real(dp),allocatable :: oldSolCoeffs(:,:,:)
     
 contains
     
 Subroutine initializeSolution(PhiW,numGaussPts,invMassMats,initPrev,readNumModes)    
    use getIC_module
    use projection_module
    implicit none
    
    integer(i4),intent(in) :: numGaussPts,readNumModes
    real(dp),   intent(in) :: PhiW(:,:)
    real(dp),   intent(in) :: invMassMats(:,:,:)
    logical ,   intent(in) :: initPrev
    !dummy
    integer(i4)     :: elem_id,j,k,n1,n2,n3
    real(dp)        :: ax(3),ay(3),x(numGaussPts),y(numGaussPts)        !change ax ay sizes later
    real(dp)        :: vec(totalModes),f(numGaussPts)
    character(80)   :: IC_type
    
    Call allocate_solutionCoeffs
    solCoeffs(:,:,:) = 0.0_dp
    
    if(initPrev)then
        Call read_solution(readNumModes)  
    else
        !------------------------------------- Initial Data ------------------------------------------!
        IC_type = 'none'

!        do j = 1, numEulerVars 
!            do k = 1, numGaussPts
!                f(k) = getIConstant(j)
!            end do
!
!            do elem_id = 1, numTri
!                call integrateCell(elem_id,PhiW,f,vec)    !vec(tm)
!                solCoeffs(j,:,elem_id) = matmul(invMassMats(:,:,elem_id),vec)            
!            end do
!        end do
!        
        do elem_id = 1, numTri
            n1 = triList(elem_id)%vtxList(1)
            n2 = triList(elem_id)%vtxList(2)
            n3 = triList(elem_id)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y
            
            Call projectMapCell(ax,ay,x,y)
            do j = 1, numEulerVars 
                do k = 1, numGaussPts
                    f(k) = getIConstant(j)
                    !f(k) = getIC_MMS(j,x(k),y(k))
                    !f(k) = getIC_MMS_sin(j,x(k),y(k))
                    !f(k) = getIC(x(k),y(k))
                    !f(k) = smooth_IC_Airfoil(j,x(k),y(k))
                    !f(k) = isentropicVortex(j,x(k),y(k))
                end do

                call integrateCell(elem_id,PhiW,f,vec)    !vec(tm)
                solCoeffs(j,:,elem_id) = matmul(invMassMats(:,:,elem_id),vec)            
            end do
        end do
    end if
 end subroutine  initializeSolution
 
 subroutine read_solution(readNumModes)   
    use inputs_module,                  only: solution_name,totalModes
    implicit none
    integer(i4),intent(in)  :: readNumModes
    integer(i4)             :: funit,elem_id,numCells,tm,totModes
    character(LEN=20)       :: comments

     funit=2
     OPEN(UNIT=funit,FILE=solution_name,STATUS='old')
     
     read(funit,*) comments
     read(funit,*) numCells,totModes
     if((numCells .ne. numTri).or.(totModes .ne. readNumModes))then
         print*,'The input solution file does not contain same number of elements or total modes: stopping'
         stop
     end if
     solCoeffs(:,:,:) = 0.0_dp
     ! read the solution coefficients
     do elem_id = 1, numCells
        do tm = 1,readNumModes
            read(funit,*) solCoeffs(1,tm,elem_id),solCoeffs(2,tm,elem_id),&
                          solCoeffs(3,tm,elem_id),solCoeffs(4,tm,elem_id)
        end do
        do tm = readNumModes + 1,totModes
            solCoeffs(:,tm,elem_id) = 0.0_dp
        end do
    end do
    close(funit)
    print*,'The initial solution was commenced from a restart'
end subroutine read_solution

 Subroutine allocate_solutionCoeffs
     Allocate(solCoeffs(numEulerVars,totalModes,numTri))
     Allocate(spRes(numEulerVars,totalModes,numTri))
     Allocate(oldSolCoeffs(numEulerVars,totalModes,numTri))
 end subroutine allocate_solutionCoeffs
 
end module initializeSolution_module
