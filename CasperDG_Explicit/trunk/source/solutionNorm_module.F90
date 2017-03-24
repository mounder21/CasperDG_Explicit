Module solutionNorm_module
    use my_kinddefs
    implicit none
    
contains

Subroutine Lp_solution_norm(solCoeffs,numGaussPts,p,norm)
    use getIC_module
    use projection_module
    use integrate_module
    use norm_module
    use Globals_module, only: numTri,triList,numFields
    use initializeBasis_module, only: Phi
    use inputs_module, only: totalModes
    implicit none
    
    real(dp),intent(in) :: solCoeffs(:,:,:),p
    integer(i4),intent(in) :: numGaussPts
    real(dp),intent(out)   :: norm(numFields)
    
    real(dp) :: ax(3),ay(3),x(numGaussPts),y(numGaussPts),res(numGaussPts),f,tempNorm(numFields,numTri)
    real(dp) :: q(numFields,numGaussPts),vec,oneOverp
    integer(i4) :: elem_id,j,k,n1,n2,n3,index1
    
    oneOverp = 1.0_dp / p
    
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
        Call projectCell(solCoeffs(:,:,elem_id),Phi,q) !q(numFields,numGaussPts)
        do j = 1, numFields 
            do k = 1, numGaussPts
                !f(k) = getIC_MMS(j,x(k),y(k))
                f = getIC_MMS_sin(j,x(k),y(k))
                res(k) = (q(j,k) - f)**p
            end do
            
            Call integrate(elem_id,res,vec)

            !--------------- Calculate norm ---------------!
            tempNorm(j,elem_id) = vec   
        end do
    end do

    do j = 1, numFields
        norm(j) = 0.0_dp
        do elem_id = 1, numTri
            norm(j) = norm(j) + tempNorm(j,elem_id)
        end do
        norm(j) = norm(j)**oneOverp
        print*,'Error of solution, euler var:', norm(j),j
    end do

end subroutine

end module
