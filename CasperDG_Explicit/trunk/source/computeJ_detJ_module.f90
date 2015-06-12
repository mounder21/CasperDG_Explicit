Module computeJ_detJ_module 
    use my_kinddefs
    use inputs_module
    use Globals_module, only: nodeList,triList,numTri,numEulerVars,edgeList,numEdges
    implicit none
    
    real(dp),pointer:: JInv(:,:,:,:),detJ(:,:)
contains

Subroutine computeJ_detJ(dmapPhi,numGaussPts)
    real(dp)    ,intent(in)  :: dmapPhi(:,:,:) !(tm,ngp,xi/eta)
    integer(i4) ,intent(in)  :: numGaussPts
    
    ! Computes the normal vectors for each face and gauss point
    ! Only handles mapBasisDegree = 1 (not possible for p = 0) 

    !dummy
    real(dp)         :: a_x(3),a_y(3)
    integer(i4)      :: n1,n2,n3
    integer(i4)      :: elem_id,numGP
    real(dp)         :: Jtemp(2,2)

    
    Call allocate_J_detJ(numGaussPts)
    
    select case (mapBasisType)
        case('H1')
            do elem_id = 1, numTri 
                if (triList(elem_id)%straight .eqv. .true.) then  
                    n1 = triList(elem_id)%vtxList(1)
                    n2 = triList(elem_id)%vtxList(2)
                    n3 = triList(elem_id)%vtxList(3)

                    a_x(1) = nodeList(n1)%x
                    a_x(2) = nodeList(n2)%x
                    a_x(3) = nodeList(n3)%x
                    a_y(1) = nodeList(n1)%y
                    a_y(2) = nodeList(n2)%y
                    a_y(3) = nodeList(n3)%y
       
                else !curved elements
                    select case (mapBasisDegreeCurved)
                        case default
                            print*, 'No curved elements yet,computeJ_detJ'
                    end select
                end if
                do numGp = 1, numGaussPts
                    !This is correct (3/8/13)
                    Jtemp(1,1) = dot_product(dmapPhi(:,numGp,1),a_x(:))
                    Jtemp(1,2) = dot_product(dmapPhi(:,numGp,1),a_y(:))
                    Jtemp(2,1) = dot_product(dmapPhi(:,numGp,2),a_x(:))
                    Jtemp(2,2) = dot_product(dmapPhi(:,numGp,2),a_y(:))

                    JInv(1,1,numGp,elem_id) =  Jtemp(2,2)
                    JInv(2,1,numGp,elem_id) = -Jtemp(2,1)
                    JInv(1,2,numGp,elem_id) = -Jtemp(1,2)
                    JInv(2,2,numGp,elem_id) =  Jtemp(1,1)
                    detJ(numGP,elem_id) = &
                        Jtemp(1,1)*Jtemp(2,2) - &
                        Jtemp(2,1)*Jtemp(1,2)

                    JInv(:,:,numGp,elem_id) = (1.0_dp/detJ(numGP,elem_id)) * JInv(:,:,numGp,elem_id) 
                end do             
            end do
            
        case('nodal')
            do elem_id = 1, numTri 
                if (triList(elem_id)%straight .eqv. .true.) then
                    n1 = triList(elem_id)%vtxList(1)
                    n2 = triList(elem_id)%vtxList(2)
                    n3 = triList(elem_id)%vtxList(3)

                    a_x(1) = nodeList(n1)%x
                    a_x(2) = nodeList(n2)%x
                    a_x(3) = nodeList(n3)%x
                    a_y(1) = nodeList(n1)%y
                    a_y(2) = nodeList(n2)%y
                    a_y(3) = nodeList(n3)%y

                else !curved elements
                    select case (mapBasisDegreeCurved)
                        case default
                            print*, 'No curved elements yet,computeJ_detJ'
                    end select
                end if
                do numGp = 1, numGaussPts
                    Jtemp(1,1) = dot_product(dmapPhi(:,numGp,1),a_x(:))
                    Jtemp(1,2) = dot_product(dmapPhi(:,numGp,1),a_y(:))
                    Jtemp(2,1) = dot_product(dmapPhi(:,numGp,2),a_x(:))
                    Jtemp(2,2) = dot_product(dmapPhi(:,numGp,2),a_y(:))

                    JInv(1,1,numGp,elem_id) =  Jtemp(2,2)
                    JInv(2,1,numGp,elem_id) = -Jtemp(2,1)
                    JInv(1,2,numGp,elem_id) = -Jtemp(1,2)
                    JInv(2,2,numGp,elem_id) =  Jtemp(1,1)
                    detJ(numGP,elem_id) = &
                        Jtemp(1,1)*Jtemp(2,2) - &
                        Jtemp(2,1)*Jtemp(1,2)

                    JInv(:,:,numGp,elem_id) = (1.0_dp/detJ(numGP,elem_id)) * JInv(:,:,numGp,elem_id) 
                end do
            end do
            
        case('monomial')
            print*,'Cant use monomial'
            stop
        case default
            print*, 'Error- choose basisType: H1, nodal'
            stop
    end select
    
    
end subroutine computeJ_detJ

Subroutine allocate_J_detJ(numGaussPts)
    integer(i4) ,intent(in)  :: numGaussPts
    
    Allocate(JInv(2,2,numGaussPts,numTri))
    Allocate(detJ(numGaussPts,numTri))
end subroutine allocate_J_detJ

end module computeJ_detJ_module 