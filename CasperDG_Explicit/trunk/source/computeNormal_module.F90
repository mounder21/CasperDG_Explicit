Module computeNormal_module
    use my_kinddefs
    use inputs_module
    use Globals_module
    implicit none

    real(dp),pointer :: normals(:,:,:)       
    real(dp),pointer :: NormalLengths(:,:,:)  
contains

Subroutine computeNormals(edgeMap_dPhi,numEdgeGaussPts,tm)
    implicit none
    
    integer(i4),intent(in)  :: tm
    integer(i4),intent(in)  :: numEdgeGaussPts
    real(dp)   ,intent(in)  :: edgeMap_dPhi(tm,numEdgeGaussPts,2,6)        !(tm,numGPs,xi/eta,face*orints)
    
   
    
    integer(i4) :: edge_id,face,numGP,leftTri,n1,n2,n3
    real(dp)    :: normVec1(2),normVec2(2),normVec3(2),vec(2),a_x(3),a_y(3)
    real(dp)    :: J_inv(2,2)
    
    Call allocate_normals(numEdgeGaussPts)
    
    normVec1(1) =  1.0_dp
    normVec1(2) =  1.0_dp
    
    normVec2(1) = -1.0_dp
    normVec2(2) =  0.0_dp
    
    normVec3(1) =  0.0_dp 
    normVec3(2) = -1.0_dp
    
    do edge_id = 1,numEdges 
        leftTri = edgeList(edge_id)%e1          ! left triangle of edge
        n1 = triList(leftTri)%vtxList(1)
        n2 = triList(leftTri)%vtxList(2)
        n3 = triList(leftTri)%vtxList(3)
        if (triList(leftTri)%straight .eqv. .true.) then
            a_x(1) = nodeList(n1)%x
            a_x(2) = nodeList(n2)%x
            a_x(3) = nodeList(n3)%x
            a_y(1) = nodeList(n1)%y
            a_y(2) = nodeList(n2)%y
            a_y(3) = nodeList(n3)%y
        else
            print*,'cant compute curved yet,need a_x and a_y,computeNormal' 
        end if 
        face = edgeList(edge_id)%locEdge1       ! local edge of left triangle 
        do numGP = 1, numEdgeGaussPts
            !n = transpose(J^(-1)) 
            
            J_Inv(1,1) =  dot_product(edgeMap_dPhi(:,numGP,2,face),a_y)      ! dPhidEta * a_y
            J_Inv(2,2) =  dot_product(edgeMap_dPhi(:,numGP,1,face),a_x)      ! dPhidxi  * a_x
            J_Inv(1,2) = -dot_product(edgeMap_dPhi(:,numGP,2,face),a_x)      !-dPhidEta * a_x
            J_Inv(2,1) = -dot_product(edgeMap_dPhi(:,numGP,1,face),a_y)      !-dPhidxi  * a_y
            select case (face)
                case(1)
                    vec(:) = matmul(transpose(J_Inv(:,:)),normVec1(:)) 
                case(2)
                    vec(:) = matmul(transpose(J_Inv(:,:)),normVec2(:))       
                case(3)
                    vec(:) = matmul(transpose(J_Inv(:,:)),normVec3(:))     
            end select
             normals(numGP,:,edge_id) = vec(:)
             NormalLengths(face,numGP,edge_id) = sqrt(vec(1)**2 + vec(2)**2)
        end do   
    end do
end subroutine

Subroutine allocate_normals(numEdgeGaussPts)
    integer(i4),intent(in) :: numEdgeGaussPts
    
    allocate(NormalLengths(4,numEdgeGaussPts,numEdges))
    allocate(normals(numEdgeGaussPts,2,numEdges))
end subroutine allocate_normals

end module computeNormal_module