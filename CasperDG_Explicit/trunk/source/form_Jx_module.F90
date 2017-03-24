Module form_Jx_module
    use my_kinddefs
    use Globals_module, only: numTri,numFields,triList,interiorEdges,edgeList,numInterior
    use inputs_module,  only: totalModes
    implicit none
    
contains

Subroutine offDiagJX(edge_id,jacobOffDiagRL,jacobOffDiagLR,x,Jx)
    integer(i4),intent(in)  :: edge_id
    real(dp), intent(in)    :: jacobOffDiagRL(numFields*totalModes,numFields*totalModes)
    real(dp), intent(in)    :: jacobOffDiagLR(numFields*totalModes,numFields*totalModes)
    real(dp),intent(in)     :: x(numFields*totalModes,numTri)
    real(dp),intent(inout)    :: Jx(numFields*totalModes,numTri)
    
    integer(i4) :: edgeNum,leftTri,rightTri
    
    
    edgeNum = interiorEdges(edge_id)
    leftTri  = edgeList(edgeNum)%e1    
    rightTri = edgeList(edgeNum)%e2
    
    Jx(:,leftTri) =  Jx(:,leftTri) &
                               + matmul(jacobOffDiagLR,x(:,rightTri))
                              
    Jx(:,rightTri) =  Jx(:,rightTri) &
                               + matmul(jacobOffDiagRL,x(:,leftTri))                                                  
    
end subroutine

Subroutine DiagJX(jacDiag,x,Jx)
    use LU_module
    
    real(dp), intent(in)    :: jacDiag(numFields*totalModes,numFields*totalModes,numTri)
    real(dp),intent(in)     :: x(numFields*totalModes,numTri)
    real(dp),intent(inout)  :: Jx(numFields*totalModes,numTri)
    real(dp)                :: tempX(numFields*totalModes)
    integer(i4) :: tri_id   
    
    do tri_id = 1,numTri
        tempX = x(:,tri_id)
        Jx(:,tri_id) = matmul(jacDiag(:,:,tri_id),tempX)
    end do
    
end subroutine

Subroutine DiagJXLU(jacDiag,x,Jx)
    use LU_module
    
    real(dp), intent(in)    :: jacDiag(numFields*totalModes,numFields*totalModes,numTri)
    real(dp),intent(in)     :: x(numFields*totalModes,numTri)
    real(dp),intent(inout)  :: Jx(numFields*totalModes,numTri)
    real(dp)                :: tempX(numFields*totalModes)
    integer(i4) :: tri_id   
    
    do tri_id = 1,numTri
        tempX = x(:,tri_id)
        Call lu_matvec(jacDiag(:,:,tri_id),tempX,Jx(:,tri_id))
    end do
    
end subroutine

Subroutine DiagJXLUTri(jacDiag,x,tri_id,Jx)
    use LU_module
    
    real(dp), intent(in)    :: jacDiag(numFields*totalModes,numFields*totalModes,numTri)
    real(dp), intent(in)    :: x(numFields*totalModes,numTri)
    integer(i4), intent(in) :: tri_id   
    real(dp),intent(inout)  :: Jx(numFields*totalModes,numTri)
    real(dp)                :: tempX(numFields*totalModes)

    tempX = x(:,tri_id)
    Call lu_matvec(jacDiag(:,:,tri_id),tempX,Jx(:,tri_id))
    
end subroutine

Subroutine form_Jx(jacDiagLU,jacobDiag,jacobOffDiagRL,jacobOffDiagLR,x,Jx)
    logical, intent(in) :: jacDiagLU
    real(dp),intent(in) :: jacobDiag(numFields*totalModes,numFields*totalModes,numTri)
    real(dp),intent(in) :: jacobOffDiagRL(numFields*totalModes,numFields*totalModes,numInterior)
    real(dp),intent(in) :: jacobOffDiagLR(numFields*totalModes,numFields*totalModes,numInterior)
    real(dp),intent(in) :: x(numFields*totalModes,numTri)
    real(dp),intent(inout):: Jx(numFields*totalModes,numTri)
    
    integer(i4) :: edge_id
    
    if(jacDiagLU)then
        Call DiagJXLU(jacobDiag,x,Jx)
    else
        Call DiagJX(jacobDiag,x,Jx)
    end if
    
    do edge_id = 1,numInterior
        Call offDiagJX(edge_id,jacobOffDiagRL(:,:,edge_id),jacobOffDiagLR(:,:,edge_id),x,Jx)
    end do
    
end subroutine


end module form_Jx_module

