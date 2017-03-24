Module integrate_module
    use my_kinddefs
    use Globals_module
    use inputs_module, only: totalModes
    use getElementJacobian_module

    implicit none

contains
    
Subroutine integrateCell(elem_id,phiW,f,b)
    !integrates phi * f over each element
    integer(i4),intent(in) :: elem_id
    real(dp),intent(in) :: phiW(:,:)    ! (tm,np)
    real(dp),intent(in) :: f(:)         ! (np)
    real(dp),intent(out):: b(:)         ! (tm)
    integer(i4) :: tm,GP
 
    b(:) = 0
    do tm = 1,totalModes
        do GP = 1, size(f)
            b(tm) = b(tm) + phiW(tm,GP)*f(GP)*detJ(GP,elem_id)
        end do
    end do

end subroutine integrateCell


Subroutine integrateEdge(edgePhiW,f,b)
    !integrates phi * f over edge
    real(dp),intent(in) :: edgePhiW(:,:)    ! premultiplied phi * weights for Gauss Quadrature 
    real(dp),intent(in) :: f(:)
    real(dp),intent(out):: b(:)
    
    b(:) = 0.0_dp
    b(:) = matmul(edgePhiW(:,:),f)
end subroutine integrateEdge

Subroutine integrate(elem_id,f,b)
    !integrates f over cell
    use initializeBasis_module, only: GaussWeights,numGaussPts
    
    integer(i4),intent(in) :: elem_id
    real(dp),intent(in) :: f(:)         ! (np)
    real(dp),intent(out):: b            ! (np)
    integer(i4) :: GP
    
    b = 0._dp
    do GP = 1, numGaussPts
        b = b + GaussWeights(GP)*f(GP)*detJ(GP,elem_id)
    end do
end subroutine

Subroutine integrateEdgeJacDiag(edgePhiW,edgePhi,flux_q,JD)
    use initializeBasis_module, only: numEdgeGaussPts
    real(dp),intent(in) :: edgephiW(:,:)    ! (tm,np)
    real(dp),intent(in) :: edgephi(:,:)     ! (tm,np)
    real(dp),intent(in) :: flux_q(:,:,:)    ! (4,4,ngp)
    real(dp),intent(out):: JD(numFields,totalModes,numFields,totalModes)
    
    integer(i4) :: GP,m,r,n,s
    
    JD = 0._dp
    do GP = 1, numEdgeGaussPts
        do s = 1,totalModes
            do n = 1,numFields
                do r = 1,totalModes
                    do m = 1,numFields      
                        JD(m,r,n,s) = JD(m,r,n,s) + edgePhiW(r,GP)*flux_q(m,n,GP)*edgePhi(s,GP)
                    end do
                end do
            end do
        end do
    end do
    
end subroutine

Subroutine integrateVolJacDiag(elem_id,dPhiW,Phi,flux_q,JD)
    use initializeBasis_module, only: numGaussPts
    integer(i4),intent(in):: elem_id
    real(dp),intent(in) :: dPhiW(:,:)    ! (tm,np)
    real(dp),intent(in) :: phi(:,:)     ! (tm,np)
    real(dp),intent(in) :: flux_q(:,:,:)    ! (4,4))
    real(dp),intent(out):: JD(numFields,totalModes,numFields,totalModes)
    
    integer(i4) :: GP,m,r,n,s
    
    JD = 0._dp
    do GP = 1, numGaussPts
        do s = 1,totalModes
            do n = 1,numFields
                do r = 1,totalModes
                    do m = 1,numFields      
                        JD(m,r,n,s) = JD(m,r,n,s) + dPhiW(r,GP)*flux_q(m,n,GP)*Phi(s,GP)*detJ(GP,elem_id)
                    end do
                end do
            end do
        end do
    end do
    
end subroutine

       
end module integrate_module
