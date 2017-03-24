Module integrateComplex_module
    use my_kinddefs
    use Globals_module
    use inputs_module, only: totalModes
    use getElementJacobian_module

    implicit none

contains
    
Subroutine integrateCellComplex(elem_id,phiW,f,b)
    !integrates phi * f over each element
    integer(i4),intent(in) :: elem_id
    real(dp),intent(in) :: phiW(:,:)    ! (tm,np)
    
#ifdef complx
    complex(dp),intent(in) :: f(:)         ! (np)
    complex(dp),intent(out):: b(:)         ! (tm)
#else
    real(dp),intent(in) :: f(:)         ! (np)
    real(dp),intent(out):: b(:)         ! (tm)
#endif
    
    integer(i4) :: tm,GP
 
    b(:) = 0
    do tm = 1,totalModes
        do GP = 1, size(f)
            b(tm) = b(tm) + phiW(tm,GP)*f(GP)*detJ(GP,elem_id)
        end do
    end do

end subroutine integrateCellComplex


Subroutine integrateEdgeComplex(edgePhiW,f,b)
    !integrates phi * f over edge
    real(dp),intent(in) :: edgePhiW(:,:)    ! premultiplied phi * weights for Gauss Quadrature 
#ifdef complx
    complex(dp),intent(in) :: f(:)         ! (np)
    complex(dp),intent(out):: b(:)         ! (tm)
#else
    real(dp),intent(in) :: f(:)         ! (np)
    real(dp),intent(out):: b(:)         ! (tm)
#endif

    
    b(:) = 0.0_dp
    b(:) = matmul(edgePhiW(:,:),f)
end subroutine integrateEdgeComplex

Subroutine integrateComplex(elem_id,f,b)
    !integrates f over cell
    use initializeBasis_module, only: GaussWeights,numGaussPts
    
    integer(i4),intent(in) :: elem_id
    
#ifdef complx
    complex(dp),intent(in) :: f(:)         ! (np)
    complex(dp),intent(out):: b
#else
    real(dp),intent(in) :: f(:)         ! (np)
    real(dp),intent(out):: b
#endif

    integer(i4) :: GP
    
    b = 0._dp
    do GP = 1, numGaussPts
        b = b + GaussWeights(GP)*f(GP)*detJ(GP,elem_id)
    end do
end subroutine

    
end module integrateComplex_module
