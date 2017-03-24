Module norm_module
    use my_kinddefs
    implicit none
contains

Subroutine Lp_norm(vector,p,norm)
    use Globals_module, only: numTri
    implicit none
    
    real(dp),intent(in) :: vector(:,:,:)
    real(dp),intent(in) :: p
    real(dp),intent(out) :: norm

    
    integer(i4) index1,index2,index3
    real(dp) oneOverp
    
    norm = 0.0_dp
    oneOverp = 1.0_dp / p
    do index1 = 1,size(vector,1)
        do index2 = 1,size(vector,2)
            do index3 = 1,size(vector,3)
                norm = norm + vector(index1,index2,index3)**p
            end do
        end do
    end do
    norm = (norm/numTri)**(oneOverp)
    if (isnan(norm)) then
        print*, 'there are Nans'
        stop
    end if
end subroutine

Subroutine Lp_norm2D(vector,p,norm)
    real(dp),intent(in) :: vector(:,:)
    real(dp),intent(in) :: p
    real(dp),intent(out) :: norm
    
    integer(i4) index1,index2
    real(dp) oneOverp
    
    norm = 0.0_dp
    oneOverp = 1.0_dp / p
    do index1 = 1,size(vector,1)
        do index2 = 1,size(vector,2)
            norm = norm + vector(index1,index2)**p
        end do
    end do
    norm = norm**(oneOverp)
    if (isnan(norm)) then
        print*, 'there are Nans'
        stop
    end if
end subroutine
    
end module norm_module