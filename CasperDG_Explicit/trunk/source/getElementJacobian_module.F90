Module getElementJacobian_module
    use my_kinddefs
    use Globals_module
    use computeJ_detJ_module, only: detJ
    implicit none 
    
contains
Subroutine getElementJacobianStr8(elementNum,elementType,eleJacobian)
    implicit none
    
    integer(i4), intent(in) :: elementNum
    character(*),intent(in) :: elementType
    real(dp), intent(out)   :: eleJacobian ! = det(Jacobean)
    !dummy

    if (elementType .eq. 'triangle') then
        eleJacobian = detJ(1,elementNum) !take first one since the are the same
    else
        print*, 'No implementation for elements other than triangle'
    end if 
end subroutine getElementJacobianStr8

Subroutine getElementJacobianCurved(elementNum,elementType,eleJacobian)
 implicit none
    
    integer(i4), intent(in) :: elementNum
    character(*),intent(in) :: elementType
    real(dp), intent(out)   :: eleJacobian 
    
     if (elementType .eq. 'triangle') then
        eleJacobian = 0.0_dp
        print*,'no curved yet: getElementJacoCurved'
        stop
     end if
end subroutine getElementJacobianCurved

end module getElementJacobian_module
