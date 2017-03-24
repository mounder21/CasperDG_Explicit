Module assembleLocMassMatrices_module
    use my_kinddefs
    use Globals_module
    use luMatrixInverse_module
    use getElementJacobian_module
    use inputs_module, only: totalModes
    implicit none 
    
    real(dp), dimension(:,:),allocatable   :: localMassMatrixPre,locMassTemp,localMassMatrix,localInvMassMatrix
    real(dp), dimension(:,:,:), pointer    :: MassMatrices,invMassMatrices
contains

Subroutine assembleLocMassMatrix(cellPhi,cellPhiW)
    real(dp), intent(in):: cellPhi(:,:),cellPhiW(:,:)

    real(dp)                    :: elemJacobian
    integer(i4)                 :: elem_id

    Call allocate_massMatrices
    localMassMatrixPre(:,:) = matmul(cellPhi(:,:),transpose(cellPhiW(:,:)))

    do elem_id = 1, numTri
      if (triList(elem_id)%straight .eqv. .true.) then
        Call getElementJacobianStr8(elem_id,triList(elem_id)%elemType,elemJacobian)
        localMassMatrix(:,:) = localMassMatrixPre(:,:) * elemJacobian
        !need a temp since the matrix gets destroyed
        locMassTemp(:,:) = localMassMatrixPre(:,:) * elemJacobian
        Call inverseMat(locMassTemp,localInvMassMatrix,size(locMassTemp,2))
       
      else
          !Call getElementJacobianCurved(elem_id,triList(elem_id)%elemType,elemJacobian)
          !localMassMatrix(:,:) = localMassMatrixPre(:,:)* elemJacobian 
          !need a temp since the matrix gets destroyed
          !locMassTemp(:,:) = localMassMatrixPre(:,:) * elemJacobian
          !Call inverseMat(locMassTemp,localInvMassMatrix,size(locMassTemp,2))
          print*, 'We cant handle curved elements yet,assembleMass'
      end if
      MassMatrices(:,:,elem_id)     = localMassMatrix(:,:)
      invMassMatrices(:,:,elem_id)  = localInvMassMatrix(:,:)

    end do
end subroutine assembleLocMassMatrix

Subroutine allocate_massMatrices
    use inputs_module, only: totalModes
    use Globals_module, only: numTri
    implicit none
    
    Allocate(localMassMatrixPre(totalModes,totalModes))
    Allocate(locMassTemp(totalModes,totalModes))
    Allocate(localMassMatrix(totalModes,totalModes))
    Allocate(localInvMassMatrix(totalModes,totalModes))
    Allocate(MassMatrices(totalModes,totalModes,numTri))
    Allocate(invMassMatrices(totalModes,totalModes,numTri))
end subroutine

end module assembleLocMassMatrices_module