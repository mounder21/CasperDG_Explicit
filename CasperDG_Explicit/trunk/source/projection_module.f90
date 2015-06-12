Module projection_module
     use my_kinddefs
     use inputs_module, only: totalModes
     use Globals_module, only: numEulerVars
     implicit none
 contains
 
 Subroutine projectEdge(solcoeffs,locEdge,q)
     use initializeBasis_module, only: edgePhi
     implicit none
     
     real(dp),      intent(in) :: solcoeffs(:,:)      !solcoeffs(numEulerVars,totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:,:)             !q(numEulerVars,numEdgeGaussPts)
     
     !dummy
     integer(i4) :: j

     q(:,:) = 0.0_dp
     do j = 1, numEulerVars
         q(j,:) = matmul(solcoeffs(j,:),edgePhi(:,:,locEdge))
     end do 
      
 end subroutine projectEdge
 
 subroutine projectEdgeMapStr8(ax,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiStr8
     implicit none
     
     real(dp),      intent(in) :: ax(:)               !ax(totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:)             !q(numEdgeGaussPts)
     
     q(:) = matmul(ax,edgeMapPhiStr8(:,:,locEdge))
     
 end subroutine projectEdgeMapStr8
 
 subroutine projectEdgeMapCurved(ax,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiCurved
     implicit none
     
     real(dp),      intent(in) :: ax(:)               !ax(totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:)             !q(numEdgeGaussPts)
     
     q(:) = matmul(ax,edgeMapPhiCurved(:,:,locEdge))
     
 end subroutine projectEdgeMapCurved
 
 Subroutine projectCell(solcoeff,thisPhi,q)
     real(dp),      intent(in) :: solcoeff(:,:)      !solcoeffs(numEulerVars,tm)
     real(dp),      intent(in) :: thisPhi(:,:)       !thisPhi(tm,numPts)
     real(dp),      intent(out) :: q(:,:)            !q(numEulerVars,numPts)
     
     !dummy
     integer(i4) :: j
     
     !Phi(totalModes,numGaussPts)
     q(:,:) = 0.0_dp

     do j = 1, numEulerVars
         q(j,:) = matmul(solcoeff(j,:),thisPhi(:,:))
     end do 
      
 end subroutine projectCell
 
  Subroutine projectMapCell(ax,ay,x,y)
     use initializeBasis_module, only: mapPhiStr8       !size(mapTotalModesStr8 = 3,numGaussPts)
     real(dp),      intent(in) :: ax(:),ay(:)      !(tm)
     real(dp),      intent(out) :: x(:),y(:)       !(numPts)      
     
     x(:) = matmul(ax(:),mapPhiStr8(:,:))
     y(:) = matmul(ay(:),mapPhiStr8(:,:))
      
 end subroutine projectMapCell
 
 
end module projection_module