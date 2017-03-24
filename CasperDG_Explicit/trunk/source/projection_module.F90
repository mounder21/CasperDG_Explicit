Module projection_module
     use my_kinddefs
     use inputs_module, only: totalModes
     use Globals_module, only: numFields
     implicit none
 contains
 
 Subroutine projectEdge(solcoeffs,locEdge,q)
     use initializeBasis_module, only: edgePhi
     implicit none
     
     real(dp),      intent(in) :: solcoeffs(:,:)      !solcoeffs(numFields,totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:,:)             !q(numFields,numEdgeGaussPts)
     
     !dummy
     integer(i4) :: j

     q(:,:) = 0.0_dp
     do j = 1, numFields
         q(j,:) = matmul(solcoeffs(j,:),edgePhi(:,:,locEdge))
     end do 
      
 end subroutine projectEdge
 
 
 subroutine projectEdgeMapStr8(a,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiStr8
     implicit none
     
     real(dp),      intent(in) :: a(:)               !ax(totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:)             !q(numEdgeGaussPts)
     
     q(:) = matmul(a,edgeMapPhiStr8(:,:,locEdge))
     
 end subroutine projectEdgeMapStr8
 
 
 subroutine projectEdgeMapCurved(a,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiCurved
     implicit none
     
     real(dp),      intent(in) :: a(:)               !a(totalModes)
     integer(i4),   intent(in) :: locEdge
     real(dp),      intent(out) :: q(:)             !q(numEdgeGaussPts)
     
     q(:) = matmul(a,edgeMapPhiCurved(:,:,locEdge))
     
 end subroutine projectEdgeMapCurved
 
 Subroutine projectCell(solcoeff,thisPhi,q)
     real(dp),      intent(in) :: solcoeff(:,:)      !solcoeffs(numFields,tm)
     real(dp),      intent(in) :: thisPhi(:,:)       !thisPhi(tm,numPts)
     real(dp),      intent(out) :: q(:,:)            !q(numFields,numPts)
     
     !dummy
     integer(i4) :: j
     
     !Phi(totalModes,numGaussPts)
     q(:,:) = 0.0_dp

     do j = 1, numFields
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