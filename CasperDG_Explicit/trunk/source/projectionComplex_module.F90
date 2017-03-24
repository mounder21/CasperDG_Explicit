Module projectionComplex_module
     use my_kinddefs
     use inputs_module, only: totalModes
     use Globals_module, only: numFields
     implicit none
 contains
 
 
 Subroutine projectEdgeComplex(solcoeffs,locEdge,q)
     use initializeBasis_module, only: edgePhi
     implicit none
     
#ifdef complx
     complex(dp),      intent(in) :: solcoeffs(:,:)   !solcoeffs(numFields,totalModes)
     complex(dp),      intent(out):: q(:,:)           !q(numFields,numEdgeGaussPts)
#else
     real(dp),      intent(in) :: solcoeffs(:,:)      !solcoeffs(numFields,totalModes)
     real(dp),      intent(out):: q(:,:)              !q(numFields,numEdgeGaussPts)
#endif
     integer(i4),   intent(in) :: locEdge
     
     
     !dummy
     integer(i4) :: j

     q(:,:) = 0.0_dp
     do j = 1, numFields
         q(j,:) = matmul(solcoeffs(j,:),edgePhi(:,:,locEdge))
     end do 
      
 end subroutine projectEdgeComplex

 
 subroutine projectEdgeMapStr8Complex(a,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiStr8
     implicit none
#ifdef complx
     complex(dp),   intent(in) :: a(:)             !ax(totalModes)
     complex(dp),   intent(out):: q(:)             !q(numEdgeGaussPts)
#else
     real(dp),      intent(in) :: a(:)               !ax(totalModes)
     real(dp),      intent(out):: q(:)             !q(numEdgeGaussPts)
#endif   
     integer(i4),   intent(in) :: locEdge

     
     q(:) = matmul(a,edgeMapPhiStr8(:,:,locEdge))
     
 end subroutine projectEdgeMapStr8Complex
 
 
 subroutine projectEdgeMapCurvedComplex(a,locEdge,q)
     use initializeBasis_module, only: edgeMapPhiCurved
     implicit none
#ifdef complx
     complex(dp),   intent(in) :: a(:)             !a(totalModes)
     complex(dp),   intent(out):: q(:)             !q(numEdgeGaussPts)
#else
     real(dp),      intent(in) :: a(:)               !ax(totalModes)
     real(dp),      intent(out):: q(:)             !q(numEdgeGaussPts)
#endif 
     integer(i4),   intent(in) :: locEdge

     
     q(:) = matmul(a,edgeMapPhiCurved(:,:,locEdge))
     
 end subroutine projectEdgeMapCurvedComplex

 
 Subroutine projectCellComplex(solcoeff,thisPhi,q)
#ifdef complx
     complex(dp),   intent(in) :: solcoeff(:,:)             !a(totalModes)
     complex(dp),   intent(out):: q(:,:)             !q(numEdgeGaussPts)
#else
     real(dp),      intent(in) :: solcoeff(:,:)               !ax(totalModes)
     real(dp),      intent(out):: q(:,:)             !q(numEdgeGaussPts)
#endif 
     real(dp),      intent(in) :: thisPhi(:,:)       !thisPhi(tm,numPts)

     !dummy
     integer(i4) :: j
     
     !Phi(totalModes,numGaussPts)
     q(:,:) = 0.0_dp

     do j = 1, numFields
         q(j,:) = matmul(solcoeff(j,:),thisPhi(:,:))
     end do 
      
 end subroutine projectCellComplex

 
 Subroutine projectMapCellComplex(ax,ay,x,y)
     use initializeBasis_module, only: mapPhiStr8       !size(mapTotalModesStr8 = 3,numGaussPts)
     implicit none
#ifdef complx
     complex(dp),   intent(in) :: ax(:),ay(:)      !(tm)
     complex(dp),   intent(out) :: x(:),y(:)       !(numPts)      
#else
     real(dp),      intent(in) :: ax(:),ay(:)      !(tm)
     real(dp),      intent(out) :: x(:),y(:)       !(numPts)      
#endif 
     
     
     x(:) = matmul(ax(:),mapPhiStr8(:,:))
     y(:) = matmul(ay(:),mapPhiStr8(:,:))
      
 end subroutine projectMapCellComplex
 
 
end module projectionComplex_module
