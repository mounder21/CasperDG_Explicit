 module basis_module
    use my_kinddefs
    implicit none 
 contains
 
 subroutine getPhi(basisType,basisDegree,totalModes,numGaussPts,pointsMatrix,Phi)
     character(*), intent(in)                   :: basisType
     integer(i4), intent(in)                    :: basisDegree,totalModes
     integer(i4), intent(in)                    :: numGaussPts
     real(dp), intent(in)                       :: pointsMatrix(:,:)
     real(dp), intent(out)                      :: Phi(:,:)
     
     select case (basisType)
         case('H1')  
             select case (basisDegree)
                 case(0)       
                     Phi(1,:) = 1.0_dp
                 case(1)
                     Phi(1,:) = -0.5_dp*(pointsMatrix(1,:) + pointsMatrix(2,:))
                     Phi(2,:) =  0.5_dp*(1._dp + pointsMatrix(1,:))
                     Phi(3,:) =  0.5_dp*(1._dp + pointsMatrix(2,:)) 
                 case(2)
                     Phi(1,:) = -0.5_dp * (pointsMatrix(1,:) + pointsMatrix(2,:))
                     Phi(2,:) =  0.5_dp * (1._dp + pointsMatrix(1,:))
                     Phi(3,:) =  0.5_dp * (1._dp + pointsMatrix(2,:))
                     Phi(4,:) =  Phi(2,:)*Phi(3,:)
                     Phi(5,:) =  Phi(1,:)*Phi(3,:) 
                     Phi(6,:) =  Phi(1,:)*Phi(2,:)
                 case default      
                    print*,'we dont have a better basis right now,getPhi H1'
             end select
         case ('nodal')
             select case (basisDegree)
                 case(0)       
                     Phi(1,:) = 1.0_dp
                 case(1)
                     Phi(1,:) = -0.5_dp*(pointsMatrix(1,:) + pointsMatrix(2,:)) 
                     Phi(2,:) =  0.5_dp*(1._dp + pointsMatrix(1,:))
                     Phi(3,:) =  0.5_dp*(1._dp + pointsMatrix(2,:)) 
                 case(2)
                     Phi(1,:) = -0.5_dp*(pointsMatrix(1,:) + pointsMatrix(2,:)) * &
                                     (pointsMatrix(1,:) + pointsMatrix(2,:) + 1._dp)
                     Phi(2,:) =  0.5_dp * (1._dp + pointsMatrix(1,:)) * pointsMatrix(1,:)
                     Phi(3,:) =  0.5_dp * (1._dp + pointsMatrix(2,:)) * pointsMatrix(2,:)
                     Phi(4,:) =  4.0_dp * Phi(2,:)*Phi(3,:)
                     Phi(5,:) =  4.0_dp * Phi(1,:)*Phi(3,:)
                     Phi(6,:) =  4.0_dp * Phi(1,:)*Phi(2,:)
                 case default      
                    print*,'we dont have a better basis right now,getPhi nodal'
             end select
         case ('monomial')
              select case (basisDegree)
                 case(0)       
                     Phi(1,:) = 1.0_dp
                 case(1)
                     Phi(1,:) = 1.0_dp
                     Phi(2,:) = pointsMatrix(1,:)
                     Phi(3,:) = pointsMatrix(2,:)
                 case(2)
                     Phi(1,:) = 1.0_dp
                     Phi(2,:) = pointsMatrix(1,:)
                     Phi(3,:) = pointsMatrix(2,:)
                     Phi(4,:) = pointsMatrix(1,:) ** 2
                     Phi(5,:) = pointsMatrix(2,:) ** 2
                     Phi(6,:) = pointsMatrix(1,:) * pointsMatrix(2,:)
                 case default      
                    print*,'we dont have a better basis right now,getPhi monomial'
             end select
     end select
     
 end subroutine getPhi

 
 subroutine getdPhi(basisType,basisDegree,totalModes,numGaussPts,pointsMatrix,dPhi)
     character(*), intent(in)                   :: basisType    !use Solution_moduleBasisType
     integer(i4), intent(in)                    :: basisDegree,totalModes
     integer(i4), intent(in)                    :: numGaussPts
     real(dp), intent(in)                       :: pointsMatrix(:,:)
     real(dp), intent(out)                      :: dPhi(:,:,:) !(totalModes,numGaussPts,xi/eta)
     
     select case (basisType)
         case('H1')  
             select case (basisDegree)
                 case(0)       
                     dPhi(1,:,:) = 0.0_dp
                 case(1)
                     dPhi(1,:,1) = -0.5_dp   !dPhi(1)/dXi
                     dPhi(2,:,1) =  0.5_dp   !dPhi(2)/dXi
                     dPhi(3,:,1) =  0.0_dp   !dPhi(3)/dXi
                     
                     dPhi(1,:,2) = -0.5_dp   !dPhi(1)/dEta
                     dPhi(2,:,2) =  0.0_dp   !dPhi(2)/dEta
                     dPhi(3,:,2) =  0.5_dp   !dPhi(3)/dEta
                 case(2)
                     dPhi(1,:,1) = -0.5_dp   !dPhi(1)/dXi
                     dPhi(2,:,1) =  0.5_dp   !dPhi(2)/dXi
                     dPhi(3,:,1) =  0.0_dp   !dPhi(3)/dXi
                     
                     dPhi(1,:,2) = -0.5_dp   !dPhi(1)/dEta
                     dPhi(2,:,2) =  0.0_dp   !dPhi(2)/dEta
                     dPhi(3,:,2) =  0.5_dp   !dPhi(3)/dEta
                     
                     dPhi(4,:,1) =   0.25_dp*(1.0_dp + pointsMatrix(2,:))
                     dPhi(5,:,1) =  -0.25_dp*(1.0_dp + pointsMatrix(2,:))
                     dPhi(6,:,1) =  -0.25_dp*(1.0_dp + 2.0_dp*pointsMatrix(1,:) + pointsMatrix(2,:))
                     
                     dPhi(4,:,2) =   0.25_dp*(1.0_dp + pointsMatrix(1,:))
                     dPhi(5,:,2) =  -0.25_dp*(1.0_dp + pointsMatrix(1,:) + 2.0_dp*pointsMatrix(2,:))
                     dPhi(6,:,2) =  -0.25_dp*(1.0_dp + pointsMatrix(1,:))
                 case default      
                    print*,'we dont have a better basis right now,getdPhi h1'
             end select
         case ('nodal')
             select case (basisDegree)
                 case(0)       
                     dPhi(1,:,:) = 0.0_dp
                 case(1)
                     dPhi(1,:,1) = -0.5_dp   !dPhi(1)/dXi
                     dPhi(2,:,1) =  0.5_dp   !dPhi(2)/dXi
                     dPhi(3,:,1) =  0.0_dp   !dPhi(3)/dXi
                     
                     dPhi(1,:,2) = -0.5_dp   !dPhi(1)/dEta
                     dPhi(2,:,2) =  0.0_dp   !dPhi(2)/dEta
                     dPhi(3,:,2) =  0.5_dp   !dPhi(3)/dEta
                 case(2)
                     dPhi(1,:,1) = -0.5_dp * (1.0_dp + 2.0_dp*pointsMatrix(1,:) + 2.0_dp*pointsMatrix(2,:))
                     dPhi(2,:,1) =  0.5_dp * (1.0_dp + 2.0_dp*pointsMatrix(1,:))
                     dPhi(3,:,1) =  0.0_dp
                     dPhi(4,:,1) =(pointsMatrix(2,:)+pointsMatrix(2,:)**2)*(1.0_dp+2.0_dp*pointsMatrix(1,:))
                     dPhi(5,:,1) =-(pointsMatrix(2,:)+pointsMatrix(2,:)**2)*&
                                      (1.0_dp + 2.0_dp*pointsMatrix(2,:) + 2.0_dp*pointsMatrix(1,:))
                     dPhi(6,:,1) =-(4._dp*pointsMatrix(1,:)**3 +6._dp*pointsMatrix(1,:)**2 +&
                                       2._dp*pointsMatrix(1,:) +6._dp*pointsMatrix(2,:)*pointsMatrix(1,:) +&
                                       6._dp*pointsMatrix(2,:)*pointsMatrix(1,:)**2 +pointsMatrix(2,:) + &
                                       pointsMatrix(2,:)**2 + 2._dp*pointsMatrix(1,:)*pointsMatrix(2,:)**2)
                     
                     dPhi(1,:,2) = -0.5_dp *(1.0_dp + 2.0_dp*pointsMatrix(2,:) + 2.0_dp*pointsMatrix(1,:))
                     dPhi(2,:,2) =  0.0_dp
                     dPhi(3,:,2) =  0.5_dp * (1._dp + 2.0_dp*pointsMatrix(2,:)) 
                     dPhi(4,:,2) =(pointsMatrix(1,:)+pointsMatrix(1,:)**2)*(1.0_dp+2.0_dp*pointsMatrix(2,:))
                     dPhi(5,:,2) =-(4._dp*pointsMatrix(2,:)**3 +6._dp*pointsMatrix(2,:)**2 +&
                                       2._dp*pointsMatrix(2,:) +6._dp*pointsMatrix(1,:)*pointsMatrix(2,:) +&
                                       6._dp*pointsMatrix(1,:)*pointsMatrix(2,:)**2 +pointsMatrix(1,:) + &
                                       pointsMatrix(1,:)**2 + 2._dp*pointsMatrix(2,:)*pointsMatrix(1,:)**2) 
                     dPhi(6,:,2) =  -(pointsMatrix(1,:)+pointsMatrix(1,:)**2)*&
                                      (1.0_dp + 2.0_dp*pointsMatrix(2,:) + 2.0_dp*pointsMatrix(1,:))
                 case default      
                    print*,'we dont have a better basis right now,getdPhi nodal'
             end select
         case ('monomial')
              select case (basisDegree)
                 case(0)       
                     dPhi(1,:,:) = 0.0_dp
                 case(1)
                     dPhi(1,:,1) = 0.0_dp
                     dPhi(2,:,1) = 1.0_dp
                     dPhi(3,:,1) = 0.0_dp
                     
                     dPhi(1,:,2) = 0.0_dp
                     dPhi(2,:,2) = 0.0_dp
                     dPhi(3,:,2) = 1.0_dp
                 case(2)
                     dPhi(1,:,1) = 0.0_dp
                     dPhi(2,:,1) = 1.0_dp
                     dPhi(3,:,1) = 0.0_dp
                     dPhi(4,:,1) = 2.0_dp*pointsMatrix(1,:)
                     dPhi(5,:,1) = 0.0_dp
                     dPhi(6,:,1) = pointsMatrix(2,:)
                     
                     dPhi(1,:,2) = 1.0_dp
                     dPhi(2,:,2) = 0.0_dp
                     dPhi(3,:,2) = 1.0_dp
                     dPhi(4,:,2) = 0.0_dp
                     dPhi(5,:,2) = 2.0_dp*pointsMatrix(2,:)
                     dPhi(6,:,2) = pointsMatrix(1,:) 
                 case default      
                    print*,'we dont have a better basis right now, getdPhi monomial'
             end select
     end select
     
 end subroutine getdPhi

 
 
 subroutine getEdgePhi(basisType,basisDegree,totalModes,numEdgeGaussPts,edgeGaussPts,edgePhi)
     character(*), intent(in)                   :: basisType
     integer(i4), intent(in)                    :: basisDegree,totalModes
     integer(i4), intent(in)                    :: numEdgeGaussPts
     real(dp), intent(in)                       :: edgeGaussPts(:)
     real(dp), pointer, intent(out)             :: edgePhi(:,:,:)
     
     ! dummy
     integer(i4) :: faceIndex,n
     
     real(dp) :: temp(totalModes, numEdgeGaussPts),ptsMatrix(2,numEdgeGaussPts)
     real(dp) :: revEdgeGaussPts(numEdgeGaussPts)
     real(dp) :: xiV(numEdgeGaussPts),etaV(numEdgeGaussPts)
     real(dp) :: sqrtTwo
     sqrtTwo = 1.4142135623730950488016887242096980785696718753769480731766797379907324784621_dp
     
     !edgePhi(totalModes {phi1,phi2,phi3,phi_n} , numGaussPts, faceIndex* orientation)
     
     temp(:,:) = 0.0_dp

     edgePhi(:,:,:) = 0.0_dp
     n = size(edgeGaussPts)
     revEdgeGaussPts = edgeGaussPts(n:1:-1)

     do faceIndex = 1, 3
         select case (faceIndex)
             !hypotenuse 
             case(1)
                !phi1
                 xiV =  edgeGaussPts
                 etaV = -xiV
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)

                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,1) = temp(:,:)
                
                 xiV =  revEdgeGaussPts
                 etaV = -xiV
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,2) = temp(:,:)
                 
             !vertical edge (xi = -1, eta)   
             case(2)
                 xiV(:) = -1   
                 etaV = edgeGaussPts
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,3) = temp(:,:)
                 
                 etaV = revEdgeGaussPts                 
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,4) = temp(:,:)

             !horizontal edge (xi, eta = -1)   
             case(3)
                 xiV = revEdgeGaussPts
                 etaV(:) = -1
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,5) = temp(:,:)

                 xiV = edgeGaussPts
                 ptsmatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgePhi(:,:,6) = temp(:,:)        
         end select
     end do
     
 end subroutine getEdgePhi
 
 subroutine getEdge_dPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,edgeGaussPts,edgedPhi)
 character(*), intent(in)                   :: basisType
     integer(i4), intent(in)                :: basisDegree,totalModes
     integer(i4), intent(in)                :: numEdgeGaussPts
     real(dp), intent(in)                   :: edgeGaussPts(:)
     real(dp), intent(out)                  :: edgedPhi(:,:,:,:)
     
     ! dummy
     integer(i4) :: faceIndex,n
     
     
     real(dp) :: temp(totalModes, numEdgeGaussPts,2),ptsMatrix(2,numEdgeGaussPts)
     real(dp) :: xiV(numEdgeGaussPts),etaV(numEdgeGaussPts),revEdgeGaussPts(numEdgeGaussPts)
     real(dp) :: sqrtTwo
     sqrtTwo = 1.4142135623730950488016887242096980785696718753769480731766797379907324784621_dp
     
     !edgePhi(totalModes {phi1,phi2,phi3,phi_n} , numGaussPts, faceIndex* orientation) 
     temp(:,:,:) = 0.0_dp
     edgedPhi(:,:,:,:) = 0.0_dp
     n = size(edgeGaussPts)
     revEdgeGaussPts = edgeGaussPts(n:1:-1)
     
     do faceIndex = 1, totalModes
         select case (faceIndex)
             !hypotenuse 
             case(1)
                !phi1
                 xiV =  edgeGaussPts 
                 etaV = -xiV
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,1) = temp(:,:,:)

                 xiV = revEdgeGaussPts 
                 etaV = -xiV
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,2) = temp(:,:,:)

             !vertical edge (xi = -1, eta)   
             case(2)
                 xiV = -1
                 etaV = edgeGaussPts
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,3) = temp(:,:,:)

                 etaV = revEdgeGaussPts
                 ptsmatrix(2,:) = etaV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,4) = temp(:,:,:)

             !horizontal edge (xi, eta = -1)   
             case(3)
                 xiV = revEdgeGaussPts
                 etaV = -1
                 ptsMatrix(1,:) = xiV(:)
                 ptsmatrix(2,:) = etaV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,5) = temp(:,:,:)

                 xiV = edgeGaussPts
                 ptsmatrix(1,:) = xiV(:)
                 Call getdPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,ptsMatrix,temp)
                 edgedPhi(:,:,:,6) = temp(:,:,:)

         end select
     end do
 end subroutine getEdge_dPhi
 end module basis_module