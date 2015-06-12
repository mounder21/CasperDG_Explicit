Module residual_module
    use my_kinddefs
    use projection_module
    use flux_module
    use Globals_module, only: edgeList,numEdges,interiorEdges,boundEdges,numInterior,numTri,bcFlag
    use inputs_module,  only: totalModes,basisDegree,sourceTerm 
    use initializeBasis_module, only: PhiW,edgePhiW,dPhiW,Phi,edgePhi,numEdgeGaussPts,numGaussPts
    use integrate_module
    use bc_module
    use norm_module
    use computeJ_detJ_module, only: JInv,detJ
    use getIC_module
    use computeNormal_module, only: normals
    implicit none
    
contains

subroutine Residual(solCoeffs,spRes)
    real(dp),   intent(in) :: solCoeffs(:,:,:)
    real(dp),   intent(out) :: spRes(:,:,:)
    
    integer(i4) :: edge_id,leftTri,rightTri
    integer(i4) :: loc1,loc2,edgeNum,n,bc_type,face1,face2,GP
    real(dp)    :: qL(numEulerVars,numEdgeGaussPts)
    real(dp)    :: qR(numEulerVars,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numEulerVars,numEdgeGaussPts)

                  
    spRes(:,:,:)                = 0.0_dp

    ! Interior Edge Fluxes
    !print*,'before              :'
    do edge_id = 1, numInterior
        edgeNum = interiorEdges(edge_id)
        
        leftTri  = edgeList(edgeNum)%e1    
        rightTri = edgeList(edgeNum)%e2
        
        face1 = edgeList(edgeNum)%locEdge1
        face2 = edgeList(edgeNum)%locEdge2

        loc1 = 2*face1  - 1   !get correct edge and orientation (positive orient)
        loc2 = 2*face2       !get correct edge and orientation (negative orient)
        
        call projectEdge(solCoeffs(:,:,leftTri),loc1,qL)
        call projectEdge(solCoeffs(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,'LaxF',flux) !flux(4,#EdgeGPs)

        do n = 1,numEulerVars     
            call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
            call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
            spRes(n,:,leftTri)  = spRes(n,:,leftTri)  - bL(:)       !(4,totalModes,numTri)
            spRes(n,:,rightTri) = spRes(n,:,rightTri) + bR(:)
        end do

    end do
    !print*,'after           :'
    
    ! Boundary Edge Fluxes
    triList(:)%used = .false.
    do edge_id = 1, (numEdges-numInterior)
        edgeNum = boundEdges(edge_id)
        bc_type = bcFlag(edgeNum)
        leftTri = edgeList(edgeNum)%e1
        
        face1 = edgeList(edgeNum)%locEdge1
        loc1 = 2*face1 -1     !get correct end and orientation (positive orient)

        call projectEdge(solCoeffs(:,:,leftTri),loc1,qL)
        
        !form qR
        select case(bc_type)
            case (12)  
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                    call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,'LaxF',flux)
                    do n =1,numEulerVars     
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                        call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                        spRes(n,:,leftTri)  = spRes(n,:,leftTri)  - bL(:)       !(4,totalModes,numTri)
                        spRes(n,:,rightTri) = spRes(n,:,rightTri) + bR(:)
                    end do
                end if
            case(14)    !inflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
  
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n =1,numEulerVars     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
                end do
            case(15)    !outflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                
                do n = 1,numEulerVars     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
                end do
            case default 
                do GP = 1,numEdgeGaussPts
                    call get_bc(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP))
                end do
                
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                 
                do n = 1,numEulerVars     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
                end do
        end select
    end do
    
    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxRes(solCoeffs,spRes)
    end if
    if (sourceTerm) then
        Call sourceTermRes(spRes)
    end if
end subroutine Residual

Subroutine integralFluxRes(solCoef,spRes)
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numEulerVars,tm,numTri)
    real(dp),   intent(inout) :: spRes(numEulerVars,totalModes,numTri)
    
    real(dp) :: Res(numEulerVars,totalModes),F_xi(numEulerVars,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numEulerVars,numGaussPts),F_eta(numEulerVars,numGaussPts),flux(numEulerVars,2)
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)
    integer(i4) :: n,tri_id,GP

    ! dPhiW(totalModes,numGaussPts,2)
    Res(:,:) = 0.0_dp
    do tri_id = 1,numTri
        !This is correct (3/8/13)
        dXidx(:)    = JInv(1,1,:,tri_id) ! (numGPs)
        dXidy(:)    = JInv(2,1,:,tri_id) 
        dEtadx(:)   = JInv(1,2,:,tri_id) 
        dEtady(:)   = JInv(2,2,:,tri_id)
        
        call projectCell(solCoef(:,:,tri_id),Phi,q) !q(numEulerVars,numGaussPts)
        
        do n = 1,numEulerVars  
            do GP = 1,numGaussPts
                call getNativeFlux(q(:,GP),flux) !-- flux(4,2)
                F_xi(n,GP)  = flux(n,1)*dXidx(GP)  + flux(n,2)*dXidy(GP)
                F_eta(n,GP) = flux(n,1)*dEtadx(GP) + flux(n,2)*dEtady(GP)
            end do
            
            call integrateCell(tri_id,dPhiW(:,:,1),F_xi(n,:), b1)
            call integrateCell(tri_id,dPhiW(:,:,2),F_eta(n,:),c1)   
            Res(n,:) = b1(:)  +  c1(:) 
            spRes(n,:,tri_id) = spRes(n,:,tri_id) + Res(n,:)             !Res(4,totalModes)
        end do
    end do
    
end subroutine integralFluxRes

Subroutine sourceTermRes(spRes)
    use getSourceTerm_module
    real(dp), intent(inout) :: spRes(numEulerVars,totalModes,numTri)
    
    real(dp) :: Res(numEulerVars,totalModes)
    real(dp) :: source(numEulerVars,numGaussPts)
    integer(i4) :: tri_id, n
    
    do tri_id = 1,numTri  
        Call getSourceTerm(tri_id,numGaussPts,source) !source(numEulerVars,numGausspts)
        do n = 1,numEulerVars  
            call integrateCell(tri_id,PhiW(:,:),source(n,:),Res(n,:))    
            spRes(n,:,tri_id) = spRes(n,:,tri_id) + Res(n,:)             !Res(4,totalModes)
        end do  
    end do
    
end subroutine sourceTermRes
    

end module residual_module
