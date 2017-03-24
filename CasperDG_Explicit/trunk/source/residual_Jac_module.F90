Module residual_Jac_module
    use my_kinddefs
    use projection_module
    use flux_module
    use Globals_module, only: edgeList,numEdges,interiorEdges,boundEdges,numInterior,numTri,bcFlag
    use inputs_module,  only: totalModes,basisDegree,sourceTerm,fluxScheme 
    use initializeBasis_module, only: PhiW,edgePhiW,dPhiW,Phi,edgePhi,numEdgeGaussPts,numGaussPts
    use integrate_module
    use bc_module
    use norm_module
    use computeJ_detJ_module, only: JInv,detJ
    use getIC_module
    use computeNormal_module, only: normals
    implicit none
    
    real(dp),pointer :: spRes(:,:,:),jacobDiag(:,:,:),jacobOffDiagRL(:,:,:),jacobOffDiagLR(:,:,:)
    real(dp),pointer :: spResTri(:,:,:),jacobDiagTri(:,:,:),jacobOffDiagRLTri(:,:,:),jacobOffDiagLRTri(:,:,:)
contains

subroutine Residual_Jac(solCoeffs)
    real(dp),   intent(in) :: solCoeffs(:,:,:)
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)
    
    real(dp)    :: flux_ql(numFields,numFields,numEdgeGaussPts)
    real(dp)    :: flux_qr(numFields,numFields,numEdgeGaussPts)

    real(dp)    :: dqbdq(numEdgeGaussPts,numFields,numFields)
    real(dp)    :: JD(numFields*totalModes,numFields*totalModes)

                  
    spRes(:,:,:)                = 0.0_dp
    jacobDiag(:,:,:)            = 0.0_dp
    jacobOffDiagLR(:,:,:)       = 0.0_dp
    jacobOffDiagRL(:,:,:)       = 0.0_dp
    ! Interior Edge Fluxes

    do edge_id = 1, numInterior
        edgeNum = interiorEdges(edge_id)
        
        leftTri  = edgeList(edgeNum)%e1    
        rightTri = edgeList(edgeNum)%e2
        
        face1 = edgeList(edgeNum)%locEdge1
        face2 = edgeList(edgeNum)%locEdge2

        loc1 = 2*face1  - 1   !get correct edge and orientation (positive orient)
        loc2 = 2*face2        !get correct edge and orientation (negative orient)
        
        call projectEdge(solCoeffs(:,:,leftTri),loc1,qL)
        call projectEdge(solCoeffs(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFlux_q(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux,flux_ql,flux_qr) 

        do n = 1,numFields     
            call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
            call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
            spRes(n,:,leftTri)  = spRes(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
            spRes(n,:,rightTri) = spRes(n,:,rightTri) - bR(:)
        end do
        
        !Diagonal Jacobian
        Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
        jacobDiag(:,:,leftTri) = jacobDiag(:,:,leftTri)     + JD
        
        Call integrateEdgeJacDiag(edgePhiW(:,:,loc2),edgePhi(:,:,loc2),flux_qr,JD)
        jacobDiag(:,:,rightTri) = jacobDiag(:,:,rightTri)   - JD
        
        !Off Diagonal Jacobian
        Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc2),flux_qr,JD)
        jacobOffDiagLR(:,:,edge_id) = jacobOffDiagLR(:,:,edge_id) + JD
        
        Call integrateEdgeJacDiag(edgePhiW(:,:,loc2),edgePhi(:,:,loc1),flux_ql,JD)
        jacobOffDiagRL(:,:,edge_id) = jacobOffDiagRL(:,:,edge_id) - JD

    end do
    
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
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
                    call getFlux_q(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme ,flux,flux_ql,flux_qr)
                    do n =1,numFields     
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                        call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                        spRes(n,:,leftTri)  = spRes(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                        spRes(n,:,rightTri) = spRes(n,:,rightTri) - bR(:)
                    end do
                    
                    !Diagonal Jacobian
                    Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                    jacobDiag(:,:,leftTri) = jacobDiag(:,:,leftTri) + JD
                end if
            case(14)    !inflow
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
  
                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                end do
                
                do n =1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                end do
                
                !Diagonal Jacobian
                Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                jacobDiag(:,:,leftTri) = jacobDiag(:,:,leftTri)     + JD
            case(15)    !outflow
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
                
                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                end do
                
                do n = 1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
                
                !Diagonal Jacobian
                Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                jacobDiag(:,:,leftTri) = jacobDiag(:,:,leftTri) + JD
            case default 
                do GP = 1,numEdgeGaussPts
                    call get_bc_q(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP),dqbdq(GP,:,:))
                end do

                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                     
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))  
                end do
                 
                do n = 1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spRes(n,:,leftTri) = spRes(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
                
                !Diagonal Jacobian
                Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                jacobDiag(:,:,leftTri) = jacobDiag(:,:,leftTri) + JD
        end select
    end do
    
    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxRes_Jac(solCoeffs)
    end if
    if (sourceTerm) then
        Call sourceTermRes()
    end if
end subroutine Residual_Jac

Subroutine integralFluxRes_Jac(solCoef)
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    
    real(dp) :: Res(numFields,totalModes),F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)
    real(dp) :: flux_q(numFields,2,numFields)
    real(dp) :: flux_qXi(numFields,numFields,numGaussPts),flux_qEta(numFields,numFields,numGaussPts)
    real(dp) :: JD(numFields*totalModes,numFields*totalModes)
    integer(i4) :: n,tri_id,GP

    ! dPhiW(totalModes,numGaussPts,2)
    Res(:,:) = 0.0_dp
    do tri_id = 1,numTri
        !This is correct (3/8/13)
        dXidx(:)    = JInv(1,1,:,tri_id) ! (numGPs)
        dXidy(:)    = JInv(2,1,:,tri_id) 
        dEtadx(:)   = JInv(1,2,:,tri_id) 
        dEtady(:)   = JInv(2,2,:,tri_id)
        
        call projectCell(solCoef(:,:,tri_id),Phi,q) !q(numFields,numGaussPts)
        
        do n = 1,numFields  
            do GP = 1,numGaussPts
                call getNativeFlux_q(q(:,GP),flux,flux_q) !-- flux(4,2)
                F_xi(n,GP)  = flux(n,1)*dXidx(GP)  + flux(n,2)*dXidy(GP)
                F_eta(n,GP) = flux(n,1)*dEtadx(GP) + flux(n,2)*dEtady(GP)
                flux_qXi(:,:,GP)  = flux_q(:,1,:)*dXidx(GP)   + flux_q(:,2,:)*dXidy(GP)
                flux_qEta(:,:,GP) = flux_q(:,1,:)*dEtadx(GP)  + flux_q(:,2,:)*dEtady(GP)
            end do
            
            call integrateCell(tri_id,dPhiW(:,:,1),F_xi(n,:), b1)
            call integrateCell(tri_id,dPhiW(:,:,2),F_eta(n,:),c1)   
            Res(n,:) = b1(:)  +  c1(:) 
            spRes(n,:,tri_id) = spRes(n,:,tri_id) - Res(n,:)             !Res(4,totalModes)
        end do
        
        Call integrateVolJacDiag(tri_id,dPhiW(:,:,1),Phi,flux_qXi,JD) 
        jacobDiag(:,:,tri_id) = jacobDiag(:,:,tri_id) - JD
        Call integrateVolJacDiag(tri_id,dPhiW(:,:,2),Phi,flux_qEta,JD)
        jacobDiag(:,:,tri_id) = jacobDiag(:,:,tri_id) - JD
    end do
    
end subroutine integralFluxRes_Jac

Subroutine sourceTermRes()
    use getSourceTerm_module
    
    
    real(dp) :: Res(numFields,totalModes)
    real(dp) :: source(numFields,numGaussPts)
    integer(i4) :: tri_id, n
    
    do tri_id = 1,numTri  
        Call getSourceTerm(tri_id,numGaussPts,source) !source(numFields,numGausspts)
        do n = 1,numFields  
            call integrateCell(tri_id,PhiW(:,:),source(n,:),Res(n,:))    
            spRes(n,:,tri_id) = spRes(n,:,tri_id) - Res(n,:)             !Res(4,totalModes)
        end do  
    end do
    
end subroutine sourceTermRes


!_________________________ ALLOCATE   _____________________________________!

Subroutine allocate_spRes
    Allocate(spRes(numFields,totalModes,numTri))
end subroutine allocate_spRes

Subroutine allocate_spResTRI
    Allocate(spResTri(numFields,totalModes,numTri))
end subroutine allocate_spResTri

Subroutine allocate_Jacobian
    Allocate(jacobDiag(numFields*totalModes,numFields*totalModes,numTri))
    Allocate(jacobOffDiagRL(numFields*totalModes,numFields*totalModes,numInterior))
    Allocate(jacobOffDiagLR(numFields*totalModes,numFields*totalModes,numInterior))
end subroutine allocate_Jacobian

Subroutine allocate_JacobianTri
    Allocate(jacobDiagTri(numFields*totalModes,numFields*totalModes,numTri))
    Allocate(jacobOffDiagRLTri(numFields*totalModes,numFields*totalModes,numInterior))
    Allocate(jacobOffDiagLRTri(numFields*totalModes,numFields*totalModes,numInterior))
end subroutine allocate_JacobianTri

!_________________________ RES per triangle _______________________________!

subroutine Residual_JacTriangle(solCoeffs,triangleNum)!,perturbTri,epsil,k,tm,edgeN)
    real(dp),   intent(in) :: solCoeffs(:,:,:)!,epsil
    integer(i4),intent(in) :: triangleNum!,k,tm,edgeN
    !logical, intent(inout) :: perturbTri
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP,numSides
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)
    real(dp)    :: flux_ql(numFields,numFields,numEdgeGaussPts)
    real(dp)    :: flux_qr(numFields,numFields,numEdgeGaussPts)
    real(dp)    :: dqbdq(numEdgeGaussPts,numFields,numFields)
    real(dp)    :: JD(4*totalModes,4*totalModes)
    integer(i4) :: intLocEdges(3),boundLocEdges(2),boundCount
    real(dp)    :: solCoeffsTemp(numFields,totalModes,numTri)
    

    solCoeffsTemp = solCoeffs 
    spResTri(:,:,:)          = 0.0_dp
    jacobDiagTri(:,:,:)      = 0.0_dp
    jacobOffDiagLRTri(:,:,:) = 0.0_dp
    jacobOffDiagRLTri(:,:,:) = 0.0_dp
    
    intLocEdges(:)      = 0
    boundLocEdges(:)    = 0
    boundCount          = 0
    
    numSides = 0
    
    do n = 1,3           
        edgeNum = triList(triangleNum)%edgeLocList(n)
        if(edgeList(edgeNum)%e2 .ne. -1)then
            numSides = numsides + 1
            intLocEdges(numSides) = triList(triangleNum)%edgeLocList(n)
        else
            boundCount = boundCount + 1
            boundLocEdges(boundCount) = triList(triangleNum)%edgeLocList(n)
        end if
    end do
    
    
    ! Interior Edge Flux
    do edge_id = 1, numSides
        edgeNum  =  intLocEdges(edge_id)
        
        leftTri  = edgeList(edgeNum)%e1    
        rightTri = edgeList(edgeNum)%e2

        face1 = edgeList(edgeNum)%locEdge1
        face2 = edgeList(edgeNum)%locEdge2

        loc1 = 2*face1  - 1   !correct edge and orientation (positive orient)
        loc2 = 2*face2        !correct edge and orientation (negative orient)

        solCoeffsTemp = solCoeffs 
!        if(perturbTri .and. (edge_id == edgeN))then
!            if(leftTri == triangleNum)then
!                solCoeffsTemp(k,tm,rightTri) = solCoeffsTemp(k,tm,rightTri) + epsil
!            else
!                solCoeffsTemp(k,tm,leftTri)  = solCoeffsTemp(k,tm,leftTri) + epsil
!            end if
!        end if
        
        
        call projectEdge(solCoeffsTemp(:,:,leftTri),loc1,qL)
        call projectEdge(solCoeffsTemp(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFlux_q(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux,flux_ql,flux_qr) 
        do n =1,numFields     
            if(leftTri == triangleNum)then
                call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                spResTri(n,:,leftTri)  = spResTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
            else if(rightTri == triangleNum)then
                call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                spResTri(n,:,rightTri) = spResTri(n,:,rightTri) - bR(:)
            end if
        end do

        !Diagonal Jacobian
        if(leftTri == triangleNum)then
            Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
            jacobDiagTri(:,:,leftTri) = jacobDiagTri(:,:,leftTri)     + JD
        else if(rightTri == triangleNum)then
            Call integrateEdgeJacDiag(edgePhiW(:,:,loc2),edgePhi(:,:,loc2),flux_qr,JD)
            jacobDiagTri(:,:,rightTri) = jacobDiagTri(:,:,rightTri)   - JD
        end if
        
        !Off Diagonal Jacobian
        if(leftTri == triangleNum)then
            Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc2),flux_qr,JD)
            jacobOffDiagLRTri(:,:,edge_id) = jacobOffDiagLRTri(:,:,edge_id) + JD
        else if(rightTri == triangleNum)then
            Call integrateEdgeJacDiag(edgePhiW(:,:,loc2),edgePhi(:,:,loc1),flux_ql,JD)
            jacobOffDiagRLTri(:,:,edge_id) = jacobOffDiagRLTri(:,:,edge_id) - JD
        end if
    end do

    ! Boundary Edge Fluxes
    triList(:)%used = .false.
    do edge_id = 1, boundCount
        edgeNum = boundLocEdges(edge_id)
        bc_type = bcFlag(edgeNum)
        leftTri = edgeList(edgeNum)%e1

        face1 = edgeList(edgeNum)%locEdge1
        loc1 = 2*face1 -1     !get correct end and orientation (positive orient)

        call projectEdge(solCoeffs(:,:,leftTri),loc1,qL)

        select case(bc_type)
            case (12)  
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
                    call getFlux_q(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux,flux_ql,flux_qr)
                    do GP = 1,numEdgeGaussPts
                        flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                    end do
                    do n =1,numFields 
                        if(leftTri == triangleNum)then
                            call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                            spResTri(n,:,leftTri)  = spResTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                        else if(rightTri == triangleNum)then
                            call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                            spResTri(n,:,rightTri) = spResTri(n,:,rightTri) - bR(:)
                        end if
                    end do

                    !Diagonal Jacobian
                    if(leftTri == triangleNum)then
                        Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                        jacobDiagTri(:,:,leftTri) = jacobDiagTri(:,:,leftTri) + JD
                    end if
                end if
            case(14)    !inflow
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                end do
                do n =1,numFields   
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResTri(n,:,leftTri) = spResTri(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                    end if
                end do

                !Diagonal Jacobian
                if(leftTri == triangleNum)then
                    Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                    jacobDiagTri(:,:,leftTri) = jacobDiagTri(:,:,leftTri) + JD
                end if
            case(15)    !outflow
                Call specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                end do
                
                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResTri(n,:,leftTri) = spResTri(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                    end if
                end do

                !Diagonal Jacobian
                if(leftTri == triangleNum)then
                    Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                    jacobDiagTri(:,:,leftTri) = jacobDiagTri(:,:,leftTri) + JD
                end if
             case default 
                 
                do GP = 1,numEdgeGaussPts
                    call get_bc_q(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP),dqbdq(GP,:,:))             
                end do
                
                Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr)
                
                do GP = 1,numEdgeGaussPts
                    flux_ql(:,:,GP) = matmul(flux_qr(:,:,GP),dqbdq(GP,:,:))
                end do

                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResTri(n,:,leftTri) = spResTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                    end if
                end do
                
                !Diagonal Jacobian
                if(leftTri == triangleNum)then
                    Call integrateEdgeJacDiag(edgePhiW(:,:,loc1),edgePhi(:,:,loc1),flux_ql,JD)
                    jacobDiagTri(:,:,leftTri) = jacobDiagTri(:,:,leftTri) + JD
                end if
        end select
    end do

    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxRes_JacTriangle(solCoeffs,numGaussPts,triangleNum)
    end if
end subroutine Residual_JacTriangle


Subroutine integralFluxRes_JacTriangle(solCoef,numGaussPts,triangleNum)
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    integer(i4),intent(in) :: numGaussPts,triangleNum
    
    real(dp) :: F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)
    real(dp) :: flux_q(numFields,2,numFields)
    real(dp) :: flux_qXi(numFields,numFields,numGaussPts),flux_qEta(numFields,numFields,numGaussPts)
    real(dp) :: JD(4*totalModes,4*totalModes)
    integer(i4) :: n,GP

    ! dPhiW(totalModes,numGaussPts,2)

    !This is correct (3/8/13)
    dXidx(:)    = JInv(1,1,:,triangleNum) ! (numGPs)
    dXidy(:)    = JInv(2,1,:,triangleNum) 
    dEtadx(:)   = JInv(1,2,:,triangleNum) 
    dEtady(:)   = JInv(2,2,:,triangleNum)

    call projectCell(solCoef(:,:,triangleNum),Phi,q) !q(numFields,numGaussPts)

    do GP = 1,numGaussPts
        call getNativeFlux_q(q(:,GP),flux,flux_q) !-- flux(4,2)
        F_xi(:,GP)  = flux(:,1)*dXidx(GP)  + flux(:,2)*dXidy(GP)
        F_eta(:,GP) = flux(:,1)*dEtadx(GP) + flux(:,2)*dEtady(GP)
        flux_qXi(:,:,GP)  = flux_q(:,1,:)*dXidx(GP)   + flux_q(:,2,:)*dXidy(GP)
        flux_qEta(:,:,GP) = flux_q(:,1,:)*dEtadx(GP)  + flux_q(:,2,:)*dEtady(GP)
    end do

    do n = 1,numFields
        call integrateCell(triangleNum,dPhiW(:,:,1),f_xi(n,:), b1)
        call integrateCell(triangleNum,dPhiW(:,:,2),f_eta(n,:),c1)   
        spResTri(n,:,triangleNum) = spResTri(n,:,triangleNum) - b1(:)  -  c1(:)         !Res(4,totalModes)
    end do 
    

    Call integrateVolJacDiag(triangleNum,dPhiW(:,:,1),Phi, flux_qXi,JD) 
    jacobDiagTri(:,:,triangleNum) = jacobDiagTri(:,:,triangleNum) - JD
    Call integrateVolJacDiag(triangleNum,dPhiW(:,:,2),Phi,flux_qEta,JD)
    jacobDiagTri(:,:,triangleNum) = jacobDiagTri(:,:,triangleNum) - JD

    
end subroutine integralFluxRes_JacTriangle

end module residual_Jac_module