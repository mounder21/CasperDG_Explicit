Module residual_module
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
    
    real(dp),pointer :: spacialResidual(:,:,:)
    real(dp),pointer :: spacialResidualTri(:,:,:)
contains

subroutine Residual(solCoeffs)
    real(dp),   intent(in) :: solCoeffs(:,:,:)
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)
  
    
    
    spacialResidual(:,:,:) = 0.0_dp

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

        call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux) !flux(4,#EdgeGPs)

        do n = 1,numFields     
            call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
            call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
            spacialResidual(n,:,leftTri)  = spacialResidual(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
            spacialResidual(n,:,rightTri) = spacialResidual(n,:,rightTri) - bR(:)
        end do

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
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                    call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux)
                    do n =1,numFields     
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                        call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                        spacialResidual(n,:,leftTri)  = spacialResidual(n,:,leftTri)  + bL(:) !(4,totalModes,numTri)
                        spacialResidual(n,:,rightTri) = spacialResidual(n,:,rightTri) - bR(:)
                    end do

                end if
            case(14)    !inflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
  
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n =1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spacialResidual(n,:,leftTri) = spacialResidual(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                end do
                
            case(15)    !outflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)

                do n = 1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spacialResidual(n,:,leftTri) = spacialResidual(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
                
            case default 
                do GP = 1,numEdgeGaussPts
                    call get_bc(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP))
                end do
                
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n = 1,numFields     
                    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spacialResidual(n,:,leftTri) = spacialResidual(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
                
        end select
    end do
    
    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxRes(solCoeffs)
    end if
    if (sourceTerm) then
        Call sourceTermRes()
    end if
end subroutine Residual

Subroutine integralFluxRes(solCoef)
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    
    real(dp) :: F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)

    integer(i4) :: n,tri_id,GP

    ! dPhiW(totalModes,numGaussPts,2)

    do tri_id = 1,numTri
        !This is correct (3/8/13)
        dXidx(:)    = JInv(1,1,:,tri_id) ! (numGPs)
        dXidy(:)    = JInv(2,1,:,tri_id) 
        dEtadx(:)   = JInv(1,2,:,tri_id) 
        dEtady(:)   = JInv(2,2,:,tri_id)
        
        call projectCell(solCoef(:,:,tri_id),Phi,q) !q(numFields,numGaussPts)
        
        do n = 1,numFields  
            do GP = 1,numGaussPts
                call getNativeFlux(q(:,GP),flux) !-- flux(4,2)
                F_xi(n,GP)  = flux(n,1)*dXidx(GP)  + flux(n,2)*dXidy(GP)
                F_eta(n,GP) = flux(n,1)*dEtadx(GP) + flux(n,2)*dEtady(GP)
            end do
            
            call integrateCell(tri_id,dPhiW(:,:,1),F_xi(n,:), b1)
            call integrateCell(tri_id,dPhiW(:,:,2),F_eta(n,:),c1)   

            spacialResidual(n,:,tri_id) = spacialResidual(n,:,tri_id) - b1(:)  -  c1(:)       
        end do
        
    end do
    
end subroutine integralFluxRes

Subroutine sourceTermRes()
    use getSourceTerm_module
    
    real(dp) :: Res(numFields,totalModes)
    real(dp) :: source(numFields,numGaussPts)
    integer(i4) :: tri_id, n
    
    do tri_id = 1,numTri  
        Call getSourceTerm(tri_id,numGaussPts,source) !source(numFields,numGausspts)
        do n = 1,numFields  
            call integrateCell(tri_id,PhiW(:,:),source(n,:),Res(n,:))    
            spacialResidual(n,:,tri_id) = spacialResidual(n,:,tri_id) + Res(n,:)             !Res(4,totalModes)
        end do  
    end do
    
end subroutine sourceTermRes


!_________________________ ALLOCATE   _____________________________________!

Subroutine allocate_spacialResidual
    use Globals_module, only: numEdges 
    Allocate(spacialResidual(numFields,totalModes,numTri))
end subroutine allocate_spacialResidual

Subroutine allocate_spacialResidualTRI
    use Globals_module, only: numEdges 
    Allocate(spacialResidualTri(numFields,totalModes,numTri))
end subroutine allocate_spacialResidualTri


!_________________________ RES per triangle _______________________________!

subroutine Residual_Triangle(solCoeffs,triangleNum,perturbTri,epsil,k,tm,edgeN)
    real(dp),   intent(in) :: solCoeffs(:,:,:),epsil
    integer(i4),intent(in) :: triangleNum,k,tm,edgeN
    logical, intent(inout) :: perturbTri
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP,numSides
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)

    integer(i4) :: intLocEdges(3),boundLocEdges(2),boundCount
    real(dp)    :: solCoeffsTemp(numFields,totalModes,numTri)
    

    solCoeffsTemp = solCoeffs 
    spacialResidualTri(:,:,:)          = 0.0_dp

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
        if(perturbTri .and. (edge_id == edgeN))then
            if(leftTri == triangleNum)then
                solCoeffsTemp(k,tm,rightTri) = solCoeffsTemp(k,tm,rightTri) + epsil
            else
                solCoeffsTemp(k,tm,leftTri)  = solCoeffsTemp(k,tm,leftTri) + epsil
            end if
        end if
        
        
        call projectEdge(solCoeffsTemp(:,:,leftTri),loc1,qL)
        call projectEdge(solCoeffsTemp(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux) !flux(4,#EdgeGPs)
        do n =1,numFields     
            if(leftTri == triangleNum)then
                call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                spacialResidualTri(n,:,leftTri)  = spacialResidualTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
            else if(rightTri == triangleNum)then
                call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                spacialResidualTri(n,:,rightTri) = spacialResidualTri(n,:,rightTri) - bR(:)
            end if
        end do


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
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                    call getFlux(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux)

                    do n =1,numFields 
                        if(leftTri == triangleNum)then
                            call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                            spacialResidualTri(n,:,leftTri)  = spacialResidualTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                        else if(rightTri == triangleNum)then
                            call integrateEdge(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                            spacialResidualTri(n,:,rightTri) = spacialResidualTri(n,:,rightTri) - bR(:)
                        end if
                    end do

                end if
            case(14)    !inflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)

                do n = 1,numFields   
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spacialResidualTri(n,:,leftTri) = spacialResidualTri(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                    end if
                end do


            case(15)    !outflow
                Call specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)

                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spacialResidualTri(n,:,leftTri) = spacialResidualTri(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                    end if
                end do
             case default 
                do GP = 1,numEdgeGaussPts
                    call get_bc(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP))         
                end do
                
                Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)

                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spacialResidualTri(n,:,leftTri) = spacialResidualTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                    end if
                end do
        end select
    end do

    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxRes_Triangle(solCoeffs,numGaussPts,triangleNum)
    end if
end subroutine Residual_Triangle


Subroutine integralFluxRes_Triangle(solCoef,numGaussPts,triangleNum)
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    integer(i4),intent(in) :: numGaussPts,triangleNum
    
    real(dp) :: F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)

    integer(i4) :: n,GP

    ! dPhiW(totalModes,numGaussPts,2)

    !This is correct (3/8/13)
    dXidx(:)    = JInv(1,1,:,triangleNum) ! (numGPs)
    dXidy(:)    = JInv(2,1,:,triangleNum) 
    dEtadx(:)   = JInv(1,2,:,triangleNum) 
    dEtady(:)   = JInv(2,2,:,triangleNum)

    call projectCell(solCoef(:,:,triangleNum),Phi,q) !q(numFields,numGaussPts)

    do GP = 1,numGaussPts
        call getNativeFlux(q(:,GP),flux) !-- flux(4,2)
        F_xi(:,GP)  = flux(:,1)*dXidx(GP)  + flux(:,2)*dXidy(GP)
        F_eta(:,GP) = flux(:,1)*dEtadx(GP) + flux(:,2)*dEtady(GP)
    end do

    do n = 1,numFields
        call integrateCell(triangleNum,dPhiW(:,:,1),f_xi(n,:), b1)
        call integrateCell(triangleNum,dPhiW(:,:,2),f_eta(n,:),c1)   
        spacialResidualTri(n,:,triangleNum) = spacialResidualTri(n,:,triangleNum) - b1(:) - c1(:)         !Res(4,totalModes)
    end do 
    
end subroutine integralFluxRes_Triangle

end module residual_module
