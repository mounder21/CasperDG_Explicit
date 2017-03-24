Module residualComplex_module
    use my_kinddefs
    use projectionComplex_module
    use flux_module
    use Globals_module, only: edgeList,numEdges,interiorEdges,boundEdges,numInterior,numTri,bcFlag
    use inputs_module,  only: totalModes,basisDegree,sourceTerm,fluxScheme 
    use initializeBasis_module, only: PhiW,edgePhiW,dPhiW,Phi,edgePhi,numEdgeGaussPts,numGaussPts
    use integrateComplex_module
    use bc_module
    use norm_module
    use computeJ_detJ_module, only: JInv,detJ
    use getIC_module
    use computeNormal_module, only: normals
    implicit none
#ifdef complx
    complex(dp),pointer :: spResRC(:,:,:)
    complex(dp),pointer :: spResRCTri(:,:,:)
#else
    real(dp),pointer :: spResRC(:,:,:)
    real(dp),pointer :: spResRCTri(:,:,:)
#endif
    
contains


subroutine ResidualComplex(solCoeffs)
#ifdef complx
    complex(dp),   intent(in) :: solCoeffs(:,:,:)
    complex(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    complex(dp)    :: bL(totalModes),bR(totalModes)
    complex(dp)    :: flux(numFields,numEdgeGaussPts)
#else
    real(dp),   intent(in) :: solCoeffs(:,:,:)
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)
#endif    
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP
          
    spResRC(:,:,:)                = 0.0_dp

    ! Interior Edge Fluxes

    do edge_id = 1, numInterior
        edgeNum = interiorEdges(edge_id)
        
        leftTri  = edgeList(edgeNum)%e1    
        rightTri = edgeList(edgeNum)%e2
        
        face1 = edgeList(edgeNum)%locEdge1
        face2 = edgeList(edgeNum)%locEdge2

        loc1 = 2*face1  - 1   !get correct edge and orientation (positive orient)
        loc2 = 2*face2       !get correct edge and orientation (negative orient)
        
        call projectEdgeComplex(solCoeffs(:,:,leftTri),loc1,qL)
        call projectEdgeComplex(solCoeffs(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFluxComplex(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux) !flux(4,#EdgeGPs)

        do n = 1,numFields     
            call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
            call integrateEdgeComplex(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
            spResRC(n,:,leftTri)  = spResRC(n,:,leftTri)  + bL(:) !(4,totalModes,numTri)
            spResRC(n,:,rightTri) = spResRC(n,:,rightTri) - bR(:)
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

        call projectEdgeComplex(solCoeffs(:,:,leftTri),loc1,qL)
        
        !form qR
        select case(bc_type)
            case (12)  
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
                    call getFluxComplex(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux)
                    do n =1,numFields     
                        call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                        call integrateEdgeComplex(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                        spResRC(n,:,leftTri)  =  spResRC(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                        spResRC(n,:,rightTri) =  spResRC(n,:,rightTri) - bR(:)
                    end do
                end if
            case(14)    !inflow
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n =1,numFields     
                    call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                    spResRC(n,:,leftTri) = spResRC(n,:,leftTri)  + bL(:)             !(4,totalModes,numTri)
                end do

            case(15)    !outflow
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs) 
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
    
                do n = 1,numFields     
                    call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)      !resL , bL(totalModes)
                    spResRC(n,:,leftTri) = spResRC(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
                
            case default 
                do GP = 1,numEdgeGaussPts
                    call get_bcComplex(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP)) 
                end do
                
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)

                do n = 1,numFields     
                    call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)      !resL , bL(totalModes)
                    spResRC(n,:,leftTri) = spResRC(n,:,leftTri) + bL(:)             !(4,totalModes,numTri)
                end do
        end select
    end do
    
    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxResComplex(solCoeffs)
    end if

end subroutine ResidualComplex


Subroutine integralFluxResComplex(solCoef)
#ifdef complx
    complex(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    complex(dp) :: F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes) 
    complex(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
#else
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    real(dp) :: F_xi(numFields,numGaussPts),b1(totalModes),c1(totalModes)
    real(dp) :: q(numFields,numGaussPts),F_eta(numFields,numGaussPts),flux(numFields,2)
#endif

    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)
    integer(i4) :: n,tri_id,GP

    ! dPhiW(totalModes,numGaussPts,2)

    do tri_id = 1,numTri
        !This is correct (3/8/13)
        dXidx(:)    = JInv(1,1,:,tri_id) ! (numGPs)
        dXidy(:)    = JInv(2,1,:,tri_id) 
        dEtadx(:)   = JInv(1,2,:,tri_id) 
        dEtady(:)   = JInv(2,2,:,tri_id)
        
        call projectCellComplex(solCoef(:,:,tri_id),Phi,q) !q(numFields,numGaussPts)
        
        do n = 1,numFields  
            do GP = 1,numGaussPts
                call getNativeFluxComplex(q(:,GP),flux) !-- flux(4,2)
                F_xi(n,GP)  = flux(n,1)*dXidx(GP)  + flux(n,2)*dXidy(GP)
                F_eta(n,GP) = flux(n,1)*dEtadx(GP) + flux(n,2)*dEtady(GP)
            end do
            
            call integrateCellComplex(tri_id,dPhiW(:,:,1),F_xi(n,:), b1)
            call integrateCellComplex(tri_id,dPhiW(:,:,2),F_eta(n,:),c1)   

           spResRC(n,:,tri_id) = spResRC(n,:,tri_id) - b1(:)  -  c1(:)          
        end do
        
    end do
    
end subroutine integralFluxResComplex
    

!_________________________ Allocate spResComplex _______________________________!


Subroutine allocate_spResRC
    Allocate(spResRC(numFields,totalModes,numTri))
end subroutine allocate_spResRC

Subroutine allocate_spResRCTri
    Allocate(spResRCTri(numFields,totalModes,numTri))
end subroutine allocate_spResRCTri


!_________________________ RES per triangle _______________________________!


subroutine ResidualComplexTriangle(solCoeffs,triangleNum,perturbTri,epsil,k,tm,edgeN)
#ifdef complx    
    complex(dp),   intent(in) :: solCoeffs(:,:,:),epsil
    complex(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    complex(dp)    :: bL(totalModes),bR(totalModes)
    complex(dp)    :: flux(numFields,numEdgeGaussPts)
    complex(dp)    :: solCoeffsTemp(numFields,totalModes,numTri)
#else
    real(dp),   intent(in) :: solCoeffs(:,:,:),epsil
    real(dp)    :: qL(numFields,numEdgeGaussPts),qR(numFields,numEdgeGaussPts)
    real(dp)    :: bL(totalModes),bR(totalModes)
    real(dp)    :: flux(numFields,numEdgeGaussPts)
    real(dp)    :: solCoeffsTemp(numFields,totalModes,numTri)
#endif
  
    integer(i4),intent(in) :: triangleNum,k,tm,edgeN
    logical, intent(inout) :: perturbTri
    
    integer(i4) :: edge_id,leftTri,rightTri,loc1,loc2,edgeNum,n,bc_type,face1,face2,GP,numSides
    integer(i4) :: intLocEdges(3),boundLocEdges(2),boundCount

    
    solCoeffsTemp = solCoeffs 
    spResRCTri(:,:,:)          = 0.0_dp
    
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
        
        
        call projectEdgeComplex(solCoeffsTemp(:,:,leftTri),loc1,qL)
        call projectEdgeComplex(solCoeffsTemp(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        call getFluxComplex(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux) !flux(4,#EdgeGPs)
        do n =1,numFields     
            if(leftTri == triangleNum)then
                call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                spResRCTri(n,:,leftTri)  = spResRCTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
            else if(rightTri == triangleNum)then
                call integrateEdgeComplex(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                spResRCTri(n,:,rightTri) = spResRCTri(n,:,rightTri) - bR(:)
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

        call projectEdgeComplex(solCoeffs(:,:,leftTri),loc1,qL)

        select case(bc_type)
            case (12)  
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
            case (13)   !Periodic BCs
                if(triList(leftTri)%used .eqv. .false.)then
                    Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
                    call getFluxComplex(qL,qR,normals(:,:,edgeNum),numEdgeGaussPts,fluxScheme,flux)
  
                    do n =1,numFields 
                        if(leftTri == triangleNum)then
                            call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)     !resL , bL(totalModes)
                            spResRCTri(n,:,leftTri)  = spResRCTri(n,:,leftTri)  + bL(:) !(4,totalModes,numTri)
                        else if(rightTri == triangleNum)then
                            call integrateEdgeComplex(edgePhiW(:,:,loc2),flux(n,:),bR)     !resR
                            spResRCTri(n,:,rightTri) = spResRCTri(n,:,rightTri) - bR(:)
                        end if
                    end do
                end if
            case(14)    !inflow
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n =1,numFields   
                    if(leftTri == triangleNum)then
                        call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResRCTri(n,:,leftTri) = spResRCTri(n,:,leftTri)  + bL(:)  !(4,totalModes,numTri)
                    end if
                end do

            case(15)    !outflow
                Call specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffs)
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResRCTri(n,:,leftTri) = spResRCTri(n,:,leftTri)  + bL(:)  !(4,totalModes,numTri)
                    end if
                end do

            case default 
                do GP = 1,numEdgeGaussPts
                    call get_bcComplex(bc_type, normals(GP,:,edgeNum), qL(:,GP), qR(:,GP))   
                end do
                
                Call getBCFluxComplex(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)
                
                do n = 1,numFields     
                    if(leftTri == triangleNum)then
                        call integrateEdgeComplex(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
                        spResRCTri(n,:,leftTri) = spResRCTri(n,:,leftTri)  + bL(:)       !(4,totalModes,numTri)
                    end if
                end do
        end select
    end do

    ! Integral Fluxes
    if (basisDegree .ne. 0) then
        Call integralFluxResComplexTriangle(solCoeffs,numGaussPts,triangleNum)
    end if
end subroutine ResidualComplexTriangle


Subroutine integralFluxResComplexTriangle(solCoef,numGaussPts,triangleNum)
#ifdef complx
    complex(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    complex(dp) :: q(numFields,numGaussPts),b1(totalModes),c1(totalModes),flux(numFields,2)
    complex(dp) :: F_xi(numFields,numGaussPts),F_eta(numFields,numGaussPts)
#else
    real(dp),   intent(in) :: solCoef(:,:,:)        !(numFields,tm,numTri)
    real(dp) :: q(numFields,numGaussPts),b1(totalModes),c1(totalModes),flux(numFields,2)
    real(dp) :: F_xi(numFields,numGaussPts),F_eta(numFields,numGaussPts)
#endif
    
    integer(i4),intent(in) :: numGaussPts,triangleNum 
    real(dp) :: dXidx(numGaussPts),dXidy(numGaussPts),dEtadx(numGaussPts),dEtady(numGaussPts)
    integer(i4) :: n,GP

    ! dPhiW(totalModes,numGaussPts,2)

    !This is correct (3/8/13)
    dXidx(:)    = JInv(1,1,:,triangleNum) ! (numGPs)
    dXidy(:)    = JInv(2,1,:,triangleNum) 
    dEtadx(:)   = JInv(1,2,:,triangleNum) 
    dEtady(:)   = JInv(2,2,:,triangleNum)

    call projectCellComplex(solCoef(:,:,triangleNum),Phi,q) !q(numFields,numGaussPts)

    do GP = 1,numGaussPts
        call getNativeFluxComplex(q(:,GP),flux) !-- flux(4,2)
        F_xi(:,GP)  = flux(:,1)*dXidx(GP)  + flux(:,2)*dXidy(GP)
        F_eta(:,GP) = flux(:,1)*dEtadx(GP) + flux(:,2)*dEtady(GP)
    end do

    do n = 1,numFields
        call integrateCellComplex(triangleNum,dPhiW(:,:,1),f_xi(n,:), b1)
        call integrateCellComplex(triangleNum,dPhiW(:,:,2),f_eta(n,:),c1)   
        spResRCTri(n,:,triangleNum) = spResRCTri(n,:,triangleNum) - b1(:)  -  c1(:)         !Res(4,totalModes)
    end do 
    
end subroutine integralFluxResComplexTriangle

end module
