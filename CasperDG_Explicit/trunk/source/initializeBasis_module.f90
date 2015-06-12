Module initializeBasis_module
    use my_kinddefs
    implicit none 
    
    real(dp), pointer :: Phi(:,:),dPhi(:,:,:),                      PhiW(:,:),dPhiW(:,:,:)
    real(dp), pointer :: edgePhi(:,:,:),edge_dPhi(:,:,:,:),         edgePhiW(:,:,:),edge_dPhiW(:,:,:,:)
    
    real(dp), pointer :: mapPhiStr8(:,:),map_dPhiStr8(:,:,:)!,       mapPhiWStr8(:,:)
    real(dp), pointer :: mapPhiCurved(:,:),map_dPhiCurved(:,:,:)!,   mapPhiWCurved(:,:)
    
    real(dp), pointer :: edgeMapPhiStr8(:,:,:),edgeMap_dPhiStr8(:,:,:,:)
    real(dp), pointer :: edgeMapPhiCurved(:,:,:),edgeMap_dPhiCurved(:,:,:,:)
    
    real(dp), pointer :: GaussPoints(:,:),GaussWeights(:)
    real(dp), pointer :: edgeGaussPoints(:), edgeGaussWeights(:)
    real(dp), pointer :: solutionPhi(:,:) 
    
    integer(i4)       :: numGaussPts,numEdgeGaussPts
contains
    
Subroutine initializeBasis(numGaussPts,numEdgeGaussPts,numPlotPts)
    use gauss_points 
    use inputs_module, only: basisType,mapBasisType,basisDegree,mapBasisDegreeStr8,totalModes,mapTotalModesStr8
    use inputs_module, only: mapBasisDegreeCurved,mapTotalModesCurved
    use Globals_module
    use basis_module
    implicit none 
    
    integer(i4), intent(in) :: numPlotPts
    integer(i4),intent(out) :: numGaussPts,numEdgeGaussPts

    ! dummy
    integer(i4),dimension(0:20)     :: allGaussPts
    integer(i4),dimension(0:21)     :: allBarGaussPts

    integer(i4):: mode, GP
    real(dp) :: triPts(2,numPlotPts)

   !------------------------------------- Set up Gauss Points/Weights ----------------------------------!    
    Call setup_p2nqptsArea(allGaussPts)
    Call setup_p2nqptsBar(allBarGaussPts)
        
    numGaussPts = allGaussPts(2*basisDegree) 
    numEdgeGaussPts = allBarGaussPts(2*basisDegree + 1)

    
    !------------------------------------- Allocate Basis Pointers -------------------------------------! 
    Call allocate_basis
    
    !------------------------------------- Get GuassPts and Weights ------------------------------------!
    Call tri_points(GaussPoints,GaussWeights,numGaussPts)
    Call bar_points(edgeGaussPoints,edgeGaussWeights,numEdgeGaussPts)
    
    !---------------------------------------- Set up Phi arrays ----------------------------------------! 
    Call getPhi(basisType,basisDegree,totalModes,numGaussPts,GaussPoints,Phi)
    Call getdPhi(basisType,basisDegree,totalModes,numGaussPts,GaussPoints,dPhi)
    Call getEdgePhi(basisType,basisDegree,totalModes,numEdgeGaussPts,edgeGaussPoints,edgePhi)

    Call getEdge_dPhi(basisType,basisDegree,totalModes,numEdgeGaussPts,edgeGaussPoints,edge_dPhi)

    
    !Straight elements
    Call getPhi (mapBasisType,mapBasisDegreeStr8,mapTotalModesStr8,numGaussPts,GaussPoints,mapPhiStr8)
    Call getdPhi(mapBasisType,mapBasisDegreeStr8,mapTotalModesStr8,numGaussPts,GaussPoints,map_dPhiStr8)
    Call getEdgePhi(mapBasisType,mapBasisDegreeStr8,mapTotalModesStr8,numEdgeGaussPts,edgeGaussPoints,&
                    edgeMapPhiStr8)
    Call getEdge_dPhi(mapBasisType,mapBasisDegreeStr8,mapTotalModesStr8,&
                    numEdgeGaussPts,edgeGaussPoints,edgeMap_dPhiStr8)  
                    
    !Curved Elements
    !Call getPhi (mapBasisType,mapBasisDegreeCurved,mapTotalModesCurved,numGaussPts,GaussPoints,mapPhiCurved)
    !Call getdPhi(mapBasisType,mapBasisDegreeCurved,mapTotalModesCurved,numGaussPts,GaussPoints,map_dPhiCurved)
    !Call getEdgePhi(mapBasisType,mapBasisDegreeCurved,mapTotalModesCurved,numEdgeGaussPts,edgeGaussPoints,&
    !                edgeMapPhiCurved)
    !Call getEdge_dPhi(mapBasisType,mapBasisDegreeCurved,mapTotalModesCurved,&
    !                numEdgeGaussPts,edgeGaussPoints,edgeMap_dPhiCurved)     
                    
  !--------------------------------- Set up PhiW,edgePhiW matrices ---------------------------------! 
                   
    do GP = 1, numGaussPts
        do mode = 1,totalModes
            PhiW(mode,GP) = Phi(mode,GP) * GaussWeights(GP)
            dPhiW(mode,GP,:) = dPhi(mode,GP,:) * GaussWeights(GP)
        end do 
    end do
     
    do GP = 1, numEdgeGaussPts
        do mode = 1,totalModes
            edgePhiW(mode,GP,:) = edgePhi(mode,GP,:) * edgeGaussWeights(GP)
            edge_dPhiW(mode,GP,:,:) = edge_dPhi(mode,GP,:,:) * edgeGaussWeights(GP)
        end do 
    end do

    if (numPlotPts == 3)then
        triPts(1,1) = -1._dp
        triPts(2,1) = -1._dp
        triPts(1,2) =  1._dp
        triPts(2,2) = -1._dp
        triPts(1,3) = -1._dp
        triPts(2,3) =  1._dp
    else if(numPlotPts == 6)then       
        triPts(1,1) = -1._dp
        triPts(2,1) = -1._dp
        triPts(1,2) =  1._dp
        triPts(2,2) = -1._dp
        triPts(1,3) = -1._dp
        triPts(2,3) =  1._dp
        triPts(1,4) =  0._dp
        triPts(2,4) = -1._dp
        triPts(1,5) =  0._dp
        triPts(2,5) =  0._dp  
        triPts(1,6) = -1._dp
        triPts(2,6) =  0._dp
    end if
    Call getPhi(basisType,basisDegree,totalModes,numPlotPts,triPts,solutionPhi)
   
end subroutine initializeBasis

Subroutine allocate_basis
    use inputs_module, only: totalModes,mapTotalModesStr8,mapTotalModesCurved,numPlotPts
 
    
    implicit none
    
    Allocate(Phi(totalModes,numGaussPts))
    Allocate(dPhi(totalModes,numGaussPts,2))
    Allocate(edgePhi(totalModes,numEdgeGaussPts,6))
    Allocate(edge_dPhi(totalModes,numEdgeGaussPts,2,6))
    Allocate(solutionPhi(totalModes,numPlotPts))
    
    !Straight Structures
    Allocate(mapPhiStr8(mapTotalModesStr8,numGaussPts))
    Allocate(map_dPhiStr8(mapTotalModesStr8,numGaussPts,2))
    Allocate(edgeMapPhiStr8(mapTotalModesStr8,numEdgeGaussPts,6))
    Allocate(edgeMap_dPhiStr8(mapTotalModesStr8,numEdgeGaussPts,2,6))
    
    !Curved Structures
    !Allocate(mapPhiCurved(mapTotalModesCurved,numGaussPts))
    !Allocate(map_dPhiCurved(mapTotalModesCurved,numGaussPts))
    !Allocate(edgeMapPhiCurved(mapTotalModesCurved,numEdgeGaussPts))
    !Allocate(edgeMap_dPhiCurved(mapTotalModesCurved,numEdgeGaussPts))
    
    !Weighted Phi
    Allocate(PhiW(totalModes,numGaussPts))
    Allocate(dPhiW(totalModes,numGaussPts,2)) 
    Allocate(edgePhiW(totalModes,numEdgeGaussPts,6))
    Allocate(edge_dPhiW(totalModes,numEdgeGaussPts,2,6))
    
    Allocate(GaussPoints(2,numGaussPts))
    Allocate(GaussWeights(numGaussPts))
    Allocate(edgeGaussPoints(numEdgeGaussPts))
    Allocate(edgeGaussWeights(numEdgeGaussPts))
    
end subroutine allocate_basis
end module initializeBasis_module