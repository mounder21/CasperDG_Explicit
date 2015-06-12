
Program casperDG
    use All_modules
    implicit none 

    integer(i4)       :: timeSteps,readNumModes
    logical           :: restart,gmesh
    character(20)     :: timeScheme
    
    ! ------------------------------------- Pre-process Variables ----------------------------------!
    ! ----------- MMS Variables ------------ !
    a = 0.01_dp ; b = 0.0_dp ; s0 = 1.3_dp                 !coeffs for u = ax + by + s0
    c = 0.01_dp ; d = 0.0_dp ; t0 = 1.3_dp                 !coeffs for v = cx + dy + t0
    
    ! ----------- Initial Variables ------------ !
    rho_in  = 1._dp
    P_in    = 1.0_dp / gamma

    fmach   = 0.8_dp                                  ! Free-stream Mach number
    alpha   = 3._dp                                   ! angle of attack
    
    ! ----------- Profile Initializations ------------!200,392,578,800,968,1250,1352,1568,1800
    mesh_name   = '../meshes/naca.msh'              ! Mesh file location/name box200_aba
    gmesh	= .false.
    abaFile     = .false.                              !does the mesh_name have _aba.msh?
    restart     = .false.
    writeCoeffs = .false.
    
    sourceTerm  = .false.
    MMS_sin     = .false.
    periodicBCs = .false.
    
    ! ----------- Restart Options ------------ !                             
    solution_name = '../vtu_files/testp2_Coef_0039.txt' 
    readNumModes  = 3       !1, 3 or 6
    
    ! ----------- Plotting Options ------------ !
    numPlotPts  = 6 
    out_file    = '../vtu_files/naca'   ! Base name for output file (###.vtu)
    out_freq    = 100                   ! Output file frequency, solution visualizations
    
    ! ----------- Time Stepping Options ------------ !
    CFL         = 1.0_dp
    timeSteps   = 100000
    timeScheme  = 'RK4'

    ! ----------- Resolution Options ------------ !
    basisType       = 'H1'                          ! basis type: H1, nodal, monomial
    mapBasisType    = 'H1'                          ! basis type: H1, nodal
    basisDegree     = 1                             ! Basis Function Polynomial Order(0,1,2)--> solution order     
    mapBasisDegreeCurved = basisDegree + 1          ! map basis polynomial order curved elements
    
    totalModes              = (basisDegree  + 1)*(basisDegree  + 2)/2 
    mapTotalModesCurved     = (mapBasisDegreeCurved + 1)*(mapBasisDegreeCurved + 2)/2
    
    
    ! ------------------------------------- Pre-processing ------------------------------------------!
    if(gmesh)then
    	call read_gmesh 
    else
    	call read_mesh
    end if
    call constructGridData             ! Generate the edge list (edgeList(:)%n1, egdeList(:)%n1) 
    !--------- find delta_h----------!
    call get_h
    
    call initializeBasis(numGaussPts,numEdgeGaussPts,numPlotPts)
    
    
    ! -------------------------------------- Assemble Normals -----------------------------------------!
    call computeJ_detJ(map_dPhiStr8,numGaussPts)
    !call computeJ_detJ(dmapPhiCurved,numGaussPts)
    
    call computeNormals(edgeMap_dPhiStr8,numEdgeGaussPts)      !normals(ngp,x/y,numEdges)
    !call computeNormals(edgeMap_dPhiCurved,numEdgeGaussPts)   !normals(ngp,x/y,numEdges)

    call assembleLocMassMatrix(Phi,PhiW)
    call alloc
    
    
    !---------------------------------------- Initialize Solution -------------------------------------!
    call initializeSolution(PhiW,numGaussPts,invMassMatrices,restart,readNumModes)
    
    
    ! ------------------_----------- Plot Initial Solution and Residual -------------------------------!
    if(numPlotPts == 3)then
        call plotSolution(out_file,solCoeffs)
    else if (numPlotPts == 6)then
        call plotSolution6pt(out_file,solCoeffs)
    end if
                
    
    ! -------------------------------------- Time step solution ---------------------------------------!  
    Call evolveSolution(timeScheme,timeSteps)
    
    
End Program casperDG
