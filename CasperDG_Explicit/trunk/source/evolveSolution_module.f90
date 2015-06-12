Module evolveSolution_module
    use my_kinddefs
    use inputs_module, only: CFL,out_freq,out_file,numPlotPts,writeCoeffs,MMS_sin
    use Globals_module, only: numTri,numEulerVars,solutionNorm
    use initializeSolution_module, only: solCoeffs,spRes
    use computeNormal_module, only: normals
    use initializeBasis_module, only: numGaussPts,numEdgeGaussPts
    use assembleLocMassMatrices_module,  only: invMassMatrices,MassMatrices
    use get_delta_t_module, only: delta_t,min_delta_t
    use residual_module
    
    implicit none
    
contains

Subroutine evolveSolution(timeScheme,timeSteps)
    use vtu_output_module
    use norm_module
    
    character(*),intent(in) :: timeScheme
    integer(i4),intent(in)  :: timeSteps 
    
    select case(timeScheme)
        case('Forward_Euler')
            Call Forward_Euler(timeSteps)
        case('RK4')
            Call RK4(timeSteps)
    end select
    
end subroutine

Subroutine Forward_Euler(timeSteps)
    use get_delta_t_module
    use vtu_output_module
    use norm_module
    use residual_module
    use solutionNorm_module
    
    integer(i4),intent(in) :: timeSteps
    integer(i4) :: timeStep,elem_id,eulerVar
    real(dp)    :: norm
    
    call Residual(solCoeffs,spRes) !returns spRes
    call plotResidual(out_file,spRes)
    call Lp_norm(spRes(:,:,:),2._dp,norm)
    call write_Solution_Coefficients(out_file,solCoeffs)
    
    do timeStep = 1,timeSteps
        call get_delta_t(solCoeffs)
        !print*,min_delta_t
        print*,timeStep,norm
        
        call Residual(solCoeffs,spRes) !returns spRes
        call Lp_norm(spRes(:,:,:),2._dp,norm)               
            
        if (norm < 10e-12) then
            print*,'Solution reached steady state'
            Call callPlot(solCoeffs,spRes)
            if(MMS_sin)then
                Call Lp_solution_norm(solCoeffs,numGaussPts,2._dp,solutionNorm)
            end if
            stop 
        end if
        
        do elem_id = 1,numTri
            !------ update solution ------!
            do eulerVar = 1, numEulerVars
                solCoeffs(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id) + &
                    CFL*min_delta_t*& 
                    matmul(invMassMatrices(:,:,elem_id),spRes(eulerVar,:,elem_id))  
            end do
        end do
        
        !--------------- Print Solution and Residual ---------------!  
        if (MOD(timeStep,out_freq)==0) then
            Call callPlot(solCoeffs,spRes)
        end if
    end do
    
end subroutine

Subroutine RK4(timeSteps)       !Forth order time accurate
    use get_delta_t_module
    use vtu_output_module
    use norm_module
    use residual_module
    use solutionNorm_module
    
    integer(i4),intent(in) :: timeSteps
    integer(i4) :: timeStep,elem_id,eulerVar
    real(dp)    :: norm
    real(dp)    :: k1(numEulerVars,totalModes,numTri)
    real(dp)    :: k2(numEulerVars,totalModes,numTri)
    real(dp)    :: k3(numEulerVars,totalModes,numTri)
    real(dp)    :: k4(numEulerVars,totalModes,numTri)
    real(dp)    :: f1(numEulerVars,totalModes,numTri)
    real(dp)    :: f2(numEulerVars,totalModes,numTri)
    real(dp)    :: f3(numEulerVars,totalModes,numTri)
    real(dp)    :: f4(numEulerVars,totalModes,numTri)
    
    
    ! ---------------------- Notes --------------------- !
    ! This is the fourth order four-stage method
    
    call Residual(solCoeffs,spRes) !returns spRes
    call plotResidual(out_file,spRes)
    call Lp_norm(spRes(:,:,:),2._dp,norm)
    call write_Solution_Coefficients(out_file,solCoeffs)
    
    do timeStep = 1,timeSteps
        call get_delta_t(solCoeffs)
        !print*,min_delta_t
        print*,timeStep,norm
        
        call Residual(solCoeffs,spRes) !returns spRes
        call Lp_norm(spRes(:,:,:),2._dp,norm)
            
        if (norm < 10e-12) then
            print*,'Solution reached steady state '
            Call callPlot(solCoeffs,spRes)
            if(MMS_sin)then
                Call Lp_solution_norm(solCoeffs,numGaussPts,2._dp,solutionNorm)
            end if
            stop 
        end if
        
         !------ Stage 1 and 2 -,spRes-----!
        do elem_id = 1,numTri
            do eulerVar = 1, numEulerVars
                k1(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id)
                k2(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id) + &
                    0.5_dp*min_delta_t*& 
                    matmul(invMassMatrices(:,:,elem_id),spRes(eulerVar,:,elem_id))
                f1(eulerVar,:,elem_id) = matmul(invMassMatrices(:,:,elem_id),spRes(eulerVar,:,elem_id))
            end do
        end do
        
         !------ Stage 3 ------!
        call Residual(k2,spRes) !returns spRes
        do elem_id = 1,numTri
            do eulerVar = 1, numEulerVars
                f2(eulerVar,:,elem_id) = matmul(invMassMatrices(:,:,elem_id),spRes(eulerVar,:,elem_id))
                k3(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id) + &
                                         0.5_dp*min_delta_t*f2(eulerVar,:,elem_id)
            end do
        end do
        
         !------ Stage 4 ------!
        call Residual(k3,spRes) !returns spRes
        do elem_id = 1,numTri
            do eulerVar = 1, numEulerVars
                f3(eulerVar,:,elem_id) = matmul(invMassMatrices(:,:,elem_id),spRes(eulerVar,:,elem_id))
                k4(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id) + &
                                         min_delta_t*f3(eulerVar,:,elem_id)                                       
            end do
        end do
        
         !------ update solution ------!
        do elem_id = 1,numTri
            do eulerVar = 1, numEulerVars
                solCoeffs(eulerVar,:,elem_id) = solCoeffs(eulerVar,:,elem_id) + &
                                         (1._dp/6._dp) * CFL * min_delta_t * &
                                         (f1(eulerVar,:,elem_id)&
                                          + 2._dp*f2(eulerVar,:,elem_id) &
                                          + 2._dp*f3(eulerVar,:,elem_id) &
                                          + f4(eulerVar,:,elem_id))
            end do
        end do
        !--------------- Print Solution and Residual ---------------!  
        if (MOD(timeStep,out_freq)==0) then
            Call callPlot(solCoeffs,spRes)
        end if
        
    end do
    
end subroutine

Subroutine callPlot(solCo,spacialResidual)
    use vtu_output_module
    real(dp),intent(in):: solCo(:,:,:),spacialResidual(:,:,:)

    if(numPlotPts == 3)then
        call plotSolution(out_file,solCo)
    else if (numPlotPts == 6)then
        call plotSolution6pt(out_file,solCo)
    end if
    
    call plotResidual(out_file,spacialResidual)
    
    if(writeCoeffs)then
        call write_Solution_Coefficients(out_file,solCo)
    end if      
end subroutine

end module evolveSolution_module
    
