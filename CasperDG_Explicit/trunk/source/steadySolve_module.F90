Module steadySolve_module
    use my_kinddefs
    use inputs_module, only: out_freq,out_file,numPlotPts,writeCoeffs,MMS_sin
    use Globals_module, only: numTri,numFields,solutionNorm,triList
    use initializeSolution_module, only: solCoeffs
    use computeNormal_module, only: normals
    use initializeBasis_module, only: numGaussPts,numEdgeGaussPts
    use assembleLocMassMatrices_module,  only: invMassMatrices,MassMatrices
    use residual_Jac_module, only: spRes,jacobDiag,jacobOffDiagRL,jacobOffDiagLR
    use residual_Jac_module, only: spResTri,jacobDiagTri,jacobOffDiagRLTri,jacobOffDiagLRTri
    use solvers_module
    
    implicit none
    
contains

Subroutine steadySolve
    use get_delta_t_module
    use residual_Jac_module
    use solutionNorm_module
    use norm_module
    use vtu_output_module
    use inputs_module,  only: totalModes
    use LU_module
    use newtonRHS_module
    use form_Jx_module
    use assembleLocMassMatrices_module, only: MassMatrices
    use sparse_module 
    
    integer(i4) :: tri_id,edge_id,k,i,j,m,l,loop
    logical     :: updateFlag,CFLMod,CFLMod2,jacDiagLU
    real(dp)    :: newtonRHSvec(numFields*totalModes,numTri)
    real(dp)    :: Jx(numFields*totalModes,numTri)
    real(dp)    :: linRes(numFields*totalModes,numTri)
    real(dp)    :: deltaCoeffs(numFields*totalModes),x(numFields*totalModes,numTri)
    real(dp)    :: invDr(numFields*totalModes)
    real(dp)    :: solNorm,CFLmin,CFLmax,alp,norm0,resL2Ratio,linResNorm,xNorm,linResNorm0,linResNormOld
    real(dp)    :: locDt,implicitConvergeNorm
    
    real(dp)    :: Res(numFields,totalModes,numTri)
    
    CFLmin  = 100._dp
    CFLmax  = 10E10
    alp     = 1.25_dp
    loop = 0
    implicitConvergeNorm = 1.0e-6_dp
    
    
    call Residual_Jac(solCoeffs)        !returns spRes      
!    call plotResidual(out_file,spRes)
    call Lp_norm(spRes,2._dp,solNorm)
!    call write_Solution_Coefficients(out_file,solCoeffs)
    Call Lp_norm(spRes,2._dp,norm0)
    x(:,:) = 0._dp
    
    jacDiagLU = .false.
    if(solver .eq. 'direct')then
        Call allocate_sparse(numTri)
    end if

    
    do while(solNorm > 1.0e-12dp)
        x(:,:) = 0._dp
        Call get_delta_t(solCoeffs)
        Call Lp_norm(spRes,2._dp,solNorm)
        
        resL2Ratio = norm0/solNorm
        CFL = min(CFLmin*(resL2Ratio)**alp,CFLmax)
        print*,'CFL :',CFL
        
        delta_t = delta_t*CFL
        
        !form linear residual r = (Jx + R)
        !Add the mass matrix damping
        do tri_id = 1,numTri
            do i = 1,totalModes
                do j = 1,totalModes
                    do m = 1,numFields
                        k = m + (i-1)*numFields
                        l = m + (j-1)*numFields
                        jacobDiag(k,l,tri_id) = jacobDiag(k,l,tri_id) + (MassMatrices(i,j,tri_id)/delta_t(tri_id))
                    end do
                end do
            end do
        end do
        
        if (solver .eq. 'jacobi') then
            jacDiagLU = .true.
            do tri_id = 1,numTri
                Call my_lu(jacobDiag(:,:,tri_id))
            end do
        end if
        
        select case(solver)
            case('jacobi')
                Print*,'Using Jacobi Solver'
                Call jacobi_solver(x,implicitConvergeNorm)
            case('direct')
                Print*,'Using Direct Solver'
                Call direct_solver(x,spRes)
            case default
                print*,'Didnt choose a solver'
        end select 
        
        Call form_Jx(jacDiagLU,jacobDiag,jacobOffDiagRL,jacobOffDiagLR,x,Jx)
        print*, 'Jac difference:',sum(abs(Jx - reshape(spRes, (/numFields*totalModes ,numTri/))))
        
        do tri_id = 1,numTri
            solCoeffs(:,:,tri_id) = solCoeffs(:,:,tri_id) - 0.9_dp*reshape(x(:,tri_id),(/numFields,totalModes/))
        end do 

        call Residual_Jac(solCoeffs) !returns spRes 
        loop = loop + 1
        
        if (MOD(loop,out_freq)==0) then
            !Call plotResidual(out_file,spRes)
            Call callPlot(solCoeffs,spRes)
            Call write_Solution_Coefficients(out_file,solCoeffs)
        end if

        call Lp_norm(spRes,2._dp,solNorm)
        print*,'Norm:',solNorm, loop
        if (solNorm < 10e-12) then
            print*,'Solution reached steady state'
            Call callPlot(solCoeffs,spRes)
            stop 
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
end module
