Module solvers_module
    use my_kinddefs
    use residual_Jac_module,only: jacobDiag,jacobOffDiagRL,jacobOffDiagLR,spRes
    use Globals_module, only: numTri,numFields,solutionNorm,triList,numInterior
    use inputs_module,  only: totalModes,solver
    use LU_module
    use form_Jx_module
    use norm_module
    
    implicit none
    contains
    
    Subroutine jacobi_solver(x,convergeNorm)
        real(dp),intent(in)       :: convergeNorm
        real(dp),intent(out)      :: x(numFields*totalModes,numTri)  
        
        logical     :: jacDiagLU
        integer(i4) :: tri_id,edge_id
        real(dp)    :: linRes(numFields*totalModes,numTri)
        real(dp)    :: invDr(numFields*totalModes),linResNorm0
        real(dp)    :: Jx(numFields*totalModes,numTri),linResNorm
        
        integer(i4) :: count
        
        jacDiagLU = .true.
        Call form_Jx(jacDiagLU,jacobDiag,jacobOffDiagRL,jacobOffDiagLR,x,Jx)
        
        !linear residual r = (Jx - R)
        linRes = Jx - reshape(spRes,(/numFields*totalModes , numTri/))
        
        Call Lp_norm2D(linRes,2._dp,linResNorm)
        linResNorm0 = linResNorm
        !linResNorm = linResNorm0
        
        count = 0
        do while (linResNorm/linResNorm0 > convergeNorm)
            do tri_id = 1,numTri
                Call lu_solve(jacobDiag(:,:,tri_id), linRes(:,tri_id), invDr)
                x(:,tri_id) = x(:,tri_id) - 0.8_dp*invDr
            end do

            Call form_Jx(jacDiagLU,jacobDiag,jacobOffDiagRL,jacobOffDiagLR,x,Jx)

            linRes = Jx - reshape(spRes,(/ numFields*totalModes , numTri /))

            Call Lp_norm2D(linRes,2._dp,linResNorm)
            !print*,'linear norm:',linResNorm/linResNorm0 ,linResNorm   
            count = count + 1
        end do
        print*,"Linear Iterations: ",count
        
    end subroutine jacobi_solver
    
    
    Subroutine direct_solver(x,Res)
        use gaussian_elimination_module
        use sparse_module
        use form_jx_module
        implicit none
        
        real(dp),intent(inout) :: x(numFields*totalModes,numTri)  
        real(dp),intent(inout) :: Res(numFields,totalModes,numTri)  
        
        Call setup_sparse(jacobDiag,jacobOffDiagRL,jacobOffDiagLR)
        Call gaussian_elimination(x,Res)
         
        
    end subroutine direct_solver
    
end module solvers_module

    
