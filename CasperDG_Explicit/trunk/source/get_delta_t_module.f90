Module get_delta_t_module
    use my_kinddefs
    use get_h_module, only: h
    use Globals_module, only: numEulerVars,numTri
    use inputs_module, only: gamma
    use initializeBasis_module, only: Phi,numGaussPts
    use projection_module
    
    implicit none 
    real(dp),pointer         :: delta_t(:)
    real(dp)         :: min_delta_t = 100._dp
        
contains

Subroutine get_delta_t(solCoeffs)
    real(dp), intent(in)    :: solCoeffs(:,:,:)
    
    real(dp)    :: locQ(numEulerVars,numGaussPts),maxU,maxV,maxC,pOverRho(numGaussPts),c(numGaussPts)
    real(dp)    :: u(numGaussPts),v(numGaussPts),rho(numGaussPts),currMin,currhMin
    integer(i4) :: elem_id,GP
    
    currMin = 10._dp
    currhMin = 1000._dp
    do elem_id = 1, numTri
        if(currhMin > h(elem_id)) then
            currhMin = h(elem_id)
        end if
    end do
    
    do elem_id = 1, numTri
        Call projectCell(solCoeffs(:,:,elem_id),Phi,locQ) ! locQ(numEulerVars,npgs)
        maxU = 0._dp
        maxV = 0._dp
        maxC = 0._dp
        do GP = 1,numGaussPts
            u(GP) = locQ(2,GP)/locQ(1,GP)
            v(GP) = locQ(3,GP)/locQ(1,GP)
                   
            maxU = max(maxU,u(GP))
            maxV = max(maxV,v(GP))
        
            rho(GP) = locQ(1,GP)
            pOverRho(GP) = (locQ(4,GP)/rho(GP) - 0.5_dp * (u(GP)**2 + v(GP)**2)) * (gamma-1._dp) 

            c(GP) = sqrt(gamma* pOverRho(1))
            maxC = max(maxC,c(GP))
        end do
        delta_t(elem_id) =  h(elem_id)/(3._dp*(sqrt(maxU**2 + maxV**2) + maxC))
        min_delta_t = min(min_delta_t,delta_t(elem_id))
    end do
end subroutine get_delta_t

Subroutine allocate_delta_t
    Allocate(delta_t(numTri))
end subroutine

end module