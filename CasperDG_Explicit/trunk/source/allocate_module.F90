Module allocate_module
    use get_delta_t_module
    use residual_module
    use residualComplex_module
    use residual_Jac_module
    use integrate_module
    implicit none
contains 

Subroutine alloc
    use inputs_module, only: Newton
    call allocate_delta_t
    
    
#ifdef complx
    call allocate_spResRC
    call allocate_spRes  
    call allocate_Jacobian
#endif
    
    if(Newton) then
        call allocate_Jacobian
        call allocate_spRes   
    else    ! time solve solution
       call allocate_spacialResidual
    end if
end subroutine alloc


end module allocate_module