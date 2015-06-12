Module allocate_module
    use get_delta_t_module
    use integrate_module
    implicit none
contains 

Subroutine alloc
    call allocate_delta_t
end subroutine alloc
end module allocate_module