module sparse_module

	use my_kinddefs
        use Globals_module, only: numTri,numFields
        use inputs_module,  only: totalModes
	implicit none
	
    type sparse_struct
        !integer(i4) :: nrow
        !integer(i4) :: ncol
        integer(i4) :: nonzero = 0
        real(dp), pointer :: blockVal(:,:) 
        integer(i4), pointer :: pivot(:)
    end type sparse_struct

    type(sparse_struct), pointer :: sparse(:,:)
                                                                ! array of cell struct
	real(dp) :: sparse_memory = 0.0_dp
	real(dp) :: max_sparse_memory = 0.0_dp
	
contains
	
	
    subroutine allocate_sparse(numEle)
        implicit none
        integer(i4),intent(in)  :: numEle

        integer(i4) :: istat
        allocate(sparse(numEle,numEle),stat = istat)

        print*,'Allocating sparse'


    end subroutine allocate_sparse

    subroutine deallocate_sparse
        implicit none

        deallocate(sparse)

    end subroutine deallocate_sparse

    subroutine allocate_sparse_block(i,j,pivot)

            implicit none
            integer(i4),intent(in) :: i,j,pivot

            integer(i4) :: istat,n

            !sparse(i,j)%nrow = numFields*totalModes
            !sparse(i,j)%ncol = numFields*totalModes
            sparse(i,j)%nonzero = 1

            allocate(sparse(i,j)%blockVal(numFields*totalModes,numFields*totalModes),stat = istat)		

            sparse_memory = sparse_memory + &
                     8.0_dp*real(numFields*totalModes*numFields*totalModes,dp)/1024.0_dp/1024.0_dp

            if(pivot == 1) then
                n = numFields*totalModes*numFields
                allocate(sparse(i,j)%pivot(n),stat = istat)	
                sparse_memory = sparse_memory + 4.0_dp*real(n,dp)/1024.0_dp/1024.0_dp
            else
                allocate(sparse(i,j)%pivot(0),stat = istat)	
            end if

            if(max_sparse_memory < sparse_memory) then
                max_sparse_memory = sparse_memory
            end if
            !print*,'allocating more memory', memory, 'MB'

    end subroutine


    subroutine deallocate_sparse_block(i,j)

        implicit none

        integer(i4),intent(in) :: i,j

        ! probably won't ever deallocate diagonals which are pivoted

        sparse_memory = sparse_memory - &
            8.0_dp*real(numFields*totalModes*numFields*totalModes,dp)/1024.0_dp/1024.0_dp

        deallocate(sparse(i,j)%blockVal)
        sparse(i,j)%nonzero = 0
        !print*,'deallocating memory', memory, 'MB'

    end subroutine


	
end module sparse_module