subroutine gaussian_elimination(x,b)
	use my_kinddefs
	use sparse_module
	use solution_module
	use mesh_module
	implicit none
	
	real(dp), intent(inout) :: x(dof)
	real(dp), intent(inout) :: b(dof)
	
	integer(i4) :: istat,e,nftm,js,je,el,er,i,j,k,ni,nj,nk,find,nrow,ncol,pivot,info
        character(len=1) :: trans='N'
        real(dp) :: zerosum,begin_time,end_time
	
	do i=1,dof
            x(i) = b(i)
	end do
	
	print*,'sparse mem before ge', sparse_memory
    
        call cpu_time(begin_time)

	! row echelon form now upper triangular
	do k = 1,ncell
            nk = sparse(k,k)%nrow
            call dgetrf(nk,nk,sparse(k,k)%value,nk,sparse(k,k)%pivot,info) !LU 
            if(info /= 0) then; print*,'dgetrf error ge',info;stop;end if    

            do i = k+1,ncell
                if(sparse(i,k)%nonzero == 1) then
                    zerosum = sum(abs(sparse(i,k)%value))
                    if(zerosum < 1.0e-13) then		
                            call deallocate_sparse_block(i,k)
                            cycle
                    end if

                    ni = sparse(i,k)%nrow
                    !x(i) = x(i) - A(i,k) / A(k,k)*x(k);
                    call ge_modify_x(i,k,ni,nk,x(solni(i)),x(solni(k)),sparse(i,k)%value)

                    do j = k+1,ncell        
                        if(sparse(k,j)%nonzero == 1) then

                            zerosum = sum(abs(sparse(k,j)%value))
                            if(zerosum < 1.0e-13) then		
                                call deallocate_sparse_block(k,j)						
                                cycle
                            end if

                            nj = sparse(k,j)%ncol

                            if(sparse(i,j)%nonzero == 0) then					
                                pivot = 0
                                call allocate_sparse_block(i,j,ni,nj,pivot)
                                sparse(i,j)%value(:) = 0.0_dp						
                            end if

                            ! Aij = Aij - Aik * Akk^(-1) * Akj 

                            call ge_modify_sparse(i,j,k,ni,nj,nk)

                            zerosum = sum(abs(sparse(i,j)%value))
                            if(zerosum < 1.0e-13) then		
                                call deallocate_sparse_block(i,j)
                            end if

                        end if
                    end do
                    call deallocate_sparse_block(i,k)
                end if
            end do
	 end do
	call cpu_time(end_time)

	!	print*,'end ge'
	!print*,'beginning backward solve'
	print*,'sparse mem after ge', sparse_memory
	print*,'ge time', end_time-begin_time
	print*,'max sparse mem during ge', max_sparse_memory
	
	
!	do i=1,ncell
!		do j=1,ncell
!			if(sparse(i,j)%nonzero == 0) then
!				cycle
!			end if
!			zerosum = 0.0_dp
!			do k=1,sparse(i,j)%nrow**2
!				zerosum = zerosum + abs(sparse(i,j)%value(k))
!			end do
!			zerosum = sum(abs(sparse(i,j)%value))
!			if(zerosum < 1.0e-14) then
!				if (i == j) then
!					print*,i,j,'ahhh zero diagonal',sparse(i,j)%value
!					stop
!				end if
!				call deallocate_sparse_block(i,j)			
!			end if
!		end do
!	end do
!		print*,'sparse mem after ge', sparse_memory
!


        call cpu_time(begin_time)

	! backward solve
	! b(k) = (b(k) - A(k,k+1:ntri)*b(k+1:ntri))/A(k,k);
	do k= ncell,1,-1
            nk = sparse(k,k)%nrow

            do i=k+1,ncell
                if(sparse(k,i)%nonzero == 1) then
                    ni = sparse(k,i)%ncol
                    call backward_solve_modify_x(i,k,ni,nk,sparse(k,i)%value,x(solni(k):),x(solni(i):))
                end if
            end do

            call dgetrs(trans,nk,1,sparse(k,k)%value,nk,sparse(k,k)%pivot,x(solni(k):),nk,info)!LU solve
            if(info /= 0) then; print*,'dgetrs error ge 3',info;stop;end if     

	end do
    call cpu_time(end_time)
	print*,'bacward solve time', end_time-begin_time
		
	!	print*,'end backward solve'

end subroutine gaussian_elimination


subroutine backward_solve_modify_x(i,k,ni,nk,Aki,xk,xi)
    use my_kinddefs
    implicit none

    integer(i4), intent(in) :: i,k,ni,nk
    real(dp), intent(in) :: Aki(nk,ni)
    real(dp), intent(inout) :: xk(nk)
    real(dp), intent(in) :: xi(ni)

    xk = xk - matmul(Aki,xi)
	
end subroutine backward_solve_modify_x



subroutine ge_modify_x(i,k,ni,nk,xi,xk,Aik)
    use my_kinddefs
    use sparse_module
    implicit none

    integer(i4), intent(in) :: i,k
    real(dp), intent(inout) :: xi(ni)
    real(dp), intent(in) :: xk(nk)
    real(dp), intent(in) :: Aik(ni,nk)

    integer(i4) :: istat, ni,nj,nk
    integer(i4) :: info
    character(len=1) :: trans='N'
    real(dp) :: xtemp(nk)

    xtemp(:) = xk(:)

    ! xtemp = Akk^(-1) * xk
    call dgetrs(trans,nk,1,sparse(k,k)%value,nk,sparse(k,k)%pivot,xtemp,nk,info)!LU solve
    if(info /= 0) then; print*,'dgetrs error ge 2',info;stop;end if     

    xi(:) = xi(:) - matmul(Aik,xtemp)
end subroutine ge_modify_x


subroutine ge_modify_sparse(i,j,k,ni,nj,nk)
    use my_kinddefs
    use sparse_module
    implicit none

    integer(i4), intent(in) :: i,j,k,ni,nj,nk

    integer(i4) :: istat
    integer(i4) :: info
    character(len=1) :: trans='N'

    real(dp), allocatable :: Akj(:)

    allocate(Akj(nk*nj),stat=istat)

    Akj(:) = sparse(k,j)%value(:)

    ! Akj = Akk^(-1) * Akj
    call dgetrs(trans,nk,nj,sparse(k,k)%value,nk,sparse(k,k)%pivot,Akj,nk,info)
    if(info /= 0) then; print*,'dgetrs error ge',info;stop;end if     

    ! Aij = Aij - Aik*Akj where Akj = Akk^(-1) * Akj from above
    call ge_modify_ij(ni,nj,nk,sparse(i,j)%value,sparse(i,k)%value,Akj)

    deallocate(Akj)
	
end subroutine ge_modify_sparse


subroutine ge_modify_ij(ni,nj,nk,Aij,Aik,Akj)
    use my_kinddefs
    implicit none
    integer(i4), intent(in) :: ni,nj,nk
    real(dp), intent(inout) :: Aij(ni,nj)
    real(dp), intent(in) :: Aik(ni,nk)
    real(dp), intent(in) :: Akj(nk,nj)

    Aij = Aij - matmul(Aik,Akj)
	
end subroutine ge_modify_ij



subroutine setup_sparse()
    use my_kinddefs
    use sparse_module
    use solution_module
    use mesh_module
    implicit none

    integer(i4) :: e,nftm,js,je,el,er,i,j,find,nrow,ncol,pivot,info

    do i=1,ncell
        do j=1,i-1
            if(sparse(i,j)%nonzero == 1) then
                    call deallocate_sparse_block(i,j)
            end if
        end do
        do j=i+1,ncell
            if(sparse(i,j)%nonzero == 1) then
                    call deallocate_sparse_block(i,j)
            end if
        end do
    end do

    !print*,'moving diagonal jacobian'
    pivot = 1
    do e=1,ncell
        nftm = numfield*cell(e)%sol_nmode
        js = jacobiani(e)
        je = js + nftm*nftm-1

        if(sparse(e,e)%nonzero == 0) then
            call allocate_sparse_block(e,e,nftm,nftm,pivot)
        end if

        sparse(e,e)%value(:) = jacobian(js:je)

    end do
            !print*,'moving off diagonal jacobian'

    pivot = 0	
    do i=1,nface_interior
        find = face_interior(i)
        el = face(find)%cell(1)
        er = face(find)%cell(2)
        nrow = numfield*cell(er)%sol_nmode
        ncol = numfield*cell(el)%sol_nmode

        ! off diagonal jacobian R/L
        if(sparse(er,el)%nonzero == 0) then
            call allocate_sparse_block(er,el,nrow,ncol,pivot)
        end if

        js = jacobianoi(1,find)
        je = js + nrow*ncol-1
        sparse(er,el)%value(:) = jacobiano(js:je)

        ! off diagonal jacobian L/R
        if(sparse(el,er)%nonzero == 0) then
            call allocate_sparse_block(el,er,ncol,nrow,pivot)
        end if
        js = jacobianoi(2,find)
        je = js + nrow*ncol-1
        sparse(el,er)%value(:) = jacobiano(js:je)

    end do
	
	
	
end subroutine setup_sparse


