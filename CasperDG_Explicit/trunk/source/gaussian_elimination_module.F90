Module gaussian_elimination_module
    use my_kinddefs
    use LU_module
    implicit none
    
contains
    
    
subroutine gaussian_elimination(x,b)
	use my_kinddefs
	use sparse_module
        use Globals_module, only: numTri,numFields
        use inputs_module,  only: totalModes
        
	implicit none
	
	real(dp), intent(inout) :: x(numFields*totalModes,numTri)
	real(dp), intent(inout) :: b(numFields*totalModes,numTri)
	
	integer(i4) :: istat,e,nftm,js,je,el,er,i,j,k,find,pivot,info,nk
        character(len=1) :: trans='N'
        real(dp) :: zerosum,begin_time,end_time,temp(numFields*totalModes)

        x = b
	
	print*,'sparse mem before ge', sparse_memory
    
        call cpu_time(begin_time)
        nk = numFields*totalmodes
	! row echelon form now upper triangular
	do k = 1,numTri
            print*,'element:', k
#ifdef USE_LAPACK
            Call dgetrf(nk,nk,sparse(k,k)%blockVal,nk,sparse(k,k)%pivot,info)
            if(info /= 0) then; print*,'dgetrs error ge 3',info;stop;end if
#else   
            Call my_lu(sparse(k,k)%blockVal)    !LU the diagonal
#endif

            do i = k+1,numTri
                if(sparse(i,k)%nonzero == 1) then
                    zerosum = sum(abs(sparse(i,k)%blockVal))
                    if(zerosum < 1.0e-8) then		
                            call deallocate_sparse_block(i,k)
                            cycle
                    end if

                    !x(i) = x(i) - A(i,k) / A(k,k)*x(k);
                    call ge_modify_x(i,k,x(:,i),x(:,k),sparse(i,k)%blockVal)

                    do j = k+1,numTri        
                        if(sparse(k,j)%nonzero == 1) then

                            zerosum = sum(abs(sparse(k,j)%blockVal))
                            if(zerosum < 1.0e-13) then		
                                call deallocate_sparse_block(k,j)						
                                cycle
                            end if

                            if(sparse(i,j)%nonzero == 0) then					
                                pivot = 0
                                call allocate_sparse_block(i,j,pivot)
                                sparse(i,j)%blockVal(:,:) = 0.0_dp						
                            end if

                            ! Aij = Aij - Aik * Akk^(-1) * Akj 

                            call ge_modify_sparse(i,j,k)

                            zerosum = sum(abs(sparse(i,j)%blockVal))
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
	do k = numTri,1,-1
            do i = k + 1,numTri
                if(sparse(k,i)%nonzero == 1) then
                    x(:,k) = x(:,k) - matmul(sparse(k,i)%blockVal,x(:,i))
                end if
            end do
#ifdef USE_LAPACK
            call dgetrs(trans,nk,1,sparse(k,k)%blockVal,nk,sparse(k,k)%pivot,x(:,k),nk,info)
            if(info /= 0) then; print*,'dgetrs error ge 3',info;stop;end if
#else
            Call lu_solve(sparse(k,k)%blockVal,x(:,k),temp)      !this could be over writing 
            x(:,k) = temp
#endif
	end do
    call cpu_time(end_time)
	print*,'backward solve time', end_time-begin_time
		
	!	print*,'end backward solve'

end subroutine gaussian_elimination




subroutine ge_modify_x(i,k,xi,xk,Aik)
    use my_kinddefs
    use sparse_module
    implicit none

    integer(i4), intent(in) :: i,k
    real(dp), intent(inout) :: xi(numFields*totalModes)
    real(dp), intent(in) :: xk(numFields*totalModes)
    real(dp), intent(in) :: Aik(numFields*totalModes,numFields*totalModes)

    integer(i4) :: istat, ni,nj,nk
    integer(i4) :: info
    character(len=1) :: trans='N'
    real(dp) :: xtemp(numFields*totalModes)
    
    nk = numFields*totalModes
    xtemp(:) = xk(:)
    !x(i) = x(i) - A(i,k) / A(k,k)*x(k);
    ! xtemp = Akk^(-1) * xk 

    
    
    ! xtemp = Akk^(-1) * xk
#ifdef USE_LAPACK
    call dgetrs(trans,nk,1,sparse(k,k)%blockVal,nk,sparse(k,k)%pivot,xtemp,nk,info)
    if(info /= 0) then; print*,'dgetrs error ge 2',info;stop;end if     
#else
    Call  lu_solve(sparse(k,k)%blockVal,xk,xtemp)
#endif
    xi(:) = xi(:) - matmul(Aik,xtemp)
end subroutine ge_modify_x


subroutine ge_modify_sparse(i,j,k)
    use my_kinddefs
    use sparse_module
    implicit none

    integer(i4), intent(in) :: i,j,k

    integer(i4) :: istat,z,nk,nj
    integer(i4) :: info
    character(len=1) :: trans='N'

    real(dp), allocatable :: Akj(:,:)
    real(dp)    :: solVec(numFields*totalModes)

    allocate(Akj(numFields*totalModes,numFields*totalModes),stat=istat)
    ! Aij = Aij - Aik * Akk^(-1) * Akj 
    nk = numFields*totalModes
    nj = nk
    Akj(:,:) = sparse(k,j)%blockVal(:,:)

    ! Akj = Akk^(-1) * Akj
#ifdef USE_LAPACK
    call dgetrs(trans,nk,nj,sparse(k,k)%blockVal,nk,sparse(k,k)%pivot,Akj,nk,info)
    if(info /= 0) then; print*,'dgetrs error ge',info;stop;end if
#else
    do z = 1,numFields*totalModes
        Call lu_solve(sparse(k,k)%blockVal,Akj(:,z), solVec)       !this might overwrite so check here
        Akj(:,z) =  solVec
    end do
#endif
    ! Aij = Aij - Aik*Akj where Akj = Akk^(-1) * Akj from above
    sparse(i,j)%blockVal = sparse(i,j)%blockVal - matmul(sparse(i,k)%blockVal,Akj)
    
    deallocate(Akj)
	
end subroutine ge_modify_sparse



subroutine setup_sparse(jacobDiag,jacobOffDiagRL,jacobOffDiagLR)
    use my_kinddefs
    use sparse_module
    use Globals_module, only: numTri,numFields,numInterior,interiorEdges,edgeList
    implicit none
    
    real(dp),intent(in) :: jacobDiag(numFields*totalModes,numFields*totalModes,numTri)
    real(dp),intent(in) :: jacobOffDiagRL(numFields*totalModes,numFields*totalModes,numInterior)
    real(dp),intent(in) :: jacobOffDiagLR(numFields*totalModes,numFields*totalModes,numInterior)

    integer(i4) :: e,el,er,i,j,find,pivot

    do i = 1,numTri
        do j = 1, i - 1
            if(sparse(i,j)%nonzero == 1) then
                    call deallocate_sparse_block(i,j)
            end if
        end do
        do j = i + 1,numTri
            if(sparse(i,j)%nonzero == 1) then
                call deallocate_sparse_block(i,j)
            end if
        end do
    end do

    !print*,'moving diagonal jacobian'
    pivot = 1
    do e = 1,numTri
        if(sparse(e,e)%nonzero == 0) then
            call allocate_sparse_block(e,e,pivot)
        end if
        sparse(e,e)%blockVal(:,:) = jacobDiag(:,:,e)
    end do
            !print*,'moving off diagonal jacobian'

    pivot = 0	
    do i = 1,numInterior  
        find = interiorEdges(i)
        el = edgeList(find)%e1 
        er = edgeList(find)%e2

        ! off diagonal jacobian R/L
        if(sparse(er,el)%nonzero == 0) then
            call allocate_sparse_block(er,el,pivot)
        end if

        sparse(er,el)%blockVal(:,:) = jacobOffDiagRL(:,:,i)      !Check if this is RL or LR

        ! off diagonal jacobian L/R
        if(sparse(el,er)%nonzero == 0) then
            call allocate_sparse_block(el,er,pivot)
        end if
        
        sparse(el,er)%blockVal(:,:) = jacobOffDiagLR(:,:,i)      !Check if this is LR or RL

    end do
	
end subroutine setup_sparse

end Module gaussian_elimination_module


