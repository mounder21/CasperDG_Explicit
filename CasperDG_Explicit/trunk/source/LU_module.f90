Module LU_module
    use my_kinddefs
    contains
    
Subroutine my_lu(A) !In place LU decomposition
    real(dp),intent(inout)  :: A(:,:)
    
    integer(i4) :: i,j,k,m
    m = size(A,1)
    
    do i = 1, m - 1
        do j = i + 1, m
            A(j,i) = A(j,i)/A(i,i)
            A(j,i+1:m) = A(j,i+1:m) - A(j,i)*A(i,i+1:m)
            if(isnan(A(j,i))) then
                print*, 'LU decomposition produced Nans'
                stop
            end if
        end do
    end do
end subroutine


!*******************************************************
!      LU_solve : Takes L and U and b and finds x  *
!*******************************************************
subroutine lu_solve(A, b, sol_vec)
  use my_kinddefs
  implicit none
  
  real(dp),intent(in) :: A(:,:)
  real(dp),intent(in)   :: b(:)
  real(dp),intent(out)  :: sol_vec(:)
  integer(i4)           :: n ! matrix size
  integer(i4)           :: m ! index to run over
  
  !---> Local Variables
  integer(i4) i,j
  
  n = size(A,1)
  m = size(A,1)
  
  sol_vec(1:m) = b(1:m)
  do i = 2,m
     do j = 1,i-1
        sol_vec(i) = sol_vec(i) - A(i,j)*sol_vec(j)
     end do
  end do
  
  ! Now do the backwards substitution
  sol_vec(m) = sol_vec(m)/A(m,m)
  do i = m - 1, 1, -1
     do j = m, i + 1, -1
        sol_vec(i) = sol_vec(i)-A(i,j)*sol_vec(j)
     end do
     sol_vec(i) = sol_vec(i)/A(i,i)
  end do
  
end subroutine lu_solve

!*****************************************************************************80
!        lu_matvec: multiplies A*x,when A is stored in LU form.
!*****************************************************************************80
subroutine lu_matvec(A, x, v)
  use my_kinddefs
  implicit none
  
  real(dp),intent(in) :: A(:,:)
  real(dp),intent(inout):: x(:)
  real(dp),intent(out)  :: v(:)
  integer(i4)             :: n ! matrix size
  integer(i4)             :: m ! index to run over
  !---> Local Varaiables
  integer(i4) i, j
 
  n = size(A,1)
  m = size(A,1)

  !---> NEW and HOPEFULLY BETTER WAY
  v(1:m) = 0._dp
  do j = 1,m
     ! v(i) = 0._dp
     do i = 1,j
        v(i) = v(i) + A(i,j)*x(j)
     end do
  end do
  x(1:m) = v(1:m)
  do j = 1,m
     ! v(i) = 0._dp
     do i = j+1,m
        v(i) = v(i) + A(i,j)*x(j)
     end do
     ! v(i) = v(i) + x(i)
  end do
end subroutine lu_matvec


  Subroutine test_lu
      
      real(dp) :: A(4,4),Atemp(4,4),x(4),v(4)
      
      A(1,1) = 2
      A(1,2) = 1
      A(1,3) = 0
      A(1,4) = 0
      A(2,1) = 1
      A(2,2) = 2
      A(2,3) = 1
      A(2,4) = 0
      A(3,1) = 0
      A(3,2) = 1
      A(3,3) = 2
      A(3,4) = 1
      A(4,1) = 0
      A(4,2) = 0
      A(4,3) = 1
      A(4,4) = 2
      
      x(1:3) = 1
      x(4)   = 2
       
      Atemp = A
      print*, matmul(A,x)
     
      Call my_lu(Atemp)
      Call lu_matvec(Atemp,x,v)
      print*,v
      
  end subroutine


        
end module
