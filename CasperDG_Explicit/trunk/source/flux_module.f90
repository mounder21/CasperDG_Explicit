Module flux_module
    use my_kinddefs
    use Globals_module, only: numEulerVars
    implicit none 
    
contains 

subroutine getFlux(qL,qR,n,numEdgeGaussPts,fluxFunction,flux)
    real(dp),       intent(in) :: qL(:,:),qR(:,:),n(:,:) !q(ngp,numEulerVars,),n(ngp,x/y)
    character(*),   intent(in) :: fluxFunction
    integer(i4),    intent(in) :: numEdgeGaussPts
    real(dp),       intent(out):: flux(numEulerVars,numEdgeGaussPts)
   
    !dummy 
    integer(i4) :: GP
    

    flux(:,:) = 0.0_dp
    select case(fluxFunction)
        case('Roe')
            do GP = 1,numEdgeGaussPts
                Call Roe(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) 
            end do
        case('LaxF')
            do GP = 1,numEdgeGaussPts
                 Call lax_friedrichs(qL(:,GP),qR(:,GP),n(GP,:),flux(:,GP))
            end do
        case default 
            print*,'Choose either FaxL, Roe or RHLL for flux'
    end select
end subroutine getFlux

subroutine getBCFlux(q,n,numEdgeGaussPts,flux)
    integer(i4),intent(in) :: numEdgeGaussPts
    real(dp),   intent(in) :: q(numEulerVars,numEdgeGaussPts),n(numEdgeGaussPts,2) !q(numEulerVars,ngp),n(ngp,x/y)
    real(dp),   intent(out):: flux(numEulerVars,numEdgeGaussPts)
   
    !dummy 
    integer(i4) :: GP,i
    real(dp)    :: fluxLoc(numEulerVars,2)
    

    do GP = 1, numEdgeGaussPts
        fluxLoc(:,:) = 0._dp
        Call getNativeFlux(q(:,GP),fluxLoc)
        
        flux(:,GP) = fluxLoc(:,1)*n(GP,1) +  fluxLoc(:,2)*n(GP,2)
    end do
    
    
end subroutine getBCFlux

subroutine getNativeFlux(q,flux)
    real(dp),intent(in) :: q(numEulerVars) 
    real(dp),intent(out):: flux(numEulerVars,2)
    
    real(dp) :: gamma,rho,u,v,p,E
    gamma = 1.4_dp

    rho = q(1)
    u   =  q(2)/q(1)
    v   =  q(3)/q(1)
    E   = q(4)    
    p   = (gamma-1.0_dp)*(q(4) - 0.5_dp*rho*(u**2 + v**2))

    

    flux(1,1) = rho*u
    flux(2,1) = rho*u*u + p
    flux(3,1) = rho*u*v
    flux(4,1) = u*(q(4) + p)

    flux(1,2) = rho*v
    flux(2,2) = rho*u*v
    flux(3,2) = rho*v*v + p
    flux(4,2) = v*(q(4) + p) 

end subroutine getNativeFlux

!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 subroutine Roe(uL, uR, nx, ny,flux)
     real(dp),intent(in) :: uL(:), uR(:) !  Input: conservative variables rho*[1, u, v, E]
     real(dp),intent(in) :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
     real(dp),intent(out) :: flux(:)       ! Output: Roe flux function (upwind)

    !Local constants
     real(dp) :: gamma                          ! Ratio of specific heat.

    !Local variables
     real(dp) :: tx, ty       ! Tangent vector (perpendicular to the face normal)
     real(dp) :: vxL, vxR, vyL, vyR             ! Velocity components.
     real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
     real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
     real(dp) :: aL, aR, HL, HR                 ! Speeds of sound.
     real(dp) :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
     real(dp) :: drho,dvn,dvt,dpres,dV(4)  ! Wave strenghs
     real(dp) :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
     real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
     real(dp) :: mag,nxT,nyT
     integer(i4) :: i, j

     mag = sqrt( nx**2 + ny**2)  
     nxT = nx/mag
     nyT = ny/mag

!Constants.
    gamma = 1.4_dp

  tx = -nyT
  ty = nxT

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nxT+vyL*nyT
     vtL = vxL*tx+vyL*ty
      pL = (gamma-1.0_dp)*( uL(4) - 0.5_dp*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
      
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nxT+vyR*nyT
     vtR = vxR*tx+vyR*ty
      pR = (gamma-1.0_dp)*( uR(4) - 0.5_dp*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(1.0_dp+RT)
    vy = (vyL+RT*vyR)/(1.0_dp+RT)
     H = (HL+RT * HR)/(1.0_dp +RT)
     a = sqrt((1.4_dp-1.0_dp)*(H-0.5_dp*(vx*vx+vy*vy)) )
    vn = vx*nxT+vy*nyT
    vt = vx*tx+vy*ty
    
!Wave Strengths
    drho = rhoR - rhoL 
    dpres =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL
    
    
    dV(:) = 0.0_dp
    dV(1) = (dpres - rho*a*dvn )/(2.0_dp*a*a)
    dV(2) = rho*dvt/a
    dV(3) =  drho - dp/(a*a)
    dV(4) = (dpres + rho*a*dvn )/(2.0_dp*a*a)

  ws(:) = 0.0_dp
!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = 0.2_dp
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = 0.2_dp
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = 1.0_dp    
  Rv(2,1) = vx - a*nxT
  Rv(3,1) = vy - a*nyT
  Rv(4,1) =  H - vn*a

  Rv(1,2) = 0.0_dp
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = 1.0_dp
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = 0.5_dp*(vx*vx+vy*vy)

  Rv(1,4) = 1.0_dp
  Rv(2,4) = vx + a*nxT
  Rv(3,4) = vy + a*nyT
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = 0.0_dp
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nxT
  fL(3) = rhoL*vnL * vyL + pL*nyT
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nxT
  fR(3) = rhoR*vnR * vyR + pR*nyT
  fR(4) = rhoR*vnR *  HR

  flux = half * (fL + fR - diss)
  flux = flux*mag
 end subroutine Roe
  
  subroutine lax_friedrichs(ql, qr, norm, flux)
    use my_kinddefs
    
    real(dp),intent(in)  :: ql(:)
    real(dp),intent(in)  :: qr(:)
    real(dp),intent(in)  :: norm(2)
    real(dp),intent(out) :: flux(:)

    !---> Local Variables
    real(dp) :: eigl, eigr,eigl_ql(4),eigr_qr(4)  
    real(dp) :: nx, ny
    real(dp) :: al, ar,al_ql(4),ar_qr(4)
    real(dp) :: ul, vl, pl, El
    real(dp) :: ur, vr, pr, Er
    real(dp) :: flux_l
    real(dp) :: flux_r
    real(dp) :: alpha
    real(dp) :: mag
    real(dp) :: gamma
    integer(i4) :: k
    
    gamma = 1.4_dp
    mag = sqrt( norm(1)**2 + norm(2)**2 )    

    nx = norm(1)/mag
    ny = norm(2)/mag
      
    ul = ql(2)/ql(1) 
    ur = qr(2)/qr(1)       
   
    
    vl = ql(3)/ql(1)
    vr = qr(3)/qr(1)
    

    El = ql(4)/ql(1)
    Er = qr(4)/qr(1)
    

    pl = ql(1)*(gamma - 1._dp)*(El - half*( ul**2 + vl**2) ) 
    pr = qr(1)*(gamma - 1._dp)*(Er - half*( ur**2 + vr**2) )  
    
   
    al = sqrt( gamma*pl/ql(1))
    ar = sqrt( gamma*pr/qr(1))    
    
    eigl = al + abs(ul*nx + vl*ny)
    eigr = ar + abs(ur*nx + vr*ny)    
    
    !alpha = max( eigl,eigr )
    if( eigl > eigr) then
       alpha = eigl
    else if( eigl <= eigr ) then
       alpha = eigr
    end if


    !---> Now we do the flux's 1 at a time
    !---> Flux 1 
    flux_l = ql(2)*nx + ql(3)*ny
    flux_r = qr(2)*nx + qr(3)*ny
   
    flux(1) = half * ( ( flux_l + flux_r) + alpha*( ql(1) - qr(1) ))
   
    !---> Flux 2 

    flux_l = (ql(2)*ul + pl)*nx + (ql(3)*ul)*ny
    flux_r = (qr(2)*ur + pr)*nx + (qr(3)*ur)*ny

    flux(2) = half*( (flux_l + flux_r) + alpha*( ql(2) - qr(2) ))
    
    !---> Flux 3 

    flux_l = (ql(2)*vl)*nx + (ql(3)*vl + pl)*ny
    flux_r = (qr(2)*vr)*nx + (qr(3)*vr + pr)*ny

    flux(3) = half*( (flux_l + flux_r) + alpha*( ql(3) - qr(3) ))

    !---> Flux 4 

    flux_l = (ql(4)*ul + pl*ul)*nx + (ql(4)*vl + pl*vl)*ny
    flux_r = (qr(4)*ur + pr*ur)*nx + (qr(4)*vr + pr*vr)*ny 

    flux(4) = half * ( (flux_l + flux_r) + alpha*( ql(4) - qr(4) ))
    
    flux    = flux*mag

  end subroutine lax_friedrichs
  
 real(dp) function my_smooth_dabs(x)
      real(dp),intent(in)   :: x
      real(dp)              :: my_epsilon
      
      my_epsilon = 1e-9
      
      my_smooth_dabs = x*(abs(x) + 2._dp*my_epsilon)/(abs(x) + my_epsilon)**2
      
  end function

end module flux_module