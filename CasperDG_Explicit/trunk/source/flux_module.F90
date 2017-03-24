Module flux_module
    use my_kinddefs
    use Globals_module, only: numFields
    implicit none 
    
contains 

subroutine getFlux_q(qL,qR,n,numEdgeGaussPts,fluxFunction,flux,flux_ql,flux_qr)
    real(dp),       intent(in) :: qL(:,:),qR(:,:),n(:,:) !q(numFields,ngp),n(ngp,x/y)
    character(*),   intent(in) :: fluxFunction
    integer(i4),    intent(in) :: numEdgeGaussPts
    real(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
    real(dp),       intent(out):: flux_ql(numFields,numFields,numEdgeGaussPts)
    real(dp),       intent(out):: flux_qr(numFields,numFields,numEdgeGaussPts)
   
    !dummy 
    integer(i4) :: GP
    

    flux(:,:) = 0.0_dp
    select case(fluxFunction)
        case('Roe')
            do GP = 1,numEdgeGaussPts
                Call Roe(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) !Fix this
            end do
        case('RHLL')
            do GP = 1,numEdgeGaussPts
                 Call Rotated_RHLL_q(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP),&
                                    flux_ql(:,:,GP),flux_qr(:,:,GP)) 
            end do
        case('LaxF')
            do GP = 1,numEdgeGaussPts
                 Call lax_friedrichs_q(qL(:,GP),qR(:,GP),n(GP,:),flux(:,GP),&
                                    flux_ql(:,:,GP),flux_qr(:,:,GP))
            end do
        case default 
            print*,'Choose either LaxF, Roe or RHLL for flux'
    end select
end subroutine getFlux_q

subroutine getFlux(qL,qR,n,numEdgeGaussPts,fluxFunction,flux)

    real(dp),       intent(in) :: qL(:,:),qR(:,:)
    real(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
    
    real(dp),intent(in) :: n(:,:) !n(ngp,x/y)
    character(*),   intent(in) :: fluxFunction
    integer(i4),    intent(in) :: numEdgeGaussPts

   
    !dummy 
    integer(i4) :: GP
    

    flux(:,:) = 0.0_dp
    select case(fluxFunction)
        case('Roe')
            do GP = 1,numEdgeGaussPts
                !Call Roe(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) 
            end do
        case('RHLL')
            do GP = 1,numEdgeGaussPts
                 Call Rotated_RHLL(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) 
            end do
        case('LaxF')
            do GP = 1,numEdgeGaussPts
                 Call lax_friedrichs(qL(:,GP),qR(:,GP),n(GP,:),flux(:,GP))
            end do
        case default 
            print*,'Choose either FaxL, Roe or RHLL for flux'
    end select
end subroutine getFlux

subroutine getFluxComplex(qL,qR,n,numEdgeGaussPts,fluxFunction,flux)
#ifdef complx
    complex(dp),       intent(in) :: qL(:,:),qR(:,:)
    complex(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
#else
    real(dp),       intent(in) :: qL(:,:),qR(:,:)
    real(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
#endif
    
    real(dp),intent(in) :: n(:,:) !n(ngp,x/y)
    character(*),   intent(in) :: fluxFunction
    integer(i4),    intent(in) :: numEdgeGaussPts

   
    !dummy 
    integer(i4) :: GP
    

    flux(:,:) = 0.0_dp
    select case(fluxFunction)
        case('Roe')
            do GP = 1,numEdgeGaussPts
                !Call Roe(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) 
            end do
        case('RHLL')
            do GP = 1,numEdgeGaussPts
                 Call Rotated_RHLLComplex(qL(:,GP),qR(:,GP),n(GP,1),n(GP,2),flux(:,GP)) 
            end do
        case('LaxF')
            do GP = 1,numEdgeGaussPts
                 Call lax_friedrichsComplex(qL(:,GP),qR(:,GP),n(GP,:),flux(:,GP))
            end do
        case default 
            print*,'Choose either FaxL, Roe or RHLL for flux'
    end select
end subroutine getFluxComplex

subroutine getBCFlux_q(q,n,numEdgeGaussPts,flux,flux_q)
    integer(i4),intent(in) :: numEdgeGaussPts
    real(dp),   intent(in) :: q(numFields,numEdgeGaussPts),n(numEdgeGaussPts,2) !q(numFields,ngp),n(ngp,x/y)
    real(dp),   intent(out):: flux(numFields,numEdgeGaussPts)
    real(dp),   intent(out):: flux_q(numFields,numFields,numEdgeGaussPts)
   
    !dummy 
    integer(i4) :: GP,i
    real(dp)    :: fluxLoc(numFields,2),fluxLoc_q(numFields,2,numFields)
    

    do GP = 1, numEdgeGaussPts
        fluxLoc(:,:) = 0._dp
        Call getNativeFlux_q(q(:,GP),fluxLoc,fluxLoc_q)
        
        flux(:,GP) = fluxLoc(:,1)*n(GP,1) +  fluxLoc(:,2)*n(GP,2)
        
        do i = 1,numFields
            flux_q(:,i,GP) = fluxLoc_q(:,1,i)*n(GP,1) +  fluxLoc_q(:,2,i)*n(GP,2)
        end do
    end do
    
    
end subroutine getBCFlux_q

subroutine getBCFlux(q,n,numEdgeGaussPts,flux)
    real(dp),       intent(in) :: q(numFields,numEdgeGaussPts)
    real(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
    real(dp)                   :: fluxLoc(numFields,2)

    integer(i4),intent(in) :: numEdgeGaussPts
    real(dp),   intent(in) :: n(numEdgeGaussPts,2) !q(numFields,ngp),n(ngp,x/y)

    !dummy 
    integer(i4) :: GP,i


    do GP = 1, numEdgeGaussPts
        fluxLoc(:,:) = 0._dp
        Call getNativeFlux(q(:,GP),fluxLoc)
        flux(:,GP) = fluxLoc(:,1)*n(GP,1) +  fluxLoc(:,2)*n(GP,2)
    end do
      
end subroutine getBCFlux

subroutine getBCFluxComplex(q,n,numEdgeGaussPts,flux)
#ifdef complx
    complex(dp),       intent(in) :: q(numFields,numEdgeGaussPts)
    complex(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
    complex(dp)                   :: fluxLoc(numFields,2)
#else
    real(dp),       intent(in) :: q(numFields,numEdgeGaussPts)
    real(dp),       intent(out):: flux(numFields,numEdgeGaussPts)
    real(dp)                   :: fluxLoc(numFields,2)
#endif    

    integer(i4),intent(in) :: numEdgeGaussPts
    real(dp),   intent(in) :: n(numEdgeGaussPts,2) !q(numFields,ngp),n(ngp,x/y)

    !dummy 
    integer(i4) :: GP,i


    do GP = 1, numEdgeGaussPts
        fluxLoc(:,:) = 0._dp
        Call getNativeFluxComplex(q(:,GP),fluxLoc)
        flux(:,GP) = fluxLoc(:,1)*n(GP,1) + fluxLoc(:,2)*n(GP,2)
    end do
      
end subroutine getBCFluxComplex

subroutine getNativeFlux_q(q,flux,flux_q)
    real(dp),intent(in) :: q(numFields) 
    real(dp),intent(out):: flux(numFields,2),flux_q(numFields,2,numFields)
    
    real(dp) :: gamma,rho,u,v,p,E
    real(dp) :: rho_q(numFields),u_q(numFields),v_q(numFields),E_q(numFields),p_q(numFields)
    gamma = 1.4_dp

    rho = q(1)
    rho_q(1)    = 1._dp
    rho_q(2:4)  = 0._dp 
    
    u       =  q(2)/q(1)
    u_q(1)  = -q(2)/q(1)**2
    u_q(2)  = 1._dp/q(1)
    u_q(3:4)= 0._dp
    
    v       =  q(3)/q(1)
    v_q(1)  = -q(3)/q(1)**2
    v_q(2)  = 0._dp
    v_q(3)  = 1._dp/q(1)
    v_q(4)  = 0._dp
    
    E = q(4)
    E_q(1:3)= 0._dp
    E_q(4)  = 1._dp
    
    p = (gamma-1.0_dp)*(q(4) - 0.5_dp*rho*(u**2 + v**2))
    
    p_q(1)    = half*(gamma - 1._dp)*((q(2)/q(1))**2 + (q(3)/q(1))**2)        !dp/drho 
    p_q(2:3)  = -(gamma - 1._dp)*(q(2:3)/q(1))                                  !dp/drhou  
    p_q(4)    =  (gamma - 1._dp)                                                  !dp/drhoE 
    

    flux(1,1) = rho*u
    flux(2,1) = rho*u*u + p
    flux(3,1) = rho*u*v
    flux(4,1) = u*(q(4) + p)
    
    flux_q(1,1,:) = rho_q*u + rho*u_q
    flux_q(2,1,:) = rho_q*u*u + rho*2._dp*u*u_q + p_q
    flux_q(3,1,:) = rho_q*u*v + rho*(u_q*v + u*v_q)
    flux_q(4,1,:) = u_q*(q(4) + p) + u*(E_q + p_q)

    flux(1,2) = rho*v
    flux(2,2) = rho*u*v
    flux(3,2) = rho*v*v + p
    flux(4,2) = v*(q(4) + p)
    
    flux_q(1,2,:) = rho_q*v + rho*v_q
    flux_q(2,2,:) = rho_q*u*v + rho*(u_q*v + u*v_q)
    flux_q(3,2,:) = rho_q*v*v + rho*2._dp*v*v_q + p_q
    flux_q(4,2,:) = v_q*(q(4) + p) + v*(E_q + p_q)
    

end subroutine getNativeFlux_q

subroutine getNativeFlux(q,flux)
    real(dp),intent(in) :: q(numFields) 
    real(dp),intent(out):: flux(numFields,2)
    real(dp) :: rho,u,v,p,E

    real(dp) :: gamma
    gamma = 1.4_dp

    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
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


subroutine getNativeFluxComplex(q,flux)
#ifdef complx
    complex(dp),intent(in) :: q(numFields) 
    complex(dp),intent(out):: flux(numFields,2)
    complex(dp) :: rho,u,v,p,E
#else
    real(dp),intent(in) :: q(numFields) 
    real(dp),intent(out):: flux(numFields,2)
    real(dp) :: rho,u,v,p,E
#endif

    real(dp) :: gamma
    gamma = 1.4_dp

    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
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
    
end subroutine getNativeFluxComplex

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

!*****************************************************************************
!* -- Rotated-Roe-HLL Flux Function ---
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* The most robust Riemann solver known to the author (in terms of nonlinear
!* instability such as carbuncle).
!*
!* NB: At a boundary face, need to switch to a geometric normal vector:
!*               (nx2,ny2)=(nx, ny), (nx1,ny1)=(-ny,nx).
!*     This is not implemented here. It requires information on whether
!*     the geometric normal, (nx,ny), is on a boundary face or not.
!*     It shouldn't be difficult for you to implement it.
!*
!* Katate Masatsuka, February 2010. http://www.cfdbooks.com
!*****************************************************************************
 Subroutine Rotated_RHLL_q(qL, qR, nx, ny,flux,flux_ql,flux_qr)!,ul2,ul_ql2)
     use my_kinddefs
 real(dp),intent(in) :: qL(:), qR(:)    !  Input: conservative variables rho*[1, u, v, E]
 real(dp),intent(in) :: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(dp),intent(out) :: flux(:),flux_ql(:,:),flux_qr(:,:)!,ul2,ul_ql2(:) ! Output: Rotated_RHLL flux function.
!Local constants
 real(dp) :: gamma                          ! Ratio of specific heat.

 real(dp) :: eps                            ! 
!Local variables
 real(dp) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 real(dp) :: tx, ty                         ! Tangent vector (taken as n1)
 real(dp) :: alpha1, alpha2                 ! Projections of the new normals
 real(dp) :: uL, uR, vL, vR                 ! Velocity components.
 real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real(dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
 real(dp) :: RT,rho,u,v,H,a                 ! Roe-averages
 real(dp) :: vn, vt                         ! Normal and tangent velocities(Roe-average)
 real(dp) :: drho,dvn,dvt,dpres,dV(4)       ! Wave strengths
 real(dp) :: abs_dq                         ! Magnitude of the velocity difference
 real(dp) :: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigenvectors
 real(dp) :: SRp,SLm                        ! Wave speeds for the HLL part
 real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 real(dp) :: temp,mag,nxT,nyT
 real(dp) :: wsOld(4)
 
 real(dp) :: rhoL_ql(4),rhoR_qr(4)
 real(dp) :: ul_ql(4),ur_qr(4),vl_ql(4),vr_qr(4),pl_ql(4),pr_qr(4),al_ql(4),ar_qr(4)
 real(dp) :: Hl_ql(4),Hr_qr(4),vnL_ql(4),vnR_qr(4),fL_ql(4,4),fR_qr(4,4)
 real(dp) :: nx1_ql(4),nx1_qr(4),ny1_ql(4),ny1_qr(4),abs_dq_ql(4),abs_dq_qr(4)
 real(dp) :: alpha1_ql(4),alpha1_qr(4),nx2_ql(4),nx2_qr(4),ny2_ql(4),ny2_qr(4)
 real(dp) :: alpha2_ql(4),alpha2_qr(4),RT_ql(4),RT_qr(4),rho_ql(4),rho_qR(4)
 real(dp) :: u_ql(4),u_qr(4),v_ql(4),v_qr(4),H_ql(4),H_qr(4),a_ql(4),a_qr(4)
 real(dp) :: vn_ql(4),vn_qr(4),vt_ql(4),vt_qr(4), vtL_ql(4), vtR_qr(4)
 real(dp) :: drho_ql(4), drho_qr(4), dpres_ql(4), dpres_qr(4)
 real(dp) :: dvn_ql(4),dvn_qr(4),dV_ql(4,4),dV_qr(4,4),dvt_ql(4),dvt_qr(4)
 real(dp) :: ws_ql(4,4),ws_qr(4,4),abs_ws_ql(4,4),abs_ws_qr(4,4)
 real(dp) :: SRp_ql(4),SRp_qr(4),SLm_ql(4),SLm_qr(4),topLPrime_ql(4,4),topLPrime_qr(4,4),topL(4),bottom
 real(dp) :: bottomPrime_ql(4),bottomPrime_qr(4),Rv_ql(4,4,4),Rv_qr(4,4,4)
 real(dp) :: diss_ql(4,4),diss_qr(4,4),dws_ql(4,4),dws_qr(4,4),val1,val2
 real(dp) :: vnL_qr(4),vnR_ql(4),vtL_qr(4),vtR_ql(4),tx_ql(4),tx_qr(4),ty_ql(4),ty_qr(4)
 real(dp) :: topRPrime_ql(4),topRPrime_qr(4),topR
 

 integer(i4) :: i, j

    mag = sqrt( nx**2 + ny**2)  
    nxT = nx/mag
    nyT = ny/mag
!Constants.
    gamma = 1.4_dp
    eps = 1.0_dp*10._dp**(-5._dp) ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = qL(1)
    rhoL_ql(1) = 1._dp
    rhoL_ql(2) = 0._dp
    rhoL_ql(3) = 0._dp
    rhoL_ql(4) = 0._dp
    
    uL          = qL(2)/qL(1)  
    ul_ql(1)    = -ql(2)/(ql(1)**2)         !du/drho  = -u/rho
    ul_ql(2)    = 1._dp/ql(1)               !du/drhou = 1/rho
    ul_ql(3:4)  = 0._dp            

    
    vL = qL(3)/qL(1) 
    vl_ql(1) = -ql(3)/(ql(1)**2)        !dv/drho  = -v/rho
    vl_ql(2) = 0._dp                    !dv/drhou = 0
    vl_ql(3) = 1._dp/ql(1)              !dv/drhov = 1/rho
    vl_ql(4) = 0._dp                    !dv/drhoE = 0
      

    
    pL = (gamma-1.0_dp)*( qL(4) - 0.5_dp*rhoL*(uL*uL+vL*vL) )
    pl_ql(1)    = half*(gamma - 1._dp)*((ql(2)/ql(1))**2 + (ql(3)/ql(1))**2)        !dp/drho 
    pl_ql(2:3)  = -(gamma - 1._dp)*(ql(2:3)/ql(1))                                  !dp/drhou  
    pl_ql(4)    =  (gamma - 1._dp)                                                  !dp/drhoE 
    
  
    aL = sqrt(gamma*pL/rhoL)
    al_ql(1)    = 0.5_dp*(gamma*pl/ql(1))**(-0.5_dp) * gamma*(ql(1)*pl_ql(1) - pl)/(ql(1))**2   !da/drho 
    al_ql(2:4)  = 0.5_dp*(gamma*pl/ql(1))**(-0.5_dp) * gamma*(ql(1)*pl_ql(2:4))   /(ql(1))**2    
    
   
    HL = ( qL(4) + pL ) / rhoL
    Hl_ql(1)    = (ql(1)*pl_ql(1) - (ql(4) +pL))/(ql(1)**2) 
    Hl_ql(2:3)  = pl_ql(2:3)/ql(1)
    Hl_ql(4)    = (1._dp + pl_ql(4))/ql(1)
    

!  Right state
    rhoR = qR(1)
    rhoR_qr(1) = 1._dp
    rhoR_qr(2) = 0._dp
    rhoR_qr(3) = 0._dp
    rhoR_qr(4) = 0._dp
    
    uR = qR(2)/qR(1)
    ur_qr(1)    = -qr(2)/(qr(1)**2)         !du/drho  = -u/rho
    ur_qr(2)    = 1._dp/qr(1)               !du/drhou = 1/rho
    ur_qr(3:4)  = 0._dp            
    
    vR = qR(3)/qR(1)     
    vr_qr(1) = -qr(3)/(qr(1)**2)        !dv/drho  = -v/rho
    vr_qr(2) = 0._dp                    !dv/drhou = 0
    vr_qr(3) = 1._dp/qr(1)              !dv/drhov = 1/rho
    vr_qr(4) = 0._dp                    !dv/drhoE = 0
    
    
    pR = (gamma-1.0_dp)*( qR(4) - 0.5_dp*rhoR*(uR*uR + vR*vR) )
    pr_qr(1)    = half*(gamma - 1._dp)*((qr(2)/qr(1))**2 + (qr(3)/qr(1))**2)        !dp/drho  
    pr_qr(2:3)  = -(gamma - 1._dp)*(qr(2:3)/qr(1))                                    
    pr_qr(4)    =  (gamma - 1._dp)                                                  !dE/drhoE
    
    
    aR = sqrt(gamma*pR/rhoR)
    ar_qr(1)    = 0.5_dp*(gamma*pr/qr(1))**(-0.5_dp) * gamma*(qr(1)*pr_qr(1) - pr)/(qr(1))**2   !da/drho 
    ar_qr(2:4)  = 0.5_dp*(gamma*pr/qr(1))**(-0.5_dp) * gamma*(qr(1)*pr_qr(2:4))   /(qr(1))**2
    
    
    HR = ( qR(4) + pR ) / rhoR
    Hr_qr(1)    = (qr(1)*pr_qr(1) - (qr(4) + pr))/(qr(1)**2) 
    Hr_qr(2:3)  = pr_qr(2:3)/qr(1)
    Hr_qr(4)    = (1._dp + pr_qr(4))/qr(1)

    
    vnL = uL*nxT + vL*nyT
    vnL_ql(1) = -ql(2)*nxT/(ql(1)**2) -ql(3)*nyT/(ql(1)**2)
    vnL_ql(2) = nxT/ql(1)
    vnL_ql(3) = nyT/ql(1)
    vnL_ql(4) = 0._dp
    
    vnR = uR*nxT + vR*nyT
    vnR_qr(1) = -qr(2)*nxT/(qr(1)**2) -qr(3)*nyT/(qr(1)**2)
    vnR_qr(2) = nxT/qr(1)
    vnR_qr(3) = nyT/qr(1)
    vnR_qr(4) = 0._dp
    
        

!Compute the flux.
    fL(1) = rhoL*vnL
    fL_ql(1,1)   = rhoL_ql(1)*vnL + rhoL*vnL_ql(1)
    fL_ql(1,2:4) = ql(1)*vnL_ql(2:4)
    

    fL(2) = rhoL*vnL * uL + pL*nxT
    fL_ql(2,1:4) = ql(2)*vnL_ql(1:4) + pl_ql(1:4)*nxT
    fL_ql(2,2)   = fL_ql(2,2) + vnL
    

    fL(3) = rhoL*vnL * vL + pL*nyT
    fL_ql(3,1:4) = ql(3)*vnL_ql(1:4) + pl_ql(1:4)*nyT
    fL_ql(3,3)   = fL_ql(3,3) + vnL


    fL(4) = rhoL*vnL *  HL
    fL_ql(4,1:4) = vnL_ql(1:4)*(ql(4) + pl) + vnL*pl_ql(1:4)
    fL_ql(4,4)   = fL_ql(4,4) + vnL

    
    fR(1) = rhoR*vnR
    fR_qr(1,1)   = vnR + qr(1)*vnR_qr(1)
    fR_qr(1,2:4) = qr(1)*vnR_qr(2:4)
    

    fR(2) = rhoR*vnR * uR + pR*nxT
    fR_qr(2,1:4) = qr(2)*vnR_qr(1:4) + pr_qr(1:4)*nxT
    fR_qr(2,2)   = fR_qr(2,2) + vnR

    fR(3) = rhoR*vnR * vR + pR*nyT
    fR_qr(3,1:4) = qr(3)*vnR_qr(1:4) + pr_qr(1:4)*nyT
    fR_qr(3,3)   = fR_qr(3,3) + vnR

    fR(4) = rhoR*vnR *  HR
    fR_qr(4,1:4) = vnR_qr(1:4)*(qr(4) + pr) + vnR*pr_qr(1:4)
    fR_qr(4,4)   = fR_qr(4,4) + vnR
    

   
!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (uR - uL)**2 + (vR - vL)**2 )
    abs_dq_ql(1)    = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)* &
                      ((uR - uL)*(ql(2)/ql(1)**2) + (vR - vL)*(ql(3)/ql(1)**2))
                      
    abs_dq_ql(2)  = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)*(uR - uL)*(-1._dp/ql(1))
    abs_dq_ql(3)  = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)*(vR - vL)*(-1._dp/ql(1))
    abs_dq_ql(4)    = 0._dp
    
    abs_dq_qr(1)    = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)* &
                      ((uR - uL)*(-qr(2)/qr(1)**2) + (vR - vL)*(- qr(3)/qr(1)**2))
                      
    abs_dq_qr(2)  = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)*(uR - uL)*(1._dp/qR(1))
    abs_dq_qr(3)  = ((uR - uL)**2 + (vR - vL)**2)**(-0.5_dp)*(vR - vL)*(1._dp/qR(1))
    abs_dq_qr(4)    = 0._dp
    

    
  if ( abs_dq > eps) then
    nx1 = (uR - uL)/abs_dq

    nx1_ql(1) = (abs_dq * ql(2)/ql(1)**2 - ((uR - uL)* abs_dq_ql(1)))/abs_dq**2
    nx1_ql(2) = (-abs_dq * (1._dp/ql(1)) - ((uR - uL)* abs_dq_ql(2)))/abs_dq**2
    nx1_ql(3) = -(uR - uL)* abs_dq_ql(3)/abs_dq**2
    nx1_ql(4) = 0._dp

    nx1_qr(1) = (-abs_dq * qr(2)/qr(1)**2 - ((uR - uL)* abs_dq_qr(1)))/abs_dq**2
    nx1_qr(2) = (abs_dq * (1._dp/qr(1))   - ((uR - uL)* abs_dq_qr(2)))/abs_dq**2
    nx1_qr(3) = -(uR - uL)* abs_dq_qr(3)/abs_dq**2
    nx1_qr(4) = 0._dp

    ny1 = (vR - vL)/abs_dq

    ny1_ql(1) = (abs_dq * ql(3)/ql(1)**2 - ((vR - vL)* abs_dq_ql(1)))/abs_dq**2
    ny1_ql(2) = -(vR - vL)* abs_dq_ql(2)/abs_dq**2
    ny1_ql(3) =  (-abs_dq * (1._dp/ql(1)) - ((vR - vL)* abs_dq_ql(3)))/abs_dq**2
    ny1_ql(4) = 0._dp

    ny1_qr(1) = (-abs_dq * qr(3)/qr(1)**2 - ((vR - vL)* abs_dq_qr(1)))/abs_dq**2
    ny1_qr(2) = -(vR - vL)* abs_dq_qr(2)/abs_dq**2
    ny1_qr(3) = (abs_dq * (1._dp/qr(1))   - ((vR - vL)* abs_dq_qr(3)))/abs_dq**2
    ny1_qr(4) = 0._dp

  else
    nx1 = -nyT
    nx1_ql(1:4) = 0._dp
    nx1_qr(1:4) = 0._dp
    
    ny1 =  nxT
    ny1_ql(1:4) = 0._dp
    ny1_qr(1:4) = 0._dp
  endif

    alpha1 = nxT * nx1 + nyT * ny1 
    alpha1_ql(1:4) = nxT * nx1_ql(1:4) + nyT * ny1_ql(1:4)
    alpha1_qr(1:4) = nxT * nx1_qr(1:4) + nyT * ny1_qr(1:4)


    !   To make alpha1 always positive.
    temp      = sign(1.0_dp,alpha1)

    nx1      = temp * nx1
    nx1_ql   = temp * nx1_ql
    nx1_qr   = temp * nx1_qr

    ny1      = temp * ny1
    ny1_ql   = ny1_ql * temp
    ny1_qr   = ny1_qr * temp
    
    
    alpha1 = temp * alpha1
    alpha1_ql = alpha1_ql * temp
    alpha1_qr = alpha1_qr * temp
    

! Take n2 as perpendicular to n1.
    nx2 = -ny1
    nx2_ql = -ny1_ql
    nx2_qr = -ny1_qr

    ny2 =  nx1
    ny2_ql = nx1_ql
    ny2_qr = nx1_qr
    
 
    alpha2 = nxT * nx2 + nyT * ny2
    alpha2_ql(1:4) = nxT * nx2_ql(1:4) + nyT * ny2_ql(1:4)
    alpha2_qr(1:4) = nxT * nx2_qr(1:4) + nyT * ny2_qr(1:4)
    

!   To make alpha2 always positive.
    temp = sign(1.0_dp,alpha2)
    nx2 = temp * nx2
    nx2_ql   = temp * nx2_ql
    nx2_qr   = temp * nx2_qr
       
    ny2 = temp * ny2
    ny2_ql   = temp * ny2_ql
    ny2_qr   = temp * ny2_qr
    
    alpha2 = temp * alpha2
    alpha2_ql = alpha2_ql * temp
    alpha2_qr = alpha2_qr * temp
    

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tangent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
    RT_ql(1)   = 0.5_dp*(rhoR/rhoL)**(-0.5_dp) * (-rhoR/rhoL**2)
    RT_ql(2:4) = 0._dp
    RT_qr(1)   = 0.5_dp*(rhoR/rhoL)**(-0.5_dp) * (1._dp/rhoL)
    RT_qr(2:4) = 0._dp
    

    rho = RT * rhoL
    rho_ql(1:4) = RT_ql(1:4)*rhoL
    rho_ql(1) = rho_ql(1) + RT
    rho_qr(1:4) = RT_qr(1:4)*rhoL
          
    
    u = (uL + RT*uR)/(1.0_dp+RT)
    u_ql(1)   = ((1.0_dp+RT)*(-ql(2)/ql(1)**2 + RT_ql(1)*uR) - ((uL + RT*uR)*RT_ql(1)))/((1.0_dp+RT)**2)
    u_ql(2)   = ((1.0_dp+RT)*( 1._dp/ql(1) + RT_ql(2)*uR) - ((uL + RT*uR)*RT_ql(2)))/((1.0_dp+RT)**2)
    u_ql(3:4) = ((1.0_dp+RT)*(  RT_ql(3:4)*uR) - ((uL + RT*uR)*RT_ql(3:4)))/((1.0_dp+RT)**2)
    
    u_qr(1)   = ((1.0_dp+RT)*(RT_qr(1)*uR -RT*qr(2)/(qr(1)**2)) - ((uL + RT*uR)*RT_qr(1)))/((1.0_dp+RT)**2)
    u_qr(2)   = ((1.0_dp+RT)*(RT_qr(2)*uR + RT/qr(1)) - ((uL + RT*uR)*RT_qr(2)))/((1.0_dp+RT)**2)
    u_qr(3:4) = ((1.0_dp+RT)*(RT_qr(3:4)*uR) - ((uL + RT*uR)*RT_qr(3:4)))/((1.0_dp+RT)**2)
    
    
    v = (vL + RT*vR)/(1.0_dp+RT)
    v_ql(1) = ((1.0_dp+RT)*(-ql(3)/ql(1)**2 + RT_ql(1)*vR) - ((vL + RT*vR)*RT_ql(1)))/((1.0_dp+RT)**2)
    v_ql(2) = ((1.0_dp+RT)*(  RT_ql(2)*vR) - ((vL + RT*vR)*RT_ql(2)))/((1.0_dp+RT)**2)
    v_ql(3) = ((1.0_dp+RT)*(1._dp/ql(1) + RT_ql(3)*vR) - ((vL + RT*vR)*RT_ql(3)))/((1.0_dp+RT)**2)
    v_ql(4) = ((1.0_dp+RT)*(  RT_ql(4)*vR) - ((vL + RT*vR)*RT_ql(4)))/((1.0_dp+RT)**2)
     
    v_qr(1) = ((1.0_dp+RT)*(RT_qr(1)*vR -RT*qr(3)/(qr(1)**2)) - ((vL + RT*vR)*RT_qr(1)))/((1.0_dp+RT)**2)
    v_qr(2) = ((1.0_dp+RT)*(RT_qr(2)*uR) - ((uL + RT*uR)*RT_qr(2)))/((1.0_dp+RT)**2)
    v_qr(3) = ((1.0_dp+RT)*(RT_qr(3)*uR + RT/qr(1)) - ((uL + RT*uR)*RT_qr(3)))/((1.0_dp+RT)**2)
    v_qr(4) = ((1.0_dp+RT)*(RT_qr(4)*uR) - ((uL + RT*uR)*RT_qr(4)))/((1.0_dp+RT)**2)
   
    
    H = ( HL + RT*HR)/(1.0_dp+RT)
    H_ql(1:4) = ((1.0_dp+RT)*( HL_ql(1:4) + RT_ql(1:4)*HR) - (( HL + RT*HR)*(RT_ql(1:4))))/((1.0_dp+RT)**2)
    H_qr(1:4) = ((1.0_dp+RT)*(RT_qr(1:4)*HR + RT*HR_qr(1:4))-( HL + RT*HR)*(RT_qr(1:4)))/ &
                ((1.0_dp+RT)**2)
                        
                
    a = sqrt( (gamma-1.0_dp)*(H-0.5_dp*(u*u + v*v)) )
    a_ql(1:4) = 0.5_dp*((gamma-1.0_dp)*(H-0.5_dp*(u*u + v*v)))**(-0.5_dp) * &
              ((gamma-1.0_dp)*(H_ql(1:4) - (u*u_ql(1:4) + v*v_ql(1:4))))

    a_qr(1:4) = 0.5_dp*((gamma-1.0_dp)*(H-0.5_dp*(u*u + v*v)))**(-0.5_dp) * &
              ((gamma-1.0_dp)*(H_qr(1:4) - (u*u_qr(1:4) + v*v_qr(1:4))))
               
    
    vn = u*nx2 + v*ny2
    vn_ql(1:4) = u_ql(1:4)*nx2 + u*nx2_ql(1:4) + v_ql(1:4)*ny2 + v*ny2_ql(1:4)  
    vn_qr(1:4) = u_qr(1:4)*nx2 + u*nx2_qr(1:4) + v_qr(1:4)*ny2 + v*ny2_qr(1:4)
    
    vt = u*nx1 + v*ny1
    vt_ql(1:4) = u_ql(1:4)*nx1 + u*nx1_ql(1:4) + v_ql(1:4)*ny1 + v*ny1_ql(1:4) 
    vt_qr(1:4) = u_qr(1:4)*nx1 + u*nx1_qr(1:4) + v_qr(1:4)*ny1 + v*ny1_qr(1:4)
    

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = uL*nx2 + vL*ny2
    vnL_ql(1:4) = ul_ql(1:4)*nx2 + uL*nx2_ql(1:4) + vL_ql(1:4)*ny2 + vL*ny2_ql(1:4)
    vnL_qr(1:4) = uL*nx2_qr(1:4) + vL*ny2_qr(1:4)
    
    vnR = uR*nx2 + vR*ny2
    vnR_ql(1:4) = uR*nx2_ql(1:4) + vR*ny2_ql(1:4)
    vnR_qr(1:4) = ur_qr(1:4)*nx2 + uR*nx2_qr(1:4) + vR_qr(1:4)*ny2 + vR*ny2_qr(1:4)
    
    vtL = uL*nx1 + vL*ny1
    vtL_ql(1:4) = ul_ql(1:4)*nx1 + uL*nx1_ql(1:4) + vL_ql(1:4)*ny1 + vL*ny1_ql(1:4)
    vtL_qr(1:4) = uL*nx1_qr(1:4) + vL*ny1_qr(1:4)
    
    vtR = uR*nx1 + vR*ny1
    vtR_ql(1:4) = uR*nx1_ql(1:4) + vR*ny1_ql(1:4)
    vtR_qr(1:4) = ur_qr(1:4)*nx1 + uR*nx1_qr(1:4) + vR_qr(1:4)*ny1 + vR*ny1_qr(1:4) 
    
    
    drho         = rhoR - rhoL 
    drho_ql(1:4) = -rhoL_ql(1:4)
    drho_qr(1:4) =  rhoR_qr(1:4)

    dpres         =  pR - pL
    dpres_ql(1:4) = -pL_ql(1:4)
    dpres_qr(1:4) =  pR_qr(1:4)
 
    
    dvn         = vnR  - vnL
    dvn_ql(1:4) = vnR_ql(1:4) - vnL_ql(1:4)
    dvn_qr(1:4) = vnR_qr(1:4) - vnL_qr(1:4)
    
    dvt         = vtR  - vtL
    dvt_ql(1:4) = vtR_ql(1:4) - vtL_ql(1:4)
    dvt_qr(1:4) = vtR_qr(1:4) - vtL_qr(1:4)

     
    dV(1) = (dpres - rho*a*dvn )/(2.0_dp*a*a)
    dV_ql(1,1:4) = ((2.0_dp*a*a)*dpres_ql(1:4) - (rho_ql(1:4)*a*dvn + rho*(a_ql(1:4)*dvn + &
                    a*dvn_ql(1:4))) - (dpres - rho*a*dvn)*(4._dp*a*a_ql(1:4)))/((2.0_dp*a*a)**2)
    dV_qr(1,1:4) = ((2.0_dp*a*a)*dpres_qr(1:4) - (rho_qr(1:4)*a*dvn + rho*(a_qr(1:4)*dvn + &
                    a*dvn_qr(1:4))) - (dpres - rho*a*dvn)*(4._dp*a*a_qr(1:4)))/((2.0_dp*a*a)**2)
                                                 
    dV(2) =  rho*dvt/a
    dV_ql(2,1:4) = (a*(rho_ql(1:4)*dvt + rho*dvt_ql(1:4)) -(rho*dvt)*(a_ql(1:4)))/a**2
    dV_qr(2,1:4) = (a*(rho_qr(1:4)*dvt + rho*dvt_qr(1:4)) -(rho*dvt)*(a_qr(1:4)))/a**2
    
    dV(3) =  drho - dpres/(a*a)
    dV_ql(3,1:4) = drho_ql(1:4) - ((a*a)*dpres_ql(1:4) - dpres*2._dp*a*a_ql(1:4))/((a*a)**2)
    dV_qr(3,1:4) = drho_qr(1:4) - ((a*a)*dpres_qr(1:4) - dpres*2._dp*a*a_qr(1:4))/((a*a)**2)
    
    dV(4) = (dpres + rho*a*dvn )/(2.0_dp*a*a)
    dV_ql(4,1:4) = ((2.0_dp*a*a)*(dpres_ql(1:4) + rho_ql(1:4)*(a*dvn) + rho*(a_ql(1:4)*dvn + &
                    a*dvn_ql(1:4))) - (dpres + rho*a*dvn)*(4._dp*a*a_ql(1:4))) / ((2._dp*a*a)**2)
    dV_qr(4,1:4) = ((2.0_dp*a*a)*(dpres_qr(1:4) + rho_qr(1:4)*(a*dvn) + rho*(a_qr(1:4)*dvn + &
                    a*dvn_qr(1:4))) - (dpres + rho*a*dvn)*(4._dp*a*a_qr(1:4))) / ((2._dp*a*a)**2)
                                    

!Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws_ql(1,1:4) = vn_ql(1:4) - a_ql(1:4)
    ws_qr(1,1:4) = vn_qr(1:4) - a_qr(1:4)
    
    ws(2) = vn
    ws_ql(2,1:4) = vn_ql(1:4)
    ws_qr(2,1:4) = vn_qr(1:4)
    
    ws(3) = vn
    ws_ql(3,1:4) = vn_ql(1:4)
    ws_qr(3,1:4) = vn_qr(1:4)
    
    ws(4) = vn+a
    ws_ql(4,1:4) = vn_ql(1:4) + a_ql(1:4)
    ws_qr(4,1:4) = vn_qr(1:4) + a_qr(1:4)
    
    
    abs_ws  = abs(ws)
    do i = 1,4
        abs_ws_ql(i,:) =  my_smooth_dabs(ws(i)) * ws_ql(i,:)
        abs_ws_qr(i,:) =  my_smooth_dabs(ws(i)) * ws_qr(i,:)
    end do
    
    
!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
    dws(1) = 0.2_dp
    dws_ql(1,:) = 0._dp
    dws_qr(1,:) = 0._dp
    
        
   if (abs_ws(1)<dws(1)) then
       abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1) + dws(1))
       abs_ws_ql(1,:) = half*((dws(1)*(2._dp*abs(ws(1))*my_smooth_dabs(ws(1))*ws_ql(1,:))-&
                        (abs_ws(1)*abs_ws(1)*dws_ql(1,:)))&
                             /(dws(1)**2) + dws_ql(1,:))
       abs_ws_qr(1,:) = half*((dws(1)*(2._dp*abs(ws(1))*my_smooth_dabs(ws(1))*ws_qr(1,:))-&
                        (abs_ws(1)*abs_ws(1)*dws_qr(1,:)))&
                             /(dws(1)**2) + dws_qr(1,:))
   end if
          

   
    dws(4) = 0.2_dp
    dws_ql(4,:) = 0._dp
    dws_qr(4,:) = 0._dp
    
   if (abs_ws(4)<dws(4)) then
       abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4) + dws(4))
       abs_ws_ql(4,:) = half*((dws(4)*(2._dp*abs(ws(4))*my_smooth_dabs(ws(4))*ws_ql(4,:))-&
                        (abs_ws(4)*abs_ws(4)*dws_ql(4,:)))&
                             /(dws(4)**2) + dws_ql(4,:))
       abs_ws_qr(4,:) = half*((dws(4)*(2._dp*abs(ws(4))*my_smooth_dabs(ws(4))*ws_qr(4,:))-&
                        (abs_ws(4)*abs_ws(4)*dws_qr(4,:)))&
                             /(dws(4)**2) + dws_qr(4,:))
   end if


!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    SRp = max( 0.0_dp, vtR + aR, vt + a)

    val1 = vtR + aR
    val2 = vt + a
   
    if (SRp == 0.0_dp)then
        SRp_ql(:) = 0._dp
        SRp_qr(:) = 0._dp
    else if(SRp == val1)then
        SRp_ql(:) = vtR_ql(:) 
        SRp_qr(:) = vtR_qr(:) + aR_qr(:)
    else if(Srp == val2)then
        SRp_ql(:) = vt_ql(:) + a_ql(:)
        SRp_qr(:) = vt_qr(:) + a_qr(:)
    end if
    
   
    SLm = min( 0.0_dp, vtL - aL, vt - a)

    val1 = vtL - aL
    val2 = vt - a
    
    if (SLm == 0.0_dp)then
        SLm_ql(:) = 0._dp
        SLm_qr(:) = 0._dp
    else if(SLm == val1)then
        SLm_ql(:) = vtL_ql(:) - aL_ql(:)
        SLm_qr(:) = vtL_qr(:)
    else if(SLm == val2)then
        SLm_ql(:) = vt_ql(:) - a_ql(:)
        SLm_qr(:) = vt_qr(:) - a_qr(:)
    end if
    
       

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   wsOld(:) = ws(:)
   ws  =  alpha2*abs_ws - ( alpha2*(SRp+SLm)*wsOld + 2.0_dp*alpha1*SRp*SLm )/ (SRp-SLm)
   
   topL =  alpha2*(SRp+SLm)*wsOld
   topR =  2.0_dp*alpha1*SRp*SLm
   
   
   do i = 1,4
    topLPrime_ql(i,:) = alpha2_ql*(SRp+SLm)*wsOld(i)+alpha2*((SRp_ql + SLm_ql)*wsOld(i)+(SRp+SLm)*ws_ql(i,:)) 
    topLPrime_qr(i,:) = alpha2_qr*(SRp+SLm)*wsOld(i)+alpha2*((SRp_qr + SLm_qr)*wsOld(i)+(SRp+SLm)*ws_qr(i,:))
   end do
   topRPrime_ql(:) = 2.0_dp*(alpha1_ql*SRp*SLm + alpha1*(SRp_ql*SLm + SRp*SLm_ql))
   topRPrime_qr(:) = 2.0_dp*(alpha1_qr*SRp*SLm + alpha1*(SRp_qr*SLm + SRp*SLm_qr))


   bottom = (SRp-SLm)
   bottomPrime_ql(:) = SRp_ql(:) - SLm_ql(:)
   bottomPrime_qr(:) = SRp_qr(:) - SLm_qr(:)
   
   do i = 1,4
    ws_ql(i,:) = alpha2_ql*abs_ws(i) + alpha2*abs_ws_ql(i,:) - &
                ( bottom *(topLPrime_ql(i,:) + topRPrime_ql) - (topL(i) + topR)*bottomPrime_ql)/bottom**2
              
    ws_qr(i,:) = alpha2_qr*abs_ws(i) + alpha2*abs_ws_qr(i,:) - &
                ( bottom *(topLPrime_qr(i,:) + topRPrime_qr) - (topL(i) + topR)*bottomPrime_qr)/bottom**2
   end do
   

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1
  tx_ql = nx1_ql 
  tx_qr = nx1_qr
  ty_ql = ny1_ql
  ty_qr = ny1_qr
  
  Rv(1,1) = 1.0_dp
  Rv_ql(1,1,:) = 0._dp
  Rv_qr(1,1,:) = 0._dp
  
  Rv(2,1) = u - a*nx2
  Rv_ql(2,1,:) = u_ql(:) - a_ql(:)*nx2 - a*nx2_ql(:)
  Rv_qr(2,1,:) = u_qr(:) - a_qr(:)*nx2 - a*nx2_qr(:)
  
  Rv(3,1) = v - a*ny2
  Rv_ql(3,1,:) = v_ql(:) - a_ql(:)*ny2 - a*ny2_ql(:)
  Rv_qr(3,1,:) = v_qr(:) - a_qr(:)*ny2 - a*ny2_qr(:)
  
  Rv(4,1) = H - vn*a
  Rv_ql(4,1,:) = H_ql(:) - a_ql(:)*vn - a*vn_ql(:)
  Rv_qr(4,1,:) = H_qr(:) - a_qr(:)*vn - a*vn_qr(:)
  
  Rv(1,2) = 0.0_dp
  Rv_ql(1,2,:) = 0._dp
  Rv_qr(1,2,:) = 0._dp
  
  Rv(2,2) = a*tx
  Rv_ql(2,2,:) = a_ql*tx + a*tx_ql
  Rv_qr(2,2,:) = a_qr*tx + a*tx_qr
  
  Rv(3,2) = a*ty
  Rv_ql(3,2,:) = a_ql*ty + a*ty_ql
  Rv_qr(3,2,:) = a_qr*ty + a*ty_qr
  
  Rv(4,2) = a*vt
  Rv_ql(4,2,:) = a_ql*vt + a*vt_ql
  Rv_qr(4,2,:) = a_qr*vt + a*vt_qr

  Rv(1,3) = 1.0_dp
  Rv_ql(1,3,:) = 0._dp
  Rv_qr(1,3,:) = 0._dp
  
  Rv(2,3) = u
  Rv_ql(2,3,:) = u_ql(:)
  Rv_qr(2,3,:) = u_qr(:)
  
  Rv(3,3) = v 
  Rv_ql(3,3,:) = v_ql(:)
  Rv_qr(3,3,:) = v_qr(:)
  
  Rv(4,3) = 0.5_dp*(u*u + v*v)
  Rv_ql(4,3,:) = u*u_ql(:) + v*v_ql(:)
  Rv_qr(4,3,:) = u*u_qr(:) + v*v_qr(:)
  
  Rv(1,4) = 1.0_dp
  Rv_ql(1,4,:) = 0._dp
  Rv_qr(1,4,:) = 0._dp
  
  Rv(2,4) = u + a*nx2
  Rv_ql(2,4,:) = u_ql(:) + a_ql(:)*nx2 + a*nx2_ql(:)
  Rv_qr(2,4,:) = u_qr(:) + a_qr(:)*nx2 + a*nx2_qr(:)
  
  Rv(3,4) = v + a*ny2
  Rv_ql(3,4,:) = v_ql(:) + a_ql(:)*ny2 + a*ny2_ql(:)
  Rv_qr(3,4,:) = v_qr(:) + a_qr(:)*ny2 + a*ny2_qr(:)
  
  Rv(4,4) = H + vn*a
  Rv_ql(4,4,:) = H_ql(:) + a_ql(:)*vn + a*vn_ql(:)
  Rv_qr(4,4,:) = H_qr(:) + a_qr(:)*vn + a*vn_qr(:)
  

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = 0.0_dp
  diss_ql(:,:) = 0._dp
  diss_qr(:,:) = 0._dp
  
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
    diss_ql(i,:) = diss_ql(i,:) + ws_ql(j,:)*(dV(j)*Rv(i,j)) + ws(j)*(dV_ql(j,:)*Rv(i,j) + dV(j)*Rv_ql(i,j,:)) 
    diss_qr(i,:) = diss_qr(i,:) + ws_qr(j,:)*(dV(j)*Rv(i,j)) + ws(j)*(dV_qr(j,:)*Rv(i,j) + dV(j)*Rv_qr(i,j,:))
   end do
  end do
  
    

!Compute the Rotated-RHLL flux.
  flux(:) = (SRp*fL - SLm*fR)/(SRp-SLm) - 0.5_dp*diss
  do i = 1,4
    flux_ql(i,:) = (((SRp-SLm)*(SRp_ql*fL(i) + SRp*fL_ql(i,:) - (SLm_ql*fR(i) ))) - &
                   ((SRp*fL(i) - SLm*fR(i))*(SRp_ql - SLm_ql)))/(SRp-SLm)**2 - 0.5_dp*diss_ql(i,:)
    flux_qr(i,:) = (((SRp-SLm)*(SRp_qr*fL(i) - (SLm_qr(:)*fR(i) + SLm*fR_qr(i,:)))) - &
                   ((SRp*fL(i) - SLm*fR(i) )*(SRp_qr - SLm_qr)))/(SRp-SLm)**2 - 0.5_dp*diss_qr(i,:)
  end do
                   
  flux      = flux*mag
  flux_ql   = flux_ql*mag
  flux_qr   = flux_qr*mag
  
  !ul2    =  flux(1)
  !ul_ql2 =  flux_qr(1,:)
  
  end subroutine Rotated_RHLL_q
  
 Subroutine Rotated_RHLL(qL, qR, nx, ny,flux)
     use my_kinddefs

     real(dp),intent(in)  :: qL(:), qR(:)    !  Input: conservative variables rho*[1, u, v, E]
     real(dp),intent(in) :: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
     real(dp),intent(out) :: flux(:)  ! Output: Rotated_RHLL flux function.
     real(dp) :: uL, uR, vL, vR                 ! Velocity components.
     real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
     real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
     real(dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
     real(dp) :: RT,rho,u,v,H,a                 ! Roe-averages
     real(dp) :: vn, vt                         ! Normal and tangent velocities(Roe-average)
     real(dp) :: drho,dvn,dvt,dpres,dV(4)       ! Wave strengths
     real(dp) :: abs_dq                         ! Magnitude of the velocity difference
     real(dp) :: ws(4), Rv(4,4)! Wave speeds and right-eigenvectors
     real(dp) :: SRp,SLm                        ! Wave speeds for the HLL part
     real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
     real(dp) :: temp,mag,nxT,nyT
     real(dp) :: wsOld(4)
    

    
    !Local constants
     real(dp) :: gamma                          ! Ratio of specific heat.
     real(dp) :: eps,abs_ws(4),dws(4)       
     ! 
    !Local variables
     real(dp) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
     real(dp) :: tx, ty                         ! Tangent vector (taken as n1)
     real(dp) :: alpha1, alpha2                 ! Projections of the new normals
   
     
     


     integer(i4) :: i, j

    mag = sqrt( nx**2 + ny**2)  
    nxT = nx/mag
    nyT = ny/mag
!Constants.
    gamma = 1.4_dp
    eps = 1.0_dp*10._dp**(-5._dp) ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = qL(1)   
    uL          = qL(2)/qL(1)    
    vL = qL(3)/qL(1) 
    pL = (gamma-1.0_dp)*( qL(4) - 0.5_dp*rhoL*(uL*uL+vL*vL) )

    aL = sqrt(gamma*pL/rhoL)

    HL = ( qL(4) + pL ) / rhoL


!  Right state
    rhoR = qR(1)
    uR = qR(2)/qR(1)
    vR = qR(3)/qR(1)        
    pR = (gamma-1.0_dp)*( qR(4) - 0.5_dp*rhoR*(uR*uR + vR*vR) )

    aR = sqrt(gamma*pR/rhoR)

    HR = ( qR(4) + pR ) / rhoR

    vnL = uL*nxT + vL*nyT
 
    vnR = uR*nxT + vR*nyT
    

!Compute the flux.
    fL(1) = rhoL*vnL
    fL(2) = rhoL*vnL * uL + pL*nxT
    fL(3) = rhoL*vnL * vL + pL*nyT
    fL(4) = rhoL*vnL *  HL

    fR(1) = rhoR*vnR
    fR(2) = rhoR*vnR * uR + pR*nxT
    fR(3) = rhoR*vnR * vR + pR*nyT
    fR(4) = rhoR*vnR *  HR

   
!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (uR - uL)**2 + (vR - vL)**2 )
    
  if ( real(abs_dq) > eps) then
    nx1 = (uR - uL)/abs_dq
    ny1 = (vR - vL)/abs_dq
  else
    nx1 = -nyT  
    ny1 =  nxT
  endif

    alpha1 = nxT * nx1 + nyT * ny1 

    !   To make alpha1 always positive.
    temp      = sign(1.0_dp,alpha1)

    nx1      = temp * nx1
    ny1      = temp * ny1   
    alpha1 = temp * alpha1
    

! Take n2 as perpendicular to n1.
    nx2 = -ny1
    ny2 =  nx1
 
    alpha2 = nxT * nx2 + nyT * ny2

!   To make alpha2 always positive.
    temp = sign(1.0_dp,alpha2)
    nx2 = temp * nx2   
    ny2 = temp * ny2
    
    alpha2 = temp * alpha2
    

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tangent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
    
    rho = RT * rhoL          
    u = (uL + RT*uR)/(1.0_dp+RT)  
    v = (vL + RT*vR)/(1.0_dp+RT) 
    H = ( HL + RT*HR)/(1.0_dp+RT)
            
    a = sqrt( (gamma-1.0_dp)*(H-0.5_dp*(u*u + v*v)) )             
  
    vn = u*nx2 + v*ny2   
    vt = u*nx1 + v*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = uL*nx2 + vL*ny2
    
    vnR = uR*nx2 + vR*ny2  
    vtL = uL*nx1 + vL*ny1
    vtR = uR*nx1 + vR*ny1
    
    drho         = rhoR - rhoL 

    dpres         =  pR - pL
    dvn         = vnR  - vnL
    dvt         = vtR  - vtL

    dV(1) = (dpres - rho*a*dvn )/(2.0_dp*a*a)                                           
    dV(2) =  rho*dvt/a   
    dV(3) =  drho - dpres/(a*a)   
    dV(4) = (dpres + rho*a*dvn )/(2.0_dp*a*a)
                                    
!Wave Speeds for Roe flux part.
    ws(1) = vn-a  
    ws(2) = vn   
    ws(3) = vn    
    ws(4) = vn+a
  
    abs_ws  = abs(real(ws))
    
    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    !only for the nonlinear fields.
   dws(1) = 0.2_dp
           
   if (abs_ws(1)<dws(1)) then
       abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1) + dws(1))
   end if
 
   dws(4) = 0.2_dp
    
   if (abs_ws(4)<dws(4)) then
       abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4) + dws(4))
   end if


    !HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    SRp = max( 0.0_dp, real(vtR) + real(aR), real(vt) + real(a))
   
    SLm = min( 0.0_dp, real(vtL) - real(aL), real(vt) - real(a))
    
    !Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
    wsOld(:) = ws(:)
    ws  =  alpha2*abs_ws - ( alpha2*(SRp+SLm)*wsOld + 2.0_dp*alpha1*SRp*SLm )/ (SRp-SLm)

    !Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1
    ty = ny1

    Rv(1,1) = 1.0_dp 
    Rv(2,1) = u - a*nx2 
    Rv(3,1) = v - a*ny2
    Rv(4,1) = H - vn*a
    Rv(1,2) = 0.0_dp
    Rv(2,2) = a*tx
    Rv(3,2) = a*ty
    Rv(4,2) = a*vt
    Rv(1,3) = 1.0_dp
    Rv(2,3) = u
    Rv(3,3) = v 
    Rv(4,3) = 0.5_dp*(u*u + v*v)
    Rv(1,4) = 1.0_dp
    Rv(2,4) = u + a*nx2
    Rv(3,4) = v + a*ny2 
    Rv(4,4) = H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = 0.0_dp
  
  do i=1,4
    do j=1,4
        diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
    end do
  end do
  
!Compute the Rotated-RHLL flux.
  flux  = (SRp*fL - SLm*fR)/(SRp-SLm) - 0.5_dp*diss              
  flux  = flux*mag

  end subroutine Rotated_RHLL
  
  Subroutine Rotated_RHLLComplex(qL, qR, nx, ny,flux)
     use my_kinddefs
#ifdef complx
     complex(dp),intent(in)  :: qL(:), qR(:)    !  Input: conservative variables rho*[1, u, v, E]
     complex(dp),intent(out) :: flux(:)  ! Output: Rotated_RHLL flux function.
     complex(dp) :: uL, uR, vL, vR                 ! Velocity components.
     complex(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
     complex(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
     complex(dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
     complex(dp) :: RT,rho,u,v,H,a                 ! Roe-averages
     complex(dp) :: vn, vt                         ! Normal and tangent velocities(Roe-average)
     complex(dp) :: drho,dvn,dvt,dpres,dV(4)       ! Wave strengths
     complex(dp) :: abs_dq                         ! Magnitude of the velocity difference
     complex(dp) :: ws(4), Rv(4,4)! Wave speeds and right-eigenvectors
     complex(dp) :: SRp,SLm                        ! Wave speeds for the HLL part
     complex(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
     complex(dp) :: temp,mag,nxT,nyT
     complex(dp) :: wsOld(4)
#else
     real(dp),intent(in)  :: qL(:), qR(:)    !  Input: conservative variables rho*[1, u, v, E]
     real(dp),intent(out) :: flux(:)  ! Output: Rotated_RHLL flux function.
     real(dp) :: uL, uR, vL, vR                 ! Velocity components.
     real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
     real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
     real(dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
     real(dp) :: RT,rho,u,v,H,a                 ! Roe-averages
     real(dp) :: vn, vt                         ! Normal and tangent velocities(Roe-average)
     real(dp) :: drho,dvn,dvt,dpres,dV(4)       ! Wave strengths
     real(dp) :: abs_dq                         ! Magnitude of the velocity difference
     real(dp) :: ws(4), Rv(4,4)! Wave speeds and right-eigenvectors
     real(dp) :: SRp,SLm                        ! Wave speeds for the HLL part
     real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
     real(dp) :: temp,mag,nxT,nyT
     real(dp) :: wsOld(4)
#endif
    
     real(dp),intent(in) :: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
    
    !Local constants
     real(dp) :: gamma                          ! Ratio of specific heat.
     real(dp) :: eps,abs_ws(4),dws(4)       
     ! 
    !Local variables
     real(dp) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
     real(dp) :: tx, ty                         ! Tangent vector (taken as n1)
     real(dp) :: alpha1, alpha2                 ! Projections of the new normals
   
     
     


     integer(i4) :: i, j

    mag = sqrt( nx**2 + ny**2)  
    nxT = nx/mag
    nyT = ny/mag
!Constants.
    gamma = 1.4_dp
    eps = 1.0_dp*10._dp**(-5._dp) ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = qL(1)   
    uL          = qL(2)/qL(1)    
    vL = qL(3)/qL(1) 
    pL = (gamma-1.0_dp)*( qL(4) - 0.5_dp*rhoL*(uL*uL+vL*vL) )

    aL = sqrt(gamma*pL/rhoL)

    HL = ( qL(4) + pL ) / rhoL


!  Right state
    rhoR = qR(1)
    uR = qR(2)/qR(1)
    vR = qR(3)/qR(1)        
    pR = (gamma-1.0_dp)*( qR(4) - 0.5_dp*rhoR*(uR*uR + vR*vR) )

    aR = sqrt(gamma*pR/rhoR)

    HR = ( qR(4) + pR ) / rhoR

    vnL = uL*nxT + vL*nyT
 
    vnR = uR*nxT + vR*nyT
    

!Compute the flux.
    fL(1) = rhoL*vnL
    fL(2) = rhoL*vnL * uL + pL*nxT
    fL(3) = rhoL*vnL * vL + pL*nyT
    fL(4) = rhoL*vnL *  HL

    fR(1) = rhoR*vnR
    fR(2) = rhoR*vnR * uR + pR*nxT
    fR(3) = rhoR*vnR * vR + pR*nyT
    fR(4) = rhoR*vnR *  HR

   
!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (uR - uL)**2 + (vR - vL)**2 )
    
  if ( real(abs_dq) > eps) then
    nx1 = (uR - uL)/abs_dq
    ny1 = (vR - vL)/abs_dq
  else
    nx1 = -nyT  
    ny1 =  nxT
  endif

    alpha1 = nxT * nx1 + nyT * ny1 

    !   To make alpha1 always positive.
    temp      = sign(1.0_dp,alpha1)

    nx1      = temp * nx1
    ny1      = temp * ny1   
    alpha1 = temp * alpha1
    

! Take n2 as perpendicular to n1.
    nx2 = -ny1
    ny2 =  nx1
 
    alpha2 = nxT * nx2 + nyT * ny2

!   To make alpha2 always positive.
    temp = sign(1.0_dp,alpha2)
    nx2 = temp * nx2   
    ny2 = temp * ny2
    
    alpha2 = temp * alpha2
    

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tangent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
    
    rho = RT * rhoL          
    u = (uL + RT*uR)/(1.0_dp+RT)  
    v = (vL + RT*vR)/(1.0_dp+RT) 
    H = ( HL + RT*HR)/(1.0_dp+RT)
            
    a = sqrt( (gamma-1.0_dp)*(H-0.5_dp*(u*u + v*v)) )             
  
    vn = u*nx2 + v*ny2   
    vt = u*nx1 + v*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = uL*nx2 + vL*ny2
    
    vnR = uR*nx2 + vR*ny2  
    vtL = uL*nx1 + vL*ny1
    vtR = uR*nx1 + vR*ny1
    
    drho         = rhoR - rhoL 

    dpres         =  pR - pL
    dvn         = vnR  - vnL
    dvt         = vtR  - vtL

    dV(1) = (dpres - rho*a*dvn )/(2.0_dp*a*a)                                           
    dV(2) =  rho*dvt/a   
    dV(3) =  drho - dpres/(a*a)   
    dV(4) = (dpres + rho*a*dvn )/(2.0_dp*a*a)
                                    
!Wave Speeds for Roe flux part.
    ws(1) = vn-a  
    ws(2) = vn   
    ws(3) = vn    
    ws(4) = vn+a
  
    abs_ws  = abs(real(ws))
    
    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    !only for the nonlinear fields.
   dws(1) = 0.2_dp
           
   if (abs_ws(1)<dws(1)) then
       abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1) + dws(1))
   end if
 
   dws(4) = 0.2_dp
    
   if (abs_ws(4)<dws(4)) then
       abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4) + dws(4))
   end if


    !HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    SRp = max( 0.0_dp, real(vtR) + real(aR), real(vt) + real(a))
   
    SLm = min( 0.0_dp, real(vtL) - real(aL), real(vt) - real(a))
    
    !Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
    wsOld(:) = ws(:)
    ws  =  alpha2*abs_ws - ( alpha2*(SRp+SLm)*wsOld + 2.0_dp*alpha1*SRp*SLm )/ (SRp-SLm)

    !Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1
    ty = ny1

    Rv(1,1) = 1.0_dp 
    Rv(2,1) = u - a*nx2 
    Rv(3,1) = v - a*ny2
    Rv(4,1) = H - vn*a
    Rv(1,2) = 0.0_dp
    Rv(2,2) = a*tx
    Rv(3,2) = a*ty
    Rv(4,2) = a*vt
    Rv(1,3) = 1.0_dp
    Rv(2,3) = u
    Rv(3,3) = v 
    Rv(4,3) = 0.5_dp*(u*u + v*v)
    Rv(1,4) = 1.0_dp
    Rv(2,4) = u + a*nx2
    Rv(3,4) = v + a*ny2 
    Rv(4,4) = H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = 0.0_dp
  
  do i=1,4
    do j=1,4
        diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
    end do
  end do
  
!Compute the Rotated-RHLL flux.
  flux  = (SRp*fL - SLm*fR)/(SRp-SLm) - 0.5_dp*diss              
  flux  = flux*mag

  end subroutine Rotated_RHLLComplex
  
  subroutine lax_friedrichs_q(ql, qr, norm, flux,flux_ql,flux_qr)
    use my_kinddefs
    
    real(dp),intent(in)  :: ql(:)
    real(dp),intent(in)  :: qr(:)
    real(dp),intent(in)  :: norm(:)
    real(dp),intent(out) :: flux(:),flux_ql(:,:),flux_qr(:,:)

    !---> Local Variables
    real(dp) :: eigl, eigr,eigl_ql(4),eigr_qr(4)  
    real(dp) :: nx, ny
    real(dp) :: al, ar,al_ql(4),ar_qr(4)
    real(dp) :: ul, vl, pl, El
    real(dp) :: ul_ql(4), vl_ql(4), pl_ql(4), El_ql(4)
    real(dp) :: ur, vr, pr, Er
    real(dp) :: ur_qr(4), vr_qr(4), pr_qr(4), Er_qr(4)
    real(dp) :: flux_l, flux_l_ql(4)
    real(dp) :: flux_r, flux_r_qr(4)
    real(dp) :: alpha, alpha_ql(4),alpha_qr(4)
    real(dp) :: mag
    real(dp) :: gamma
    integer(i4) :: k
    
    gamma = 1.4_dp
    mag = sqrt( norm(1)**2 + norm(2)**2 )    

    nx = norm(1)/mag
    ny = norm(2)/mag
      
    ul = ql(2)/ql(1) 
    ur = qr(2)/qr(1)
    
    ul_ql(1)    = -ql(2)/(ql(1)**2)         !du/drho  = -u/rho
    ul_ql(2)    = 1._dp/ql(1)               !du/drhou = 1/rho
    ul_ql(3:4)  = 0._dp            
    
    ur_qr(1)    = -qr(2)/(qr(1)**2)         !du/drho  = -u/rho
    ur_qr(2)    = 1._dp/qr(1)               !du/drhou = 1/rho
    ur_qr(3:4)  = 0._dp            
   
    
    vl = ql(3)/ql(1)
    vr = qr(3)/qr(1)
    
    vl_ql(1) = -ql(3)/(ql(1)**2)        !dv/drho  = -v/rho
    vl_ql(2) = 0._dp                    !dv/drhou = 0
    vl_ql(3) = 1._dp/ql(1)              !dv/drhov = 1/rho
    vl_ql(4) = 0._dp                    !dv/drhoE = 0
    
    vr_qr(1) = -qr(3)/(qr(1)**2)        !dv/drho  = -v/rho
    vr_qr(2) = 0._dp                    !dv/drhou = 0
    vr_qr(3) = 1._dp/qr(1)              !dv/drhov = 1/rho
    vr_qr(4) = 0._dp                    !dv/drhoE = 0

    El = ql(4)/ql(1)
    Er = qr(4)/qr(1)
    
    El_ql(1)    = -ql(4)/(ql(1)**2)         !dE/drho 
    El_ql(2:3)  = 0._dp          
    El_ql(4)    = 1._dp/ql(1)               !dE/drhoE 
    
    Er_qr(1)    = -qr(4)/(qr(1)**2)         !dE/drho  
    Er_qr(2:3)  = 0._dp                     !dE/drhou 
    Er_qr(4)    = 1._dp/qr(1)               !dE/drhoE 
    

    pl = ql(1)*(gamma - 1._dp)*(El - half*( ul**2 + vl**2) ) 
    pr = qr(1)*(gamma - 1._dp)*(Er - half*( ur**2 + vr**2) ) 
    
    pl_ql(1)    = half*(gamma - 1._dp)*((ql(2)/ql(1))**2 + (ql(3)/ql(1))**2)        !dp/drho 
    pl_ql(2:3)  = -(gamma - 1._dp)*(ql(2:3)/ql(1))                                  !dp/drhou  
    pl_ql(4)    =  (gamma - 1._dp)                                                  !dp/drhoE 
    
    pr_qr(1)    = half*(gamma - 1._dp)*((qr(2)/qr(1))**2 + (qr(3)/qr(1))**2)        !dp/drho  
    pr_qr(2:3)  = -(gamma - 1._dp)*(qr(2:3)/qr(1))                                    
    pr_qr(4)    =  (gamma - 1._dp)                                                  !dE/drhoE 
    
   
    al = sqrt( gamma*pl/ql(1))
    ar = sqrt( gamma*pr/qr(1))

    al_ql(1)    = 0.5_dp*(gamma*pl/ql(1))**(-0.5_dp) * gamma*(ql(1)*pl_ql(1) - pl)/(ql(1))**2   !da/drho 
    al_ql(2:4)  = 0.5_dp*(gamma*pl/ql(1))**(-0.5_dp) * gamma*(ql(1)*pl_ql(2:4))   /(ql(1))**2    
    
    ar_qr(1)    = 0.5_dp*(gamma*pr/qr(1))**(-0.5_dp) * gamma*(qr(1)*pr_qr(1) - pr)/(qr(1))**2   !da/drho 
    ar_qr(2:4)  = 0.5_dp*(gamma*pr/qr(1))**(-0.5_dp) * gamma*(qr(1)*pr_qr(2:4))   /(qr(1))**2
    
    
    eigl = al + my_smooth_abs(ul*nx + vl*ny)
    eigr = ar + my_smooth_abs(ur*nx + vr*ny)
    
    eigl_ql(1:4) = al_ql(1:4) + my_smooth_dabs(ul*nx + vl*ny) * (ul_ql(1:4)*nx + vl_ql(1:4)*ny)
    eigr_qr(1:4) = ar_qr(1:4) + my_smooth_dabs(ur*nx + vr*ny) * (ur_qr(1:4)*nx + vr_qr(1:4)*ny)
    
    
    !alpha = max( eigl,eigr )
    if( eigl > eigr) then
       alpha = eigl
       alpha_ql(1:4) = eigl_ql(1:4)
       alpha_qr(1:4) = 0._dp
    else if( eigl <= eigr ) then
       alpha = eigr
       alpha_qr(1:4) = eigr_qr(1:4)
       alpha_ql(1:4) = 0._dp
    end if


    !---> Now we do the flux's 1 at a time
    !---> Flux 1 
    flux_l = ql(2)*nx + ql(3)*ny
    flux_r = qr(2)*nx + qr(3)*ny
    
    flux_l_ql(1) = 0._dp 
    flux_l_ql(2) = nx
    flux_l_ql(3) = ny
    flux_l_ql(4) = 0._dp
    
    flux_r_qr(1) = 0._dp 
    flux_r_qr(2) = nx
    flux_r_qr(3) = ny
    flux_r_qr(4) = 0._dp
   
    flux(1) = half * ( ( flux_l + flux_r) + alpha*( ql(1) - qr(1) ))
    
    flux_ql(1,1)    = half * (flux_l_ql(1)   + alpha_ql(1)  *( ql(1) - qr(1) ) + alpha)
    flux_ql(1,2:4)  = half * (flux_l_ql(2:4) + alpha_ql(2:4)*( ql(1) - qr(1) ))
    
    flux_qr(1,1)    = half * (flux_r_qr(1)   + alpha_qr(1)  *( ql(1) - qr(1) ) - alpha)
    flux_qr(1,2:4)  = half * (flux_r_qr(2:4) + alpha_qr(2:4)*( ql(1) - qr(1) ))
   
    !---> Flux 2 

    flux_l = (ql(2)*ul + pl)*nx + (ql(3)*ul)*ny
    flux_r = (qr(2)*ur + pr)*nx + (qr(3)*ur)*ny
    
    
    flux_l_ql(1) = (-ul**2 + pl_ql(1))*nx - (ul*vl)*ny 
    flux_l_ql(2) = (2._dp*ul + pl_ql(2))*nx + vl*ny
    flux_l_ql(3) = pl_ql(3)*nx + ul*ny
    flux_l_ql(4) = pl_ql(4)*nx
    
    flux_r_qr(1) = (-ur**2 + pr_qr(1))*nx - (ur*vr)*ny 
    flux_r_qr(2) = (2._dp*ur + pr_qr(2))*nx + vr*ny 
    flux_r_qr(3) = pr_qr(3)*nx + ur*ny
    flux_r_qr(4) = pr_qr(4)*nx

    flux(2) = half*( (flux_l + flux_r) + alpha*( ql(2) - qr(2) ))
   
    flux_ql(2,1)    = half * (flux_l_ql(1)  + alpha_ql(1)  *( ql(2) - qr(2) ))
    flux_ql(2,2)    = half * (flux_l_ql(2)  + alpha_ql(2)  *( ql(2) - qr(2) ) + alpha)
    flux_ql(2,3:4)  = half * (flux_l_ql(3:4)+ alpha_ql(3:4)*( ql(2) - qr(2) ))
    
    flux_qr(2,1)    = half * (flux_r_qr(1)  + alpha_qr(1)  *( ql(2) - qr(2) ))
    flux_qr(2,2)    = half * (flux_r_qr(2)  + alpha_qr(2)  *( ql(2) - qr(2) ) - alpha)
    flux_qr(2,3:4)  = half * (flux_r_qr(3:4)+ alpha_qr(3:4)*( ql(2) - qr(2) ))
    
    !---> Flux 3 

    flux_l = (ql(2)*vl)*nx + (ql(3)*vl + pl)*ny
    flux_r = (qr(2)*vr)*nx + (qr(3)*vr + pr)*ny
    
    flux_l_ql(1) = -ul*vl*nx + (-vl**2 + pl_ql(1))*ny
    flux_l_ql(2) = vl*nx + pl_ql(2)*ny
    flux_l_ql(3) = ul*nx + (2._dp*vl + pl_ql(3))*ny
    flux_l_ql(4) = pl_ql(4)*ny
    
    flux_r_qr(1) = -ur*vr*nx + (-vr**2 + pr_qr(1))*ny
    flux_r_qr(2) = vr*nx + pr_qr(2)*ny
    flux_r_qr(3) = ur*nx + (2._dp*vr + pr_qr(3))*ny
    flux_r_qr(4) = pr_qr(4)*ny

    flux(3) = half*( (flux_l + flux_r) + alpha*( ql(3) - qr(3) ))
    
    flux_ql(3,1:2)  = half * (flux_l_ql(1:2)+ alpha_ql(1:2)*( ql(3) - qr(3) ))
    flux_ql(3,3)    = half * (flux_l_ql(3)  + alpha_ql(3)  *( ql(3) - qr(3) ) + alpha)
    flux_ql(3,4)    = half * (flux_l_ql(4)  + alpha_ql(4)  *( ql(3) - qr(3) ))
    
    flux_qr(3,1:2)  = half * (flux_r_qr(1:2)+ alpha_qr(1:2)*( ql(3) - qr(3) ))
    flux_qr(3,3)    = half * (flux_r_qr(3)  + alpha_qr(3)  *( ql(3) - qr(3) ) - alpha)
    flux_qr(3,4)    = half * (flux_r_qr(4)  + alpha_qr(4)  *( ql(3) - qr(3) ))

    !---> Flux 4 

    flux_l = (ql(4)*ul + pl*ul)*nx + (ql(4)*vl + pl*vl)*ny
    flux_r = (qr(4)*ur + pr*ur)*nx + (qr(4)*vr + pr*vr)*ny 
    
    flux_l_ql(1) = (-El*ul+pl_ql(1)*ul-pl*ql(2)/(ql(1)**2))*nx+(-El*vl+pl_ql(1)*vl-pl*ql(3)/(ql(1)**2))*ny
    flux_l_ql(2) = (El + pl_ql(2)*ul + pl/ql(1))*nx + pl_ql(2)*vl*ny
    flux_l_ql(3) = pl_ql(3)*ul*nx + (El + pl_ql(3)*vl + pl/ql(1))*ny
    flux_l_ql(4) = (ul + pl_ql(4)*ul)*nx + (vl + pl_ql(4)*vl)*ny
    
    flux_r_qr(1) = (-Er*ur+pr_qr(1)*ur-pr*qr(2)/(qr(1)**2))*nx+(-Er*vr+pr_qr(1)*vr-pr*qr(3)/(qr(1)**2))*ny
    flux_r_qr(2) = (Er + pr_qr(2)*ur + pr/qr(1))*nx + pr_qr(2)*vr*ny
    flux_r_qr(3) = pr_qr(3)*ur*nx + (Er + pr_qr(3)*vr + pr/qr(1))*ny
    flux_r_qr(4) = (ur + pr_qr(4)*ur)*nx + (vr + pr_qr(4)*vr)*ny

    flux(4) = half * ( (flux_l + flux_r) + alpha*( ql(4) - qr(4) ))
    
    flux_ql(4,4)    = half * (flux_l_ql(4)  + alpha_ql(4)  *( ql(4) - qr(4) ) + alpha)
    flux_ql(4,1:3)  = half * (flux_l_ql(1:3)+ alpha_ql(1:3)*( ql(4) - qr(4) ))
    
    flux_qr(4,4)    = half * (flux_r_qr(4)  + alpha_qr(4)  *( ql(4) - qr(4) ) - alpha)
    flux_qr(4,1:3)  = half * (flux_r_qr(1:3)+ alpha_qr(1:3)*( ql(4) - qr(4) ))
    
    flux    = flux*mag
    flux_ql = flux_ql*mag
    flux_qr = flux_qr*mag

  end subroutine lax_friedrichs_q
  


  
  subroutine lax_friedrichs(ql, qr, norm, flux)
    use my_kinddefs

    real(dp),intent(in)  :: ql(:)
    real(dp),intent(in)  :: qr(:)
    real(dp),intent(out) :: flux(:)
    real(dp) :: eigl, eigr
    real(dp) :: nx, ny
    real(dp) :: al, ar
    real(dp) :: ul, vl, pl, El
    real(dp) :: ur, vr, pr, Er
    real(dp) :: flux_l
    real(dp) :: flux_r
    real(dp) :: alpha

    real(dp) :: mag
    real(dp) :: gamma
    real(dp),intent(in)  :: norm(:)


    !---> Local Variables
   
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
  
    eigl = al + my_smooth_abs(ul*nx + vl*ny)
    eigr = ar + my_smooth_abs(ur*nx + vr*ny)

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
    
    flux    = flux * mag

  end subroutine lax_friedrichs
  
  subroutine lax_friedrichsComplex(ql, qr, norm, flux)
    use my_kinddefs
#ifdef complx
    complex(dp),intent(in)  :: ql(:)
    complex(dp),intent(in)  :: qr(:)
    complex(dp),intent(out) :: flux(:)
    complex(dp) :: eigl, eigr
    complex(dp) :: nx, ny
    complex(dp) :: al, ar
    complex(dp) :: ul, vl, pl, El
    complex(dp) :: ur, vr, pr, Er
    complex(dp) :: flux_l
    complex(dp) :: flux_r
    complex(dp) :: alpha
#else
    real(dp),intent(in)  :: ql(:)
    real(dp),intent(in)  :: qr(:)
    real(dp),intent(out) :: flux(:)
    real(dp) :: eigl, eigr
    real(dp) :: nx, ny
    real(dp) :: al, ar
    real(dp) :: ul, vl, pl, El
    real(dp) :: ur, vr, pr, Er
    real(dp) :: flux_l
    real(dp) :: flux_r
    real(dp) :: alpha
#endif
    real(dp) :: mag
    real(dp) :: gamma
    real(dp),intent(in)  :: norm(:)


    !---> Local Variables
   
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
  
    eigl = al + my_smooth_absComplex(ul*nx + vl*ny)
    eigr = ar + my_smooth_absComplex(ur*nx + vr*ny)

    !alpha = max( eigl,eigr )
    if( real(eigl) > real(eigr)) then
       alpha = eigl
    else if( real(eigl) <= real(eigr) ) then
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
    
    flux    = flux * mag

  end subroutine lax_friedrichsComplex
  
 real(dp) function my_smooth_dabs(x)
      real(dp),intent(in)   :: x
      real(dp)              :: my_epsilon
      
      my_epsilon = 1e-9
      
      my_smooth_dabs = x*(abs(x) + 2._dp*my_epsilon)/(abs(x) + my_epsilon)**2
      
  end function
  
  real(dp) function my_smooth_abs(x)
      real(dp),intent(in)   :: x
      real(dp)              :: my_epsilon
      
      my_epsilon = 1e-9
      
      my_smooth_abs = x*x/(sign(1._dp,x)*x + my_epsilon)
      
  end function
  

  complex(dp) function my_smooth_absComplex(x)
#ifdef complx
      complex(dp),intent(in)   :: x
#else
      real(dp),intent(in)   :: x
#endif
      real(dp)              :: my_epsilon
      
      my_epsilon = 1e-9
      
      my_smooth_absComplex = x*x/(sign(1._dp,real(x,dp))*x + my_epsilon)
      
  end function
  

end module flux_module