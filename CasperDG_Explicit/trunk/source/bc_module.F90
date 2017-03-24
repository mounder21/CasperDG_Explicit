module bc_module
  implicit none
contains

  subroutine get_bc(bc_type, norm, q, qb)
    use my_kinddefs
    use inputs_module

    integer(i4),intent(in)            :: bc_type
    real(dp),dimension(:),intent(in)  :: norm
    
    real(dp),dimension(:),intent(in)  :: q
    real(dp),dimension(:),intent(out) :: qb
    real(dp) :: u, v, p, aT, ub, vb, rhob, E, un, Mn, Rp, Rm, ab, s, unb, ut, pb,rho
    real(dp) :: CvT, Eb, term


    !---> Local Variables
    real(dp) mag, nx, ny, gm1
    
    
    gm1 = gamma - 1._dp
    
    !---> Extract the primatives
    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
    E   = q(4)/q(1)
    
    !---> Get normals
    mag = sqrt(norm(1)*norm(1) + norm(2)*norm(2))
    nx = norm(1)/mag
    ny = norm(2)/mag
    
    select case(bc_type)
    case(1)
       !---> Freestream
       qb(1) = 1._dp
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp/(gamma*(gamma - 1._dp)) + half*(fmach**2)
       
    case(2)
       !---> Inviscid Wall   
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       
       !---> Get the boundary velocity
       ub = u - (u*nx + v*ny)*nx
       vb = v - (u*nx + v*ny)*ny

       !---> Now set the conservative vector
       qb(1) = q(1)
       qb(2) = ub*qb(1)
       qb(3) = vb*qb(1)
       qb(4) = q(4)

    case(3)
       !---> Viscous Wall
       
       CvT = E - half*(u**2 + v**2)
       
       qb(1) = q(1)
       qb(2) = 0._dp
       qb(3) = 0._dp
       qb(4) = CvT*q(1)
       
    case(5)
       !---> Normal velocity
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       
       u = q(2)/q(1)
       v = q(3)/q(1)
       E = q(4)/q(1)
      
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       aT = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/aT
       
       !---> Compute the Riemann invariants
       Rp =  un + 2._dp*aT/(gamma - 1._dp)
       
       Rm = fmach*cos(alpha*pi/180._dp)*nx + fmach*sin(alpha*pi/180._dp)*ny &
            - 2._dp/(gamma - 1._dp)

       unb = half*(Rp + Rm)
          
       ab = (gamma - 1._dp)/4._dp*(Rp - Rm)
       
       if( real(unb) <= 0._dp) then
          !--->Inflow  Entropy and tangential Velocity
          s = 1._dp/gamma
        
          ut = -fmach*cos(alpha*pi/180._dp)*ny + fmach*sin(alpha*pi/180._dp)*nx

       else if( real(unb) > 0._dp) then
          !---> Outflow entropy and tangential Velocity
          s = p/(q(1)**gamma)
         
          ut = -u*ny + v*nx
       end if

       !---> Compute the new density and pressure
       rhob = (s*gamma/(ab**2))**(1._dp/(1._dp - gamma))
       pb = s*rhob**gamma

       ub = nx*unb - ny*ut
       vb = nx*ut + ny*unb

       if( real(Mn) < -1._dp) then
           
          !---> Super sonic inflow, everything imposed from outside
          qb(1) = rho_in 
          qb(2) = rho_in*fmach*cos(alpha*pi/180._dp)
          qb(3) = rho_in*fmach*sin(alpha*pi/180._dp)
          qb(4) = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0)

       else if( real(Mn) > 1._dp) then
          !---> Super sonic outflow, everything comes from inside
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)

       else
          !---> Subsonic in-outflow , use the stuff computed above
          qb(1) = rhob
          qb(2) = rhob*ub
          qb(3) = rhob*vb
          qb(4) = pb/(gamma - 1._dp) + half*rhob*(ub**2 + vb**2)
       end if
       
  
    case(7)
       
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       a = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/a
       
       cvT = E - half*(u**2 + v**2)
       
       pb = 1._dp/gamma
       rhob = pb/(cvT*(gamma - 1._dp))
       
       Eb = E !pb/(rhob*(gamma - 1._dp)) + half*(u**2 + v**2)
       if( real(Mn) > 1._dp) then
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)
       else
          qb(1) = rhob
          qb(2) = u*rhob
          qb(3) = v*rhob
          qb(4) = rhob*Eb
       end if
       
       
    case(8)
       !---> Super sonic outflow
       qb(:) = q(:)
       !---> Super sonic outflow
    
    case(9)
       qb(1) = 1._dp
       qb(2) = fmach
       qb(3) = 0._dp
       pb = 1._dp/gamma
       qb(4) = pb/(gamma - 1._dp) + half*(fmach*fmach)
       
    case(10)
       term = ((gamma + 1._dp)*fmach*fmach)/(2._dp + &
            (gamma - 1._dp)*fmach*fmach)
       
       rhob = term
       ub = fmach/term
       pb = 1._dp/gamma*(1._dp + &
            2._dp*gamma/(gamma + 1._dp)*(fmach*fmach - 1._dp))
       qb(1) = rhob
       qb(2) = ub*rhob 
       qb(3) = 0._dp
       qb(4) = pb/(gamma - 1._dp) + rhob*half*ub*ub

    case(11)
       pb  = (gamma - 1._dp)* (E - 0.5_dp*rho_in*(u**2 + v**2))
       qb(1) = q(1)
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp !q(4)

   end select
    
  end subroutine get_bc
  
  subroutine get_bcComplex(bc_type, norm, q, qb)
    use my_kinddefs
    use inputs_module

    integer(i4),intent(in)            :: bc_type
    real(dp),dimension(:),intent(in)  :: norm
    
#ifdef complx
    complex(dp),dimension(:),intent(in)  :: q
    complex(dp),dimension(:),intent(out) :: qb
    complex(dp) :: u, v, p, aT, ub, vb, rhob, E, un, Mn, Rp, Rm, ab, s, unb, ut, pb,rho
    complex(dp) :: CvT, Eb, term
#else
    real(dp),dimension(:),intent(in)  :: q
    real(dp),dimension(:),intent(out) :: qb
    real(dp) :: u, v, p, aT, ub, vb, rhob, E, un, Mn, Rp, Rm, ab, s, unb, ut, pb,rho
    real(dp) :: CvT, Eb, term
#endif

    !---> Local Variables
    real(dp) mag, nx, ny, gm1
    
    
    gm1 = gamma - 1._dp
    
    !---> Extract the primatives
    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
    E   = q(4)/q(1)
    
    !---> Get normals
    mag = sqrt(norm(1)*norm(1) + norm(2)*norm(2))
    nx = norm(1)/mag
    ny = norm(2)/mag
    
    select case(bc_type)
    case(1)
       !---> Freestream
       qb(1) = 1._dp
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp/(gamma*(gamma - 1._dp)) + half*(fmach**2)
       
    case(2)
       !---> Inviscid Wall   
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       
       !---> Get the boundary velocity
       ub = u - (u*nx + v*ny)*nx
       vb = v - (u*nx + v*ny)*ny

       !---> Now set the conservative vector
       qb(1) = q(1)
       qb(2) = ub*qb(1)
       qb(3) = vb*qb(1)
       qb(4) = q(4)

    case(3)
       !---> Viscous Wall
       
       CvT = E - half*(u**2 + v**2)
       
       qb(1) = q(1)
       qb(2) = 0._dp
       qb(3) = 0._dp
       qb(4) = CvT*q(1)
       
    case(5)
       !---> Normal velocity
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       u = q(2)/q(1)
       v = q(3)/q(1)
       E = q(4)/q(1)
      
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       aT = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/aT
       
       !---> Compute the Riemann invariants
       Rp =  un + 2._dp*aT/(gamma - 1._dp)
       
       Rm = fmach*cos(alpha*pi/180._dp)*nx + fmach*sin(alpha*pi/180._dp)*ny &
            - 2._dp/(gamma - 1._dp)

       unb = half*(Rp + Rm)
          
       ab = (gamma - 1._dp)/4._dp*(Rp - Rm)
       
       if( real(un) <= 0._dp) then
          !--->Inflow  Entropy and tangential Velocity
          s = 1._dp/gamma
        
          ut = -fmach*cos(alpha*pi/180._dp)*ny + fmach*sin(alpha*pi/180._dp)*nx

       else if( real(un) > 0._dp) then
          !---> Outflow entropy and tangential Velocity
          s = p/(q(1)**gamma)
         
          ut = -u*ny + v*nx
       end if

       !---> Compute the new density and pressure
       rhob = (s*gamma/(ab**2))**(1._dp/(1._dp - gamma))
       pb = s*rhob**gamma

       ub = nx*unb - ny*ut
       vb = nx*ut + ny*unb

       if( real(Mn) < -1._dp) then
           
          !---> Super sonic inflow, everything imposed from outside
          qb(1) = rho_in 
          qb(2) = rho_in*fmach*cos(alpha*pi/180._dp)
          qb(3) = rho_in*fmach*sin(alpha*pi/180._dp)
          qb(4) = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0)

       else if( real(Mn) > 1._dp) then
          !---> Super sonic outflow, everything comes from inside
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)
       else
          !---> Subsonic in-outflow , use the stuff computed above
          qb(1) = rhob
          qb(2) = rhob*ub
          qb(3) = rhob*vb
          qb(4) = pb/(gamma - 1._dp) + half*rhob*(ub**2 + vb**2)
       end if
       
    case(7)
       
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       a = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/a
       
       cvT = E - half*(u**2 + v**2)
       
       pb = 1._dp/gamma
       rhob = pb/(cvT*(gamma - 1._dp))
       
       Eb = E !pb/(rhob*(gamma - 1._dp)) + half*(u**2 + v**2)
       if( real(Mn) > 1._dp) then
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)
       else
          qb(1) = rhob
          qb(2) = u*rhob
          qb(3) = v*rhob
          qb(4) = rhob*Eb
       end if
       
       
    case(8)
       !---> Super sonic outflow
       qb(:) = q(:)
       !---> Super sonic outflow
    
    case(9)
       qb(1) = 1._dp
       qb(2) = fmach
       qb(3) = 0._dp
       pb = 1._dp/gamma
       qb(4) = pb/(gamma - 1._dp) + half*(fmach*fmach)
       
    case(10)
       term = ((gamma + 1._dp)*fmach*fmach)/(2._dp + &
            (gamma - 1._dp)*fmach*fmach)
       
       rhob = term
       ub = fmach/term
       pb = 1._dp/gamma*(1._dp + &
            2._dp*gamma/(gamma + 1._dp)*(fmach*fmach - 1._dp))
       qb(1) = rhob
       qb(2) = ub*rhob 
       qb(3) = 0._dp
       qb(4) = pb/(gamma - 1._dp) + rhob*half*ub*ub

    case(11)
       pb  = (gamma - 1._dp)* (E - 0.5_dp*rho_in*(u**2 + v**2))
       qb(1) = q(1)
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp !q(4)

   end select
    
  end subroutine get_bcComplex

!*****************************************************************************80
!>\brief get_bc_q: computes the boundary q based on the interior state and the 
!             boundary condition
!> \author Nick Burgess
!!  CFD Lab Dept. of Mechanical Engineering
!>  Univ. of Wyoming
!>  \version 3.2
!> \date 9/3/2010 
!   ARGUMENT LIST :
!> \param[in] bc_type The type of boundary condition it is
!> \param[in] d       Normal distance from wall to cell center for k-omega b.c
!> \param[in] norm    Normal vector components
!> \param[in] q       Solution vector
!> \param[out] qb     Boundary solution vector 
!*****************************************************************************80
  subroutine get_bc_q(bc_type, norm, q, qb, dqbdq)
    use my_kinddefs
    use inputs_module

    integer(i4),intent(in)            :: bc_type
    real(dp),dimension(:),intent(in)  :: norm
    real(dp),dimension(:),intent(in)  :: q
    real(dp),dimension(:),intent(out) :: qb
    real(dp),dimension(:,:),intent(out) :: dqbdq

    !---> Local Variables
    real(dp) mag, nx, ny
    real(dp) u, v, p, aT, ub, vb, rhob, E, un, Mn, Rp, Rm, ab, s, unb, ut, pb,rho
    real(dp) CvT, Eb, term, gm1
    
    gm1 = gamma - 1._dp
    
    !---> Extract the primatives
    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
    E   = q(4)/q(1)
    
    !---> Get normals
    mag = sqrt(norm(1)*norm(1) + norm(2)*norm(2))
    nx = norm(1)/mag
    ny = norm(2)/mag
    
    dqbdq(:,:) = 0._dp
  
    select case(bc_type)
    case(1)
       !---> Freestream
       qb(1) = 1._dp
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp/(gamma*(gamma - 1._dp)) + half*(fmach**2)
       
    case(2)
       !---> Inviscid Wall   
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       
       
       !---> Get the boundary velocity
       ub = u - (u*nx + v*ny)*nx
       vb = v - (u*nx + v*ny)*ny

       !---> Now set the conservative vector
       qb(1) = q(1)
       qb(2) = ub*qb(1)
       qb(3) = vb*qb(1)
       qb(4) = q(4)
       
       dqbdq(1,1) = 1._dp
       dqbdq(2,2) = 1._dp - nx*nx
       dqbdq(2,3) = -nx*ny
       dqbdq(3,2) = -nx*ny
       dqbdq(3,3) = 1._dp - ny*ny
       dqbdq(4,4) = 1._dp 
      
    case(3)
       !---> Viscous Wall
       
       CvT = E - half*(u**2 + v**2)
       
       qb(1) = q(1)
       qb(2) = 0._dp
       qb(3) = 0._dp
       qb(4) = CvT*q(1)
       
       !---> No slip wall
       dqbdq(1,1) = 1._dp
       dqbdq(4,1) = half*(u*u + v*v)
       dqbdq(4,2) = -u
       dqbdq(4,3) = -v
       dqbdq(4,4) = 1._dp
    case(5)
       !---> Normal velocity
       mag = sqrt(norm(1)**2 + norm(2)**2)
       nx = norm(1)/mag
       ny = norm(2)/mag
       u = q(2)/q(1)
       v = q(3)/q(1)
       E = q(4)/q(1)
      
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       aT = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/aT
       
       !---> Compute the Riemann invariants
       Rp =  un + 2._dp*aT/(gamma - 1._dp)
       
       Rm = fmach*cos(alpha*pi/180._dp)*nx + fmach*sin(alpha*pi/180._dp)*ny &
            - 2._dp/(gamma - 1._dp)

       unb = half*(Rp + Rm)
          
       ab = (gamma - 1._dp)/4._dp*(Rp - Rm)
       
       if( un <= 0._dp) then
          !--->Inflow  Entropy and tangential Velocity
          s = 1._dp/gamma
        
          ut = -fmach*cos(alpha*pi/180._dp)*ny + fmach*sin(alpha*pi/180._dp)*nx

       else if( un > 0._dp) then
          !---> Outflow entropy and tangential Velocity
          s = p/(q(1)**gamma)
         
          ut = -u*ny + v*nx
       end if

       !---> Compute the new density and pressure
       rhob = (s*gamma/(ab**2))**(1._dp/(1._dp - gamma))
       pb = s*rhob**gamma

       ub = nx*unb - ny*ut
       vb = nx*ut + ny*unb

       if( Mn < -1._dp) then
           
          !---> Super sonic inflow, everything imposed from outside
          qb(1) = rho_in 
          qb(2) = rho_in*fmach*cos(alpha*pi/180._dp)
          qb(3) = rho_in*fmach*sin(alpha*pi/180._dp)
          qb(4) = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0)

       else if( Mn > 1._dp) then
          !---> Super sonic outflow, everything comes from inside
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)
          
       else
          !---> Subsonic in-outflow , use the stuff computed above
          qb(1) = rhob
          qb(2) = rhob*ub
          qb(3) = rhob*vb
          qb(4) = pb/(gamma - 1._dp) + half*rhob*(ub**2 + vb**2)
       end if
       
       !---> Linearization of Characteristic In/Outflow Boundary
       call lin_charbc(q, norm, dqbdq(1:4,1:4))
  
    case(7)
       
       p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
       a = sqrt(p/q(1)*gamma)
       
       un = u*nx + v*ny
       Mn = un/a
       
       cvT = E - half*(u**2 + v**2)
       
       pb = 1._dp/gamma
       rhob = pb/(cvT*(gamma - 1._dp))
       
       Eb = E !pb/(rhob*(gamma - 1._dp)) + half*(u**2 + v**2)
       if( Mn > 1._dp) then
          qb(1) = q(1)
          qb(2) = q(2)
          qb(3) = q(3)
          qb(4) = q(4)
       else
          qb(1) = rhob
          qb(2) = u*rhob
          qb(3) = v*rhob
          qb(4) = rhob*Eb
       end if
       
       call lin_outflowbc(q, norm, dqbdq(1:4,1:4))
       
    case(8)
       !---> Super sonic outflow
       qb(:) = q(:)
       !---> Super sonic outflow
       dqbdq(1,1) = 1._dp
       dqbdq(2,2) = 1._dp
       dqbdq(3,3) = 1._dp
       dqbdq(4,4) = 1._dp
    
    case(9)
       qb(1) = 1._dp
       qb(2) = fmach
       qb(3) = 0._dp
       pb = 1._dp/gamma
       qb(4) = pb/(gamma - 1._dp) + half*(fmach*fmach)
       
    case(10)
       term = ((gamma + 1._dp)*fmach*fmach)/(2._dp + &
            (gamma - 1._dp)*fmach*fmach)
       
       rhob = term
       ub = fmach/term
       pb = 1._dp/gamma*(1._dp + &
            2._dp*gamma/(gamma + 1._dp)*(fmach*fmach - 1._dp))
       qb(1) = rhob
       qb(2) = ub*rhob 
       qb(3) = 0._dp
       qb(4) = pb/(gamma - 1._dp) + rhob*half*ub*ub

    case(11)
       pb  = (gamma - 1._dp)* (E - 0.5_dp*rho_in*(u**2 + v**2))
       qb(1) = q(1)
       qb(2) = fmach*cos(alpha*pi/180._dp)
       qb(3) = fmach*sin(alpha*pi/180._dp)
       qb(4) = 1._dp !q(4)
       
       dqbdq(1,1) = 1._dp
       !dqbdq(4,4) = 1._dp

   end select
    
  end subroutine get_bc_q
  
  subroutine get_bclin(bc_type, norm, q, dqbdq)
    use my_kinddefs
    use inputs_module, only: fmach,gamma,alpha
    
    integer(i4),intent(in)              :: bc_type
    real(dp),dimension(:),intent(in)    :: norm
    real(dp),dimension(:),intent(in)    :: q
    real(dp),dimension(:,:),intent(out) :: dqbdq

    !---> Local Variables
    real(dp) u, v, rho, mag, nx, ny, E, CvT,gm1
      

    gm1 = gamma - 1._dp
    
    dqbdq(:,:) = 0._dp
    
    !---> Extract some primatives
    mag = sqrt(norm(1)**2 + norm(2)**2)
    nx = norm(1)/mag
    ny = norm(2)/mag
    
    rho = q(1)
    u = q(2)/q(1)
    v = q(3)/q(1)
    E = q(4)/q(1)
    
    select case(bc_type)
    case(1)
       !---> Freestream Do nothing already set = 0
    case(2)
      !---> Inviscid Wall
       dqbdq(1,1) = 1._dp
       dqbdq(2,2) = 1._dp - nx*nx
       dqbdq(2,3) = -nx*ny
       dqbdq(3,2) = -nx*ny
       dqbdq(3,3) = 1._dp - ny*ny
       dqbdq(4,4) = 1._dp  

    case(3)
       !---> No slip wall
       dqbdq(1,1) = 1._dp
       dqbdq(4,1) = half*(u*u + v*v)
       dqbdq(4,2) = -u
       dqbdq(4,3) = -v
       dqbdq(4,4) = 1._dp
             
    case(5)      
       !---> Linearization of Characteristic In/Outflow Boundary
       call lin_charbc(q, norm, dqbdq(1:4,1:4))
    case(7)
       call lin_outflowbc(q, norm, dqbdq(1:4,1:4))
    case(8) !---> Super sonic outflow
       dqbdq(1,1) = 1._dp
       dqbdq(2,2) = 1._dp
       dqbdq(3,3) = 1._dp
       dqbdq(4,4) = 1._dp
    case(9)
       dqbdq(:,:) = 0._dp
    case(11)
       dqbdq(1,1) = 1._dp
       !dqbdq(4,4) = 1._dp
    end select

  end subroutine get_bclin
  
  subroutine lin_charbc(q, norm, dqbdq)
    use my_kinddefs
    use inputs_module, only: fmach,gamma,alpha,rho_in,P_in
    
    real(dp),dimension(:),intent(in)    :: q
    real(dp),dimension(:),intent(in)    :: norm
    real(dp),dimension(:,:),intent(out) :: dqbdq 

    !--->Local Variables
    real(dp) u, v, E, p, c, un, ut, nx, ny, mag, Mn, rho, Rp, Rm,uref
    real(dp) :: rho_0,rhou_0, rhov_0,rhoE_0,u_0,v_0,p_0,c_0
    real(dp) unb, cb, rhob, pb, ub, vb, sb, term,gm1
    real(dp),dimension(4) :: dundq, dpdq, dcdq, drhodq, dRpdq, dutdq,dudq,dvdq
    real(dp),dimension(4) :: dpbdq, dsdq, dcbdq, dunbdq, dtermdq,qb,dubdq,dvbdq


    dqbdq = 0.0_dp
    
    gm1 = gamma - 1._dp
    
    !---> Extract the primatives
    rho = q(1)
    u   = q(2)/q(1)
    v   = q(3)/q(1)
    E   = q(4)/q(1)
    
    rho_0  = rho_in
    rhou_0 = rho_in*fmach*cos(alpha*pi/180._dp)
    rhov_0 = rho_in*fmach*sin(alpha*pi/180._dp)
    rhoE_0 = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0)

    u_0 = rhou_0/rho_0
    v_0 = rhov_0/rho_0
    p_0 = gm1*(rhoE_0 - 0.5_dp * rho_0 * (u_0*u_0+v_0*v_0))
    c_0 = sqrt(gamma*p_0/rho_0)

    mag = sqrt(norm(1)**2 + norm(2)**2)
    nx = norm(1)/mag
    ny = norm(2)/mag

    p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
    c = sqrt(p/rho*gamma)

    un = u*nx + v*ny
    uref = u_0*nx + v_0*ny    

    Mn = un/c

    !---> Compute the Riemann invariants
    Rp =  un + 2._dp*c/(gm1)
    
    Rm = fmach*cos(alpha*pi/180._dp)*nx + fmach*sin(alpha*pi/180._dp)*ny &
         - 2._dp/(gm1)

    
    unb = half*(Rp + Rm)
    
    cb = 0.25_dp*gm1*(Rp - Rm)
    
    !---------------------------- Primative Derivatives ------------------------
    !---> drhodq
    drhodq(:) = 0._dp
    drhodq(1) = 1._dp
    
    dudq(:) = 0.0_dp
    dudq(1) = -u/rho
    dudq(2) = 1.0_dp/rho

    dvdq(:) = 0.0_dp
    dvdq(1) = -v/rho
    dvdq(3) = 1.0_dp/rho

    !---> dundq(:)
    dundq(:) = 0._dp
    dundq(1) = -u/rho*nx - v/rho*ny
    dundq(2) = 1._dp/rho*nx
    dundq(3) = 1._dp/rho*ny
    
    !---> dpdq(:)
    dpdq(1) = half*gm1*(u**2 + v**2) 
    dpdq(2) = -gm1*u
    dpdq(3) = -gm1*v
    dpdq(4) = gm1
    

    !---> dqdq(:) 
    dcdq(:) = half*gamma/(c)*(dpdq(:)/rho - drhodq(:)*p/rho**2)

    !---> dRpdq(:)
    dRpdq(:) = dundq(:) + 2._dp/(gm1)*dcdq(:)
    
    !---> d(unb)/dq
    dunbdq(:) = half*(dRpdq(:))
                
    !---> d(ab)/dq
    dcbdq(:) = (gm1)/4._dp*(dRpdq(:))
    
    if( abs(unb) >= abs(cb)) then
        ! supersonic

        if(unb < 0.0_dp) then
            !  Super sonic inflow, everything imposed from outside
            qb(1) = rho_0
            qb(2) = rhou_0
            qb(3) = rhov_0
            qb(4) = rhoE_0

            dqbdq(1:4,1:4) = 0.0_dp
        else
            !  Super sonic outflow, everything comes from inside
            qb(1:4) = q(1:4)

            dqbdq(1,1) = 1.0_dp
            dqbdq(2,2) = 1.0_dp
            dqbdq(3,3) = 1.0_dp
            dqbdq(4,4) = 1.0_dp

        end if

    else
        ! subsonic

        if(unb < 0.0_dp) then
            ! Subsonic inflow
            sb = p_0/rho_0**gamma
            dsdq(:) = 0.0_dp

            ub = u_0 + (unb-uref)*nx
            vb = v_0 + (unb-uref)*ny

            dubdq(:) = dunbdq(:)*nx
            dvbdq(:) = dunbdq(:)*ny

        else
            ! Subsonic outflow
            sb = p/rho**gamma
            dsdq(:) = dpdq(:)/rho**gamma - gamma*p*drhodq(:)/(rho*rho**gamma)

            ub = u + (unb-un)*nx
            vb = v + (unb-un)*ny

            dubdq(:) = dudq(:) +(dunbdq(:)-dundq(:))*nx
            dvbdq(:) = dvdq(:) +(dunbdq(:)-dundq(:))*ny

        end if

        term = sb*gamma/cb**2
        rhob = term**(-1.0_dp/gm1)
        pb = sb*rhob**gamma
        dtermdq(:) = gamma*dsdq(:)/cb**2 - 2._dp*sb*gamma*dcbdq(:)/cb**3


        qb(1) = rhob
        qb(2) = rhob*ub
        qb(3) = rhob*vb
        qb(4) = pb/gm1 + 0.5_dp*rhob*(ub*ub+vb*vb)


        dqbdq(1,1:4) = 1._dp/(1._dp - gamma)*term**(gamma/(1._dp - gamma))*dtermdq(:)

        !---> We need to do the pressure here
        dpbdq(:) = gamma*rhob**(gm1)*sb*dqbdq(1,1:4) + (rhob**gamma)*dsdq(:)

        !---> d(rhob*ub)/dq
        dqbdq(2,1:4) = dqbdq(1,1:4)*ub + dubdq(:)*rhob

        !---> d(rhob*vb)/dq
        dqbdq(3,1:4) = dqbdq(1,1:4)*vb + dvbdq(:)*rhob

        !---> d(rhob*Eb)/dq
        dqbdq(4,1:4) = dpbdq(:)/gm1 + 0.5_dp*dqbdq(1,1:4)*(ub**2 + &
                        vb**2) + rhob*( dubdq(:)*ub + dvbdq(:)*vb)

    end if
  end subroutine lin_charbc

!*****************************************************************************80
!    lin_outflowbc: linearizes the outflow boundary condition relative to left
!                   state
!*****************************************************************************80
  subroutine lin_outflowbc(q, norm, dqbdqe)
    use my_kinddefs
    use inputs_module, only: fmach,gamma
    
    real(dp),dimension(:),intent(in)    :: q
    real(dp),dimension(:),intent(in)    :: norm
    real(dp),dimension(:,:),intent(out) :: dqbdqe

    !---> Local Variables
    real(dp) u, v, E, rho, rhob, pb, denom, cvT, a, nx, ny, mag, Mn, p,gm1
    real(dp),dimension(4) :: drhobdq, dudq, dvdq, dEdq
   
    gm1 = gamma - 1._dp
    
    !---> Extract the primatives
    rho = q(1)
    u = q(2)/q(1)
    v = q(3)/q(1)
    E = q(4)/q(1)
    
    mag = sqrt(norm(1)**2 + norm(2)**2)
    nx = norm(1)/mag
    ny = norm(2)/mag

    p = q(1)*(gamma - 1._dp)*(E - half*(u**2 + v**2))
    a = sqrt(p/rho*gamma)

    Mn = (u*nx + v*ny)/a
    pb = 1._dp/gamma
    cvT = E - half*(u**2 + v**2)
    rhob = pb/(cvT*(gamma - 1._dp))
    !---> First compute important denominator term
    denom = (gamma - 1._dp)*(cvT)**2
    
    drhobdq(1) = -pb*(-E/rho + u**2/rho  + v**2/rho)/denom
    drhobdq(2) = pb*u/(denom*rho)
    drhobdq(3) = pb*v/(denom*rho)
    drhobdq(4) = -pb/(denom*rho)
    
    dudq(:) = 0._dp
    dudq(1) = -u/rho
    dudq(2) = 1._dp/rho
    
    dvdq(:) = 0._dp
    dvdq(1) = -v/rho
    dvdq(3) = 1._dp/rho
    
    dEdq(:) = 0._dp
    dEdq(1) = -E/rho
    dEdq(4) = 1._dp/rho
 
    if( Mn > 1._dp) then
       dqbdqe(1,1) = 1._dp
       dqbdqe(2,2) = 1._dp
       dqbdqe(3,3) = 1._dp
       dqbdqe(4,4) = 1._dp
    else
       dqbdqe(1,:) = drhobdq(:)
       dqbdqe(2,:) = drhobdq(:)*u + rhob*dudq(:)
       dqbdqe(3,:) = drhobdq(:)*v + rhob*dvdq(:)
       dqbdqe(4,:) = drhobdq(:)*E + rhob*dEdq(:)
    end if
    
  end subroutine lin_outflowbc

  
  Subroutine getFixedBC(tri,edgeNum,loc1,numEdgeGaussPts,qR)
    use my_kinddefs
    use Globals_module, only: triList,edgeList,nodeList,numFields
    use getIC_module
    use projection_module
    
    integer(i4),intent(in) :: edgeNum,numEdgeGaussPts,tri,loc1
    real(dp),intent(out)  :: qR(:,:)
    integer(i4) :: n1,n2,n3,j,k
    

    real(dp)    ::ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts)
    
    n1 = triList(tri)%vtxList(1)
    n2 = triList(tri)%vtxList(2)
    n3 = triList(tri)%vtxList(3)
    ax(1) = nodeList(n1)%x
    ay(1) = nodeList(n1)%y
    ax(2) = nodeList(n2)%x
    ay(2) = nodeList(n2)%y
    ax(3) = nodeList(n3)%x
    ay(3) = nodeList(n3)%y

    Call projectEdgeMapStr8(ax,loc1,x) 
    Call projectEdgeMapStr8(ay,loc1,y) 
    
    do j = 1, numFields 
        do k = 1, numEdgeGaussPts
            qR(j,k) = getIC_MMS(j,x(k),y(k))
        end do
    end do
  end subroutine getFixedBC
  
  subroutine specializedBC_q(bc_type,edgeNum,leftTri,loc1,qL,qR,dqbdq)
    use my_kinddefs
    use inputs_module, only: gamma
    use Globals_module, only: triList,nodeList
    use projection_module
    use flux_module
    use initializeSolution_module, only: solCoeffs
    use getIC_module
    use initializeBasis_module, only: numEdgeGaussPts
    
    integer(i4),intent(in)  :: bc_type,edgeNum,leftTri,loc1
    real(dp),intent(in)     :: ql(:,:)
    real(dp),intent(out)    :: qR(:,:), dqbdq(:,:,:)
    
    integer(i4) :: face2,loc2,n1,n2,n3,GP,rightTri,n
    real(dp)    :: ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts),p,p_0,u_0,v_0,u,v
    real(dp)    :: u_ql(4),v_ql(4),p_ql(4),M_ql(4),E_ql(4),aT_ql(4),aT
    real(dp)    :: E,rhoE_0,M,M_0,aT_0,rho_0,RT 

    select case(bc_type)
        case (12)   !fixed bc--- never used
            !Call getFixedBC(leftTri,edgeNum,loc1,numEdgeGaussPts,qR)
            !Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_q)  !flux(4,#EdgeGPs)
            !do n =1,numFields     
            !    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
            !    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
            !end do
            print*,'This bc type (12) is not in use '
            stop
         case (13)   !Periodic BCs
            rightTri = triList(leftTri)%periodic_Tri_Twin
            triList(leftTri)%used  = .true.
            triList(rightTri)%used = .true.

            !get the twin face of the right triangle
            face2 = triList(rightTri)%twin_Face
            loc2 = 2*face2
            call projectEdge(solCoeffs(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)
            dqbdq(:,:,:) = 0._dp    !qR does not depend on the left state
        case(14)    !inflow
            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                u_0 = getIC_MMS_sin(2,x(GP),y(GP))
                v_0 = getIC_MMS_sin(3,x(GP),y(GP))

                p   = (gamma - 1._dp)* (qL(4,GP) - 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2))
                qR(1,GP) = qL(1,GP)
                qR(2,GP) = qL(1,GP)*u_0
                qR(3,GP) = qL(1,GP)*v_0
                qR(4,GP) = p/(gamma - 1._dp) + 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2)
                qR(1,GP) = getIC_MMS_sin(1,x(GP),y(GP))
                qR(2,GP) = getIC_MMS_sin(2,x(GP),y(GP))
                qR(3,GP) = getIC_MMS_sin(3,x(GP),y(GP))
                qR(4,GP) = getIC_MMS_sin(4,x(GP),y(GP))
            end do
            dqbdq(:,:,:) = 0._dp
            
        case(15)    !outflow

            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                rho_0   = getIC_MMS_sin(1,x(GP),y(GP))
                u_0     = getIC_MMS_sin(2,x(GP),y(GP))
                v_0     = getIC_MMS_sin(3,x(GP),y(GP))
                rhoE_0  = getIC_MMS_sin(4,x(GP),y(GP))
                p_0     = rho_0*(gamma - 1._dp)*(rhoE_0/rho_0 - 0.5_dp*(u_0**2 + v_0**2))
                aT_0    = sqrt(p_0/rho_0*gamma)

                u           = qL(2,GP)/qL(1,GP)
                u_ql(1)     = -qL(2,GP)/qL(1,GP)**2
                u_ql(2)     = 1._dp/qL(1,GP)
                u_ql(3:4)   = 0._dp

                v       = qL(3,GP)/qL(1,GP)
                v_ql(1) = -qL(3,GP)/qL(1,GP)**2
                v_ql(2) = 0._dp
                v_ql(3) = 1._dp/qL(1,GP)
                v_ql(4) = 0._dp

                E         = qL(4,GP)/qL(1,GP)
                E_ql(1)   = -qL(4,GP)/qL(1,GP)**2
                E_ql(2:3) = 0._dp
                E_ql(4)   = 1._dp/qL(1,GP)

                p          = qL(1,GP)*(gamma - 1._dp)*(E - 0.5_dp*(u**2 + v**2))
                p_ql(1)    = half*(gamma - 1._dp)*((ql(2,GP)/ql(1,GP))**2 + (ql(3,GP)/ql(1,GP))**2)        
                p_ql(2:3)  = -(gamma - 1._dp)*(ql(2:3,GP)/ql(1,GP))                                 
                p_ql(4)    =  (gamma - 1._dp)                                                

                aT      = sqrt(p/qL(1,GP)*gamma)
                aT_ql(1)  = 0.5_dp*(gamma*p/ql(1,GP))**(-0.5_dp)*gamma*(ql(1,GP)*p_ql(1)-p)/(ql(1,GP))**2  
                aT_ql(2:4)= 0.5_dp*(gamma*p/ql(1,GP))**(-0.5_dp)*gamma*(ql(1,GP)*p_ql(2:4))/(ql(1,GP))**2

                M   = sqrt(u**2 + v**2)/aT
                M_ql= (aT*(0.5_dp*(u**2 + v**2)**(-0.5_dp)*(2._dp*u*u_ql + 2._dp*v*v_ql)) - &
                       at_ql*(sqrt(u**2 + v**2)))/aT**2

                M_0 = sqrt(u_0**2 + v_0**2)/aT_0

                qR(:,GP)  = qL(:,GP)*(M/M_0) 
                dqbdq(GP,1,1) = (M/M_0) 
                dqbdq(GP,2,2) = (M/M_0)
                dqbdq(GP,3,3) = (M/M_0)
                dqbdq(GP,4,4) = (M/M_0)
                do n = 1,4
                    dqbdq(GP,n,:) =   dqbdq(GP,n,:) + qL(:,GP)*(M_ql(n)/M_0)
                end do 
            end do
            
    end select
  end subroutine specializedBC_q
  
  
 subroutine specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
    use my_kinddefs
    use inputs_module, only: gamma
    use Globals_module, only: triList,nodeList
    use projection_module
    use projection_module
    use initializeSolution_module, only: solCoeffs
    use flux_module
    use getIC_module
    use initializeBasis_module, only: numEdgeGaussPts
    
    integer(i4),intent(in)  :: bc_type,edgeNum,leftTri,loc1
    real(dp),intent(in)     :: ql(:,:)
    real(dp),intent(out)    :: qR(:,:)


    real(dp)    :: ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts)
    integer(i4) :: face2,loc2,n1,n2,n3,GP,rightTri,n
    
    real(dp)    :: p,p_0,u_0,v_0,u,v
    real(dp)    :: aT
    real(dp)    :: E,rhoE_0,M,M_0,aT_0,rho_0,RT 

    select case(bc_type)
        case (12)   !fixed bc--- never used
            !Call getFixedBC(leftTri,edgeNum,loc1,numEdgeGaussPts,qR)
            !Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)  !flux(4,#EdgeGPs)
            !do n =1,numFields     
            !    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
            !    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
            !end do
            print*,'This bc type (12) is not in use '
            stop
         case (13)   !Periodic BCs
            rightTri = triList(leftTri)%periodic_Tri_Twin
            triList(leftTri)%used  = .true.
            triList(rightTri)%used = .true.

            !get the twin face of the right triangle
            face2 = triList(rightTri)%twin_Face
            loc2 = 2*face2
            call projectEdge(solCoeffs(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        case(14)    !inflow
            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                u_0 = getIC_MMS_sin(2,x(GP),y(GP))
                v_0 = getIC_MMS_sin(3,x(GP),y(GP))

                p   = (gamma - 1._dp)* (qL(4,GP) - 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2))
                qR(1,GP) = qL(1,GP)
                qR(2,GP) = qL(1,GP)*u_0
                qR(3,GP) = qL(1,GP)*v_0
                qR(4,GP) = p/(gamma - 1._dp) + 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2)
                qR(1,GP) = getIC_MMS_sin(1,x(GP),y(GP))
                qR(2,GP) = getIC_MMS_sin(2,x(GP),y(GP))
                qR(3,GP) = getIC_MMS_sin(3,x(GP),y(GP))
                qR(4,GP) = getIC_MMS_sin(4,x(GP),y(GP))
            end do

            
        case(15)    !outflow

            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                rho_0   = getIC_MMS_sin(1,x(GP),y(GP))
                u_0     = getIC_MMS_sin(2,x(GP),y(GP))
                v_0     = getIC_MMS_sin(3,x(GP),y(GP))
                rhoE_0  = getIC_MMS_sin(4,x(GP),y(GP))
                p_0     = rho_0*(gamma - 1._dp)*(rhoE_0/rho_0 - 0.5_dp*(u_0**2 + v_0**2))
                aT_0    = sqrt(p_0/rho_0*gamma)

                u = qL(2,GP)/qL(1,GP)
                v = qL(3,GP)/qL(1,GP)
                E = qL(4,GP)/qL(1,GP)
                p = qL(1,GP)*(gamma - 1._dp)*(E - 0.5_dp*(u**2 + v**2))                                            
                aT = sqrt(p/qL(1,GP)*gamma)
                M = sqrt(u**2 + v**2)/aT
                M_0 = sqrt(u_0**2 + v_0**2)/aT_0
                qR(:,GP)  = qL(:,GP)*(M/M_0) 
   
            end do
            
    end select
  end subroutine specializedBC
  
  subroutine specializedBC_Complex(bc_type,edgeNum,leftTri,loc1,qL,qR,solCoeffsRC)
    use my_kinddefs
    use inputs_module, only: gamma
    use Globals_module, only: triList,nodeList
    use projection_module
    use projectionComplex_module
    use flux_module
    use getIC_module
    use initializeBasis_module, only: numEdgeGaussPts
    
    integer(i4),intent(in)  :: bc_type,edgeNum,leftTri,loc1
#ifdef complx
    complex(dp),intent(in)     :: ql(:,:),solCoeffsRC(:,:,:)
    complex(dp),intent(out)    :: qR(:,:)

#else
    real(dp),intent(in)     :: ql(:,:),solCoeffsRC(:,:,:)
    real(dp),intent(out)    :: qR(:,:)

#endif
    real(dp)                :: ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts)
    integer(i4) :: face2,loc2,n1,n2,n3,GP,rightTri,n
    
    real(dp)    :: p,p_0,u_0,v_0,u,v
    real(dp)    :: aT
    real(dp)    :: E,rhoE_0,M,M_0,aT_0,rho_0,RT 

    select case(bc_type)
        case (12)   !fixed bc--- never used
            !Call getFixedBC(leftTri,edgeNum,loc1,numEdgeGaussPts,qR)
            !Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux)  !flux(4,#EdgeGPs)
            !do n =1,numFields     
            !    call integrateEdge(edgePhiW(:,:,loc1),flux(n,:),bL)          !resL , bL(totalModes)
            !    spRes(n,:,leftTri) = spRes(n,:,leftTri)  - bL(:)             !(4,totalModes,numTri)
            !end do
            print*,'This bc type (12) is not in use '
            stop
         case (13)   !Periodic BCs
            rightTri = triList(leftTri)%periodic_Tri_Twin
            triList(leftTri)%used  = .true.
            triList(rightTri)%used = .true.

            !get the twin face of the right triangle
            face2 = triList(rightTri)%twin_Face
            loc2 = 2*face2
            call projectEdgeComplex(solCoeffsRC(:,:,rightTri),loc2,qR)       !sol(numEulerV,tm,numTri)

        case(14)    !inflow
            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                u_0 = getIC_MMS_sin(2,x(GP),y(GP))
                v_0 = getIC_MMS_sin(3,x(GP),y(GP))

                p   = (gamma - 1._dp)* (qL(4,GP) - 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2))
                qR(1,GP) = qL(1,GP)
                qR(2,GP) = qL(1,GP)*u_0
                qR(3,GP) = qL(1,GP)*v_0
                qR(4,GP) = p/(gamma - 1._dp) + 0.5_dp*qL(1,GP)*(u_0**2 + v_0**2)
                qR(1,GP) = getIC_MMS_sin(1,x(GP),y(GP))
                qR(2,GP) = getIC_MMS_sin(2,x(GP),y(GP))
                qR(3,GP) = getIC_MMS_sin(3,x(GP),y(GP))
                qR(4,GP) = getIC_MMS_sin(4,x(GP),y(GP))
            end do

            
        case(15)    !outflow

            n1 = triList(leftTri)%vtxList(1)
            n2 = triList(leftTri)%vtxList(2)
            n3 = triList(leftTri)%vtxList(3)
            ax(1) = nodeList(n1)%x
            ay(1) = nodeList(n1)%y
            ax(2) = nodeList(n2)%x
            ay(2) = nodeList(n2)%y
            ax(3) = nodeList(n3)%x
            ay(3) = nodeList(n3)%y

            Call projectEdgeMapStr8(ax,loc1,x) 
            Call projectEdgeMapStr8(ay,loc1,y) 

            do GP = 1, numEdgeGaussPts
                rho_0   = getIC_MMS_sin(1,x(GP),y(GP))
                u_0     = getIC_MMS_sin(2,x(GP),y(GP))
                v_0     = getIC_MMS_sin(3,x(GP),y(GP))
                rhoE_0  = getIC_MMS_sin(4,x(GP),y(GP))
                p_0     = rho_0*(gamma - 1._dp)*(rhoE_0/rho_0 - 0.5_dp*(u_0**2 + v_0**2))
                aT_0    = sqrt(p_0/rho_0*gamma)

                u = qL(2,GP)/qL(1,GP)
                v = qL(3,GP)/qL(1,GP)
                E = qL(4,GP)/qL(1,GP)
                p = qL(1,GP)*(gamma - 1._dp)*(E - 0.5_dp*(u**2 + v**2))                                            
                aT = sqrt(p/qL(1,GP)*gamma)
                M = sqrt(u**2 + v**2)/aT
                M_0 = sqrt(u_0**2 + v_0**2)/aT_0
                qR(:,GP)  = qL(:,GP)*(M/M_0) 
   
            end do
            

    end select
  end subroutine specializedBC_Complex
  

end module bc_module