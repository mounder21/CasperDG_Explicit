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
       
    case(8)
       !---> Super sonic outflow
       qb(:) = q(:)
    
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
  
  Subroutine getFixedBC(tri,edgeNum,loc1,numEdgeGaussPts,qR)
    use my_kinddefs
    use Globals_module, only: triList,edgeList,nodeList,numEulerVars
    use getIC_module
    use projection_module
    
    integer(i4),intent(in) :: edgeNum,numEdgeGaussPts,tri,loc1
    real(dp),intent(out)  :: qR(:,:)
    integer(i4) :: n1,n2,n3,j,k
    real(dp) ::ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts)
    
    n1 = triList(tri)%vtxList(1)
    n2 = triList(tri)%vtxList(2)
    n3 = triList(tri)%vtxList(3)
    ax(1) = nodeList(n1)%x
    ay(1) = nodeList(n1)%y
    ax(2) = nodeList(n2)%x
    ay(2) = nodeList(n2)%y
    ax(3) = nodeList(n3)%x
    ay(3) = nodeList(n3)%y
    !print*,edgeList(edgeNum)%n1,edgeList(edgeNum)%n2
    !print*,n1,n2,n3,loc1

    Call projectEdgeMapStr8(ax,loc1,x) 
    Call projectEdgeMapStr8(ay,loc1,y) 
    
    do j = 1, numEulerVars 
        do k = 1, numEdgeGaussPts
            qR(j,k) = getIC_MMS(j,x(k),y(k))
        end do
    end do
  end subroutine getFixedBC
  
  subroutine specializedBC(bc_type,edgeNum,leftTri,loc1,qL,qR)
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
    real(dp),intent(out)    :: qR(:,:)
    
    integer(i4) :: face2,loc2,n1,n2,n3,GP,rightTri,n
    real(dp)    :: ax(3),ay(3),x(numEdgeGaussPts),y(numEdgeGaussPts),p,p_0,u_0,v_0,u,v
    real(dp)    :: u_ql(4),v_ql(4),p_ql(4),M_ql(4),E_ql(4),aT_ql(4),aT
    real(dp)    :: E,rhoE_0,M,M_0,aT_0,rho_0,RT 

    select case(bc_type)
        case (12)   !fixed bc--- never used
            !Call getFixedBC(leftTri,edgeNum,loc1,numEdgeGaussPts,qR)
            !Call getBCFlux(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_q)  !flux(4,#EdgeGPs)
            !do n =1,numEulerVars     
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

                u           = qL(2,GP)/qL(1,GP)

                v       = qL(3,GP)/qL(1,GP)

                E         = qL(4,GP)/qL(1,GP)

                p          = qL(1,GP)*(gamma - 1._dp)*(E - 0.5_dp*(u**2 + v**2))
                aT      = sqrt(p/qL(1,GP)*gamma)
                M   = sqrt(u**2 + v**2)/aT

                M_0 = sqrt(u_0**2 + v_0**2)/aT_0

                qR(:,GP)  = qL(:,GP)*(M/M_0) 

            end do
            
!            case(15)    !outflow
!                n1 = triList(leftTri)%vtxList(1)
!                n2 = triList(leftTri)%vtxList(2)
!                n3 = triList(leftTri)%vtxList(3)
!                ax(1) = nodeList(n1)%x
!                ay(1) = nodeList(n1)%y
!                ax(2) = nodeList(n2)%x
!                ay(2) = nodeList(n2)%y
!                ax(3) = nodeList(n3)%x
!                ay(3) = nodeList(n3)%y
!
!                Call projectEdgeMap(ax,loc1,x) 
!                Call projectEdgeMap(ay,loc1,y) 
!                
!                do GP = 1, numEdgeGaussPts
!                    rho_0 = getIC_MMS_sin(1,x(GP),y(GP))
!                    u_0 = getIC_MMS_sin(2,x(GP),y(GP))
!                    v_0 = getIC_MMS_sin(3,x(GP),y(GP))
!                    E_0 = getIC_MMS_sin(4,x(GP),y(GP))
!                    p_0 = (gamma - 1._dp)*(E_0 - 0.5_dp*rho_0*(u_0**2 + v_0**2))
!                    aT_0 = sqrt(p_0/rho_0*gamma)
!
!                    u = (qL(2,GP)/qL(1,GP))
!                    v = (qL(3,GP)/qL(1,GP))
!                    E = qL(4,GP)
!                    p = (gamma - 1._dp)*(E - 0.5_dp*qL(1,GP)*(u**2 + v**2))
!                    aT = sqrt(p/qL(1,GP)*gamma)
!                    RT = p/qL(1,GP)
!                    
!                    
!                    M = sqrt(u**2 + v**2)/aT
!                    M_0 = sqrt(u_0**2 + v_0**2)/aT_0
!
!                    qR(1,GP) = p_0/RT
!                    qR(2,GP) = p_0/RT * u
!                    qR(3,GP) = p_0/RT * v
!                    qR(4,GP) = p_0/(gamma-1._dp) + 0.5_dp*p_0/RT*(u**2+v**2)
!                    !qR(:,GP) = qL(:,GP)
!                    
!                end do  
    end select
  end subroutine
  

end module bc_module