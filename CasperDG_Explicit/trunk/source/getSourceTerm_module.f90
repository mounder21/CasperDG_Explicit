Module getSourceTerm_module
    use my_kinddefs
    implicit none
    
    contains
Subroutine getSourceTerm(elem_id,numGaussPts,source)
    use inputs_module, only: MMS_sin
    use Globals_module, only: numEulerVars,triList,nodeList 
    use projection_module
    
    integer(i4),intent(in) :: elem_id
    integer(i4),intent(in) :: numGaussPts
    real(dp),intent(out) :: source(numEulerVars,numGaussPts)
    
    integer(i4) :: n1,n2,n3,GP,k
    real(dp)    :: ax(3),ay(3)
    real(dp)    :: x(numGaussPts),y(numGaussPts)
    
    n1 = triList(elem_id)%vtxList(1)
    n2 = triList(elem_id)%vtxList(2)
    n3 = triList(elem_id)%vtxList(3)
    ax(1) = nodeList(n1)%x
    ay(1) = nodeList(n1)%y
    ax(2) = nodeList(n2)%x
    ay(2) = nodeList(n2)%y
    ax(3) = nodeList(n3)%x
    ay(3) = nodeList(n3)%y

    Call projectMapCell(ax,ay,x,y)
    if(MMS_sin)then
        do k = 1,numEulerVars
            do GP = 1,numGaussPts
                source(k,GP) = getSourceSin(k,x(GP),y(GP))
            end do
        end do
    else 
        do k = 1,numEulerVars
            do GP = 1,numGaussPts
                source(k,GP) = getSource(k,x(GP),y(GP))
            end do
        end do
    end if
    
    
end subroutine

real(dp) function getSource(numEulerVar,x,y)
    use inputs_module, only: a,b,c,d,s0,t0
    use inputs_module, only: rho_in,p_in,gamma
    
    integer(i4),intent(in) :: numEulerVar
    real(dp),intent(in) :: x,y
    real(dp) :: u,v,ux,uy,vx,vy,E,Ex,Ey,p,px,py
    ! rho = rho_in
    
    u  = a*x + b*y + s0
    v  = c*x + d*y + t0
    E  = (p_in / (gamma - 1._dp)) 
    p  = (gamma - 1._dp)* (E - 0.5_dp*rho_in*(u**2 + v**2))
    ux = a
    uy = b
    vx = c
    vy = d
    Ex = 0._dp
    Ey = 0._dp
    px = (gamma - 1._dp)*(Ex - 0.5_dp*rho_in*(2._dp*u*ux + 2._dp*v*vx))
    py = (gamma - 1._dp)*(Ey - 0.5_dp*rho_in*(2._dp*u*uy + 2._dp*v*vy))
    select case(numEulerVar)
        case(1)
            getSource = rho_in * (ux + vy)
        case(2)
            getSource = rho_in * (2._dp*u*ux + uy*v + u*vy) + px
        case(3)
            getSource = rho_in * (ux*v + u*vx + 2._dp*v*vy) + py
        case(4)
            getSource = ux*E + u*Ex + ux*p +u*px + vy*E + v*Ey + vy*p + v*py
    end select

    return
end function  getSource

real(dp) function getSourceSin(numEulerVar,x,y)
    use inputs_module, only: a,b,c,d,s0,t0
    use inputs_module, only: rho_in,p_in,gamma
    
    integer(i4),intent(in) :: numEulerVar
    real(dp),intent(in) :: x,y
    real(dp) :: u,v,ux,uy,vx,vy,rhoE,rhoEx,rhoEy,p,px,py
    ! rho = rho_in

    u  = a*sin(0.1_dp*pi*x) + b*sin(0.1_dp*pi*y) + s0
    v  = c*sin(0.1_dp*pi*x) + d*sin(0.1_dp*pi*y) + t0
    rhoE  = (p_in / (gamma - 1._dp)) + 0.5_dp*rho_in*(u**2 + v**2) 
    p  = (gamma - 1._dp)*(rhoE - 0.5_dp*rho_in*(u**2 + v**2))
    
    ux = 0.1_dp*pi*a*cos(0.1_dp*pi*x)
    uy = 0.1_dp*pi*b*cos(0.1_dp*pi*y)
    vx = 0.1_dp*pi*c*cos(0.1_dp*pi*x)
    vy = 0.1_dp*pi*d*cos(0.1_dp*pi*y)
    rhoEx = 0.5_dp*rho_in*(2._dp*u*ux + 2._dp*v*vx)
    rhoEy = 0.5_dp*rho_in*(2._dp*u*uy + 2._dp*v*vy)
    px = (gamma - 1._dp)*(rhoEx - 0.5_dp*rho_in*(2._dp*u*ux + 2._dp*v*vx))
    py = (gamma - 1._dp)*(rhoEy - 0.5_dp*rho_in*(2._dp*u*uy + 2._dp*v*vy))
    
    select case(numEulerVar)
        case(1)
            getSourceSin = rho_in * (ux + vy)
        case(2)
            getSourceSin = rho_in * (2._dp*u*ux + uy*v + u*vy) + px
        case(3)
            getSourceSin = rho_in * (ux*v + u*vx + 2._dp*v*vy) + py
        case(4)
            getSourceSin = ux*rhoE + u*rhoEx + ux*p +u*px + vy*rhoE + v*rhoEy + vy*p + v*py
    end select

    return
end function  getSourceSin

end module getSourceTerm_module
    