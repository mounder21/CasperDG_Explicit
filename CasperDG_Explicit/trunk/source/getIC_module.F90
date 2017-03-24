Module getIC_module
    use my_kinddefs
    use inputs_module
    implicit none
    
contains

real(dp) function getIC(x,y)
    real(dp),intent(in) :: x,y
    
    getIC = 5._dp*x + 3._dp*y

    return
end function getIC

real(dp) function getIC_MMS(eulerVar,x,y)
    use inputs_module, only: a,b,c,d,s0,t0
    
    integer(i4), intent(in) :: eulerVar
    real(dp),intent(in) :: x,y
    real(dp) :: u,v
    u = a*x + b*y + s0
    v = c*x + d*y + t0

    select case(eulerVar)
        case(1)
            getIC_MMS = rho_in 
        case(2)
            getIC_MMS = rho_in*u
        case(3)
            getIC_MMS = rho_in*v
        case(4)
            getIC_MMS = (p_in / (gamma - 1.0_dp))  
    end select
    return
end function getIC_MMS

real(dp) function getIC_MMS_sin(eulerVar,x,y)
    use inputs_module, only: a,b,c,d,s0,t0
    
    integer(i4), intent(in) :: eulerVar
    real(dp),intent(in) :: x,y
    real(dp) :: u,v
    u  = a*sin(0.1_dp*pi*x) + b*sin(0.1_dp*pi*y) + s0
    v  = c*sin(0.1_dp*pi*x) + d*sin(0.1_dp*pi*y) + t0

    select case(eulerVar)
        case(1)
            getIC_MMS_sin = rho_in 
        case(2)
            getIC_MMS_sin = rho_in*u
        case(3)
            getIC_MMS_sin = rho_in*v
        case(4)
            getIC_MMS_sin = (p_in / (gamma - 1.0_dp)) + 0.5_dp*rho_in*(u**2 + v**2)
    end select
    return
end function getIC_MMS_sin

real(dp) function getIConstant(eulerVar)
    integer(i4),intent(in) :: eulerVar

    select case (eulerVar)
            case (1)
                getIConstant = rho_in 
            case (2)
                getIConstant = rho_in*fmach*cos(alpha*pi/180._dp)
            case (3)
                getIConstant = rho_in*fmach*sin(alpha*pi/180._dp)
            case (4)
                getIConstant = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0)
    end select
    return
end function getIConstant

real(dp) function smooth_IC_Airfoil(eulerVar,x,y)
    integer(i4),intent(in) :: eulerVar
    real(dp),intent(in)    :: x,y
    real(dp) :: Rad,fs
    
    Rad  = sqrt((x-5)**2 + (y-5)**2) 
    fs = 0.0
    select case (eulerVar)
        case (1)
            smooth_IC_Airfoil = rho_in
        case (2)
            fs = rho_in*fmach*cos(alpha*pi/180._dp)
            if(Rad >= 2.75_dp)then
                smooth_IC_Airfoil = fs 
            else if(Rad >= 0.25_dp)then
                smooth_IC_Airfoil = 0.4_dp*fs*Rad - 0.1_dp*fs
            else
                smooth_IC_Airfoil = 0.0_dp
            end if
        case (3)
            fs = rho_in*fmach*sin(alpha*pi/180._dp)
            if(Rad >= 2.75_dp)then
                smooth_IC_Airfoil = fs 
            else if(Rad >= 0.25_dp)then
                smooth_IC_Airfoil = 0.4_dp*fs*Rad - 0.1_dp*fs
            else
                smooth_IC_Airfoil = 0.0_dp
            end if
        case (4)
            smooth_IC_Airfoil = (P_in / (gamma - 1.0_dp))  + 0.5_dp*rho_in*(fmach**2.0) 
    end select
    
    return
end function smooth_IC_Airfoil

real(dp) function isentropicVortex(eulerVar,x,y)
    integer(i4),intent(in) :: eulerVar
    real(dp),intent(in)    :: x,y
    real(dp) :: r,beta,x0,y0,time
    real(dp) :: xmut,ymvt,u,v,rho,p
    
   
    beta = 5._dp
    x0 = 5._dp
    y0 = 5._dp
    time = 0._dp                   
    r  = sqrt((x-x0)**2 + (y-y0)**2) 
     
    rho = 1._dp
    u = 0._dp
    v = 0.1_dp
    p = 1._dp

    xmut = x-u*time
    ymvt = y-v*time
    r = sqrt((xmut-x0)**2 + (ymvt-y0)**2)


    u   = u - beta*exp(1._dp-r**2)*(ymvt-y0)/(2._dp*pi)
    v   = v + beta*exp(1._dp-r**2)*(xmut-x0)/(2._dp*pi)
    rho_in = (1._dp - ((gamma-1._dp)*beta**2*exp(2._dp*(1._dp-r**2))/(16._dp*gamma*pi*pi)))**(1._dp/(gamma-1._dp))

    P_in = rho_in**gamma
  
    select case (eulerVar)
        case (1)
            if(r < 6.1_dp)then
                isentropicVortex = rho_in
            else
                isentropicVortex = rho
            end if
        case (2)
            if(r < 6.1_dp)then
                isentropicVortex = rho_in*u
            else
                isentropicVortex = rho*u
            end if
        case (3)
            if(r < 6.1_dp)then
                isentropicVortex = rho_in*v
            else
                isentropicVortex = rho*v
            end if
        case (4)
            if(r < 6.1_dp)then
                isentropicVortex = P_in/(gamma-1._dp) + 0.5*rho_in*(u**2 + v**2);
            else
                isentropicVortex = P_in/(gamma-1._dp) + 0.5*rho*(u**2 + v**2);
            end if
    end select
    
    return
    end function isentropicVortex

end module getIC_module