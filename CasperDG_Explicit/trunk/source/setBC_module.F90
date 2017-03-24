Module setBC_module
    use my_kinddefs
    use inputs_module, only: mesh_name,alpha
    use Globals_module, only: boundEdges,edgeList,nodeList,bcFlag
    implicit none
    contains
    
Subroutine setBC(extCount)
    integer(i4),intent(in) :: extCount
    
    integer(i4) :: edgeCount,bdedge,n1,n2
    real(dp)    :: x1,x2,y1,y2,aveX,aveY,Rad
!Construct boundary types
!assign based on location: inviscid wall(2), sonic flow(5)

select case(mesh_name)
        case('meshes/test.msh')
             bcFlag(:) = 5
!             do edgeCount = 1, extCount-1
!                bdedge =  boundEdges(edgeCount)
!                if (bdedge == 14 .or. bdedge == 18) then
!                     bcFlag(bdedge) = 2
!                end if
!             end do
           
        case('meshes/naca_10.msh')
            print*,"entered mesh"
            bcFlag(:) = 5
            do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
                !assign based on location: inviscid wall(2), sonic flow(5)
                ! [-4,-4] to [5,4]
                !left and right edges are super sonic
                aveX = 0.5_dp *(x1 + x2)
                aveY = 0.5_dp *(y1 + y2)
                Rad  = sqrt(aveX**2 + aveY**2) 
                
!                if(alpha == 0.0_dp)then
!                    if(y1 < -1.9 .or. y2 < -1.9) then
!                        bcFlag(bdedge) = 2
!                    else if(y1 > 1.9 .or. y2 > 1.9 ) then
!                        bcFlag(bdedge) = 2
!                    end if
!                    if(x1 > 2.9 .or. x2 > 2.9) then
!                        bcFlag(bdedge) = 8
!                    end if
!                    if(x1 < -2.9 .or. x2 < -2.9) then
!                        bcFlag(bdedge) = 5
!                    end if
!                end if
                
                if (Rad<1.0_dp) then         
                    bcFlag(bdedge) = 2
                end if   
            end do
            
        case('meshes/naca.msh') 
            print*,"entered mesh"
            bcFlag(:) = 5
            do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
                !assign based on location: inviscid wall(2), sonic flow(5)
                ! x[-3,3] to y[-2,2]
                !left and right edges are super sonic
                aveX = 0.5_dp *(x1 + x2)
                aveY = 0.5_dp *(y1 + y2)
                Rad  = sqrt(aveX**2 + aveY**2) 
                
                if(alpha == 0.0_dp)then
                    if(y1 < -1.9 .or. y2 < -1.9) then
                        bcFlag(bdedge) = 5
                    else if(y1 > 1.9 .or. y2 > 1.9 ) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 > 2.9 .or. x2 > 2.9) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 < -2.9 .or. x2 < -2.9) then
                        bcFlag(bdedge) = 5
                    end if
                end if
                
                if (Rad<1.0_dp) then         
                    bcFlag(bdedge) = 2
                end if   
            end do
            
        case('gmeshes/naca0012.msh')   
            bcFlag(:) = 5
            do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
                !assign based on location: inviscid wall(2), sonic flow(5)
                ! x[-3,3] to y[-2,2]
                !left and right edges are super sonic
                aveX = 0.5_dp *(x1 + x2)
                aveY = 0.5_dp *(y1 + y2)
                Rad  = sqrt(aveX**2 + aveY**2) 
                
                if(alpha == 0.0_dp)then
                    if(y1 < -1.9 .or. y2 < -1.9) then
                        bcFlag(bdedge) = 5
                    else if(y1 > 1.9 .or. y2 > 1.9 ) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 > 2.9 .or. x2 > 2.9) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 < -2.9 .or. x2 < -2.9) then
                        bcFlag(bdedge) = 5
                    end if
                end if
                
                if (Rad<2.0_dp) then         
                    bcFlag(bdedge) = 2
                end if   
            end do
            
        case('meshes/bump.msh') 
            ![-5,0] to [6,10]
             bcFlag(:) = 5
             do edgeCount = 1, extCount-1
                bdedge = boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
              
                if (y1 < 1._dp .or. y2 < 1.0_dp) then    
                    if((x1 > -4.9 .or. x2 > -4.9).and.(x1 < 5.9 .or. x2 < 5.9))then
                        bcFlag(bdedge) = 2
                    end if
                end if   
            end do

        case('../meshes/bump_small.msh')
            ![-5,0] to [6,10]
             bcFlag(:) = 5
             do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
              
                if (y1 < 1._dp .or. y2 < 1.0_dp) then    
                    if((x1 > -4.9 .or. x2 > -4.9).and.(x1 < 5.9 .or. x2 < 5.9))then
                        bcFlag(bdedge) = 2
                    end if
                end if   
            end do
        case('../meshes/wall.msh')
            ![0,0] to [1.2,1.2]
             bcFlag(:) = 2
             do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
              
                if (x1 < 0.1_dp .or. x2 < 0.1_dp) then    
                    bcFlag(bdedge) = 5
                end if   
            end do
        case('../meshes/box_big.msh')
             bcFlag(:) = 2
             do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
              
                if (x1 <= 0.01_dp .and. x2 <= 0.01_dp) then    !left wall: periodic BC
                    bcFlag(bdedge) = 13    
                else if(x1 >= 9.99_dp .and. x2 >= 9.99_dp)then !right wall: characteristic BC
                    bcFlag(bdedge) = 13
                end if
                if(y1 >= 4.90_dp .and. y2 >= 4.90_dp) then     !top wall: inviscid wall
                    bcFlag(bdedge) = 2
                else if(y1 <= -4.99_dp .and. y2 <= -4.99_dp)then    !wall wall: inviscid wall
                    bcFlag(bdedge) = 2
                end if   
            end do
         case('../meshes/box.msh')
             bcFlag(:) = 2
             do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
              
                if (x1 <= 0.01_dp .and. x2 <= 0.01_dp) then    
                    bcFlag(bdedge) = 13
                else if(x1 >= 9.99_dp .and. x2 >= 9.99_dp)then
                    bcFlag(bdedge) = 13
                end if
                if(y1 >= 4.90_dp .and. y2 >= 4.90_dp) then
                    bcFlag(bdedge) = 2
                else if(y1 <= -4.99_dp .and. y2 <= -4.99_dp)then
                    bcFlag(bdedge) = 2
                end if   
            end do
        case default
            bcFlag(:) = 5
            do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount)
                n1 = edgeList(bdedge)%n1
                n2 = edgeList(bdedge)%n2
                !get node coordinates
                x1 = nodeList(n1)%x
                x2 = nodeList(n2)%x
                y1 = nodeList(n1)%y
                y2 = nodeList(n2)%y
                !assign based on location: inviscid wall(2), sonic flow(5)
                ! x[-3,3] to y[-2,2]
                !left and right edges are super sonic
                aveX = 0.5_dp *(x1 + x2)
                aveY = 0.5_dp *(y1 + y2)
                Rad  = sqrt(aveX**2 + aveY**2) 
                
                if(alpha == 0.0_dp)then
                    if(y1 < -1.9 .or. y2 < -1.9) then
                        bcFlag(bdedge) = 5
                    else if(y1 > 1.9 .or. y2 > 1.9 ) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 > 2.9 .or. x2 > 2.9) then
                        bcFlag(bdedge) = 5
                    end if
                    if(x1 < -2.9 .or. x2 < -2.9) then
                        bcFlag(bdedge) = 5
                    end if
                end if
                
                if (Rad<2.0_dp) then         
                    bcFlag(bdedge) = 2
                end if   
            end do
            
!             bcFlag(:) = 5 !inflow
!             do edgeCount = 1, extCount-1
!                bdedge =  boundEdges(edgeCount)
!                n1 = edgeList(bdedge)%n1
!                n2 = edgeList(bdedge)%n2
!                !get node coordinates
!                x1 = nodeList(n1)%x
!                x2 = nodeList(n2)%x
!                y1 = nodeList(n1)%y
!                y2 = nodeList(n2)%y
!              
!                if (x1 <= 0.01_dp .and. x2 <= 0.01_dp) then    !left wall: inflow
!                    bcFlag(bdedge) = 14    
!                else if(x1 >= 9.99_dp .and. x2 >= 9.99_dp)then !right wall: outflow
!                    bcFlag(bdedge) = 15
!                end if
!                if(y1 >= 9.90_dp .and. y2 >= 9.90_dp) then     !top wall: outflow
!                    bcFlag(bdedge) = 15
!                else if(y1 <= 0.01_dp .and. y2 <= 0.01_dp)then    !bottom wall: inflow
!                    bcFlag(bdedge) = 14
!                end if   
!            end do
            print*,'do not have mesh boundary types set up yet '
    end select
end subroutine setBC
end module setBC_module