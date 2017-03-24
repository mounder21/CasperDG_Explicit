Module constructGridData_module
    use my_kinddefs
    use Globals_module
    use inputs_module
    use setBC_module
    implicit none 
    integer(i4) :: extCount
contains 
!**************************************************************************************************************!
!* 1. Edge data:
!*    edgeList(:)%n1, n2  = End nodes of each edge (edge points n1 -> n2)
!*    edgeList(:)%e1, e2  = Left and right elements of each edge
!**************************************************************************************************************!


Subroutine constructGridData
    use Globals_module, only: numNodes,nodeList,triList,numTri,numEdges,edgeList,boundEdges,interiorEdges
    use Globals_module, only: numInterior,periodicBCs
    implicit none 
    
    integer(i4) :: minVtx_List(numNodes), edgeLinkedList(numEdges),edgeCount,inCount
    integer(i4) :: triNumber,j,v1,v2,v3,currentMinVtx,node1,node2,edge_number,currEdge,temp,tempEdge
    integer(i4) :: triNumberTemp,nodeTemp1,nodeTemp2,leftTri,rightTri
    real(dp)    :: a1,a2,b1,b2,l1,l2,l3,y1,y2,yTemp1,yTemp2
    logical :: coMatch

    ! Initialize array values 
    minVtx_List(:) = -1
    edgeLinkedList(:) = -1
    edgeList(1:numEdges)%e1 = -1    
    edgeList(1:numEdges)%e2 = -1
    edgeList(1:numEdges)%locEdge1 = -1
    edgeList(1:numEdges)%locEdge2 = -1
    
    edge_number = 1
    numInterior = 0
    
    triList(:)%periodic_Tri_Twin = 0
    triList(:)%nghbrList(1) = 0
    triList(:)%nghbrList(2) = 0
    triList(:)%nghbrList(3) = 0
    triList(:)%edgeLocList(1) = 0
    triList(:)%edgeLocList(2) = 0
    triList(:)%edgeLocList(3) = 0
    
    do triNumber = 1, numTri
        v1 = triList(triNumber)%vtxList(1)
        v2 = triList(triNumber)%vtxList(2)
        v3 = triList(triNumber)%vtxList(3)
        
        l3 = sqrt((nodeList(v2)%x - nodeList(v1)%x)**2 + (nodeList(v2)%y - nodeList(v1)%y)**2)
        l1 = sqrt((nodeList(v3)%x - nodeList(v2)%x)**2 + (nodeList(v3)%y - nodeList(v2)%y)**2)
        l2 = sqrt((nodeList(v1)%x - nodeList(v3)%x)**2 + (nodeList(v1)%y - nodeList(v3)%y)**2)
        
        !Compute cross product to ensure orientation
        a1 = nodeList(v3)%x - nodeList(v2)%x
        a2 = nodeList(v3)%y - nodeList(v2)%y
        b1 = nodeList(v1)%x - nodeList(v2)%x
        b2 = nodeList(v1)%y - nodeList(v2)%y
        
        triList(triNumber)%crossProduct = a1*b2 - a2*b1
        ! switch since it is orientated the wrong way
        if (triList(triNumber)%crossProduct < 0) then 
            print*,'cross product was negative'
            temp = triList(triNumber)%vtxList(1)
            triList(triNumber)%vtxList(1) = triList(triNumber)%vtxList(2)
            triList(triNumber)%vtxList(2) = temp
            triList(triNumber)%crossProduct = abs(triList(triNumber)%crossProduct)
        end if
        
        !area is 0.5 * cross product since were in 2D 
        triList(triNumber)%area  = 0.5_dp * (a1*b2 - a2*b1)
        triList(triNumber)%perimeter = l1 + l2 + l3
        triList(triNumber)%l1 = l1
        triList(triNumber)%l2 = l2
        triList(triNumber)%l3 = l3
        edgesPerTri : do j = 1, 3
        select case (j)
            case (1)
                node1 = v2
                node2 = v3    
            case (2)
                node1 = v3
                node2 = v1
            case (3)
                node1 = v1
                node2 = v2
        end select
        ! Get current min vertex
            currentMinVtx = min(node1,node2)
        ! Search minVtx_List
            currEdge = minVtx_List(currentMinVtx)
            
            if (currEdge .eq. -1) then !add the edge to the list, and add the tri num 
                edgeList(edge_number)%n1 = node1
                edgeList(edge_number)%n2 = node2
                edgeList(edge_number)%e1 = triNumber !left triangle   
                edgeList(edge_number)%locEdge1 = j
                minVtx_List(currentMinVtx) = edge_number
                edge_number = edge_number + 1
            else  
                call noMatchRoutine(node1,node2,currEdge,edge_number,triNumber,edgeLinkedList,j)
            end if

        end do edgesPerTri
    end do
    
    inCount = 1
    extCount = 1
    Allocate(interiorEdges(numInterior))
    Allocate(boundEdges(numEdges - numInterior))
    Allocate(bcFlag(numEdges))

    bcFlag(:) = 0
    
    do edgeCount = 1, numEdges
        if(edgeList(edgeCount)%e2 /= -1) then
            interiorEdges(inCount) = edgeCount
            inCount = inCount + 1
        else
            boundEdges(extCount) = edgeCount
            extCount = extCount + 1
        end if
        
        ! Assemble element neighbors and local edge list
        leftTri  = edgeList(edgeCount)%e1
        do j = 1,3
            if (triList(leftTri)%edgeLocList(j) == 0)then
                triList(leftTri)%edgeLocList(j) = edgeCount
                exit 
            end if
        end do   
        
        rightTri = edgeList(edgeCount)%e2
        if(rightTri /= -1)then
            do j = 1,3
                if((triList(rightTri)%edgeLocList(j) == 0))then
                    triList(rightTri)%edgeLocList(j) = edgeCount
                    exit 
                end if
            end do
        end if
        
        if(rightTri /= -1)then
            do j = 1,3
                if (triList(leftTri)%nghbrList(j) == 0)then
                    triList(leftTri)%nghbrList(j) = rightTri
                    exit 
                end if
            end do     
            do j = 1,3
                if (triList(rightTri)%nghbrList(j) == 0)then
                    triList(rightTri)%nghbrList(j) = leftTri
                    exit 
                end if
            end do 
        end if
    end do
    
    !set bc types
    Call setBC(extCount)
    
    if(periodicBCs) then
        !form periodic_Tri_twins (must have symmetric mesh)
        do edgeCount = 1,extCount-1
            currEdge = boundEdges(edgeCount) 
            !check if periodic
            if(bcFlag(currEdge) == 13)then       !inflow periodic 
                triNumber = edgeList(currEdge)%e1
                node1 = edgeList(currEdge)%n1
                node2 = edgeList(currEdge)%n2
                y1 = nodeList(node1)%y
                y2 = nodeList(node2)%y
                !now search other boundary triangles
                do inCount = 1,extCount-1
                    tempEdge = boundEdges(inCount)
                    triNumberTemp = edgeList(tempEdge)%e1
                    if((bcFlag(tempEdge) == 13).and.(triNumberTemp.ne.triNumber))then !inflow periodic 
                        nodeTemp1 = edgeList(tempEdge)%n1
                        nodeTemp2 = edgeList(tempEdge)%n2
                        yTemp1 = nodeList(nodeTemp1)%y
                        yTemp2 = nodeList(nodeTemp2)%y
                        !now check y coordinates
                        Call checkCoordinates(y1,y2,yTemp1,yTemp2,coMatch)
                        if(coMatch)then
                            triList(triNumber)%periodic_Tri_Twin = triNumberTemp
                            triList(triNumberTemp)%periodic_Tri_Twin = triNumber
                            triList(triNumberTemp)%twin_Face = edgeList(tempEdge)%locEdge1  
                            triList(triNumber)%twin_Face = edgeList(currEdge)%locEdge1
                        end if
                    end if
                end do
            end if 
        end do
    end if
    
end subroutine constructGridData


Subroutine checkNodes(node1,node2,nodeA,nodeB,match)
    integer(i4),intent(in) :: node1,node2,nodeA,nodeB
    logical, intent(out) :: match
    
    if (((node1 == nodeA).and.(node2 == nodeB)) .or. ((node1 == nodeB).and.(node2 == nodeA))) then
        match = .true.
    else
        match = .false.
    end if            
end subroutine checkNodes

Subroutine checkCoordinates(y1,y2,yTemp1,yTemp2,match)
    real(dp),intent(in) :: y1,y2,yTemp1,yTemp2
    logical, intent(out) :: match
    
    if (((y1 == yTemp1).and.(y2 == yTemp2)) .or. ((y1 == yTemp2).and.(y2 == yTemp1))) then
        match = .true.
    else
        match = .false.
    end if            
end subroutine checkCoordinates


recursive Subroutine noMatchRoutine(node1,node2,currEdge,edge_number,triNumber,edgeLinkedList,j)
    use Globals_module,only    : edgeList
    integer(i4),intent(inout) :: edgeLinkedList(numEdges),currEdge
    integer(i4),intent(inout) :: node1,node2,edge_number,triNumber
    integer(i4),intent(inout) :: j 
    logical                   :: match 
    integer(i4)               :: edgeLLVal
    !check for edge match at minVtx_List(currentMinVtx): n1,n2
    call checkNodes(node1,node2,edgeList(currEdge)%n1,edgeList(currEdge)%n2,match)
    
    if (match)  then
        ! we have a match so don't add an edge but add the right triangle 
        edgeList(currEdge)%e2 = triNumber               !right triangle
        edgeList(currEdge)%locEdge2 = j
        numInterior = numInterior + 1
        
    else ! no match
        edgeLLVal = edgeLinkedList(currEdge)
        if (edgeLLVal .eq. -1) then ! New edge, add the edge to the list, and add the tri num 
            edgeList(edge_number)%n1 = node1
            edgeList(edge_number)%n2 = node2 
            edgeLinkedList(currEdge) = edge_number
            edgeList(edge_number)%e1 = triNumber        !left triangle
            edgeList(edge_number)%locEdge1 = j
            edge_number = edge_number + 1
        else ! do not know if its a new edge, check for match
            currEdge = edgeLLVal
            call noMatchRoutine(node1,node2,currEdge,edge_number,triNumber,edgeLinkedList,j)
        end if
    end if
end subroutine noMatchRoutine

end module constructGridData_module