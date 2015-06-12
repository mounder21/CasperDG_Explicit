
Module read_mesh_module
    use my_kinddefs
    use Globals_module
    implicit none 
contains 

!**************************************************************************************************************!
Subroutine read_mesh
    use inputs_module,  only: mesh_name,abaFile      !(dataModules)
    use Globals_module, only: numNodes,nodeList,numTri,numEdges
    
    implicit none 
    
integer             :: funit,i,temp
character(LEN=30)   :: comments
real(dp) x,y

funit=2
OPEN(UNIT=funit,FILE=mesh_name,STATUS='old')

READ(funit,*) comments
READ(funit,*) numNodes, numTri, numEdges
    
call AllocateMesh   ! Located Below
if(abaFile)then   
    ! Read in the point locations
    do i=1,numNodes
        read(funit,*) temp,nodeList(i)%x , nodeList(i)%y
    end do

    ! Read in the triangle node list (3 per volume)
    do i=1,numTri
       triList(i)%numVtx = 3
       ALLOCATE(triList(i)%vtxList(3))
       read(funit,*) temp,triList(i)%vtxList(1), triList(i)%vtxList(2), triList(i)%vtxList(3)
       triList(i)%elemType = 'triangle'
       triList(i)%straight = .true.
    end do
else
    ! Read in the point locations
    do i=1,numNodes
        if((mesh_name .eq. '../meshes/box_big.msh').or.(mesh_name .eq. '../meshes/box.msh'))then
            read(funit,*)  x,y  
            nodeList(i)%x = 10._dp*x
            nodeList(i)%y = 10._dp*y -5._dp
        else
            read(funit,*) nodeList(i)%x , nodeList(i)%y
        end if
    end do

    ! Read in the triangle node list (3 per volume)
    do i=1,numTri
       triList(i)%numVtx = 3
       ALLOCATE(triList(i)%vtxList(3))
       read(funit,*) triList(i)%vtxList(1), triList(i)%vtxList(2), triList(i)%vtxList(3)
       triList(i)%elemType = 'triangle'
       triList(i)%straight = .true.
    end do
end if
end subroutine read_mesh


Subroutine read_gmesh
    use inputs_module,  only: mesh_name      !(dataModules)
    use Globals_module, only: numNodes,nodeList,numTri,numEdges
    
    implicit none 
    
integer             :: funit,i,temp
character(LEN=30)   :: comments
real(dp) x,y

funit=2
OPEN(UNIT=funit,FILE=mesh_name,STATUS='old')

READ(funit,*) comments
READ(funit,*) comments
READ(funit,*) comments
READ(funit,*) comments
READ(funit,*) numNodes
    
call AllocateMesh   ! Located Below
    ! Read in the point locations
do i=1,numNodes
    read(funit,*) temp , nodeList(i)%x , nodeList(i)%y , temp
end do

READ(funit,*) comments
READ(funit,*) numTri

    ! Read in the triangle node list (3 per volume)
    do i=1,numTri
       triList(i)%numVtx = 3
       ALLOCATE(triList(i)%vtxList(3))
       read(funit,*) temp,temp,temp, triList(i)%vtxList(1), triList(i)%vtxList(2), triList(i)%vtxList(3)
       triList(i)%elemType = 'triangle'
       triList(i)%straight = .true.
    end do

end subroutine read_gmesh

!**************************************************************************************************************!
Subroutine AllocateMesh
    use Globals_module
    
    allocate(nodeList(numNodes))
    allocate(triList(numTri))
    allocate(edgeList(numEdges))

end subroutine AllocateMesh
!**************************************************************************************************************!
end module read_mesh_module
