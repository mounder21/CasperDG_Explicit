Module Globals_module
    use my_kinddefs
    use grid_data_type_module, only : node_type, elm_type, edge_type, bgrid_type
    implicit none 
    
!    public :: numNodes, nodeList 
!    public :: numTri, numQuad, numElms, triList
!    public :: numEdges, edgeList
!    public :: numBoundTypes, boundEdges,interiorEdges
!    public :: inlet
!    public :: numFields ,printCount,periodicBCs
!    public :: solutionNorm
    
 !  Solution Data

    real(dp), dimension(4) :: inlet
    integer(i4),parameter :: numFields = 4,dim = 2
    real(dp)    :: solutionNorm(numFields)
    logical     :: periodicBCs
    integer(i4) :: printCount = 0
    integer(i4) :: printCountRes = 0
    integer(i4) :: printCountCoeffs = 0

!  Node data
   integer(i4)                             :: numNodes !total number of nodes
   type(node_type), dimension(:), pointer  :: nodeList !array of nodes

!  Element data (element=cell)
   integer(i4)                                 :: numTri    !total number of triangles
   integer(i4)                                 :: numQuad   !total number of quadrilaterals
   integer(i4)                                 :: numElms   !total number of elements
   type(elm_type),  dimension(:), pointer      :: triList   !array of elements

!  Edge data
   integer(i4)                                 :: numEdges      !total number of edges
   type(edge_type), dimension(:), pointer      :: edgeList      !array of edges
   integer(i4),dimension(:), pointer           :: boundEdges     !array of boundary edges
   integer(i4),dimension(:), pointer           :: bcFlag
   integer(i4),dimension(:), pointer           :: interiorEdges  !array of interior edges
   integer(i4)                                 :: numInterior
!  Boundary data
   integer(i4)                                 :: numBoundTypes !total number of boundary types
   
end module Globals_module