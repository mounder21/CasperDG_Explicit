Module inputs_module
    use my_kinddefs
    implicit none 

    character(LEN=80)   :: mesh_name        ! Name of the mesh file
    character(LEN=80)   :: solution_name    ! Name of solution file to initialize solution 
    character(LEN=80)   :: out_file         ! Name of the .tec file to write
    character(LEN=80)   :: solver           ! Name of implicit solver
    character(LEN=80)   :: fluxScheme       ! Name of implicit solver
    integer(i4)         :: out_freq         ! Output file frequency
    logical             :: sourceTerm,writeCoeffs,MMS_sin,abaFile,Newton
    
    real(dp),parameter  :: gamma = 1.4_dp   ! Ratio of specific heat for gas
    real(dp)            :: fmach            ! Free-stream Mach number
    real(dp)            :: alpha            ! angle of attack in degrees
    real(dp)            :: CFL              ! CFL number
    
    integer(i4)             :: basisDegree      ! basis function polynomial order
    integer(i4),parameter   :: mapBasisDegreeStr8 = 1
    integer(i4)             :: mapBasisDegreeCurved
    integer(i4)             :: totalModes
    integer(i4),parameter   :: mapTotalModesStr8 = 3
    integer(i4)             :: mapTotalModesCurved
    character(20)           :: basisType
    character(20)           :: mapBasisType
    integer(i4)             :: numPlotPts
    
    real(dp)            :: rho_in 
    real(dp)            :: p_in 
    real(dp)            :: u_in 
    real(dp)            :: v_in 
    real(dp)            :: a,b,c,d,s0,t0
    
end module inputs_module
 !------------------------------------------------------------------------------------------------------------!

module grid_data_type_module

  use my_kinddefs

  implicit none

  private

  public ::  node_type
  public ::   elm_type
  public ::  edge_type
  public :: bgrid_type
 

!Data type for nodal quantities
  type node_type
!  to be read from a grid file
   real(dp)                         :: x, y     !nodal coordinates
!  to be constructed in the code
   integer(i4)                      :: numNghbrs    !number of neighbors
   integer(i4),dimension(:),pointer :: nghbrList    !list of neighbors
   integer(i4)                      :: numElms      !number of elements
   integer(i4),dimension(:),pointer :: elmList      !list of elements
   
!  to be computed in the conservative variables (rho,rhou,rhov,rhoE)
   real(dp)                         :: Pressure     !Pressure
   real(dp),  dimension(4)          :: Q            !conservative variables (rho,rhou,rhov,rhoE)
   real(dp),  dimension(4)          :: dQ           !change in conservative variables 
   real(dp),  dimension(4)          :: res          !residual (rhs)
   real(dp),  dimension(4)          :: primQ        !primitive variables (rho,u,v,E)
   real(dp),  dimension(4,2)        :: gradPrimQ    !gradient of primU
   real(dp),  dimension(4)          :: primQ_exact  !exact solutions (primitive)
   real(dp)                         :: dt           !local time step
   real(dp)                         :: wsn          !Half the max wave speed at face

  end type node_type


!Data type for element quantities
  type elm_type
!  to be read from a grid file
   integer(i4)                          :: numVtx       !number of vertices
   integer(i4),   dimension(:), pointer :: vtxList      !list of vertices
   
   integer(i4),   dimension(3)          :: nghbrList    !list of neighbors
   integer(i4),   dimension(3)          :: edgeLocList  !list of edge
   
   real(dp)                             :: area
   real(dp)                             :: perimeter,crossProduct
   real(dp)                             :: l1,l2,l3,l4
   character(80)                        :: elemType
   logical                              :: straight
   integer(i4)                          :: periodic_Tri_Twin,twin_Face
   logical                              :: used
  end type elm_type


!Data type for edge quantities
  type edge_type
!  to be constructed in the code
   integer(i4)                      :: n1, n2 !associated nodes
   integer(i4)                      :: e1, e2 !associated elements
   integer(i4)                      :: locEdge1,locEdge2 ! left and right elements

  end type edge_type


!Data type for boundary quantities
  type bgrid_type
!  to be read from a boundary grid file
   character(80)                    :: bc_type      !type of boundary condition
   integer(i4)                      :: numbnodes    !# of boundary nodes
   integer(i4),dimension(:),pointer :: bnodeList    !list of boundary nodes

  end type bgrid_type


 end module grid_data_type_module
