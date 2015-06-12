Module vtu_output_module
    use Globals_module
    use my_kinddefs
    use projection_module
    use Globals_module, only: nodeList,numEdges,edgeList,boundEdges,interiorEdges,bcFlag
    use constructGridData_module, only: extCount
    use initializeBasis_module, only: solutionPhi
    implicit none
     
contains


!****************************************************************************80
!   write_grid: writes just a file for the grid
!****************************************************************************80
  subroutine plotSolution(project,solCoeffs)
    implicit none 
    
    character(*),intent(in)                      :: project
    real(dp),    intent(in)                      :: solCoeffs(:,:,:)
    
    !dummy
    integer(i4) :: elem_id,boundaryTri,edgeCount,bdedge,bc_type,n1,n2,n3
    real(dp)    :: locQ(numEulerVars,3)
    logical     :: isBoundElem
    ! ----------------------------------- Project Solution ---------------------------------------!

    !---> Local Variables
    character(80) :: filename
    integer(i4)   :: iunit
    integer(i4)   :: n,io
    CHARACTER(LEN=90) :: file_num
    real(dp) :: u(3),v(3),w(3)
    
    write(file_num,'(I4.4)') printCount
    
    iunit = 2
    
    print*,'# Writing to File...',printCount
    printCount = printCount + 1
    filename = adjustr(trim(project)) // '_' // &
                    adjustr(trim(file_num)) // '.vtu'
    
    open(unit = iunit, file = filename, status = 'replace', iostat = io)
    if( io /= 0) then
       print*,'ERROR: Opening soln vtu file'
    else
       
       write(iunit,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
       write(iunit,*)' <UnstructuredGrid>'
       write(iunit,*)'  <Piece NumberOfPoints="',3*numTri,'" NumberOfCells="',numTri,'">'
       
       write(iunit,*)'   <PointData Scalars="rho">'
       write(iunit,*)'    <DataArray type="Float32" Name="rho" Format="ascii">'
       
       do elem_id = 1, numTri
          
         Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)

         !plot rho
         write(iunit,*)'    ',locQ(1,1),locQ(1,2),locQ(1,3)
       end do
       
       write(iunit,*)'    </DataArray>'

!       write(iunit,*)'    <DataArray type="Float32"  NumberOfComponents="3" Name="velocities" Format="ascii">'
!       do elem_id = 1, numTri
!           
!        Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
!        !plot velocities
!        u(1) = locQ(2,1)/locQ(1,1)
!        u(2) = locQ(2,2)/locQ(1,2)
!        u(3) = locQ(2,3)/locQ(1,3)
!        v(1) = locQ(3,1)/locQ(1,1)
!        v(2) = locQ(3,2)/locQ(1,2)
!        v(3) = locQ(3,3)/locQ(1,3)
!        w(1) = 0.0_dp
!        w(2) = 0.0_dp
!        w(3) = 0.0_dp
!        !write(iunit,*)'    ',locQ(2,1),locQ(3,1),w(1),locQ(2,2),locQ(3,2),w(2),locQ(2,3),locQ(3,3),w(3)
!        write(iunit,*)'    ',u(1),v(1),w(1),u(2),v(2),w(2),u(3),v(3),w(3)
!       end do
!       write(iunit,*)'    </DataArray>'
       
       write(iunit,*)'    <DataArray type="Float32" Name="u" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
            write(iunit,*)'    ',locQ(2,1)/locQ(1,1),locQ(2,2)/locQ(1,2),locQ(2,3)/locQ(1,3)
        end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'    <DataArray type="Float32" Name="v" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
            write(iunit,*)'    ',locQ(3,1)/locQ(1,1),locQ(3,2)/locQ(1,2),locQ(3,3)/locQ(1,3)
        end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'    <DataArray type="Float32" Name="w" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
            write(iunit,*)'    ',0.0,0.0,0.0
        end do
       write(iunit,*)'    </DataArray>'


       write(iunit,*)'    <DataArray type="Float32" Name="rhoE" Format="ascii">'
       do elem_id = 1, numTri
        Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
        !plot rhoE
        write(iunit,*)'    ',locQ(4,1),locQ(4,2),locQ(4,3)  
       end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'   </PointData>'
       write(iunit,*)'   <CellData>'
       write(iunit,*)'    <DataArray type="Int32" Name="BC_assignment" Format="ascii">'
       do elem_id = 1, numTri
           isBoundElem = .false.
           do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount) 
                boundaryTri = edgeList(bdedge)%e1
                if (elem_id == boundaryTri) then
                    isBoundElem = .true.
                    bc_type = bcFlag(bdedge)
                end if
           end do
           
           select case(isBoundElem)
               case(.true.)
                    write(iunit,*)'    ',bc_type 
               case(.false.)
                    write(iunit,*)'    ',0   
           end select
       end do
       write(iunit,*)'    </DataArray>'

       write(iunit,*)'    <DataArray type="Int32" Name="Element Number" Format="ascii">'
       do elem_id = 1, numTri
           write(iunit,*)'    ',elem_id      
       end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'   </CellData>'
       write(iunit,*)'   <Points>'
       write(iunit,*)'    <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
       
       do n = 1,numTri
        n1 = triList(n)%vtxList(1)
        n2 = triList(n)%vtxList(2)
        n3 = triList(n)%vtxList(3)
        
        write(iunit,*)'    ', nodeList(n1)%x, nodeList(n1)%y, 0.0 &
                            , nodeList(n2)%x, nodeList(n2)%y, 0.0 &
                            , nodeList(n3)%x, nodeList(n3)%y, 0.0
       end do 
       
       write(iunit,*) '    </DataArray>'
       write(iunit,*) '   </Points>'

       
       write(iunit,*) '   <Cells>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
       do n = 1,numTri
           
           write(iunit,*)'    ',3*n-3,3*n-2,3*n-1
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="offsets" Format="ascii">'
       do n = 1,numTri
            write(iunit,*)'    ',3*n
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="types" Format="ascii">'
          do n = 1,numTri
            write(iunit,*)'    ',5
       end do
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '   </Cells>'
       write(iunit,*) '  </Piece>'
       write(iunit,*) ' </UnstructuredGrid>'
       write(iunit,*) '</VTKFile>'
    end if
       close(iunit) 
  end subroutine plotSolution
  
  
  
  !_______________________________________ 6 pt solution plot __________________________________________!
  
  subroutine plotSolution6pt(project,solCoeffs)
    implicit none 
    
    character(*),intent(in)                      :: project
    real(dp),    intent(in)                      :: solCoeffs(:,:,:)
    
    !dummy
    integer(i4) :: elem_id,edgeCount,n1,n2,n3
    real(dp)    :: locQ(numEulerVars,6)

    ! ----------------------------------- Project Solution ---------------------------------------!

    !---> Local Variables
    character(80) :: filename
    integer(i4)   :: iunit
    integer(i4)   :: n,io
    CHARACTER(LEN=90) :: file_num
    real(dp) :: u(6),v(6),rho(6),rhoE(6)
    real(dp) :: x1,x2,x3,x4,x5,x6
    real(dp) :: y1,y2,y3,y4,y5,y6
    
    write(file_num,'(I4.4)') printCount
    
    iunit = 2
    
    print*,'# Writing to File...',printCount
    printCount = printCount + 1
    filename = adjustr(trim(project)) // '_' // &
                    adjustr(trim(file_num)) // '.vtu'
    
    open(unit = iunit, file = filename, status = 'replace', iostat = io)
    if( io /= 0) then
       print*,'ERROR: Opening soln vtu file'
    else
       
       write(iunit,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
       write(iunit,*)' <UnstructuredGrid>'
       write(iunit,*)'  <Piece NumberOfPoints="',12*numTri,'" NumberOfCells="',4*numTri,'">'
       
       write(iunit,*)'   <PointData Scalars="rho">'
       write(iunit,*)'    <DataArray type="Float32" Name="rho" Format="ascii">'
       
       do elem_id = 1, numTri
          
         Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,6 nodes per tri)
         rho(1) = locQ(1,1)
         rho(2) = locQ(1,2)
         rho(3) = locQ(1,3)
         rho(4) = locQ(1,4)
         rho(5) = locQ(1,5)
         rho(6) = locQ(1,6)
         !plot rho
         write(iunit,*)'    ',rho(1),rho(4),rho(6)
         write(iunit,*)'    ',rho(4),rho(2),rho(5)
         write(iunit,*)'    ',rho(6),rho(5),rho(3)
         write(iunit,*)'    ',rho(4),rho(5),rho(6)
       end do
       
       write(iunit,*)'    </DataArray>'

       
       write(iunit,*)'    <DataArray type="Float32" Name="u" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,6 nodes per tri)
            u(1) = locQ(2,1)/locQ(1,1)
            u(2) = locQ(2,2)/locQ(1,2)
            u(3) = locQ(2,3)/locQ(1,3)
            u(4) = locQ(2,4)/locQ(1,4)
            u(5) = locQ(2,5)/locQ(1,5)
            u(6) = locQ(2,6)/locQ(1,6)
            write(iunit,*)'    ',u(1),u(4),u(6)
            write(iunit,*)'    ',u(4),u(2),u(5)
            write(iunit,*)'    ',u(6),u(5),u(3)
            write(iunit,*)'    ',u(4),u(5),u(6)
                                 
        end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'    <DataArray type="Float32" Name="v" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,6 nodes per tri)
            v(1) = locQ(3,1)/locQ(1,1)
            v(2) = locQ(3,2)/locQ(1,2)
            v(3) = locQ(3,3)/locQ(1,3)
            v(4) = locQ(3,4)/locQ(1,4)
            v(5) = locQ(3,5)/locQ(1,5)
            v(6) = locQ(3,6)/locQ(1,6)
            write(iunit,*)'    ',v(1),v(4),v(6)
            write(iunit,*)'    ',v(4),v(2),v(5)
            write(iunit,*)'    ',v(6),v(5),v(3)
            write(iunit,*)'    ',v(4),v(5),v(6)
           
        end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'    <DataArray type="Float32" Name="w" Format="ascii">'
        do elem_id = 1, numTri
            Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,6 nodes per tri)
            write(iunit,*)'    ',0.0,0.0,0.0
            write(iunit,*)'    ',0.0,0.0,0.0
            write(iunit,*)'    ',0.0,0.0,0.0
            write(iunit,*)'    ',0.0,0.0,0.0
        end do
       write(iunit,*)'    </DataArray>'


       write(iunit,*)'    <DataArray type="Float32" Name="rhoE" Format="ascii">'
       do elem_id = 1, numTri
        Call projectCell(solCoeffs(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,6 nodes per tri)
        rhoE(1) = locQ(4,1)
        rhoE(2) = locQ(4,2)
        rhoE(3) = locQ(4,3)
        rhoE(4) = locQ(4,4)
        rhoE(5) = locQ(4,5)
        rhoE(6) = locQ(4,6)
!         plot rho
        write(iunit,*)'    ',rhoE(1),rhoE(4),rhoE(6)
        write(iunit,*)'    ',rhoE(4),rhoE(2),rhoE(5)
        write(iunit,*)'    ',rhoE(6),rhoE(5),rhoE(3)
        write(iunit,*)'    ',rhoE(4),rhoE(5),rhoE(6)
                             
       end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'   </PointData>'

       write(iunit,*)'   <Points>'
       write(iunit,*)'    <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
       
       do n = 1,numTri
            n1 = triList(n)%vtxList(1)
            n2 = triList(n)%vtxList(2)
            n3 = triList(n)%vtxList(3)
            x1 = nodeList(n1)%x
            y1 = nodeList(n1)%y
            x2 = nodeList(n2)%x
            y2 = nodeList(n2)%y
            x3 = nodeList(n3)%x
            y3 = nodeList(n3)%y
            x4 = (x2 + x1) / 2._dp
            y4 = (y2 + y1) / 2._dp
            x5 = (x2 + x3)/ 2._dp
            y5 = (y2 + y3)/ 2._dp
            x6 = (x3 + x1)/ 2._dp
            y6 = (y3 + y1)/ 2._dp
            
        write(iunit,*) x1, y1, 0.0 &
                     , x4, y4, 0.0 & 
                     , x6, y6, 0.0
                            
        write(iunit,*) x4, y4, 0.0 &
                     , x2, y2, 0.0 &
                     , x5, y5, 0.0 
                            
        write(iunit,*) x6, y6, 0.0 &
                     , x5, y5, 0.0 &
                     , x3, y3, 0.0 
                            
        write(iunit,*) x4, y4, 0.0 &
                     , x5, y5, 0.0 &
                     , x6, y6, 0.0 
       end do 

       
       write(iunit,*) '    </DataArray>'
       write(iunit,*) '   </Points>'

       
       write(iunit,*) '   <Cells>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
       do n = 1,4*numTri  
               write(iunit,*)'    ',3*n-3,3*n-2,3*n-1  
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="offsets" Format="ascii">'
       do n = 1,4*numTri
            write(iunit,*)'    ',3*n
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="types" Format="ascii">'
          do n = 1,4*numTri
            write(iunit,*)'    ',5
       end do
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '   </Cells>'
       write(iunit,*) '  </Piece>'
       write(iunit,*) ' </UnstructuredGrid>'
       write(iunit,*) '</VTKFile>'
       
    end if
      close(iunit)
  end subroutine plotSolution6pt
  
  
  
  
  !___________________________________________ PLOT RESIDUAL _________________________________________!
  
  
  subroutine plotResidual(project,res)
    implicit none 
    
    character(*),intent(in)                      :: project
    real(dp),    intent(in)                      :: res(:,:,:) !res(numEulerVars,totalModes,numTri)
    
    !dummy
    integer(i4) :: elem_id,boundaryTri,edgeCount,bdedge,bc_type,n1,n2,n3
    real(dp) :: locQ(numEulerVars,3)
    logical  :: isBoundElem
    ! ----------------------------------- Project Solution ---------------------------------------!

    !---> Local Variables
    character(80) :: filename
    integer(i4)   :: iunit
    integer(i4)   :: n,io
    CHARACTER(LEN=90) :: file_num
    real(dp) :: eps
    eps = 1e-15

    
    write(file_num,'(I4.4)') printCountRes
    
    iunit = 2
    
    printCountRes = printCountRes + 1
    filename = adjustr(trim(project)) // '_Res_' // &
                    adjustr(trim(file_num)) // '.vtu'
    
    open(unit = iunit, file = filename, status = 'replace', iostat = io)
    if( io /= 0) then
       print*,'ERROR: Opening soln vtu file'
    else
       
       write(iunit,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
       write(iunit,*)' <UnstructuredGrid>'
       write(iunit,*)'  <Piece NumberOfPoints="',3*numTri,'" NumberOfCells="',numTri,'">'
       
       write(iunit,*)'   <PointData Scalars="rho">'
       write(iunit,*)'    <DataArray type="Float32" Name="Res_Rho" Format="ascii">'
       
       do elem_id = 1, numTri
         Call projectCell(res(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
         write(iunit,*)'    ',abs(locQ(1,1))+eps,abs(locQ(1,2))+eps,abs(locQ(1,3))+eps
       end do
       write(iunit,*)'    </DataArray>'
       
       write(iunit,*)'    <DataArray type="Float32" Name="Res_U" Format="ascii">'
       do elem_id = 1, numTri
         Call projectCell(res(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
         write(iunit,*)'    ',abs(locQ(2,1))+eps,abs(locQ(2,2))+eps,abs(locQ(2,3))+eps
       end do
       write(iunit,*)'    </DataArray>'

      write(iunit,*)'    <DataArray type="Float32" Name="Res_V" Format="ascii">'
       do elem_id = 1, numTri
         Call projectCell(res(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
         write(iunit,*)'    ',abs(locQ(3,1))+eps,abs(locQ(3,2))+eps,abs(locQ(3,3))+eps
       end do
       write(iunit,*)'    </DataArray>'
       
       write(iunit,*)'    <DataArray type="Float32" Name="Res_rhoE" Format="ascii">'
       do elem_id = 1, numTri
         Call projectCell(res(:,:,elem_id),solutionPhi,locQ) ! locQ(numEulerVars,3 nodes per tri)
         write(iunit,*)'    ',abs(locQ(4,1))+eps,abs(locQ(4,2))+eps,abs(locQ(4,3))+eps
       end do
       write(iunit,*)'    </DataArray>'
!       
       write(iunit,*)'   </PointData>'
       write(iunit,*)'   <CellData>'
       write(iunit,*)'    <DataArray type="Int32" Name="BC_assignment" Format="ascii">'
       do elem_id = 1, numTri
           isBoundElem = .false.
           do edgeCount = 1, extCount-1
                bdedge =  boundEdges(edgeCount) 
                boundaryTri = edgeList(bdedge)%e1
                if (elem_id == boundaryTri) then
                    isBoundElem = .true.
                    bc_type = bcFlag(bdedge)
                end if
           end do
           select case(isBoundElem)
               case(.true.)
                    write(iunit,*)'    ',bc_type
               case(.false.)
                    write(iunit,*)'    ',0   
           end select
       end do
       write(iunit,*)'    </DataArray>'
       
       write(iunit,*)'    <DataArray type="Int32" Name="Per Tri Twins" Format="ascii">'
       do elem_id = 1, numTri
        write(iunit,*)'    ',triList(elem_id)%periodic_Tri_Twin
       end do
       write(iunit,*)'    </DataArray>'
       
       write(iunit,*)'    <DataArray type="Int32" Name="Element Number" Format="ascii">'
       do elem_id = 1, numTri
           write(iunit,*)'    ',elem_id      
       end do
       write(iunit,*)'    </DataArray>'
       write(iunit,*)'   </CellData>'
       write(iunit,*)'   <Points>'
       write(iunit,*)'    <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
       
       do n = 1,numTri
        n1 = triList(n)%vtxList(1)
        n2 = triList(n)%vtxList(2)
        n3 = triList(n)%vtxList(3)
        
        write(iunit,*)'    ', nodeList(n1)%x, nodeList(n1)%y, 0.0 &
                            , nodeList(n2)%x, nodeList(n2)%y, 0.0 &
                            , nodeList(n3)%x, nodeList(n3)%y, 0.0
       end do 
       
       write(iunit,*) '    </DataArray>'
       write(iunit,*) '   </Points>'

       
       write(iunit,*) '   <Cells>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
       do n = 1,numTri
           
           write(iunit,*)'    ',3*n-3,3*n-2,3*n-1
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="offsets" Format="ascii">'
       do n = 1,numTri
            write(iunit,*)'    ',3*n
       end do 
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '    <DataArray type="Int32" Name="types" Format="ascii">'
          do n = 1,numTri
            write(iunit,*)'    ',5
       end do
       write(iunit,*) '    </DataArray>'
       
       write(iunit,*) '   </Cells>'
       write(iunit,*) '  </Piece>'
       write(iunit,*) ' </UnstructuredGrid>'
       write(iunit,*) '</VTKFile>'
       
    end if
    close(iunit)
  end subroutine plotResidual
  
  !_______________________________________ Write Solution Coeffs __________________________________!
  
  Subroutine write_Solution_Coefficients(project,solCoeffs)
    use inputs_module, only: totalModes
    implicit none 
    
    character(*),intent(in)                      :: project
    real(dp),    intent(in)                      :: solCoeffs(:,:,:) !res(numEulerVars,totalModes,numTri)
    
    !---> Local Variables
    character(80) :: filename
    integer(i4)   :: iunit
    integer(i4)   :: io,tm,elem_id
    CHARACTER(LEN=90) :: file_num

    
    write(file_num,'(I4.4)') printCountCoeffs
    
    iunit = 2
    
    printCountCoeffs = printCountCoeffs + 1
    filename = adjustr(trim(project)) // '_Coef_' // &
                    adjustr(trim(file_num)) // '.txt'
    
    open(unit = iunit, file = filename, status = 'replace', iostat = io)
    if( io /= 0) then
       print*,'ERROR: Opening soln vtu file'
    else
        write(iunit,*)' Solution Coefficients for file: ',filename   
        write(iunit,*) numTri,totalModes
        do elem_id = 1, numTri
            do tm = 1,totalModes
                write(iunit,*) solCoeffs(1,tm,elem_id),solCoeffs(2,tm,elem_id),&
                               solCoeffs(3,tm,elem_id),solCoeffs(4,tm,elem_id)
            end do
        end do
    end if
    close(iunit)
  end subroutine write_Solution_Coefficients
end module vtu_output_module