Module newtonRHS_module
    use my_kinddefs
    use Globals_module, only: numTri,numFields,triList,interiorEdges,edgeList
    use inputs_module,  only: totalModes
    implicit none
    
contains

Subroutine newtonRHS(edge_id,jacobOffDiagRL,jacobOffDiagLR,solCoeffs,newtonRHSvec)
    integer(i4),intent(in)  :: edge_id
    real(dp), intent(in)    :: jacobOffDiagRL(numFields*totalModes,numFields*totalModes)
    real(dp), intent(in)    :: jacobOffDiagLR(numFields*totalModes,numFields*totalModes)
    real(dp),intent(in)     :: solCoeffs(numFields*totalModes,numTri)
    real(dp),intent(out)    :: newtonRHSvec(numFields*totalModes,numTri)
    
    integer(i4) :: edgeNum,leftTri,rightTri
    
    
    edgeNum = interiorEdges(edge_id)
    leftTri  = edgeList(edgeNum)%e1    
    rightTri = edgeList(edgeNum)%e2
    
    newtonRHSvec(:,leftTri) = newtonRHSvec(:,leftTri) &
                               + matmul(jacobOffDiagLR,solCoeffs(:,rightTri))
                              
    newtonRHSvec(:,rightTri) = newtonRHSvec(:,rightTri) &
                               + matmul(jacobOffDiagRL,solCoeffs(:,leftTri))                                                  
    
end subroutine

end module newtonRHS_module
