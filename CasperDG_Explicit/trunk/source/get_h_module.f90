Module get_h_module
    use my_kinddefs
    use Globals_module, only: triList,numTri
    use inputs_module, only: basisDegree
    implicit none
    
    real(dp),pointer:: h(:)
contains

Subroutine get_h
    integer(i4) :: elem_id 
    real(dp) :: per
    
    Allocate(h(numTri))
    h = 1000._dp
    do elem_id = 1,numTri
        per = triList(elem_id)%perimeter
        h(elem_id)= triList(elem_id)%area / (per*(basisDegree+1._dp)**2)
    end do
    
end subroutine get_h

end module