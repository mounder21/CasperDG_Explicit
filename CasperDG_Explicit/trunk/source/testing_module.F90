Module testing_module
    use my_kinddefs
    use residual_Jac_module
    use residualComplex_module
    implicit none
contains
  
#ifdef  complx
subroutine test_jacob
    use inputs_module,  only: totalModes
    use form_Jx_module
    use initializeSolution_module, only: solCoeffs
    
    real(dp) :: finDiff(numFields,totalModes),finDiff2(numFields*totalModes)
    real(dp) :: epsil,maxerr,MaxOverall
    real(dp) :: resA(numFields,totalModes),resB(numFields,totalModes)
    real(dp) :: R2(numFields*totalModes),R1(numFields*totalModes)
    
    complex(dp) :: solComplex(numFields,totalModes,numTri)
    complex(dp) :: h
    
    real(dp) :: totalCalcDRes(numFields*totalModes,numFields*totalModes)
    integer(i4) :: k,n,triangleNum,edgeNum,tm,i,tri_id,edge_id
    real(dp)    :: Jx(numFields*totalModes,numTri)
    real(dp)    :: x(numFields*totalModes,numTri),xtemp(numFields,totalModes,numTri)
    
    
    MaxOverall = 0
    Call Residual_Jac(solCoeffs)

    x = reshape(solCoeffs,(/numFields*totalModes,numTri/))
    xtemp = solCoeffs
    
    Call form_Jx(.false.,jacobDiag,jacobOffDiagRL,jacobOffDiagLR,x,Jx)
    
    epsil = 1E-30_qp
    h = (0._qp,1E-30_qp)
    solComplex = solCoeffs + h*xtemp
    call ResidualComplex(solComplex)
        
    do tri_id = 1,numTri
        finDiff  = Imag(spResRC(:,:,tri_id))/epsil
        finDiff2 = reshape(finDiff,(/numFields*totalModes/))
        
        maxerr = 0
        maxerr = max(maxerr,maxval(abs(finDiff2 - Jx(:,tri_id))))
        
        !print*, finDiff2
        !print*, Jx(:,tri_id)
        print*,'error:      ', maxerr,tri_id
        MaxOverall = max(MaxOverall,maxerr)
        
    end do
    print*,'Max err:',MaxOverall

  end subroutine
  
  
    
  subroutine test_jacobian
      use bc_module
      real(dp) :: ql(4,3),qr(4,3),qrA(4),qrB(4)
      real(dp) :: flux1(4),flux2(4),norm(2),qp(4),flux1a(4,2),flux1b(4,2),normal(3,2)
      real(dp) :: fluxa(4),fluxb(4)
      real(dp) :: flux_ql(4,4),flux_qr(4,4)!,ua,u_qla(4),ub,u_qlb(4),flux_q(4,2,4)
      real(dp) :: fluxBC_q(4,4)
      real(dp) :: finDiffFlux(4)
      real(dp) :: dqbdq(4,4)
      real(dp) :: epsil
      integer(i4) :: k
      
      ql(1,1) = 1._dp
      ql(2,1) = 0.3_dp
      ql(3,1) = 0.2_dp
      ql(4,1) = ql(1,1)/0.4_dp + 0.5_dp*((ql(2,1)/ql(1,1))**2 + (ql(3,1)/ql(1,1))**2)
      
      qr(1,1) = 1.0_dp
      qr(2,1) = 0.3_dp
      qr(3,1) = 0.1_dp
      qr(4,1) = 10.1_dp
     
      ql(1,2) = 1.0_dp
      ql(2,2) = 12._dp
      ql(3,2) = 13.3_dp
      ql(4,2) = 23.1_dp
      
      qr(1,2) = 1.0_dp
      qr(2,2) = 0.3_dp
      qr(3,2) = 0.1_dp
      qr(4,2) = 10.1_dp
      
      ql(1,3) = 1.0_dp
      ql(2,3) = 12._dp
      ql(3,3) = 13.3_dp
      ql(4,3) = 23.1_dp
      
      qr(1,3) = 1.0_dp
      qr(2,3) = 0.3_dp
      qr(3,3) = 0.1_dp
      qr(4,3) = 10.1_dp
      
      norm(1) = 2._dp
      norm(2) = 1._dp
      
      normal(1,1) = 2._dp
      normal(1,2) = 1._dp
      
      normal(2,1) = 2._dp
      normal(2,2) = 1._dp
      
      normal(3,1) = 2._dp
      normal(3,2) = 1._dp
      
      epsil = 1E-7_dp
      
      !call get_bc_q(2, norm, ql(:,1), qrA,dqbdq)
      !Call getBCFlux_q(qR,normals(:,:,edgeNum),numEdgeGaussPts,flux,flux_qr) !flux(4,numEdgeGaussPts)
      call get_bc_q(5,norm, ql(:,1), qrB, dqbdq)          
      call getBCFlux_q(qrB,norm,1,fluxa,fluxBC_q)
      
      !print*, fluxa
      !call getNativeFlux_q(ql(:,2),flux1a,flux_q)
      !call Rotated_RHLL(qL, qR, norm(1), norm(2) ,flux1,flux_ql,flux_qr)!,ua,u_qla)
      !call lax_friedrichs_q(ql, qr, norm, flux1,flux_ql,flux_qr)!,ua,u_qla)
      !print*,ql(:,1)
      
     do k = 1,4
         qp(:) = ql(:,1)
         qp(k) = qp(k) + epsil
         
         call get_bc_q (3,norm, qp, qrB, dqbdq)
         !print*,qrB
         call getBCFlux_q(qrB,norm,1,fluxb,fluxBC_q)
         !print*, fluxb
         
         flux_ql(:,:) = matmul(fluxBC_q,dqbdq)
         finDiffFlux = (fluxb - fluxa)/epsil
         
         print*,'fin        :', finDiffFlux
         print*,'calculated :', flux_ql(:,k)
         print*,' ' 
         
         !print*,dqbdq
         !call getBCFlux_q(qp,normal,3,fluxb,fluxBC_q)
         !call getNativeFlux_q(qp(:,2),flux1b,flux_q)
         !call Rotated_RHLL_q(qp, qr, norm(1), norm(2) ,flux2,flux_ql,flux_qr)!,ub,u_qlb)
         !call lax_friedrichs_q(ql, qp, norm, flux2,flux_ql,flux_qr)!,ub,u_qlb)
         
         !finDiffFlux = (flux1b(:,2) - flux1a(:,2))/epsil
         
         !print*,'flux a', fluxa(k,1)
         !print*,'flux b', fluxb
         !finDiffFlux = (qrB - qrA)/epsil
         
         
         !print*,'calc:', flux_q(:,2,k)
         !finDiffFlux = (flux2 - flux1) / epsil
         !print*,'finite diff:', finDiffFlux
         !print*,'calculated :', fluxBC_q(:,1,k)
         !print*,'error:      ', sum(finDiffFlux - fluxBC_q(:,2,k))
     end do
     
  end subroutine
  
  subroutine test_Offjacobian
    use initializeSolution_module, only: solCoeffs
    use inputs_module,  only: totalModes

    
    real(dp)    :: finDiff(4,totalModes)
    real(dp)    :: epsil
    real(dp)    :: solCoeffsA(4,totalModes,numTri)
    real(dp)    :: totalCalcDRes(4*totalModes,4*totalModes)
    integer(i4) :: k,n,edgeNum,numSides,intLocEdges(3),edge_id,leftTri,rightTri,i,tm,tri_id
    real(dp)    :: resA(3,4,totalModes),resB(3,4,totalModes),maxerr,MaxOverall
    logical     :: perturbTri
    
    complex(dp) :: solComplex(numFields,totalModes,numTri)
    complex(dp) :: h
    
    solComplex = solCoeffs
    
    MaxOverall = 0
    epsil = 1E-30_dp
    h = (0._dp,1E-30_dp)
    
    do tri_id = 1,1
        totalCalcDRes = 0._dp
        numSides = 0

        call Residual_JacTriangle(solCoeffs,tri_id)
        
        do n = 1,3           
            edgeNum = triList(tri_id)%edgeLocList(n)
            if(edgeList(edgeNum)%e2 .ne. -1)then
                numSides = numsides + 1
                intLocEdges(numSides) = triList(tri_id)%edgeLocList(n)
            end if
        end do

        do edge_id = 1, numSides
            edgeNum  =  intLocEdges(edge_id)
            leftTri  = edgeList(edgeNum)%e1    
            rightTri = edgeList(edgeNum)%e2
            if(tri_id == leftTri)then
                totalCalcDRes = jacobOffDiagLRTri(:,:,edgeNum)
            elseif (tri_id == rightTri)then
                totalCalcDRes = jacobOffDiagRLTri(:,:,edgeNum)
            end if
        end do

        perturbTri = .true.
        do tm = 1,totalModes
            do k = 1,4
                maxerr = 0
                i = (tm - 1)*4 + k 
                do edge_id = 1, numSides
                    edgeNum  = intLocEdges(edge_id)
                    leftTri  = edgeList(edgeNum)%e1    
                    rightTri = edgeList(edgeNum)%e2

                    call ResidualComplexTriangle(solComplex,tri_id,perturbTri,h,k,tm,edge_id)
                    
                    if(tri_id == leftTri)then
                        finDiff =  Imag(spResRCTri(:,:,tri_id))/epsil
                    elseif (tri_id == rightTri)then
                        finDiff =  Imag(spResRCTri(:,:,tri_id))/epsil
                    end if
                    

                    !print*,'error      :', maxval(abs(reshape(finDiff,(/4*totalModes/)) - totalCalcDRes(:,i)))
                    print*,'finite diff:', reshape(finDiff,(/4*totalModes/))
                    print*,'calculated :', totalCalcDRes(:,i)
                    !print*,''

                    !print*,'error:      ', finDiff - totalCalcDRes

                    maxerr = max(maxerr,maxval(abs(reshape(finDiff,(/4*totalModes/)) - totalCalcDRes(:,i))))
                    MaxOverall = max(MaxOverall,maxerr)
                end do
                !print*, maxerr,tri_id
            end do
        end do
        !print*,MaxOverall,tri_id
    end do

  end subroutine
  
  subroutine test_both   
    use my_kinddefs
    use initializeSolution_module, only: solCoeffs
    use inputs_module,  only: totalModes
   
    implicit none 

    integer(i4) :: i,ii,j,k,ss,js,tm,m,n,ind,nftm,astat,el,er,ssl,ssr,tml,tmr,jsr,jsl,find
    real(dp) :: error,mem

    complex(dp) :: solComplex(numFields,totalModes,numTri)
    complex(dp) :: resc(numFields*totalModes)
    complex(dp) :: h
    real(dp)    :: jacobian(numFields*totalModes*numFields*totalModes)
    real(dp)    :: jacobiano(numFields*totalModes*numFields*totalModes)
  
    call Residual_Jac(solCoeffs)
   
    print*,'checking block diagonal jacobian'
   
    error = 0.0_dp
    h = (0._qp,1E-31_qp)
    solComplex= solCoeffs
    error = 0.0_dp
    do i=1,numTri    
        print*,'element num:',i
        
        ss = 1
        js = 1
        do j=1,totalModes
            do k=1,numFields
                solComplex(k,j,i) = solComplex(k,j,i) + h
                call ResidualComplex(solComplex)
                solComplex(k,j,i) = solComplex(k,j,i) - h
                
                ind = 1
                
                resc = reshape(spResRC(:,:,i),(/numFields*totalModes/))
                jacobian = reshape(jacobDiag(:,:,i),(/numFields*totalModes*numFields*totalModes/))
                do m=1,totalModes
                    do n=1,numFields
!                        print*,'approx:',imag(resc(ind))/imag(h)
!                        print*,'JacDia:' ,jacobian(js)
                        print*,'error',abs(jacobian(js)-imag(resc(ind))/imag(h))
                        error = max(error,abs(jacobian(js)-imag(resc(ind))/imag(h)))

                        js = js + 1
                        ind = ind + 1
                       
                    end do                               
                end do
            end do
        end do
        print*,error,i
    end do
!          
   
    print*,'checking block off diagonal jacobian'
    error = 0.0_dp
    do i=1,numInterior
        
        find = interiorEdges(i)
        el = edgeList(find)%e1    
        er = edgeList(find)%e2

        ssl = 1     
        jsr = 1

        ! check R/L
        do j=1,totalModes
            do k=1,numFields

                solComplex(k,j,el) = solComplex(k,j,el) + h
                call ResidualComplex(solComplex)
                solComplex(k,j,el) = solComplex(k,j,el) - h

                ind = 1
                resc = reshape(spResRC(:,:,er),(/numFields*totalModes/))

                jacobiano = reshape(jacobOffDiagRL(:,:,i),(/numFields*totalModes*numFields*totalModes/))

                do m=1,totalModes
                    do n=1,numFields
                         !print*,'JacobRL',jacobiano(jsr)
                         !print*,'approx ',imag(resc(ind))/imag(h)
!                         print*,'error', abs(jacobiano(jsr)-imag(resc(ind))/imag(h))
                        error = max(error,abs(jacobiano(jsr)-imag(resc(ind))/imag(h)))

                        jsr = jsr + 1
                        ind = ind + 1
                    end do                               
                end do              

                ssl = ssl + 1

            end do
            print*,'RL error',error,el,er
        end do

        ssr = 1
        jsl = 1   

        ! check L/R
        do j=1,totalModes
            do k=1,numFields

                solComplex(k,j,er) = solComplex(k,j,er) + h
                call ResidualComplex(solComplex)
                solComplex(k,j,er) = solComplex(k,j,er) - h

                ind = 1
                resc = reshape(spResRC(:,:,el),(/numFields*totalModes/))

                jacobiano = reshape(jacobOffDiagLR(:,:,i),(/numFields*totalModes*numFields*totalModes/))

                do m=1,totalModes
                    do n=1,numFields
!                        print*,'error', abs(jacobiano(jsl)-imag(resc(ind))/imag(h))
                        error = max(error, abs(jacobiano(jsl)-imag(resc(ind))/imag(h)))

                        jsl = jsl + 1
                        ind = ind + 1
                    end do                               
                end do              
                ssr = ssr + 1
            end do
        end do
        print*,'LR error',error,el,er
    end do
   
    print*,'total jacobian error',error
   
    end subroutine test_both
    
#endif    
  
  end module 
  
  


 