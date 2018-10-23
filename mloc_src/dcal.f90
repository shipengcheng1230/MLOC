      subroutine dircal (it)
      
      implicit none
      
      include 'mloc.inc'
      
      integer iev, it, i, j, mf, nf, kbayes, ngt(0:11), ncum, igt, ngt5
      real f
      real cv22(2,2), alpha, al, bl, kcrit2, xl, yl
      real pc, avc
      real t11, t12, t22, kcrit1, ddep, dot
      real area
      real sccv(nevmax,4,4), shcv(4,4), dccv(nevmax,4,4)
      character*100 outfil
      
      data pc/0.90/
      
      outfil = trim(outfile)//'.dcal'
      open (io_cal,file=outfil,status='new')
      
      sccv = 0.
      shcv = 0.
      dccv = 0.
      ngt = 0
      ncum = 0
      ngt5 = 0
      
      do iev = 1,nev
      
         ! Scaled covariance matrices
         if (mindx(iev,1) .ne. 0 .and. mindx(iev,2) .ne. 0) then
            mf=2
            kbayes=3
            nf=kbayes+ndatc(iev,it)-mf
            call fstat1 (mf, nf, pc, f)
            kcrit2 = mf*((real(kbayes)+eciev(iev,it))/real(nf))*f 
            cv22(1,1) = ccv(iev,1,1)
            cv22(1,2) = ccv(iev,1,2)
            cv22(2,1) = ccv(iev,1,2)
            cv22(2,2) = ccv(iev,2,2)
            call elips (cv22, alpha, al, bl)
            xl = sqrt(kcrit2/al)
            yl = sqrt(kcrit2/bl)
            avc = xl*yl*pi ! Geometric area of 90% confidence ellipse for cluster vector
            
            ! Convert 90% confidence ellipse back to matrix form under the asumption of kcrit = 1.
            call ell2cv (1.0, alpha, xl, yl, t11, t12, t22)
            sccv(iev,1,1) = t11
            sccv(iev,1,2) = t12
            sccv(iev,2,1) = t12
            sccv(iev,2,2) = t22
         end if
         
         ! Statistics for univariate parameters
         mf = 1
         kbayes = 3
         nf=kbayes+ndatc(iev,it)-mf
         call fstat1 (mf, nf, pc, f)
         kcrit1 = mf*((real(kbayes)+eciev(iev,it))/real(nf))*f 
!         write (io_cal,'(i3,3x,a,f10.3)') iev, 'kcrit1 = ', kcrit1
         
         ! Depth uncertainty
         if (mindx(iev,3) .ne. 0) then
            ddep = sqrt(kcrit1*ccv(iev,3,3))
            sccv(iev,3,3) = ddep*ddep
         end if
         
         ! Origin time uncertainty
         if (mindx(iev,4) .ne. 0) then
            dot = sqrt(kcrit1*ccv(iev,4,4))
            sccv(iev,4,4) = dot*dot
         end if
         
      end do ! End of loop over events
      
      ! Scaled hypocentroid covariance matrix
      call ell2cv (1.0, alphah, xl1h, xl2h, t11, t12, t22)
      shcv(1,1) = t11
      shcv(1,2) = t12
      shcv(2,1) = t12
      shcv(2,2) = t22
      shcv(3,3) = sdxhath(3)*sdxhath(3)
      shcv(4,4) = sdxhath(4)*sdxhath(4)
      
      write (io_cal,'(/a)') 'Cumulative uncertainty of direct-calibrated events'
      write (io_cal,'(a)') 'iev     alpha        xl        yl      area       eqr      ddep       dot'
      
      do iev = 1,nev
         do i = 1,4
            do j = 1,4
               dccv(iev,i,j) = shcv(i,j) + sccv(iev,i,j)
            end do
         end do
!         write (io_cal,'(i3,4f10.3)') iev, dccv(iev,1,1), dccv(iev,1,2), dccv(iev,2,1), dccv(iev,2,2)
         cv22(1,1) = dccv(iev,1,1)
         cv22(1,2) = dccv(iev,1,2)
         cv22(2,1) = dccv(iev,2,1)
         cv22(2,2) = dccv(iev,2,2)
         call elips (cv22, alpha, al, bl)
         xl = sqrt(1./al)
         yl = sqrt(1./bl)
         igt = nint(yl)
         if (igt .ge. 0 .and. igt .le. 10) then
            ngt(igt) = ngt(igt) + 1
         else
            ngt(11) = ngt(11) + 1 ! This holds the number of events worse than CE10
         end if
         area = xl*yl*pi ! Geometric area
         eqr(iev) = sqrt(area/pi)
         ddep = sqrt(dccv(iev,3,3))
         dot = sqrt(dccv(iev,4,4))
         write (io_cal,'(i3,7f10.3)') iev, alpha, xl, yl, area, eqr(iev), ddep, dot
         alphadc(iev) = alpha
         xl1dc(iev) = xl
         xl2dc(iev) = yl
         ddepdc(iev) = ddep
         dotdc(iev) = dot
      end do
      
      ! Summary of calibration levels
      write (*,'(/a)') 'CE    N   Cumulative'
      write (io_cal,'(/a)') 'CE    N   Cumulative'
      do i=0,11
         ncum = ncum + ngt(i)
         if (i .eq. 5) ngt5 = ncum
         write (*,'(i2,2x,i3,7x,i4)') i, ngt(i), ncum
         write (io_cal,'(i2,2x,i3,7x,i4)') i, ngt(i), ncum
      end do
      write (*,'(/a,i3,a,i3,a/)') 'Direct calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
      write (io_cal,'(/a,i3,a,i3,a/)') 'Direct calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
      
      close (io_cal)
      
      return
      end
