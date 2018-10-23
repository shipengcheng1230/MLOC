      subroutine cal_shift (mode, it, otsp)
      
!When calibration locations are available for one of more cluster events,
!calculate the shift needed to bring the cluster into best alignment
!(the "optimal" shift vector). The original method (mode = 0) was just
!take the average of all the individual shift vectors. this is still
!supported, but the preferred way is take into account the
!uncertainties of the calibration data and the cluster vectors, through
!their covariance matrices, and also to consider the possibility of
!bias that is not included in the formal covariance matrices.
!Covariance matrices for calibration locations are input by the user,
!using the command "calX" or they can be calculayed from the epicenter
!confidence ellipse of a GTXX event in an MNF input file.
!
!The calibration event CVs (rcv) are assumed to be already scaled to a
!90% confidence level. This could be done, for example, by applying a
!utility program that converts from the confidence ellipse calculated
!in some other code (e.g., a single event location code) to an
!equivalent covariance matrix. In the case of InSAR analysis or
!geologic information (fault trace) it may be seat-of-the-pants. The
!CVs for cluster vectors (ccv) are first converted to 90% confidence
!ellipses with the normal Bayesian approach, then converted back to raw
!CVs (subroutine ell2cv), after which they can be added to the
!calibration event CVs (because they are completely independent
!estimates), yielding the combined CVs (scv). These are converted (no
!scaling) directly to 90% confidence ellipses for the combination. The
!area of these combined confidence ellipses is used for inverse
!weighting of the shift vectors estimated for each calibration event, to
!get the estimate of the optimal shift vector for the entire cluster.
!
!The uncertainty of the optimal shift vector (gcv) is based on the
!weighted (same weights as for the individual shift vectors)
!combination of the combined CVs (scv), expressed as a confidence
!ellipse (90%). This is only part of the uncertainty, however. There is
!additional uncertainty, which can be thought of as associated with the
!"unmodelled" or "biased" part of the indirect calibration problem; it
!arises if any of the calibration event locations are biased. Such bias
!could arise with local network solutions if the velocity model is poorly
!chosen, or if readings from too great a distance have been used, or if
!some phases are mis-identified or outliers are not properly weighted. It
!can also arise if a serious error (e.g., an outlier reading, with poor
!azimuthal coverage) has caused bias in the estimate of a cluster
!vector in the HD analysis. Such bias is seen by comparing the
!individual estimates of shift vector to the optimal shift vector. If
!there are departures that are statistically unlikely, given the
!alleged confidence of the calibration locations and cluster vectors,
!then there is need to account for the unmodeled error.
!Graphically, it is easily seen by plotting residual shift vectors
!(subtract the optimal shift vector) and comparing to the confidence
!ellipses of each shift vector (scv).

!Instead of taking the RMS of the residual shift vectors to estimate the
!level of inconsistency (this is too sensitive to outliers and overestimates
!the uncertainty), I do a test based on the coverage statistic for the
!residual calibration shift vectors, for which the cumulative binomial
!probability distribution gives the probability of the observed number of
!"uncovered" vectors. See 'rdbt_test2' for more explanation. If the null
!hypothesis is rejected, the test is repeated after adding a small amount
!to each covariance matrix, equivalent to adding a circular uncertainty.
!When the null hypothesis cannot be rejected, this extra covariance is
!converted to a "radius of doubt", which is added to the modelled
!uncertainty (gcv) to yield an "augmented GT covariance matrix"
!(agcv). The area of the equivalent ellipse is converted to a circular
!area, and the radius of that circle is a useful measure of the calibration
!level of the cluster. The calibration levels of individual events (accv) are
!found by adding their cluster vector CVs to agcv.

! I brought back the older algorithm (rdbt_test1) for radius of doubt (v8.1) that is based
! on a test of the hypothesis that all residual shift vectors have zero length.
! I have kept the test based on coverage statistics for now, but it has the problem
! that we seldom have more than a couple calibration events and the 90% coverage
! requirement translates into 100% coverage requirement if there are fewer than 10
! calibration events. If there are 10 or more calibration events then the larger of 
! the two estimates is used. Otherwise the test based on zero-length residual vectors
! is used.
      
      implicit none
      
      include 'mloc.inc'
      include 'cal_shift.inc'
      
      integer mode, iev, it, i, j, mf, nf, kbayes, nr, ngt(0:11), ncum, igt, ngt5, nrc
      real delsumlat, delsumlon, delsumdep, delsumtim, otsp(nevmax,0:itmax1), hms2s, f, lon_test
      real cv22(2,2), alpha, al, bl, kcrit2, xl, yl, accv_1(4,4), accv_2(4,4), yl_1, yl_2
      real pc, sumw1, sumw2, sumw3, sumw4, avc, avr, gcv(4,4), agcv(4,4)
      real t11, t12, t22, kcrit1, ddep, dot
      real scale12, scale3, scale4
      real sumrvl2, rvlrms, rvl, rvlmax, dx, dy, dgkmlo
      real sumdd, sumdt, dtiev, ddiev, ddrms, dtrms, area
      real sccv(nevmax,4,4)
      real coverage_percent, covp, rdbt1, rdbt2
      character*100 outfil
      real depsn, otsn, deprsv(nevmax), otrsv(nevmax)
      character*12 caltype_pr_iev
      character(132)::msg
      
      data pc/0.90/ ! Confidence level
      
      ngt = 0
      ncum = 0
      ngt5 = 0
      
      outfil = trim(outfile)//'.cal'
      open (io_cal,file=outfil,status='new')
      
      write (io_cal,'(a,i4,a)') ' On ', ncal(3), ' calibration events: '
      write (*,'(a,i4,a)') ' On ', ncal(3), ' calibration events: '
      
      delsumlat = 0.
      delsumlon = 0.
      delsumdep = 0.
      delsumtim = 0.
      
      if (mode .eq. 0) then ! Average of individual shift vectors, no weighting
         do iev = 1,nev
            if (cal_event(iev,3)) then ! Difference between calibration and inversion solution
               if (mindx(iev,1) .ne. 0 .and. mindx(iev,2) .ne. 0) then
                  w12(iev) = 1.
               else
                  w12(iev) = 0.
               end if
               if (mindx(iev,3) .ne. 0) then
                  w3(iev) = 1.
               else
                  w3(iev) = 0.
               end if
               if (mindx(iev,4) .ne. 0) then
                  w4(iev) = 1.
               else
                  w4(iev) = 0.
               end if
               delsumlat = delsumlat + (cal_lat(iev,3) - latp(iev,it+1))*w12(iev)
               lon_test = lonp(iev,it+1)
               call set_longitude_range (lon_test, longitude_range)
               delsumlon = delsumlon + (cal_lon(iev,3) - lon_test)*w12(iev)
               delsumdep = delsumdep + (cal_dep(iev,3) - depthp(iev,it+1))*w3(iev)
               delsumtim = delsumtim + (hms2s(cal_hr(iev,3), cal_min(iev,3), cal_sec(iev,3)) - otsp(iev,it+1))*w4(iev)
            end if
         end do
         del_cal_lat = delsumlat/ncal(3) ! Average correction for latitude
         del_cal_lon = delsumlon/ncal(3) ! Average correction for longitude
         del_cal_dep = delsumdep/ncal(3) ! Average correction for depth
         del_cal_tim = delsumtim/ncal(3) ! Average correction for origin time
         write (io_cal,'(a,4f10.3)') ' Average (REF-SOL): ', del_cal_lat, del_cal_lon, del_cal_dep, del_cal_tim
         
      else if (mode .eq. 1) then ! Weighted mean, epicenters weighted by inverse of area of confidence ellipse
         sumw1 = 0.
         sumw2 = 0.
         sumw3 = 0.
         sumw4 = 0.
         w12 = 0.
         w3 = 0.
         w4 = 0.
         sccv = 0.
         scv = 0.
         sgcv = 0.
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
            ! write (io_cal,'(3x,a,f10.3)') 'kcrit1 = ', kcrit1
            
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
            
            ! Calibration events
            if (cal_event(iev,3)) then
               write (io_cal,'(/a,i3)') 'Event ', iev
               
               ! Area of 90% confidence ellipse for cluster vector + calibration location and inverse weight factor
               if (mindx(iev,1) .ne. 0 .and. mindx(iev,2) .ne. 0) then
                  write (io_cal,'(a,i3)') 'Cluster vector CV for event ', iev
                  write (io_cal,'(3x,a,f8.3)') 'kcrit2: ', kcrit2
                  write (io_cal,'(3x,a,3f8.3)') '90% confidence ellipse: ', alpha, xl, yl
                  write (io_cal,'(3x,a,f8.3,a)') 'Area: ', avc, ' km**2'
                  
                  ! Calibration location confidence ellipse and area
                  cv22(1,1) = rcv(iev,3,1,1)
                  cv22(1,2) = rcv(iev,3,1,2)
                  cv22(2,1) = rcv(iev,3,1,2)
                  cv22(2,2) = rcv(iev,3,2,2)
                  call elips (cv22, alpha, al, bl)
                  xl = sqrt(1./al)
                  yl = sqrt(1./bl)
                  avr = xl*yl*pi ! Geometric area of 90% confidence ellipse for calibration location
                  write (io_cal,'(a,i3)') 'Calibration location CV for event ', iev
                  write (io_cal,'(3x,a,3f8.3)') '90% confidence ellipse: ', alpha, xl, yl
                  write (io_cal,'(3x,a,f8.3,a)') 'Area: ', avr, ' km**2'
                  w12(iev) = 10./(avc + avr) ! Factor of ten makes the numbers closer to 1
               else
                  w12(iev) = 1.
               end if
               
               ! Depth uncertainty
               if (mindx(iev,3) .ne. 0) then
                  w3(iev) = 1./(sccv(iev,3,3)+rcv(iev,3,3,3))
               else
                  w3(iev) = 1.
               end if
               
               ! Origin time uncertainty
               if (mindx(iev,4) .ne. 0) then
                  w4(iev) = 1./(sccv(iev,4,4)+rcv(iev,3,4,4))
               else
                  w4(iev) = 1.
               end if
               
               write (io_cal,'(a,i3)') 'Weights for event ', iev
               write (io_cal,'(3x,a,f8.3)') 'w12 = ', w12(iev)
               write (io_cal,'(3x,a,f8.3)') 'w3 =  ', w3(iev)
               write (io_cal,'(3x,a,f8.3)') 'w4 =  ', w4(iev)
               
               ! Add the re-scaled cluster vector covariance matrix and that of the calibration parameters
               scv(iev,1,1) = sccv(iev,1,1) + rcv(iev,3,1,1)
               scv(iev,1,2) = sccv(iev,1,2) + rcv(iev,3,1,2)
               scv(iev,2,1) = sccv(iev,2,1) + rcv(iev,3,2,1)
               scv(iev,2,2) = sccv(iev,2,2) + rcv(iev,3,2,2)
               scv(iev,3,3) = sccv(iev,3,3) + rcv(iev,3,3,3)
               scv(iev,4,4) = sccv(iev,4,4) + rcv(iev,3,4,4)
               
               write (io_cal,'(a,i3)') 'Covariance matrices for calibration event ', iev
               write (io_cal,'(3x,a,6f10.3)') ' ccv: ', ccv(iev,1,1), ccv(iev,1,2), ccv(iev,2,1), ccv(iev,2,2),&
                ccv(iev,3,3), ccv(iev,4,4)
               write (io_cal,'(3x,a,6f10.3)') 'sccv: ', sccv(iev,1,1), sccv(iev,1,2), sccv(iev,2,1), sccv(iev,2,2),&
                sccv(iev,3,3), sccv(iev,4,4)
               write (io_cal,'(3x,a,6f10.3)') ' rcv: ', rcv(iev,3,1,1), rcv(iev,3,1,2), rcv(iev,3,2,1), rcv(iev,3,2,2),&
                rcv(iev,3,3,3), rcv(iev,3,4,4)
               write (io_cal,'(3x,a,6f10.3)') ' scv: ', scv(iev,1,1), scv(iev,1,2), scv(iev,2,1), scv(iev,2,2),&
                scv(iev,3,3), scv(iev,4,4)
               
               ! Weighted sum of shift vector components
               delsumlat = delsumlat + (cal_lat(iev,3) - latp(iev,it+1))*w12(iev)
               sumw1 = sumw1 + w12(iev)
               lon_test = lonp(iev,it+1)
               call set_longitude_range (lon_test, longitude_range)
               delsumlon = delsumlon + (cal_lon(iev,3) - lon_test)*w12(iev)
               sumw2 = sumw2 + w12(iev)
               delsumdep = delsumdep + (cal_dep(iev,3) - depthp(iev,it+1))*w3(iev)
               sumw3 = sumw3 + w3(iev)
               delsumtim = delsumtim + (hms2s(cal_hr(iev,3), cal_min(iev,3), cal_sec(iev,3)) - otsp(iev,it+1))*w4(iev)
               sumw4 = sumw4 + w4(iev)
            end if
         end do ! End of loop over events
         
         ! Weighted mean of calibration shifts for all calibration events
         del_cal_lat = delsumlat/sumw1 ! Weighted mean correction for latitude
         del_cal_lon = delsumlon/sumw2 ! Weighted mean correction for longitude
         del_cal_dep = delsumdep/sumw3 ! Weighted mean correction for depth
         del_cal_tim = delsumtim/sumw4 ! Weighted mean correction for origin time
         write (io_cal,'(/a,4f10.3)') 'Weighted mean (REF-SOL) calibration shift: '
         write (io_cal,'(3x,a,f8.3,a)') 'Latitude : ', del_cal_lat, ' deg'
         write (io_cal,'(3x,a,f8.3,a)') 'Longitude: ', del_cal_lon, ' deg'
         write (io_cal,'(3x,a,f8.3,a)') 'Depth    : ', del_cal_dep, ' km'
         write (io_cal,'(3x,a,f8.3,a)') 'OT       : ', del_cal_tim, ' sec'
         ! Screen display
         write (*,'(/a,4f10.3)') 'Weighted mean (REF-SOL) calibration shift: '
         write (*,'(3x,a,f8.3,a)') 'Latitude : ', del_cal_lat, ' deg'
         write (*,'(3x,a,f8.3,a)') 'Longitude: ', del_cal_lon, ' deg'
         write (*,'(3x,a,f8.3,a)') 'Depth    : ', del_cal_dep, ' km'
         write (*,'(3x,a,f8.3,a)') 'OT       : ', del_cal_tim, ' sec'
         
         ! Renormalize the weights for printed output (so they sum to ncal(3))
         scale12 = real(ncal(3))/sumw1
         scale3 = real(ncal(3))/sumw3
         scale4 = real(ncal(3))/sumw4
         do iev = 1,nev
            if (cal_event(iev,3)) then
               w12(iev) = w12(iev)*scale12
               w3(iev) = w3(iev)*scale3
               w4(iev) = w4(iev)*scale4
            end if
         end do
         
         ! Initialize covariance matrices for uncertainty of the calibration shift
         gcv = 0.
         agcv = 0.
         
         ! Weighted average covariance matrix for all calibration events
         do iev = 1,nev
            if (cal_event(iev,3)) then
               gcv(1,1) = gcv(1,1) + scv(iev,1,1)*w12(iev)
               gcv(1,2) = gcv(1,2) + scv(iev,1,2)*w12(iev)
               gcv(2,1) = gcv(2,1) + scv(iev,2,1)*w12(iev)
               gcv(2,2) = gcv(2,2) + scv(iev,2,2)*w12(iev)
               gcv(3,3) = gcv(3,3) + scv(iev,3,3)*w3(iev)
               gcv(4,4) = gcv(4,4) + scv(iev,4,4)*w4(iev)
            end if
         end do
         do i = 1,4
            do j = 1,4
               gcv(i,j) = gcv(i,j)/real(ncal(3))
            end do
         end do
         write (io_cal,'(/a)') 'Weighted average CV for all calibration events'
         write (io_cal,'(a,6f10.3)') 'gcv: ', gcv(1,1), gcv(1,2), gcv(2,1), gcv(2,2), gcv(3,3), gcv(4,4)
         
         ! RMS length of residual vectors (cluster vector - calibration shift).
         ! Similar for depth and OT.
         ! Used as estimates of unmodelled bias in the estimation of GT shift vectors.
         nr = 0 ! Number of calibration events
         nrc = 0 ! Number of calibration events with coverage after the calibration shift
         sumrvl2 = 0
         sumdd = 0
         sumdt = 0
         write (io_cal,'(/a)') 'Residual calibration shift vectors and COVP, based on SGCV'
         write (io_cal,'(a,5x,a,5x,a,7x,a,8x,a,8x,a,6x,a)') 'iev', 'dtiev', 'ddiev', 'rvl', 'dx', 'dy', 'covp'
         rvlmax = 0.
         do iev = 1,nev
            if (cal_event(iev,3)) then
               nr = nr + 1
               dy = (del_cal_lat - (cal_lat(iev,3)-latp(iev,it+1)))/dgkmla
               lon_test = lonp(iev,it+1)
               call set_longitude_range (lon_test, longitude_range)
               dx = (del_cal_lon - (cal_lon(iev,3)-lon_test))/dgkmlo(latp(iev,it+1))
               rsvx(iev) = dx
               rsvy(iev) = dy
               rvl = sqrt(dy*dy+dx*dx)
               rvlmax = amax1(rvl,rvlmax)
               ! sgcv is used for the coverage tests of calibration events. It includes the uncertainty of each
               ! individual calibration shift vector (scv) and the uncertainty of the calibration shift (gcv)
               sgcv(iev,1,1) = scv(iev,1,1) + gcv(1,1)
               sgcv(iev,1,2) = scv(iev,1,2) + gcv(1,2)
               sgcv(iev,2,1) = scv(iev,2,1) + gcv(2,1)
               sgcv(iev,2,2) = scv(iev,2,2) + gcv(2,2)
               call coverage_ellipse (dx, dy, sgcv(iev,1,1), sgcv(iev,1,2), sgcv(iev,2,2), covp)
               if (covp .le. 1.0) nrc = nrc + 1
               sumrvl2 = sumrvl2 + rvl*rvl
               ddiev = del_cal_dep - (cal_dep(iev,3) - depthp(iev,it+1))
               sumdd = sumdd + ddiev*ddiev
               deprsv(nr) = ddiev ! Input for subroutine croux, below
               dtiev = del_cal_tim - (hms2s(cal_hr(iev,3), cal_min(iev,3), cal_sec(iev,3)) - otsp(iev,it+1))
               sumdt = sumdt + dtiev*dtiev
               otrsv(nr) = dtiev ! Input for subroutine croux, below
               write (io_cal,'(i3,6f10.3)') iev, otrsv(nr), deprsv(nr), rvl, dx, dy, covp
            else
               rsvx(iev) = 0.
               rsvy(iev) = 0.
            end if
         end do
         if (nr .gt. 0) then
            coverage_percent = real(nrc)/real(nr)
         else
            coverage_percent = 0.
         end if
         write (io_cal,'(a,i3,a,i3,a,f6.3)') 'Coverage: ', nrc, ' / ', nr, ' = ', coverage_percent
         rvlrms = sqrt(sumrvl2/real(nr))
         ddrms = sqrt(sumdd/real(nr))
         dtrms = sqrt(sumdt/real(nr))
         
         ! Radius of doubt (inconsistency between shift vectors, unmodelled errors in calibration events and HD)
         ! If there are enough calibration events (10 or more) the coverage test can be considered and the larger
         ! of the two estimates is used. Otherwise the radius of doubt from the test of zero-length residual
         ! shift vectors is used.
         if (nr .ge. 2) then
            call rdbt_test1 (rvlmax, pc) ! Test of residual shift vectors having zero length
            rdbt1 = rdbt
            call rdbt_test2 (rvlmax, pc) ! Test based on coverage statistics
            rdbt2 = rdbt
            write (*,'(/a,f4.1,a)') 'Radius of doubt = ', rdbt1,&
            ' km, based on test of zero-length residual calibration shift vectors'
            write (*,'(a,f4.1,a)') 'Radius of doubt = ', rdbt2, ' km, based on test of coverage statistics'
            if (nr .ge. 10) then
               rdbt = amax1(rdbt1,rdbt2)
            else
               rdbt = rdbt1
            end if
            write (*,'(a,f4.1)') 'Radius of doubt = ', rdbt
         else
            rdbt = 0. ! Radius of doubt = 0 if there is only one calibration event.
         end if
         
         ! Use spread of reduced shift vectors for depth and OT for equivalent "radius of doubt" for OT and depth
         if (nr .ge. 2) then
            call croux (deprsv, nr, depsn)
            call croux (otrsv, nr, otsn)
         else
            depsn = 0.
            otsn = 0.
         end if
         
         write (io_cal,'(/a)') 'Inconsistency (between multiple calibration data) terms for augmented GT covariance: '
         write (io_cal,'(3x,a,2(f8.3,a))') 'Epicenter: ', rdbt,  ' km (rms = ', rvlrms, ' km)'
         write (io_cal,'(3x,a,2(f8.3,a))') 'Depth    : ', depsn, ' km (rms = ', ddrms, ' km)'
         write (io_cal,'(3x,a,2(f8.3,a))') 'OT       : ', otsn, ' sec (rms = ', dtrms, ' sec)'
         
         ! Augmented GT covariance matrix
         agcv(1,1) = gcv(1,1) + rdbt*rdbt
         agcv(1,2) = gcv(1,2)
         agcv(2,1) = gcv(2,1)
         agcv(2,2) = gcv(2,2) + rdbt*rdbt
         agcv(3,3) = gcv(3,3) + depsn*depsn
         agcv(4,4) = gcv(4,4) + otsn*otsn
         write (io_cal,'(/a)') 'Augmented GT covariance matrix'
         write (io_cal,'(a,6f10.3)') 'agcv: ', agcv(1,1), agcv(1,2), agcv(2,1), agcv(2,2), agcv(3,3), agcv(4,4)
                  
         ! Calibration shift uncertainties from augmented CV. Convert area to circular calibration level for cluster.
         cv22(1,1) = agcv(1,1)
         cv22(1,2) = agcv(1,2)
         cv22(2,1) = agcv(2,1)
         cv22(2,2) = agcv(2,2)
         call elips (cv22, alpha, al, bl)
         alphagt = alpha
         s1gt = sqrt(1./al)
         s2gt = sqrt(1./bl)
         area = s1gt*s2gt*pi ! Geometric area
         gtlevel = sqrt(area/pi)
         depgt = sqrt(agcv(3,3))
         otgt = sqrt(agcv(4,4))
         write (io_cal,'(/a)') 'Uncertainties from the augmented GT covariance matrix:'
         write (io_cal,'(3x,a,3f10.3)') 'Confidence ellipse: ', alphagt, s1gt, s2gt
         write (io_cal,'(3x,a,f8.3,a)') 'Area of ellipse: ', area, ' km**2'
         write (io_cal,'(3x,a,f8.3,a)') 'Equivalent circular radius (CE level) = ', gtlevel, ' km'
         write (io_cal,'(3x,a,f10.3)') 'Depth shift uncertainty: ', depgt
         write (io_cal,'(3x,a,f10.3)') 'OT shift uncertainty: ', otgt
         
         ! Calibration level for the hypocentroid
         cal_level = s2gt
                    
         ! Augmented CV for all cluster events after calibration shift is applied. Normally the augmented GT shift CV is added
         ! to the event's cluster vector CV. If an event has a calibration location, however, there are three options,
         ! controlled by the variable icaltype:
         ! 1 = Traditional: calibration events get the uncertainties of the calibration data
         ! 2 = Systematic: all events done same way, add calibration uncertainty to cluster vector uncertainty.
         ! 3 = Optimal: use traditional or systematic method, whichever has shorter semi-major axis.
         write (io_cal,'(/a)') 'Cumulative uncertainty of calibration-shifted cluster events'
         write (io_cal,'(a)') 'iev     alpha        xl        yl      area       eqr      ddep       dot      covp     caltype'
         ngt5 = 0
         nrc = 0
         do iev = 1,nev
            if (cal_event(iev,3)) then
               if (icaltype .eq. 1) then
                  caltype_pr_iev = 'traditional'
                  accv = 0.
                  accv(iev,1,1) = rcv(iev,3,1,1) + rdbt*rdbt
                  accv(iev,1,2) = rcv(iev,3,1,2)
                  accv(iev,2,1) = rcv(iev,3,1,2)
                  accv(iev,2,2) = rcv(iev,3,2,2) + rdbt*rdbt
                  accv(iev,3,3) = rcv(iev,3,3,3) + depsn*depsn
                  accv(iev,4,4) = rcv(iev,3,4,4) + otsn*otsn
               else if (icaltype .eq. 2) then
                  caltype_pr_iev = 'systematic'
                  do i = 1,4
                     do j = 1,4
                        accv(iev,i,j) = agcv(i,j) + sccv(iev,i,j)
                     end do
                  end do
               else if (icaltype .eq. 3) then
                  accv_1 = 0.
                  accv_1(1,1) = rcv(iev,3,1,1) + rdbt*rdbt
                  accv_1(1,2) = rcv(iev,3,1,2)
                  accv_1(2,1) = rcv(iev,3,1,2)
                  accv_1(2,2) = rcv(iev,3,2,2) + rdbt*rdbt
                  accv_1(3,3) = rcv(iev,3,3,3) + depsn*depsn
                  accv_1(4,4) = rcv(iev,3,4,4) + otsn*otsn
                  cv22(1,1) = accv_1(1,1)
                  cv22(1,2) = accv_1(1,2)
                  cv22(2,1) = accv_1(2,1)
                  cv22(2,2) = accv_1(2,2)
                  call elips (cv22, alpha, al, bl)
                  yl_1 = sqrt(1./bl)
                  do i = 1,4
                     do j = 1,4
                        accv_2(i,j) = agcv(i,j) + sccv(iev,i,j)
                     end do
                  end do
                  cv22(1,1) = accv_2(1,1)
                  cv22(1,2) = accv_2(1,2)
                  cv22(2,1) = accv_2(2,1)
                  cv22(2,2) = accv_2(2,2)
                  call elips (cv22, alpha, al, bl)
                  yl_2 = sqrt(1./bl)
                  accv = 0.
                  if (yl_1 .le. yl_2) then
                     caltype_pr_iev = 'traditional'
                     accv(iev,1,1) = accv_1(1,1)
                     accv(iev,1,2) = accv_1(1,2)
                     accv(iev,2,1) = accv_1(2,1)
                     accv(iev,2,2) = accv_1(2,2)
                     accv(iev,3,3) = accv_1(3,3)
                     accv(iev,4,4) = accv_1(4,4)
                  else
                     caltype_pr_iev = 'systematic'
                     accv(iev,1,1) = accv_2(1,1)
                     accv(iev,1,2) = accv_2(1,2)
                     accv(iev,2,1) = accv_2(2,1)
                     accv(iev,2,2) = accv_2(2,2)
                     accv(iev,3,3) = accv_2(3,3)
                     accv(iev,4,4) = accv_2(4,4)
                  end if
               end if
            else ! Not a calibration event
               caltype_pr_iev = '            '
               do i = 1,4
                  do j = 1,4
                     accv(iev,i,j) = agcv(i,j) + sccv(iev,i,j)
                  end do
               end do
            end if
            !write (io_cal,'(i3,4f10.3)') iev, accv(iev,1,1), accv(iev,1,2), accv(iev,2,1), accv(iev,2,2)
            cv22(1,1) = accv(iev,1,1)
            cv22(1,2) = accv(iev,1,2)
            cv22(2,1) = accv(iev,2,1)
            cv22(2,2) = accv(iev,2,2)
            call elips (cv22, alpha, al, bl)
            xl = sqrt(1./al)
            yl = sqrt(1./bl)
            igt = nint(yl)
            if (igt .ge. 0 .and. igt .le. 10) then
               ngt(igt) = ngt(igt) + 1
            else
               ngt(11) = ngt(11) + 1 ! This holds the number of events worse than CE10
            end if
            area = xl*yl*pi ! geometric area
            eqr(iev) = sqrt(area/pi) ! Radius of equivalent circle
            ddep = sqrt(accv(iev,3,3))
            dot = sqrt(accv(iev,4,4))
            covp = 0.
            call coverage_ellipse (rsvx(iev), rsvy(iev), accv(iev,1,1), accv(iev,1,2), accv(iev,2,2), covp)
            if (cal_event(iev,3) .and. covp .le. 1.0) nrc = nrc + 1
            write (io_cal,'(i3,8f10.3,2x,a)') iev, alpha, xl, yl, area, eqr(iev), ddep, dot, covp, trim(caltype_pr_iev)
            alphacg(iev) = alpha
            xl1cg(iev) = xl
            xl2cg(iev) = yl
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
         write (*,'(/a,i3,a,i3,a/)') 'Indirect calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
         write (io_cal,'(/a,i3,a,i3,a/)') 'Indirect calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
         
         ! Coverage percentage for calibration events
         if (nr .gt. 0) then
            coverage_percent = real(nrc)/real(nr)
            write (*,'(a,i3,a,i3,a,f6.3/)') 'Coverage for calibration events: ', nrc, ' / ', nr, ' = ', coverage_percent
            write (io_cal,'(a,i3,a,i3,a,f6.3/)') 'Coverage for calibration events: ', nrc, ' / ', nr, ' = ', coverage_percent
         end if
         
      else
         write (msg,'(a,i4)') 'cal_shift: bad mode - ', mode
         call oops (trim(msg))
      end if
      
      close (io_cal)
      
      return
      end


!-----------------------------------------------------------------------      
      subroutine ell2cv (kcrit, alpha, sminor, smajor, t11, t12, t22)
      
      ! Given a confidence ellipse and a critical value, returns the 
      ! equivalent 2x2 covariance matrix
      
      ! Input
      ! kcrit   critical value
      ! alpha   strike of semi-minor axis (positive clockwise from north)
      ! sminor  semi-minor axis length
      ! smajor  semi-major axis length
      
      ! Output
      ! t11     CV(1,1)
      ! t12     CV(1,2) and CV(2,1)
      ! t22     CV(2,2)
      
      implicit none
      
      real rpd
      parameter (rpd=1.7453293e-2)
      
      real lambda1, lambda2, fi, cf2, sf2, csf
      real alpha, sminor, smajor, t11, t22, t12, kcrit
      
      lambda1 = (sminor*sminor)/kcrit
      lambda2 = (smajor*smajor)/kcrit
      
      fi = rpd*alpha
      
      cf2 = cos(fi)*cos(fi)
      sf2 = sin(fi)*sin(fi)
      csf = cos(fi)*sin(fi)
      
      t11 = lambda1*cf2 + lambda2*sf2
      t22 = lambda1*sf2 + lambda2*cf2
      t12 = (lambda1-lambda2)*csf

      return
      end


!-----------------------------------------------------------------------
      subroutine rdbt_test1 (rvlmax, pc)
      
! Test of the null hypothesis that all residual vectors are zero length.
! If this can be rejected, we want to know how much the covariance matrices
! would need to be inflated to be able to reject the null hypothesis.
! Redo the test after adding a circular component, ie. on the diagonal elements,
! of uncertainty (incr) to each covariance matrix. This is repeated
! until the null hypothesis cannot be rejected or a maximum limit is reached. The
! radius of doubt is then the accumulated increase to the covariance matrices.

! This ought to be combined with a recalulation of the optimal shift, to get updated
! residual shift vectors before each test.

      implicit none
      
      include 'mloc.inc'
      include 'cal_shift.inc'
      
      logical nullhyp
      integer iev, nr, i, ii, nf, mf
      real incr, rvlmax, pc, rdbtsq, kcrit90, f
      character(132)::msg
      
      data incr /0.2/
      
      write (io_cal,'(/a)') 'Radius of doubt test based on null hypothesis that all residual cluster vectors have zero length'
      
      rdbt = 0.
      
      ! Fill arrays for calibration events
      covs = 0.
      cvs = 0.
      nr = 0
      do iev = 1,nev
        if (cal_event(iev,3)) then
           nr = nr + 1
!           write (io_cal,'(i3,2f8.2)') iev, rsvx(iev), rsvy(iev)
           cvs(nr,1) = rsvy(iev) ! dy = Latitude
           cvs(nr,2) = rsvx(iev) ! dx = Longitude
           ii = nr*2 - 1
           covs(ii,ii) = scv(iev,1,1)
           covs(ii,ii+1) = scv(iev,1,2)
           covs(ii+1,ii) = scv(iev,2,1)
           covs(ii+1,ii+1) = scv(iev,2,2)
        end if
      end do
      
      ! Critical value
      write (msg,'(a,i4)') 'rdbt_test1: nr  = ', nr
      if (debug) call debugger (trim(msg))      
      write (msg,'(a,i6)') 'rdbt_test1: ntc = ', ntc
      if (debug) call debugger (trim(msg))      
      write (msg,'(a,i5)') 'rdbt_test1: nqc = ', nqc
      if (debug) call debugger (trim(msg))      
      write (msg,'(a,f10.3)') 'rdbt_test1: shatsqc = ', shatsqc
      if (debug) call debugger (trim(msg))      
      mf = (nr - 1)*2
      nf = ntc - nqc - mf
      write (msg,'(a,i4)') 'rdbt_test1: mf = ', mf
      if (debug) call debugger (trim(msg))      
      write (msg,'(a,i6)') 'rdbt_test1: nf = ', nf
      if (debug) call debugger (trim(msg))      
      call fstat1 (mf, nf, pc, f)
      write (msg,'(a,f10.3)') 'rdbt_test1: f = ', f
      if (debug) call debugger (trim(msg))      
      kcrit90 = real(mf)*sngl(shatsqc)*f
!      kcrit90 = 1. ! I thought this was correct, since the CVs have already been scaled, but it can't be.
                    ! It does not take into account the number of calibration events. With a cluster like
                    ! Pahute with 52 calibration events, the difference is clear. kcrit90 has to scale with
                    ! number of calibration events (mf) or else you get ridiculously large rdbt values.
      write (msg,'(a,f8.1,a)') 'rdbt_test1: critical value = ', kcrit90, ' at 90%'
      if (debug) call debugger (trim(msg))      
      write (io_cal,'(a,f8.1,a)') 'Critical value = ', kcrit90, ' at 90%'
      
      call cvzero (nr, nullhyp, kcrit90)
      
      do while (.not.nullhyp)
         rdbt = rdbt + incr
         rdbtsq = rdbt*rdbt
         do i = 1,nr
            ii = i*2 - 1
            covs(ii,ii) = covs(ii,ii) + rdbtsq
            covs(ii+1,ii+1) = covs(ii+1,ii+1) + rdbtsq
         end do
         call cvzero (nr, nullhyp, kcrit90)
         if (rdbt .gt. rvlmax) then
            write (msg,'(a,f4.1,a)') 'rdbt_test1: rdbt_test exceeds rvlmax (', rvlmax,')'
            if (debug) call debugger (trim(msg))      
            write (io_cal,'(a,f4.1,a)') 'rdbt_test exceeds rvlmax (', rvlmax,')'
            exit
         end if
      end do
      
      write (msg,'(a,f6.2)') 'rdbt_test1: rdbt1 = ', rdbt
      if (debug) call debugger (trim(msg))
      write (io_cal,'(a,f6.2)') 'rdbt1 = ', rdbt
      
      return
      end
      

!************************************************************************
      subroutine cvzero (nr, nullhyp, kcrit90)
      
      ! Test whether all reduced shift vectors can have zero length.
      ! See Jordan and Sverdrup [1981].
      
      implicit none
      save
      
      include 'mloc.inc'
      include 'cal_shift.inc'
      
      double precision ac(nmax2,nmax2), unn(nmax2,nmax2), vnn(nmax2,nmax2), qn(nmax2)
      double precision wnn(nmax2,nmax2), tempay(nmax2,nmax2), y(nmax2,nmax2) 
      real ys(nmax2,nmax2), dxos(nmax2), temp(nmax2)
      integer irank(nmax2), indx(nmax2), nr
      real kocs2, kcrit90
      integer i, n, j, k, iev
!      integer ii, jj 
      real sumx, sumy, avex, avey
      character*40 nulltext
      logical nullhyp
      character(132)::msg
            
      nullhyp = .true.
      
      ac = 0.d0
      unn = 0.d0
      vnn = 0.d0
      qn = 0.d0
      wnn = 0.d0
      tempay = 0.d0
      y = 0.d0
      ys = 0.
      dxos = 0.
      temp = 0.
      
      sumy = 0.
      sumx = 0. 
      do i = 1,nr
         sumy = sumy + cvs(i,1)
         sumx = sumx + cvs(i,2)
      end do
      avey = sumy/(nr)
      avex = sumx/(nr)
      write (msg,'(a,f10.3)') 'cvzero: avey = ', avey
      if (debug) call debugger (trim(msg)) 
      write (msg,'(a,f10.3)') 'cvzero: avex = ', avex
      if (debug) call debugger (trim(msg)) 
      
      n = nr*2
      
      ! Put cluster vectors into linear array
      
      j = 0
      do i = 1,nr
         j = j + 1
         dxos(j*2-1) = cvs(i,1) - avey
         dxos(j*2) = cvs(i,2) - avex
         write (msg,'(a,i2,2f10.3)') 'cvzero: dxos ', j, dxos(j*2-1), dxos(j*2)
         if (debug) call debugger (trim(msg)) 
      end do 
      
      ! Inverse of covariance matrix by singular value decomposition.
      ! The covariance matrix always has nullity = 2 , because
      ! the relative locations are given relative to a fixed epicentroid.
      ! The generalized inverse is formed by zeroing the contribution of
      ! the two corresponding eigenfunctions.
      
!      if (debug) print *, 'ac:'
!      do i = 1,n
!         ii = i+1*2-2
!         do j = 1,n
!            jj = j+1*2-2
!            ac(i,j) = dble(covs(ii,jj))
!         end do
!         if (debug) print *, (ac(i,j),j=1,n)
!      end do
      ac = dble(covs)
      call dsvd (ac, n, n, nmax2, nmax2, unn, vnn, qn)
      if (debug) then
         write (io_log,'(a)') 'qn:'
         write (io_log,'(10e10.3)') (qn(j),j=1,n)
      end if
      call dindexx (n, qn, indx) 
      if (debug) then
         write (io_log,'(a)') 'indx:'
         write (io_log,'(10e10.3)') (indx(j),j=1,n)
      end if
      call rank (n, indx, irank)
      if (debug) then
         write (io_log,'(a)') 'irank:'
         write (io_log,'(10e10.3)') (irank(j),j=1,n)
      end if
      
      do j = 1,n
         do i = 1,n
            wnn(i,j) = 0.0d0
         end do
         if (irank(j) .gt. 2) wnn(j,j) = 1.0d0/qn(j)
      end do
      if (debug) then
         write (io_log,'(a)') 'wnn:' 
         write (io_log,'(10e10.3)') (wnn(j,j),j=1,n)
      end if
      call ddot3 (wnn, unn, tempay, nmax2, nmax2, nmax2, nmax2, nmax2, nmax2, n, n, n)
      call ddot1 (vnn, tempay, y, nmax2, nmax2, nmax2, nmax2, nmax2, nmax2, n, n, n)
      do j = 1,n
         do i = 1,n
            ys(i,j) = sngl(y(i,j))
         end do
      end do
      if (debug) then
         write (io_log,'(a)') 'ys:'
         write (io_log,'(10e10.3)') ((ys(i,j),i=1,n),j=1,n)
      end if
      
      ! Form statistic and compare with critical value.
      
      call dot1 (ys, dxos, temp, nmax2, nmax2, nmax2, 1, nmax2, 1, n, n, 1)
      if (debug) then
         write (io_log,'(a)') 'temp = '
         write (io_log,'(10e10.3)') (temp(i),i=1,n)
      end if
      kocs2 = 0.
      do i = 1,n
         kocs2 = kocs2 + temp(i)*dxos(i)
      end do
      write (io_cal,'(/a)') 'Individual event contributions to kocs2:'
      i = 0
      do iev = 1,nev
         if (cal_event(iev,3)) then
            i = i + 1
            k = i*2
            j = k-1
            write (io_cal,'(i3,2x,f10.5)') iev, temp(j)*dxos(j) + temp(k)*dxos(k)
         end if
      end do
      
      if (kocs2 .le. kcrit90) then
         nulltext = '; null hypothesis cannot be rejected'
         nullhyp = .true.
      else
         nulltext = '; null hypothesis is rejected'
         nullhyp = .false.
      end if
      
      write (io_cal,'(a,f6.2,a,f6.2,a)') 'rdbt_test = ', rdbt,'; observed value = ', kocs2, nulltext
      
      return
      end
      

!************************************************************************
      subroutine shorten (c1, c2, incr)
      
      implicit none
      
      real :: c1, c2, cl, f, incr
      logical :: debug
      character(len=132) :: msg
      
      data debug /.false./
            
      cl = sqrt(c1*c1+c2*c2)
      if (cl .lt. incr) return
      
      f = (cl-incr)/cl
      
      c1 = amax1(c1*f, 0.)
      c2 = amax1(c2*f, 0.)
      
      if (debug) then
         write (msg,'(a,2f10.4)') 'c1,c2 = ', c1, c2
         call logit (trim(msg))
      end if
      
      return
      end


!-----------------------------------------------------------------------      
      subroutine rdbt_test2 (rvlmax, pc)

! Test for inconsistency between the relative locations of calibration events as
! determined by HD and the relative locations as inferred from the provided calibration
! data. Both sets can have bias. The test is based on the calculation of coverage statistics
! for the 90% confidence ellipses (combining HD and calibration uncertainties) and the
! residual shift vectors. The cumulative binomial probability distribution is used to
! determine if the observed number of "uncovered" shift vectors is consistent with the
! stated uncertainty of the associated confidence ellipses. If the number of "uncovered"
! shift vectors is too large, an increment (rcdv) is added to the covariance matrices of all
! events and the test is repeated. When the test is satisfied, the incremental covariance is
! converted to a "radius of doubt". This test is really only appropriate if there are a large
! (> 10) number of calibration events.
      
      implicit none
      
      include 'mloc.inc'
      include 'cal_shift.inc'
      
      logical nullhyp
            
      integer :: iev, nr, nrc, k, j
      real :: incr, rvlmax, prob, betai, pc, pc1, pc2, covp, coverage, rdbtsq
      character(len=40) :: nulltext
      character(len=132) :: msg
      
      data incr /0.2/ ! Increment in uncertainty, in km
      data pc2 /0.10/ ! Threshold of probability (not coverage) to decide if null hypothesis is rejected
      
      write (io_cal,'(/a)') 'Radius of doubt based on coverage statistics'
      
      rdbt = 0.
      
      ! Fill arrays for calibration events
      cvs = 0.
      nr = 0
      nrc = 0
      do iev = 1,nev
         if (cal_event(iev,3)) then
            nr = nr + 1
!            write (io_cal,'(i3,2f8.2)') iev, rsvx(iev), rsvy(iev)
            cvs(nr,1) = rsvy(iev) ! dy = Latitude
            cvs(nr,2) = rsvx(iev) ! dx = Longitude
            covp = 0.
            call coverage_ellipse (cvs(nr,2), cvs(nr,1), sgcv(iev,1,1), sgcv(iev,1,2), sgcv(iev,2,2), covp)
            if (covp .le. 1.0) nrc = nrc + 1
         end if
      end do
      
      write (io_cal,'(/a,f4.2)') 'Threshold precentage of probability for radius of doubt test: ', pc2
      
      ! Test of whether the observed number of events with covp > 1.0 (outside the confidence
      ! ellipse) is significant (at confidence level given by pc), using cumulative binomial
      ! probability distribution (Numerical Methods, section 6.3).
      ! 
      k = nr - nrc ! number of events with covp > 1.0
      j = nr - k + 1
      pc1 = 1.0 - pc
      ! Incomplete beta function is not defined for k = 0 (all calibration events covered).
      if (k .gt. 0) then
         prob = betai(real(k),real(j),pc1)
      else if (k .eq. 0) then
         prob = 1.
      else
         write (msg,'(a,i3)') 'rdbt_test2: illegal value of k: ', k
         call oops (trim(msg))
      end if
      nullhyp = (prob .ge. pc2)
      if (nullhyp) then
         nulltext = '; null hypothesis cannot be rejected'
      else
         nulltext = '; null hypothesis is rejected'
      end if
      coverage = real(nrc)/real(nr)
      write (io_cal,'(t19,a,t27,a,t33,a,t40,a,t45,a,t56,a)') 'rdbt', 'nr', 'nrc', 'k','coverage', 'P(X  k)'
      write (io_cal,'(a,f10.3,3i6,2f10.3,a)') 'rdbt_test = ', rdbt, nr, nrc, k, coverage, prob, trim(nulltext)
      
      do while (.not.nullhyp)
         nr = 0
         nrc = 0
         rdbt = rdbt + incr
         rdbtsq = rdbt*rdbt
         do iev = 1,nev
            if (cal_event(iev,3)) then
               nr = nr + 1
               covp = 0.
               call coverage_ellipse (cvs(nr,2), cvs(nr,1), sgcv(iev,1,1)+rdbtsq, sgcv(iev,1,2), sgcv(iev,2,2)+rdbtsq, covp)
               if (covp .le. 1.0) nrc = nrc + 1
            end if
         end do
         k = nr - nrc
         j = nr - k + 1
         ! Incomplete beta function is not defined for k = 0 (all calibration events covered).
         if (k .gt. 0) then
            prob = betai(real(k),real(j),pc1)
         else if (k .eq. 0) then
            prob = 1.
         else
            write (msg,'(a,i3)') 'rdbt_test2: illegal value of k: ', k
            call oops (trim(msg))
         end if
         nullhyp = (prob .ge. pc2)
         if (nullhyp) then
            nulltext = '; null hypothesis cannot be rejected'
         else
            nulltext = '; null hypothesis is rejected'
         end if
         coverage = real(nrc)/real(nr)
         write (io_cal,'(a,f10.3,3i6,2f10.3,a)') 'rdbt_test = ', rdbt, nr, nrc, k, coverage, prob, trim(nulltext)
         if (rdbt .gt. rvlmax) then
            write (io_cal,'(a,f4.1,a)') 'rdbt_test exceeds rvlmax (', rvlmax,')'
            exit
         end if
      end do
      
      write (io_cal,'(a,f6.2)') 'rdbt2 = ', rdbt
      
      return
      end


!************************************************************************
      subroutine coverage_ellipse (dx, dy, cv11, cv12, cv22, covp)

! Determines if a point defined by (dx,dy) is inside the ellipse defined by
! the covariance matrix (cv11,cv12,cv22). (dx,dy) is relative to the origin of
! the ellipse. covp is less than or equal to 1.0 if the ellipse covers the point.

      implicit none

      real dx, dy, cv11, cv12, cv22, cv(2,2), alpha_sminor, alpha_smajor, sminor, smajor
      real fi, fir, rpd, u, v, covp, al, bl, a, b
      logical debug
      character(len=132) :: msg
      
      data debug/.false./
      
      rpd = 3.14159/180.
      
      if (debug) then
         write (msg,'(a,5f10.3)') 'coverage_ellipse: ', dx, dy, cv11, cv12, cv22
         call debugger (trim(msg))
      end if
      
      ! Convert covariance matrix into ellipse parameters
      cv(1,1) = cv11
      cv(1,2) = cv12
      cv(2,1) = cv12
      cv(2,2) = cv22
      call elips (cv, alpha_sminor, al, bl)
      sminor = sqrt(1./al)
      smajor = sqrt(1./bl)
      alpha_smajor = alpha_sminor + 90.
      if (debug) then
         write (msg,'(a,4f10.3)') 'coverage_ellipse: ', alpha_sminor, sminor, alpha_smajor, smajor
         call debugger (trim(msg))
      end if
      
      ! Rotate the point into the ellipse coordinate system (semi-major axis defines X axis)
      fi = 450. - alpha_smajor
      if (fi .gt. 360.) fi = fi - 360.
      fir = fi*rpd
      u = dx*cos(fir) + dy*sin(fir)
      v = dy*cos(fir) - dx*sin(fir)
      if (debug) then
         write (msg,'(a,4f10.3)') 'coverage_ellipse: ', fi, fir, u, v
         call debugger (trim(msg))
      end if
      
      ! Substitute point coordinates into ellipse equation
      a = (u*u) / (smajor*smajor)
      b = (v*v) / (sminor*sminor)
      covp = a + b
      if (debug) then
         write (msg,'(a,3f10.3)') 'coverage_ellipse: ', a, b, covp
         call debugger (trim(msg))
      end if
      
      return
      end
      
