!***********************************************************************************************************************************
      subroutine mlocout_summary (it, otsh, otsp)

! Always prints two output files:
!   Printed summary output (.summary)
!   HDF format (.hdf)
! Optionally:
!   Covariance info (.cv)
!   HDF format for calibrated locations (.hdf_cal)

! January 26, 1989 by eric bergman
! Modified 9/20/94 by eab
! Individual event data importances added, 11/14/2001 by eab
! New format for covariance matrix output, 2/3/05 by eab.

      implicit none
      
      include 'mloc.inc'
      
      character*3 binlbl(12)
      character latprh*3, lonprh*3, depthprh*5, timeprh*11
      character*100 outfil
      character*12 dimpcode
      integer ndath(0:itmax1), ndatdlh(0:itmax1), ndatprh(0:itmax1), ndatflh, it, i, j, iev, iazh, ndtotal, ibin
      integer iacvi, iacvc, iacvf, iazc, iazt, it1, hourhr, minhr, iazp, iazes
      real prmax(0:itmax1), sumdx(4), sumdx0(4), sumdxp(4), nsvsprd, xsigma, hdh, cv1i, cv2i, dgkmlo, dcvi, dcvc, cv1f, cv2f, dcvf
      real dlatpkm, dlonpkm, sechr, otsh(0:itmax1), otsp(nevmax,0:itmax1), hms2s, hdt, hdc
      real dcal_lat, dcal_lon, dcal_dep, dcal_ot, hdp, dlatkm, dlonkm
      real delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg
      real sol_input_s, sol_input_lat, sol_input_lon, sol_input_dep
      real cal_input_s, cal_input_lat, cal_input_lon, cal_input_dep
      real :: lon_test, lon_test1, lon_test2
      character(len=132) :: msg
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone
      integer :: values(8)
      
      data binlbl/'NNE',' NE','ENE','ESE',' SE','SSE','SSW',' SW', 'WSW','WNW',' NW','NNW'/
      data dimpcode/'            '/
      
      it1 = it + 1
      
      prmax(0) = wind2
      do i = 1,itmax1
         prmax(i) = wind2
      end do
      
      outfil = trim(outfile)//'.summary'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_summary: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      open (io_out,file=outfil,status='new')
      
      ! Various data sets and procedures used in the analysis.
      write (io_out,'(2a)') ' Program: ', trim(version)
      write (io_out,'(2a)') ' Run: ', trim(basename)
      write (io_out,'(2a)') ' Author: ', trim(mloc_author)
      call date_and_time (date, time, zone, values)
      write (io_out,'(4a)') ' Date: ', date, ' ', time
      if (locmod) then
         write (io_out,'(5a)') ' Travel times: ', trim(taup_model), ' (', trim(locmodfname), ')'
      else
         write (io_out,'(2a)') ' Travel times: ', trim(taup_model)
      end if
      write (io_out,'(a,f6.2,a,f6.2,a)') ' Lg travel times: ', lg_a, ' + ', lg_b, ' * delta'
      write (io_out,'(a,f6.2,a,f6.2,a)') ' T-phase travel times: ', tphase_a, ' + ', tphase_b, ' * delta'
      select case (tt_corr)
         case (0)
            write (io_out,'(2a)') ' Station elevation correction: off',&
             ' (focal depth referenced to surface in direct calibration)'
         case (1)
            write (io_out,'(2a)') ' Station elevation correction: on',&
             ' (focal depth referenced to geoid in direct calibration)'
      end select
      if (suppstn) then
         do i = 1,n_supp_stn_file
            write (io_out,'(2a)') ' Supplemental station file: ', trim(suppfilnam(i))
         end do
      end if
      if (n_skip .gt. 0) then
         do i = 1,n_skip
            write (io_out,'(a,1x,a8,1x,a8,1x,a8)') ' Skipping station/phase/author :', (skip_params(i,j),j=1,3)
         end do
      end if
      write (io_out,'(2a)') ' Data flags used: ', dflag_pr
      write (io_out,'(2a)') ' Residuals weighted by reading errors: ', data_weight_pr
      write (io_out,'(2a)') ' Assume perfect theoretical TTs for hypocentroid (ttsprd = 0): ', pttt_pr
      if (data_weight) write (io_out,'(2a)') ' Reading errors: ', trim(rderrfname)
      write (io_out,'(a,3(f5.2,a))') ' Minimum allowed reading errors: ', rderr_min_loc, ' (local) /', rderr_min, ' (general) /',&
       rderr_min_depth, ' (depth phases)'
      if (rels_set) then
         write (io_out,'(a,2(f5.2,a))') ' Reading errors for local phases: ', rderr_loc_p, ' (P) /', rderr_loc_s, '(S)'
      else
         write (io_out,'(a)') ' Reading errors for local phases: Not set'      
      end if
      write (io_out,'(a,f6.2)') ' Cluster vector fudge factor (non-gaussian): ', radius_cvff
      write (io_out,'(2a)') ' Travel time spread: ', trim(ttsprdfname)
      write (io_out,'(a,2f5.1)') ' Windowing: ', wind1, wind2
      if (bias_corr) then
         write (io_out,'(a)') ' Hypocentroid bias correction: on'
      else
         write (io_out,'(a)') ' Hypocentroid bias correction: off'
      end if
      if (dflag) then
         write (io_out,'(a)') ' Data flags: used'
      else
         write (io_out,'(a)') ' Data flags: not used'
      end if
      if (read_hdf) write (io_out,'(2a)') ' Starting locations: ', trim(rhdf_filnam)
      
      if (direct_cal) then
         write (io_out,'(a)') ' Direct calibration: on'
      else
         write (io_out,'(a)') ' Direct calibration: off'
      end if
      
      ! Limits used to select arrival time data.
      write (io_out,'(a)') ' Distance ranges:'
      if (ponly) then
         if (nhlim .eq. 1) write (io_out,'(4x,a,1(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
         if (nhlim .eq. 2) write (io_out,'(4x,a,2(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
         if (nhlim .eq. 3) write (io_out,'(4x,a,3(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
      else
         write (io_out,'(4x,a,3(2f6.1,3x))') 'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim)
      end if
      write (io_out,'(4x,a,3(2f6.1,3x))') 'Cluster vectors: ', (clim(i,1),clim(i,2),i=1,nclim)
      
! Hypocentroid data

      ! Location parameters.
      if (latfh) then
         latprh = 'lat'
      else
         latprh = '   '
      End if
      if (lonfh) then
         lonprh = 'lon'
      else
         lonprh = '   '
      end if
      if (depthfh) then
         depthprh = 'depth'
      else
         depthprh = '     '
      end if
      if (timefh) then
         timeprh = 'origin time'
      else
         timeprh = '           '
      end if
      write (io_out,'(a,2(a3,a),a5,a,a11,a)') ' Hypocentroid location parameters: (', latprh, ', ', lonprh, ', ', depthprh,&
       ', ', timeprh, ')'
      
      ! Starting point for relocation.
      write (io_out,'(106x,a)') 'NDATA LOST'
      write (io_out,'(25x,a,8x,a,5x,a,7x,a,11x,a,3x,a,2x,a,2x,a,1x,a,1X,a)') 'TIME', 'LATITUDE', 'LONGITUDE', 'DEPTH', 'E**2/NSTA',&
       'SHAT', 'PRMAX',  'DELT', 'PRES', 'FLAG'
      lon_test = lonh(0)
      call set_longitude_range (lon_test, longitude_range)
      write (io_out,'(a,16x,2(i2,a),f5.2,f12.3,f14.3,f12.1/)') ' START', hourh(0), ':', minh(0), ':', sech(0), lath(0), lon_test,&
       depthh(0)
      
      ! Corrections at each iteration.
      do i = 0,it
         ndath(i) = 0
         ndatdlh(i) = 0
         ndatprh(i) = 0
         ndatflh = 0
         do iev = 1,nev
            ndath(i) = ndath(i) + ndat(iev,i)
            ndatdlh(i) = ndatdlh(i) + ndatdl(iev,i)
            ndatprh(i) = ndatprh(i) + ndatpr(iev,i)
            ndatflh = ndatflh + ndatfl(iev)
         end do
         write (io_out,'(a,i2,20x,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,2x,f8.1,a,i5,2x,f5.2,f7.1,1x,3i5)') ' ITER', i, delx0(4,i),&
          ' S', -delx0(1,i), ' KM', delx0(2,i), ' KM', delx0(3,i), ' KM', ehatsqh(i), '/', ndath(i), shath(i), prmax(i),&
          ndatdlh(i), ndatprh(i), ndatflh
      end do
      
      ! Bias correction term
      if (.not. bias_corr) then
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a)') ' BIAS', bcorr(4), ' S', bcorr(1), ' KM', bcorr(2), ' KM',&
          bcorr(3), ' KM'
      end if
      
      ! Cumulative difference: solution - starting point
      do i = 1,4
         sumdx0(i) = 0.
      end do
      do i = 0,it
         sumdx0(1) = sumdx0(1) - delx0(1,i) ! Sign change for geocentric latitude
         sumdx0(2) = sumdx0(2) + delx0(2,i)
         sumdx0(3) = sumdx0(3) + delx0(3,i)
         sumdx0(4) = sumdx0(4) + delx0(4,i)
      end do
      if (.not. bias_corr) then
         sumdx0(1) = sumdx0(1) - bcorr(1) ! Sign change for geocentric latitude
         sumdx0(2) = sumdx0(2) + bcorr(2)
         sumdx0(3) = sumdx0(3) + bcorr(3)
         sumdx0(4) = sumdx0(4) + bcorr(4)
      end if
      if ((sumdx0(1)**2 + sumdx0(2)**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
         hdh = sqrt(sumdx0(1)**2 + sumdx0(2)**2)
         iazh = nint(atan2(sumdx0(2),sumdx0(1))*xdpr)
      else
         hdh = 0.
         iazh = 0
      end if  
      if (iazh .lt. 0) iazh = iazh + 360
      write (io_out,'(/a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') ' CUMULATIVE', sumdx0(4), ' S',&
       sumdx0(1), ' KM', sumdx0(2), ' KM', sumdx0(3), ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
      
      ! Solution
      lon_test = lonh(it1)
      call set_longitude_range (lon_test, longitude_range)
      write (io_out,'(/a,13x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' SOLUTION', hourh(it1), ':', minh(it1), ':', sech(it1), lath(it1),&
       lon_test, depthh(it1)
      
      ! Standard errors and confidence ellipse
      write (io_out,'(a,t28,f6.2,a,4x,f6.3,a,4x,f6.3,a,3x,f6.2,a,3(f6.1,a))') ' STANDARD ERRORS', sdxhath(4), ' S', sdxhath(1),&
       ' DEG', sdxhath(2), ' DEG', sdxhath(3), ' KM    ELLIPSE:', alphah, ' DEG', xl1h, ' AND ', xl2h, ' KM'
      
      ! Calibration shift
      if (calibration) then
         dlatkm = del_cal_lat/dgkmla
         dlonkm = del_cal_lon/dgkmlo(lath(it1))
         if ((dlatkm**2 + dlonkm**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
            hdh = sqrt(dlatkm**2 + dlonkm**2)
            iazh = nint(atan2(dlonkm,dlatkm)*xdpr)
         else
            hdh = 0.
            iazh = 0
         end if  
         if (iazh .lt. 0) iazh = iazh + 360
         write (io_out,'(/a)')&
          ' ***************************************************************************************************************'
         write (io_out,'(a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f4.1,a)') ' SHIFT', del_cal_tim, ' S',&
          del_cal_lat, ' DEG', del_cal_lon, ' DEG', del_cal_dep, ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
          
         call timecr (otsh(it+1)+del_cal_tim, hourhr, minhr, sechr) 
         lon_test = lonh(it1) + del_cal_lon
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t23,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' CALIBRATED', hourhr, ':', minhr, ':', sechr,&
          lath(it1)+del_cal_lat, lon_test, depthh(it1)+del_cal_dep
          
         ! Calibration event correction vectors
         write (io_out,'(/a/)') ' Calibration event correction vectors (weights):'
         do iev = 1,nev
            if (cal_event(iev,3)) then ! Calculate difference between calibration and inversion solution
               dcal_lat = cal_lat(iev,3) - latp(iev,it1)
               lon_test = lonp(iev,it1)
               call set_longitude_range (lon_test, longitude_range)
               dcal_lon = cal_lon(iev,3) - lon_test
               dcal_dep = cal_dep(iev,3) - depthp(iev,it1)
               dcal_ot = hms2s(cal_hr(iev,3), cal_min(iev,3), cal_sec(iev,3)) - otsp(iev,it1)
               dlatpkm = dcal_lat/dgkmla
               dlonpkm = dcal_lon/dgkmlo(cal_lat(iev,3))
               if ((dlatpkm**2 + dlonpkm**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
                  hdp = sqrt(dlatpkm**2 + dlonpkm**2)
                  iazp = nint(atan2(dlonpkm,dlatpkm)*xdpr)
               else
                  hdp = 0.
                  iazp = 0
               end if  
               if (iazp .lt. 0) iazp = iazp + 360
               write (io_out,&
                '(a7,i3,2x,a16,t30,f6.2,a,f4.2,a,t44,f6.3,a,f4.2,a,t58,f6.3,a,f4.2,a,t70,f7.2,a,f4.2,a,t84,a,i4,a,t103,a,f4.1,a)')&
                ' EVENT ', iev, evtnam(iev)(1:16), dcal_ot, '(', w4(iev), ')', dcal_lat, '(', w12(iev), ')',& 
                dcal_lon, '(', w12(iev), ')', dcal_dep, '(', w3(iev), ')', 'AZIM = ', iazp, ' DEG', 'DIST = ', hdp, ' KM'
            end if
         end do
         
      end if
     
! Cluster vectors

      if (nev .eq. 1) go to 999
      
      write (io_out,'(/a)')&
       ' *************************************************************************************************************'
      write (io_out,'(a26,t30,a4,t40,a9,t55,a5)') ' cluster vector statistics:', 'ITER', 'E**2/NSTA', 'SHATC'
      do i = 0,it
         ndtotal = 0
         do iev = 1,nev
            ndtotal = ndtotal + ndatc(iev,i)
         end do
         write (io_out,'(t32,i1,t38,f8.1,a1,i5,t55,f5.2)') i, ehatsqc(i), '/', ndtotal, shatc(i)
      end do
      call croux (shatsqci, nev, nsvsprd)
      write (io_out,'(a,f6.3)') ' Spread of event normalized sample variances: ', nsvsprd
      
      ! Initial and final cluster vectors
      write (io_out,'(/t40,a,t69,a,t98,a)') 'INITIAL', 'CHANGE', 'FINAL'
      do iev = 1,nev ! Loop over events
      
         ! Initial cluster vector
         cv1i = (latp(iev,0) - lath(0))/dgkmla
         cv2i = (lonp(iev,0) - lonh(0))/dgkmlo(lath(0))
         if ((cv1i**2 + cv2i**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvi = sqrt(cv1i**2 + cv2i**2)
            iacvi = nint(atan2(cv2i,cv1i)*xdpr)
         else
            dcvi = 0.
            iacvi = 0
         end if  
         if (iacvi .lt. 0) iacvi = iacvi + 360
         
         ! Change in cluster vector
         do i = 1,4
            sumdxp(i) = 0.
         end do
         do i = 0,it 
            sumdxp(1) = sumdxp(1) - dxp(iev,1,i) ! Sign change for geocentric latitude
            sumdxp(2) = sumdxp(2) + dxp(iev,2,i)
            sumdxp(3) = sumdxp(3) + dxp(iev,3,i)
            sumdxp(4) = sumdxp(4) + dxp(iev,4,i)
         end do
         if ((sumdxp(1)**2 + sumdxp(2)**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvc = sqrt(sumdxp(1)**2 + sumdxp(2)**2)
            iacvc = nint(atan2(sumdxp(2),sumdxp(1))*xdpr)
         else
            dcvc = 0.
            iacvc = 0
         end if 
         if (iacvc .lt. 0) iacvc = iacvc + 360
         
         ! Final cluster vector
         cv1f = (latp(iev,it1) - lath(it1))/dgkmla
         cv2f = (lonp(iev,it1) - lonh(it1))/dgkmlo(lath(it1))
         if ((cv1f**2 + cv2f**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvf = sqrt(cv1f**2 + cv2f**2)
            iacvf = nint(atan2(cv2f,cv1f)*xdpr)
         else
            dcvf = 0.
            iacvf = 0
         end if 
         if (iacvf .lt. 0) iacvf = iacvf + 360
         
         write (io_out,'(a7,i3,2x,a16,3(f8.1,a9,i3,a4,5x))') ' EVENT ', iev, evtnam(iev)(1:16), dcvi, ' KM  AT  ', iacvi, ' DEG',&
          dcvc, ' KM  AT  ', iacvc, ' DEG', dcvf, ' KM  AT  ', iacvf, ' DEG'
      end do
      
      ! Data importance distribution, by 30-degree wedges
      write (io_out,'(/a)')&
       ' *******************************************************************************************************************'
      write (io_out,'(a,t30,12(3x,a3))') ' DATA IMPORTANCE:', (binlbl(i),i=1,12)
      write (io_out,'(2x,a,t30,12(f6.3))') 'HYPOCENTROID:', (dimph(i),i=1,12)
      do iev = 1,nev
         do ibin = 1,12
            if (dimpciev(iev,ibin) .lt. 0.010) then
               dimpcode(ibin:ibin) = '0'
            else
               dimpcode(ibin:ibin) = '-'
            end if
         end do
         write (io_out,'(2x,a,i4,2x,a16,t30,12(f6.3),3x,a)') 'EVENT', iev, evtnam(iev)(1:16), (dimpciev(iev,i),i=1,12), dimpcode
      end do
      write (io_out,'(2x,a,t30,12(f6.3))') 'CLUSTER VECTORS:', (dimpc(i),i=1,12)
      
      ! Cluster vector info about each event.
      shatsqc=shatc(it)*shatc(it)
      do iev=1,nev ! Loop over events
         write (io_out,'(/2a)') ' ***********************************************************************************************',&
         '*******************************'
         write (io_out,'(a,i3,t20,a,t60,a)') ' CLUSTER EVENT ', iev, trim(evtnam(iev)), trim(infile(iev))
         write (io_out,'(a,t9,f3.1)') ' MAG', rmag(iev)
         write (io_out,'(a,2(a3,a),a5,a,a11,a)') ' Cluster vector location parameters: (', latpr(iev), ', ', lonpr(iev),&
          ', ', depthpr(iev), ', ', timepr(iev), ')'
         
         ! Input hypocenter from data file and starting point for relocation
         write (io_out,'(t13,a,9x,a,8x,a,5x,a,7x,a,9x,a,4x,a,2x,a,2x,a,1x,a,1x,a)') 'DATE', 'TIME', 'LATITUDE', 'LONGITUDE',&
          'DEPTH', 'EC**2/NSTA', 'EH**2/NSTA', 'PRMAX', 'DELT', 'PRES', 'FLAG'
         lon_test = lon_inp(iev)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t10,2(i2,a),i4,3x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' INPUT', mone(iev), '/', idye(iev), '/', iyre(iev),&
          hour_inp(iev), ':', min_inp(iev), ':', sec_inp(iev), lat_inp(iev), lon_test, depth_inp(iev,1)
         lon_test = lonp(iev,0)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t10,2(i2,a),i4,3x,2(i2,a),f5.2,f12.3,f14.3,f12.1/)') ' START', mone(iev), '/', idye(iev), '/',&
          iyre(iev), hourp(iev,0), ':', minp(iev,0), ':', secp(iev,0), latp(iev,0), lon_test, depthp(iev,0)
         
         ! Corrections at each iteration
         do i = 0,it
            write (io_out,'(a,i2,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,f6.1,a,i4,3x,f6.1,a,i4,f7.1,1x,3i5)') ' ITER', i,&
             dxp(iev,4,i), ' S', -dxp(iev,1,i), ' KM', dxp(iev,2,i), ' KM', dxp(iev,3,i), ' KM', eciev(iev,i), '/',&
             ndatc(iev,i), ehiev(iev,i), '/', ndat(iev,i), prmax(i), ndatdl(iev,i), ndatpr(iev,i), ndatfl(iev)
         end do
         
         ! Cumulative difference: solution - starting point
         do i = 1,4
            sumdxp(i) = 0.
         end do
         do i = 0,it 
            sumdxp(1) = sumdxp(1) - dxp(iev,1,i) ! Sign change for geocentric latitude
            sumdxp(2) = sumdxp(2) + dxp(iev,2,i)
            sumdxp(3) = sumdxp(3) + dxp(iev,3,i)
            sumdxp(4) = sumdxp(4) + dxp(iev,4,i)
         end do
         if ((sumdxp(1)**2 + sumdxp(2)**2) .ge. 0.01) then ! No vector calculated unless length is 100 m or more
            hdc = sqrt(sumdxp(1)**2 + sumdxp(2)**2)
            iazc = nint(atan2(sumdxp(2),sumdxp(1))*xdpr)
         else
            hdc = 0.
            iazc = 0
         end if
         if (iazc .lt. 0) iazc = iazc + 360
         
         write (io_out,'(/a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') '  CUMULATIVE', sumdxp(4), ' S',&
          sumdxp(1), ' KM', sumdxp(2), ' KM', sumdxp(3), ' KM', 'AZIM = ', iazc, ' DEG', 'DIST = ', hdc, ' KM'
         
         ! Hypocentroid change
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a)') ' +HYPOCENTROID', sumdx0(4), ' S', sumdx0(1), ' KM',&
          sumdx0(2), ' KM', sumdx0(3), ' KM'
         
         ! Total change (cluster + hypocentroid)
         do i = 1,4
            sumdx(i) = sumdx0(i) + sumdxp(i)
         end do
         if ((sumdx(1)**2 + sumdx(2)**2) .ge. 0.01) then ! No vector calculated unless length is 100 m or more
            hdt = sqrt(sumdx(1)**2 + sumdx(2)**2)
            iazt = nint(atan2(sumdx(2),sumdx(1))*xdpr)
         else
            hdt = 0.
            iazt = 0
         end if  
         if (iazt .lt. 0) iazt = iazt + 360
         
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') '  TOTAL', sumdx(4), ' S', sumdx(1),&
          ' KM', sumdx(2), ' KM', sumdx(3), ' KM', 'AZIM = ', iazt, ' DEG', 'DIST = ', hdt, ' KM'
         
         ! Solution
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(/a,2(i2,a),i4,2x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' SOLUTION ', mone(iev), '/', idye(iev), '/',&
          iyre(iev), hourp(iev,it1), ':', minp(iev,it1), ':', secp(iev,it1), latp(iev,it1), lon_test, depthp(iev,it1)
         
         ! Standard errors and confidence ellipse
         write (io_out,'(a,t28,f6.2,a,4x,f6.3,a,4x,f6.3,a,3x,f6.2,a,3(f6.1,a))') ' STANDARD ERRORS', sdxhatc(iev,4), ' S',&
          sdxhatc(iev,1), ' DEG', sdxhatc(iev,2), ' DEG', sdxhatc(iev,3), ' KM    ELLIPSE:', alphac(iev), ' DEG', xl1c(iev),&
          ' AND ', xl2c(iev), ' KM'
         
         ! Change from input hypocenter (data file)
         sol_input_s = (real(minp(iev,it1)*60) + secp(iev,it1)) - (real(min_inp(iev)*60) + sec_inp(iev))
         sol_input_lat = latp(iev,it1) - lat_inp(iev)
         lon_test1 = lonp(iev,it1)
         call set_longitude_range (lon_test1, longitude_range)
         lon_test2 = lon_inp(iev)
         call set_longitude_range (lon_test2, longitude_range)
         sol_input_lon = lon_test1 - lon_test2
         sol_input_dep = depthp(iev,it1) - depth_inp(iev,1)
         call delaz (lat_inp(iev), lon_inp(iev), latp(iev,it1), lonp(iev,it1), delta, deltdg, deltkm, azeqst, azesdg, azsteq,&
          azsedg, 0)
         if (azesdg .lt. 0) azesdg = azesdg + 360.
         iazes = nint(azesdg)
         write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f5.1,a)') ' SOLUTION-INPUT', sol_input_s,&
          ' S', sol_input_lat, ' DEG', sol_input_lon, ' DEG', sol_input_dep, ' KM', 'AZIM = ', iazes, ' DEG', 'DIST =',&
         deltkm, ' KM'
         
         ! Calibration
         if (calibration) then
            write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f4.1,a)') ' SHIFT', del_cal_tim, ' S',&
             del_cal_lat, ' DEG', del_cal_lon, ' DEG', del_cal_dep, ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
            call timecr (otsp_cal(iev), hourhr, minhr, sechr) 
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            write (io_out,'(a,t23,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' CALIBRATED', hourhr, ':', minhr, ':', sechr,&
             latp_cal(iev), lon_test, depthp_cal(iev)
            ! Change from input hypocenter (data file)
!            cal_input_s = (real(minp(iev,it1)*60) + secp(iev,it1)) - (real(min_inp(iev)*60) + sec_inp(iev))
            cal_input_s = (real(minhr*60) + sechr) - (real(min_inp(iev)*60) + sec_inp(iev))
            cal_input_lat = latp_cal(iev) - lat_inp(iev)
            lon_test1 = lonp_cal(iev)
            call set_longitude_range (lon_test1, longitude_range)
            lon_test2 = lon_inp(iev)
            call set_longitude_range (lon_test2, longitude_range)
            cal_input_lon = lon_test1 - lon_test2
            cal_input_dep = depthp_cal(iev) - depth_inp(iev,1)
            call delaz (lat_inp(iev), lon_inp(iev), latp_cal(iev), lonp_cal(iev), delta, deltdg, deltkm, azeqst, azesdg,&
             azsteq, azsedg, 0)
            if (azesdg .lt. 0) azesdg = azesdg + 360.
            iazes = nint(azesdg)
            write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f5.1,a)') ' CALIBRATED-INPUT',&
             cal_input_s, ' S', cal_input_lat, ' DEG', cal_input_lon, ' DEG', cal_input_dep, ' KM', 'AZIM = ', iazes, ' DEG',&
             'DIST =', deltkm, ' KM'
         end if
         
         ! Normalized sample error for cluster vector and how many sigma out from expected value
         xsigma = (shatsqci(iev) - 1.)/nsvsprd
         write (io_out,'(/a,f6.3,a,f7.2,a)') ' Normalized sample variance for cluster vector: ', shatsqci(iev), ' (', xsigma, ')'
         
      end do ! End loop over events
      
      close (io_out)
            
999   continue
      
      return
      end
      
