!*****************************************************************************************
subroutine mlocout_hdf (it)

! Writes one or more files summarizing final locations and uncertainties. There are several scenarios:

! (1) No calibration: writes a .hdf file, uncertainties for relative location only (cluster vectors).
! (2) Direct calibration: writes a .hdf_dcal file, uncertainties for absolute location (cluster vector plus hypocentroid).
! (3) Indirect calibration: writes both a .hdf file (as above) and a .hdf_cal file, with full uncertainties.
! (4) Direct and indirect calibration: Writes both a .hdf_dcal and a .hdf_cal file. Indirect calibration takes precedence.

! This subroutine also writes the '.focal_mech' file if any focal mechanism data has been provided via the 'mech' command.
! This subroutine also writes output to the .log file for the command 'subc', creating the main body of a command file
! for a subcluster that is especially well-suited for direct calibration.

   implicit none
   
   include 'mloc.inc'
   
   character(len=4) :: cal_code, fdp, fdm
   character(len=100) :: outfil
   character(len=132) :: msg, fmt
   character(len=64) :: hypocenter_list(nevmax)
   character(len=3) :: rmag_pr
   character(len=1) :: free_depth, cal2
   real :: avh, saaz1, saaz2, secpr, alpha, xl1, xl2, ddep, dot, openazc, azesiev(ntmax0), azdif
   real :: rstadelc, cstadelc, openaz_limit, lon_test, min_depth_uncertainty
   integer :: it, iev, it1, hourpr, minpr, indx(ntmax0), i, j, ii, jj, ilast, ndatx, calibration_type
   integer :: n_subc, n_subc_delt
   
   data openaz_limit/200./
   data min_depth_uncertainty/0.1/ ! To avoid rounding to zero on output
   
   hypocenter_list = ' '
   
   if (direct_cal) then
      outfil = trim(outfile)//'.hdf_dcal' ! Hypocentroid + cluster vector uncertainties
   else
      outfil = trim(outfile)//'.hdf' ! Cluster vector uncertainties only
   end if
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_hdf: opening ', trim(outfil), ' on unit ', io_out
      call fyi (trim(msg))
   end if
   open (io_out,file=outfil,status='new')
   
   it1 = it + 1
   fmt = '(i4,4i3,f6.2,f10.5,f11.5,f7.2,1x,a1,a1,f6.2,1x,a3,a2,1x,a10,3i5,f7.2,f6.2,2(1x,a4),3f6.1,2(i4,f6.2),f7.1,1x,a4)'
   
   ! Main body of the command file for a subcluster for direct calibration
   if (subc_set) then
      write (io_log,'(/a)') 'Subcluster for direct calibration (command "subc")'
      write (io_log,'(a,f6.2)') ' Epicentral distance limit = ', subc_delt
      write (io_log,'(a,i3)') ' Minimum # of readings within that distance = ', subc_nmin
      write (io_log,'(a,i3)') ' Minimum connectivity with the rest of the cluster = ', subc_nconnect
      n_subc = 0
      do iev = 1,nev
         n_subc_delt = 0
         do i = 1,nst(iev)
            if (fcode(iev,i) .eq. ' ' .and. delt(iev,i) .le. subc_delt) n_subc_delt = n_subc_delt + 1
         end do
         if (n_subc_delt .ge. subc_nmin .and. ndatc(iev,it) .ge. subc_nconnect) then
            n_subc = n_subc + 1
            write (io_log,'(a)') 'memb'
            write (io_log,'(2a)') 'even ', trim(evtnam(iev))
            write (io_log,'(2a)') 'inpu ', trim(infile20(iev))
         end if
      end do   
      write (msg,'(a,i3,a)') 'mlocout_hdf: ', n_subc, ' events selected by command "subc" for subcluster'
      call fyi (trim(msg))
   end if
   
   if (direct_cal) write (io_log,'(/a)') 'From direct calibration:' ! Logging cal_ commands
   
   do iev = 1,nev ! Loop over events
   
      if (direct_cal) then
         alpha = alphadc(iev)
         xl1 = xl1dc(iev)
         xl2 = xl2dc(iev)
         ddep = ddepdc(iev)
         dot = dotdc(iev)
         calibration_type = 1
      else
         alpha = alphac(iev)
         xl1 = xl1c(iev)
         xl2 = xl2c(iev)
         ddep = sdxhatc(iev,3)
         dot = sdxhatc(iev,4)
         calibration_type = 0
      end if
      
      ! Calibration code
      call calibration_code (iev, calibration_type, cal_code)
      cal2 = cal_code(2:2)
      call uctolc (cal2,-1)
      
      ! Uncertainty of focal depth
      fdp = ' '
      fdm = ' '
      if (mindx(iev,3) .gt. 0) then ! Free depth solution
         free_depth = 'f'
         if (ddep .le. 99.) then
            write (fdp,'(f4.1)') max(min_depth_uncertainty,ddep)
            write (fdm,'(f4.1)') max(min_depth_uncertainty,ddep)
         end if
      else ! Take from default or assigned uncertainties
         free_depth = ' '
         if (depthp_plus(iev) .ge. min_depth_uncertainty .and. depthp_plus(iev) .le. 99.) write (fdp,'(f4.1)') depthp_plus(iev)
         if (depthp_minus(iev) .ge. min_depth_uncertainty .and. depthp_minus(iev) .le. 99.) write (fdm,'(f4.1)') depthp_minus(iev)
         ddep = max(depthp_minus(iev),depthp_plus(iev))
      end if
      
      ! Magnitude
      if (rmag(iev) .gt. 0.) then
         write (rmag_pr,'(f3.1)') rmag(iev)
      else
         rmag_pr = '   '
      end if
      
      avh = xl1*xl2*3.1415 ! Geometric area of 90% confidence ellipse
      
      ! Make sure the semi-axis azimuths are between 0 and 360 degrees
      saaz1 = alpha
      if (saaz1 .lt. 0.) saaz1 = saaz1 + 360.
      saaz2 = saaz1 + 90.
      if (saaz2 .gt. 360.) saaz2 = saaz2 - 360.
      
      ! Open azimuth based on defining phases for cluster vector
      do i = 1,nst(iev) ! Index table on azimuth for sorting
         azesiev(i) = azes(iev,i)
      end do
      call indexx (nst(iev), azesiev, indx)
      openazc = 0.
      do ii = 1,nst(iev)-1
         i = indx(ii)
         ilast = i
         if (connected(iev,i)) then
            do jj = ii+1,nst(iev)
               j = indx(jj)
               if (connected(iev,j)) then
                  azdif = azesiev(j) - azesiev(i)
                  openazc = max(openazc,azdif)
                  ilast = j
                  exit
               end if
            end do
         end if
      end do
      ! Gap between the largest and smallest azimuth
      do jj = 1,nst(iev)-1
         j = indx(jj)
         if (connected(iev,j)) then
            azdif = azesiev(j) + 360. - azesiev(ilast)
            openazc = max(openazc,azdif)
            exit
         end if
      end do
      if (openazc .ge. openaz_limit) then
         write (msg,'(a,i3,1x,a,1x,f5.1)') 'mlocout_hdf: large open azimuth, event ', iev, trim(evtnam(iev)), openazc
         call fyi (trim(msg))
      end if
      
      ! Closest and furthest station, based on defining stations for cluster vector
      rstadelc = 180.
      cstadelc = 0.
      do i = 1,nst(iev)
         if (connected(iev,i)) then
            rstadelc = min(rstadelc,delt(iev,i))
            cstadelc = max(cstadelc,delt(iev,i))
         end if
      end do
      
      ! Number of readings flagged as outliers (fcode = 'x')
      ndatx = 0
      do i = 1,nst(iev)
         if (fcode(iev,i) .eq. 'x') ndatx = ndatx + 1
      end do
      
      lon_test = lonp(iev,it1)
      call set_longitude_range (lon_test, longitude_range)
      
      ! Log the parameters for the cal_ command in case we want to use these events later as calibration events
      if (calibration_type .gt. 0) then
         write (io_log,'(i3,1x,i4,2i3,1x,a4,2i3,f6.2,f10.5,f11.5,f7.2,i4,2f6.2,f6.1,f6.2)')&
          iev,& ! Event number
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          'cal'//cal2,& ! cal_ command
          hourp(iev,it1),& ! Origin hour.
          minp(iev,it1),& ! Origin minute.
          secp(iev,it1),& ! Origin seconds.
          latp(iev,it1),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp(iev,it1),& ! Final depth (km).
          nint(saaz1),& ! Semi-axis azimuth.
          xl1,& ! Semi-axis length, km.
          xl2,& ! Semi-axis length, km.
          ddep,& ! Uncertainty in depth (km)
          dot ! Uncertainty in origin time (sec).
      end if
      
      ! Simplified hypocentral data written to the log file for easy import into a document.
      ! To avoid interference with the 'cal_' command output, it is written to a character array
      ! and printed later
      write (hypocenter_list(iev),'(i3,1x,i4,4i3,f6.2,f9.4,f10.4,f6.1,1x,a4,1x,a3,a2)')&
       iev,& ! Event number
       iyre(iev),& ! Origin year.
       mone(iev),& ! Origin month.
       idye(iev),& ! Origin day.
       hourp(iev,it1),& ! Origin hour.
       minp(iev,it1),& ! Origin minute.
       secp(iev,it1),& ! Origin seconds.
       latp(iev,it1),& ! Geographic latitude.
       lon_test,& ! Geographic longitude.
       depthp(iev,it1),& ! Final depth (km).
       cal_code,& ! Calibration code (GTCNU)
       rmag_pr,& ! magnitude.
       mmag(iev) ! magnitude scale

      ! HDF file for uncalibrated or direct calibration      
      write (hdfline(iev),fmt)&
       iyre(iev),& ! Origin year.
       mone(iev),& ! Origin month.
       idye(iev),& ! Origin day.
       hourp(iev,it1),& ! Origin hour.
       minp(iev,it1),& ! Origin minute.
       secp(iev,it1),& ! Origin seconds.
       latp(iev,it1),& ! Geographic latitude.
       lon_test,& ! Geographic longitude.
       depthp(iev,it1),& ! Final depth (km).
       depset_pr(iev),& ! How depth was set
       free_depth,& ! Free depth flag
       depth_inp(iev,1),& ! Depth from input file.
       rmag_pr,& ! magnitude.
       mmag(iev),& ! magnitude scale
       evid(iev),& ! Event ID, right-most 10 characters
       ndat(iev,it),& ! Number of observations contributed to hypocentroid estimation.
       ndatc(iev,it),& ! Number of observations used for cluster vector.
       ndatx,& ! Number of observations flagged as outliers (fcode = 'x')
       shatsqci(iev),& ! Normalized sample variance for the cluster vector 
       dot,& ! Uncertainty in origin time (sec). Formerly SE of position.
       fdp,& ! + uncertainty in depth (deeper), in km.
       fdm,& ! - uncertainty in depth (shallower), in km.
       rstadelc,& ! Epicentral distance of nearest station for cluster vector
       cstadelc,& ! Epicentral distance of farthest station for cluster vector
       openazc,& ! Largest open azimuth for cluster vector.
       nint(saaz1),& ! Semi-axis azimuth.
       xl1,& ! Semi-axis length, km.
       nint(saaz2),& ! Semi-axis azimuth.
       xl2,& ! Semi-axis length, km.
       avh,& ! Area of confidence ellipse, km**2.
       cal_code ! Calibration code (GTCNU)
      write (io_out,'(2a)') trim(hdfline(iev)), ' '//annotation(iev)
   end do
   
   write (io_log,'(/a)') 'Basic hypocenter list, uncalibrated or direct calibration:'
   do iev = 1,nev
      write (io_log,'(a)') hypocenter_list(iev)
   end do
   
   close (io_out)
   
   ! HDF file for indirect calibrated locations
   if (calibration) then
      hypocenter_list = ' '
      outfil=trim(outfile)//'.hdf_cal'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_hdf: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      open (io_out,file=outfil,status='new')
      write (io_log,'(/a)') 'From Indirect calibration:' ! Logging cal_ commands

      do iev = 1,nev ! Loop over events
      
         ! Calibration code
         call calibration_code (iev, 2, cal_code)
         cal2 = cal_code(2:2)
         call uctolc (cal2,-1)
         
         ! Uncertainty of focal depth
         fdp = ' '
         fdm = ' '
         if (mindx(iev,3) .gt. 0) then ! Free depth solution
            free_depth = 'f'
            if (ddep .le. 99.) then
               write (fdp,'(f4.1)') max(min_depth_uncertainty,sqrt(accv(iev,3,3))) ! Uncertainty of depth, calibration-shifted
               write (fdm,'(f4.1)') max(min_depth_uncertainty,sqrt(accv(iev,3,3))) ! Uncertainty of depth, calibration-shifted
            end if
         else ! Take from default or assigned uncertainties
            free_depth = ' '
            if (depthp_plus(iev) .le. 99.) write (fdp,'(f4.1)') max(min_depth_uncertainty,depthp_plus(iev))
            if (depthp_minus(iev) .le. 99.) write (fdm,'(f4.1)') max(min_depth_uncertainty,depthp_minus(iev))
         end if
         
         ! Magnitude
         if (rmag(iev) .gt. 0.) then
            write (rmag_pr,'(f3.1)') rmag(iev)
         else
            rmag_pr = '   '
         end if
         
         avh = xl1cg(iev)*xl2cg(iev)*pi ! Geometric area of 90% confidence ellipse (including calibration shift uncertainty)
         call timecr (otsp_cal(iev), hourpr, minpr, secpr)
         
         ! Make sure the semi-axis azimuths are between 0 and 360 degrees. 
         saaz1 = alphacg(iev)
         if (saaz1 .lt. 0.) saaz1 = saaz1 + 360.
         saaz2 = saaz1 + 90.
         if (saaz2 .gt. 360.) saaz2 = saaz2 - 360.
            
         ! Open azimuth based on defining phases for cluster vector
         do i = 1,nst(iev) ! Index table on azimuth for sorting
            azesiev(i) = azes(iev,i)
         end do
         call indexx (nst(iev), azesiev, indx)
         openazc = 0.
         do ii = 1,nst(iev)-1
            i = indx(ii)
            ilast = i
            if (connected(iev,i)) then
               do jj = ii+1,nst(iev)
                  j = indx(jj)
                  if (connected(iev,j)) then
                     azdif = azesiev(j) - azesiev(i)
                     openazc = max(openazc,azdif)
                     ilast = j
                     exit
                  end if
               end do
            end if
         end do
         ! Gap between the largest and smallest azimuth
         do jj = 1,nst(iev)-1
            j = indx(jj)
            if (connected(iev,j)) then
               azdif = azesiev(j) + 360. - azesiev(ilast)
               openazc = max(openazc,azdif)
               exit
            end if
         end do
!         if (openazc .ge. openaz_limit) then
!            write (msg,'(a,i3,1x,a,1x,f5.1)') 'mlocout_hdf: large open azimuth, event ', iev, trim(evtnam(iev)), openazc
!            call fyi (trim(msg))
!         end if
         
         ! Closest and furthest station, based on defining stations for cluster vector
         rstadelc = 180.
         cstadelc = 0.
         do i = 1,nst(iev)
            if (connected(iev,i)) then
               rstadelc = min(rstadelc,delt(iev,i))
               cstadelc = max(cstadelc,delt(iev,i))
            end if
         end do
         
         ! Number of readings flagged as outliers (fcode = 'x')
         ndatx = 0
         do i = 1,nst(iev)
            if (fcode(iev,i) .eq. 'x') ndatx = ndatx + 1
         end do
         
         lon_test = lonp_cal(iev)
         call set_longitude_range (lon_test, longitude_range)
         
         ! Log the parameters for the cal_ command in case we want to use these events later as calibration events
         write (io_log,'(i3,1x,i4,2i3,1x,a4,2i3,f6.2,f10.5,f11.5,f7.2,i4,2f6.2,f6.1,f6.2)')&
          iev,& ! Event number
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          'cal'//cal2,& ! cal_ command
          hourpr,& ! Origin hour.
          minpr,& ! Origin minute.
          secpr,& ! Origin seconds.
          latp_cal(iev),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp_cal(iev),& ! Final depth (km).
          nint(saaz1),& ! Semi-axis azimuth.
          xl1cg(iev),& ! Semi-axis length, km.
          xl2cg(iev),& ! Semi-axis length, km.
          ddep,& ! Uncertainty in depth (km)
          dot ! Uncertainty in origin time (sec).

         ! Simplified hypocentral data written to the log file for easy import into a document.
         ! To avoid interference with the 'cal_' command output, it is written to a character array
         ! and printed later
         write (hypocenter_list(iev),'(i3,1x,i4,4i3,f6.2,f9.4,f10.4,f6.1,1x,a4,1x,a3,a2)')&
          iev,& ! Event number
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          hourp(iev,it1),& ! Origin hour.
          minp(iev,it1),& ! Origin minute.
          secp(iev,it1),& ! Origin seconds.
          latp(iev,it1),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp(iev,it1),& ! Final depth (km).
          cal_code,& ! Calibration code (GTCNU)
          rmag_pr,& ! magnitude.
          mmag(iev) ! magnitude scale

         write (hdfline(iev),fmt)&
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          hourpr,& ! Origin hour.
          minpr,& ! Origin minute.
          secpr,& ! Origin seconds.
          latp_cal(iev),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp_cal(iev),& ! Final depth (km).
          depset_pr(iev),& ! How depth was set
          free_depth,& ! Free depth flag
          depth_inp(iev,1),& ! Depth from input file.
          rmag_pr,& ! magnitude.
          mmag(iev),& ! magnitude scale
          evid(iev),& ! Event ID, right-most 10 characters
          ndat(iev,it),& ! Number of observations contributed to hypocentroid estimation.
          ndatc(iev,it),& ! Number of observations used for cluster vector.
          ndatx,& ! Number of observations flagged as large outliers (fcode = 'x')
          shatsqci(iev),& ! Normalized sample variance for the cluster vector ! Adopted here on 1/19/2012
          sqrt(accv(iev,4,4)),& ! Uncertainty of OT, calibration shifted.
          fdp,& ! + uncertainty in depth (deeper), in km.
          fdm,& ! - uncertainty in depth (shallower), in km.
          rstadelc,& ! Epicentral distance of nearest station for cluster vector
          cstadelc,& ! Epicentral distance of farthest station for cluster vector
          openazc,& ! Largest open azimuth for cluster vector.
          nint(saaz1),& ! Semi-axis azimuth.
          xl1cg(iev),& ! Semi-axis length, km.
          nint(saaz2),& ! Semi-axis azimuth.
          xl2cg(iev),& ! Semi-axis length, km.
          avh,& ! Area of confidence ellipse, km**2.
          cal_code ! Calibration code (GTCNU)
         write (io_out,'(2a)') trim(hdfline(iev)), ' '//annotation(iev)
      end do ! End of loop over events
         
      write (io_log,'(/a)') 'Basic hypocenter list from indirect calibration:'
      do iev = 1,nev
         write (io_log,'(a)') hypocenter_list(iev)
      end do
   
      close (io_out)
   end if
   
   ! Focal mechanism output file
   ! First part is identical to HDF format
   if (focal_mech) then
      outfil=trim(outfile)//'.focal_mech'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_hdf: opening ', trim(outfil), ' on unit ', io_focal_mech
         call fyi (trim(msg))
      end if
      open (io_focal_mech,file=outfil,status='new')
      
      if (calibration) then ! Indirect calibration
         do iev = 1,nev
            if (focal_mech_iev(iev)) then
               lon_test = lonp_cal(iev)
               call set_longitude_range (lon_test, longitude_range)
               write (io_focal_mech,'(i4,4i3,f6.2,f8.3,f9.3,f6.1,1x,a1,1x,a)')&
                iyre(iev),& ! Origin year.
                mone(iev),& ! Origin month.
                idye(iev),& ! Origin day.
                hourpr,& ! Origin hour.
                minpr,& ! Origin minute.
                secpr,& ! Origin seconds.
                latp_cal(iev),& ! Geographic latitude.
                lon_test,& ! Geographic longitude.
                depthp_cal(iev),& ! Final depth (km).
                depset_pr(iev),& ! How depth was set
                focal_mech_line(iev) ! Focal mechanism data
            end if
         end do         
      else ! No calibration or direct calibration
         do iev = 1,nev
            if (focal_mech_iev(iev)) then
               lon_test = lonp(iev,it1)
               call set_longitude_range (lon_test, longitude_range)
               write (io_focal_mech,'(i4,4i3,f6.2,f8.3,f9.3,f6.1,1x,a1,1x,a)')&
                iyre(iev),& ! Origin year.
                mone(iev),& ! Origin month.
                idye(iev),& ! Origin day.
                hourp(iev,it1),& ! Origin hour.
                minp(iev,it1),& ! Origin minute.
                secp(iev,it1),& ! Origin seconds.
                latp(iev,it1),& ! Geographic latitude.
                lon_test,& ! Geographic longitude.
                depthp(iev,it1),& ! Final depth (km).
                depset_pr(iev),& ! How depth was set
                focal_mech_line(iev) ! Focal mechanism data
            end if
         end do
      end if
      
      close (io_focal_mech)
   end if
   
   return
   
end subroutine mlocout_hdf
   
