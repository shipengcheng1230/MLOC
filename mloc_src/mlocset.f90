!***********************************************************************
      subroutine mlocset ()
      
! Set up location problem: establish free parameters and compute partial derivatives and residuals.
! Call inversion, update parameters, and iterate until convergence. Make calibration shift if
! necessary. Produce all requested output files and plots.

! January 26, 1989 by eric bergman, numerous changes since...
      
      implicit none
      
      include 'mloc.inc'
      
      integer, parameter :: max=60
      integer :: nphase
      real :: tt(max), dtdd(max), dtdh(max), dddp(max)
      character(len=8) :: phcd(max)

      integer :: i, j, k, l, iev, mt, it, itlim_flag, k0, nphreid, nr, ic, id0, ip, ielev_m,&
       ird, atsp_hr, atsp_min, io, ios, lunit, n_metadata_found
      real :: otsp(nevmax,0:itmax1), otsh(0:itmax1), x1sum, x2sum, x3sum, x4sum, xlat, sxlat,&
       cxlat, xlon, sxlon, cxlon, psr, dtpdh, hgtcr, arrtime, ttobs, ttcomp, hms2s, dgkmlogc,&
       usrc(2), dtim, dlat, dlon, ddep, rlatdg1, rlondg1, dt1, dts, p1, p2, theta, glat, glon,&
       rdp0, atsp, atsp_sec, res_test, minimum_weight
      character(len=8) :: phtest
      character(len=5) :: blank5
      character(len=18) :: nsmd_coord
      character(len=30) :: blank30
      character(len=34) :: line34
      character(len=100) :: outfil, outfildat0 
      character(len=132) :: msg, dum1, dum2, line132, filename, linein, file_folder
      character(len=160) :: command_line
      logical :: convrg, ptype, stype, op, ex, metadata_found, metadata_found_i
      
      data itlim_flag /2/ ! Iteration limit for changing filter flags
      data minimum_weight /0.01/
      blank5 = ' '
      blank30 = ' '
            
      ! Output file for .dat0 file (event data as read in), only in cluster mode
      outfildat0 = trim(outfile)//'.dat0_mnf'
      open (io_dat0,file=outfildat0,status='new')
      write (io_dat0,'(a,3x,a,t110,a)') 'B', 'original disposition of '//trim(outfile), ' '
      write (io_dat0,'(a,t121,a)') 'F   MNF v1.3', ' '
      
      ! Read station file.
      write (*,'(/a)') 'Station coordinates...'
      call redsta 
      
      ! Read travel time and ellipticity correction datasets
      write (*,'(/a)') 'Travel time tables...'
      call redtab
      
      ! Travel time spread for individual phases
      write (*,'(/a)') 'Travel times spreads...'
      call rdttsprd ()
      
      ! Reading errors from previous run or use default values.
      write (*,'(/a)') 'Reading errors...'
      call getrderr ()
      
      ! Get starting locations from an HDF file
      if (read_hdf) then
         write (*,'(/a)') 'Starting locations from an HDF file...'
         call hdf_read ()
      end if
      
      ! List of stations suspected of reporting bogus depth phases
      if (bdp_list) then
         write (*,'(/a)') 'Stations suspected of reporting bogus depth phases...'
         call read_bdp ()
      end if
            
      ! For each event, read epicenter information and phase data.
      write (*,'(/a)') 'Phase data for each event...'
      do iev = 1,nev
         write (msg,'(a,i3,2(2x,a))') 'mlocset: event ', iev, trim(evtnam(iev)), trim(infile(iev))
         call fyi (trim(msg))
         write (io_log,'(a)') trim(msg)
         open (io_in,file=infile(iev),status='old')
         call read_mnf (iev, io_in)
         close (io_in) 
      end do
      write (io_dat0,'(a,t110,a)') 'EOF', ' '
      close (io_dat0)
      
      ! Differential time data
      if (diffdat) then
            open (io_diffdat,file=diffdatfilnam,status='old')
            call read_mnf (0, io_diffdat)
            close (io_diffdat)
      end if
      
      ! For indirect calibration if calibration data have been given in both the input file and the command file.
      ! Command file values take precedence.
      if (calibration) then
         do iev = 1,nev
            if (cal_event(iev,1)) then ! Values from command file
               cal_lat(iev,3) = cal_lat(iev,1)
               cal_lon(iev,3) = cal_lon(iev,1)
               cal_dep(iev,3) = cal_dep(iev,1)
               cal_hr(iev,3) = cal_hr(iev,1)
               cal_min(iev,3) = cal_min(iev,1)
               cal_sec(iev,3) = cal_sec(iev,1)
               rcv(iev,3,1,1) = rcv(iev,1,1,1)
               rcv(iev,3,1,2) = rcv(iev,1,1,2)
               rcv(iev,3,2,1) = rcv(iev,1,2,1)
               rcv(iev,3,2,2) = rcv(iev,1,2,2)
               rcv(iev,3,3,3) = rcv(iev,1,3,3)
               rcv(iev,3,4,4) = rcv(iev,1,4,4)
               cal_event(iev,3) = .true.
               ncal(3) = ncal(3) + 1
            else if (cal_event(iev,2)) then ! Values from input file
               cal_lat(iev,3) = cal_lat(iev,2)
               cal_lon(iev,3) = cal_lon(iev,2)
               cal_dep(iev,3) = cal_dep(iev,2)
               cal_hr(iev,3) = cal_hr(iev,2)
               cal_min(iev,3) = cal_min(iev,2)
               cal_sec(iev,3) = cal_sec(iev,2)
               rcv(iev,3,1,1) = rcv(iev,2,1,1)
               rcv(iev,3,1,2) = rcv(iev,2,1,2)
               rcv(iev,3,2,1) = rcv(iev,2,2,1)
               rcv(iev,3,2,2) = rcv(iev,2,2,2)
               rcv(iev,3,3,3) = rcv(iev,2,3,3)
               rcv(iev,3,4,4) = rcv(iev,2,4,4)            
               cal_event(iev,3) = .true.
               ncal(3) = ncal(3) + 1
            end if
         end do
      end if
      
      ! Readings that failed the station file date range
      if (n_failed_date_range .gt. 0) then
         write (msg,'(a,i4,a)') 'mlocset: ', n_failed_date_range, ' readings failed the station file date range'
         call fyi (trim(msg))
         write (io_stn_log,'(a)') trim(msg)
      end if
      
      ! Summary of missing stations
      if (n_miss_sta_total .gt. 0) then
         write (msg,'(a,i4,a)') 'mlocset: ', n_miss_sta_total, ' codes missing from the station files'
         call fyi (trim(msg))
         write (io_stn_log,'(/a)') trim(msg)
         if (nsmd) then ! Check missing stations against NEIC metadata and create the body of a supplemental station file
            io = lunit()
            filename = trim(station_path)//dirsym//trim(nsmd_file)
            open (io,file=filename,status='old')
            metadata_found = .false.
            n_metadata_found = 0
            write (io_log,'(a)') '5 Supplemental station file from NEIC Metadata'
            do i = 1,n_miss_sta_total
               metadata_found_i = .false.
               nsmd_coord = ' '
               do
                  read (io,'(a)',iostat=ios) linein
                  if (ios .lt. 0) exit
                  if (linein(4:8) .eq. n_miss_sta_list(i)(1:5)) then
                     if (metadata_found_i) then
                        if (linein(40:57) .ne. nsmd_coord) then ! Different coordinates for the same station code
                           write (io_log,'(a)') trim(linein)
                           nsmd_coord = linein(40:57)
                           write (msg,'(2a)') 'mlocset: multiple coordinates found in NEIC metadata for ', n_miss_sta_list(i)
                           call warnings (trim(msg))
                        end if
                     else ! First match
                        metadata_found_i = .true.
                        nsmd_coord = linein(40:57)
                        n_metadata_found = n_metadata_found + 1
                        write (io_log,'(a)') trim(linein)
                     end if
                  end if
               end do
               if (metadata_found_i) then
                  metadata_found = .true.
                  write (io_stn_log,'(3x,a,1x,i3,a)') n_miss_sta_list(i), n_miss_sta(i), ' instances (see .log file for listing)'
               else
                  write (io_stn_log,'(3x,a,1x,i3,a)') n_miss_sta_list(i), n_miss_sta(i), ' instances'
               end if
               metadata_found_i = .false.
               rewind (io)
            end do
            close (io)
            if (metadata_found) then
               write (msg,'(a,i6,a)') 'mlocset: ', n_metadata_found, ' missing station codes were found in the NEIC metadata file'
               call fyi (trim(msg))
               call warnings ('mlocset: A supplemental station file must be created in order to use these stations')
            end if
         else
            if (n_miss_sta_total .gt. 0) then
               do i = 1,n_miss_sta_total
                  write (io_stn_log,'(3x,a,1x,i3,a)') n_miss_sta_list(i), n_miss_sta(i), ' instances'
               end do
            end if
         end if
      end if
      
      ! Stations used
      write (io_stn_log,'(/i6,a)') nkstat, ' stations in the arrival time dataset for which coordinates were found:'
      do i = 1,nkstat
         call geogra (stalat(kstat(i)), glat) ! Geographic latitude
         glon = stalon(kstat(i))
         call set_longitude_range (glon, 0)
         ielev_m = int(ahgtr(kstat(i))*1.0e3) ! Elevation in m
         write (io_stn_log,'(a,1x,a,1x,a,1x,f8.4,1x,f9.4,1x,i5,1x,2a)') nstr1(kstat(i)), sta_agency(kstat(i)),&
          sta_deployment(kstat(i)), glat, glon, ielev_m, sta_author(kstat(i)), trim(duplication(kstat(i)))
      end do
      
      ! Initialize
      call model_init (otsp, otsh, mt)
      
      ! Relative depth phases
      call rdp ()
                     
!***********************************************************************
      !  location iteration loop starts here                                  
      it = 0
      convrg = .false.
      
      do while (it .le. nstep .and. .not.convrg) ! Iteration loop
         write (*,'(/a,i1,a)') ' Beginning iteration ', it, '...'
         write (io_log,'(/a,i1,a)') ' Beginning iteration ', it, '...'
         write (io_log,'(a)') 'iev ndat    ndatdl  ndatpr'
         do iev = 1,nev ! Loop over events
            call depset (depthp(iev,it), usrc)
            rlatdg1 = latpgc(iev,it)
            rlondg1 = lonp(iev,it)
            do i = 1,nst(iev) ! Loop over readings
            
               ! Theoretical travel time, azimuth, back-azimuth, delta, ray parameter, elipticity correction, station height correction.
               ! All relative to current event position.
               call mlocsteq (iev, i, rlatdg1, rlondg1, stladg(iev,i), stlndg(iev,i), depthp(iev,it), ahgts(iev,i),&
                phase(iev,i), ttcomp, delt(iev,i), azes(iev,i), azse(iev,i), psr, dtpdh, elcr(iev,i), hgtcr, latp(iev,it))
               psd(iev,i) = psr*rpd ! Convert ray parameter to sec/degree
               if (verbose_log) write (io_log,'(a,2i4,1x,a6,1x,4f10.3,f6.1,f6.3,1x,a,f8.2,f10.2)') 'mlocsteq (ttcomp): ',&
                iev, i, stname(iev,i), rlatdg1, rlondg1, stladg(iev,i), stlndg(iev,i), depthp(iev,it), ahgts(iev,i),&
                phase(iev,i), delt(iev,i), ttcomp
                              
               ! If the named phase doesn't exist at this distance (ttcomp = 0.), and we're allowing phases to be renamed, reset the phase name
               if (ttcomp .lt. 0.1 .and. phidird(iev,i)) then
                  if (ptype(phase(iev,i))) then
                     if (phase(iev,i)(1:3) .ne. 'PPP') then
                        phase(iev,i) = 'UNKNOWNP' ! Unknown P-type
                        fcode(iev,i) = 'p'
                     end if
                  else if (stype(phase(iev,i))) then
                     if (phase(iev,i)(1:3) .ne. 'SSS') then
                        phase(iev,i) = 'UNKNOWNS' ! Unknown S-type
                        fcode(iev,i) = 'p'
                     end if
                  end if
               end if
               
               ! Spread of the travel time model for a given phase and epicentral distance
               if (read_ttsprd) then ! Use estimates from a previous run
                  call ttsig2 (phase(iev,i), delt(iev,i), ttsprd(iev,i), ttoff(iev,i))
               else ! default values of spread, all zero-offset
                  call ttsig  (phase(iev,i), delt(iev,i), ttsprd(iev,i), ttoff(iev,i))
               end if
               
               ! Partial derivatives
               a(iev,i,1) = psr*cos(azes(iev,i)*rpd)/radius ! geocentric
               a(iev,i,2) = -psr*sin(azes(iev,i)*rpd)/radius
               a(iev,i,3) = dtpdh
               if (rel_phase(iev,i)) then ! Can't use relative phases to determine OT
                  a(iev,1,4) = 0.
               else
                  a(iev,i,4) = 1.
               end if
               
               !write (*,'(i3,1x,2a,3f10.3)') iev, stname(iev,i), phase (iev,i), a(iev,i,1), a(iev,i,2), a(iev,i,3)
               
               ! Station corrections 
               select case (tt_corr)
                  case (0) ! None
                     stacs(iev,i) = 0.0
                     sdstcs(iev,i) = 1.0
                  case (1) ! Station elevation corrections
                     stacs(iev,i) = hgtcr
                     sdstcs(iev,i) = 1.0
               end select
               
               ! Observed travel time in seconds, relative to current origin time for each event
               if (rel_phase(iev,i)) then
                  arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                  ttobs = arrtime
               else
                  arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                  if (ipah(iev,i) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
                  ttobs = arrtime - otsp(iev,it)
               end if
               tto(iev,i) = ttobs
               
               s(iev,i,it) = elcr(iev,i) + stacs(iev,i) ! Path anomalies: ellipticity and station corrections
               dt(iev,i,it) = ttobs - ttcomp ! Travel time residuals (uncorrected for path anomalies)
               
               ! Check for one-minute errors, first iteration only, only for unflagged readings
               if (it .eq. 0 .and. fcode(iev,i) .eq. ' ') then
                  res_test = abs(dt(iev,i,it))
                  if (res_test .ge. 55. .and. res_test .le. 65.) then
                     write (msg,'(a, i3,1x,a,1x,a,1x,i5, 1x,f6.2)') 'mlocset: possible one-minute error for ',&
                      iev, stname(iev,i), phase(iev,i), mnf_line(iev,i), dt(iev,i,it)
                     call warnings (trim(msg))
                  end if
               end if

               if (verbose_log) then
                  atsp = otsp(iev,it) + ttcomp + s(iev,i,it)
                  call timecr (atsp, atsp_hr, atsp_min, atsp_sec)
                  write (io_log,'(a,2i3,f6.2)') 'mlocsteq (theoretical arrival time):', atsp_hr, atsp_min, atsp_sec
               end if
               
               ! Relative depth phase residuals
               if ((phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ') .and. rel_depth_phase(iev,i) .gt. -100.) then
                  call trtm (delt(iev,i), max, nphase, tt, dtdd, dtdh, dddp, phcd)
                  do j = 1,nphase
                     if (phcd(j) .eq. phase(iev,i)) exit
                  end do
                  if (j .lt. nphase) then
                     !for pwP-P we would need to add water depth correction to theoretical TT, so it needs to be saved for each case
                     rdp0 = tt(j) - tt(1)
                     rdp_res(iev,i,it) = rel_depth_phase(iev,i) - rdp0
                  else
                     write (msg,'(a,i3,2(a,f5.1))') 'mlocset: depth phase not found for event ', iev, '; h = ',&
                      depthp(iev,it), '; delta = ', delt(iev,i) 
                     call warnings (trim(msg))
                     rdp_res(iev,i,it) = 999.
                  end if
               end if
               
               ! Weights on the basis of absolute residual, i.e., relative to the window for that phase.
               ! For each phase the window is defined by a multiple (command WIND) of the spread calculated for that phase
               ! from a previous run (.ttsprd file). There is a roll-off at the edge if inverse weighting has been
               ! selected (command WEIG); otherwise it is a strict cut-off. The window is also offset by the
               ! mean residual from the .ttsprd file. 
               ! The window is expanded at local distance (delta < windloclim) for direct calibration,
               ! to help avoid losing important readings because of a poor starting location.
               weight(iev,i) = 1.
               p1 = wind1*ttsprd(iev,i)
               p2 = wind2*ttsprd(iev,i)
               if (direct_cal .and. delt(iev,i) .lt. windloclim) then
                  p1 = p1*2.
                  p2 = p2*2.
               end if
               dts = dt(iev,i,it) - s(iev,i,it) - ttoff(iev,i)
               if (data_weight) then
                  if (abs(dts) .le. p1) then
                     weight(iev,i) = 1.0
                  else if (abs(dts) .gt. p2) then
                     weight(iev,i) = 0.
                  else if (p2 .gt. p1) then
                     theta = pi*0.5*(abs(dts)-p1)/(p2-p1)
                     weight(iev,i) = 1.0 - sin(theta)
                     ! Make sure the weight is not too small. This can lead to divide by zero problem in mlocinv from round-off error.
                     if (weight(iev,i) .lt. minimum_weight) weight(iev,i) = minimum_weight
                  end if
               else
                  if (abs(dts) .le. p2) then
                     weight(iev,i) = 1.
                  else
                     weight(iev,i) = 0.
                  end if
               end if
               
            end do ! End of loop over readings
            
            ! Set flags (fltrh-, fltrc-) to determine if a station will be used to estimate an improved hypocentroid.
            ! To avoid convergence problems because of stations slipping in and out of the data set on alternate
            ! iterations, the flag is not reset after the iteration defined by "itlim_flag".
            if (it .le. itlim_flag) then
               call stflt (iev, it)
            else
               ndat(iev,it) = ndat(iev,it-1)
               ndatdl(iev,it) = ndatdl(iev,it-1)
               ndatpr(iev,it) = ndatpr(iev,it-1)
            end if
            write (io_log,'(i3,1x,3(i4,4x))') iev, ndat(iev,it), ndatdl(iev,it), ndatpr(iev,it)
            
         end do ! End of loop over events
                  
         ! Determine the number of other events in the cluster with which each station-phase is associated,
         ! counting only those cases in which fltrh(iev,i) .eq. .false. It only matters that this number
         ! is at least 1, so we cut the search off when ic=1 for any station-phase. fltrc is reset whenever
         ! fltrh is (see above).
         if (it .le. itlim_flag) then
            do iev = 1,nev ! Outer loop over events
               do j = 1,nst(iev) ! Outer loop over readings
                  phtest = phase(iev,j)
                  ic = 0
                  connected(iev,j) = .false.
                  if (.not.fltrc(iev,j)) then
                     k0 = kcode(iev,j)
                     id0 = idiff(iev,j)
                     do l = 1,nev ! Inner loop over events
                        if (l .eq. iev) cycle ! Don't consider the same event, to avoid problems with multiple readings for the same station-phase.
                        do k = 1,nst(l) ! Inner loop over readings
                           if (k0 .eq. kcode(l,k) .and. phtest .eq. phase(l,k) .and. id0 .eq. idiff(l,k)) then
                              if (.not.fltrc(l,k)) then
                                 ic = ic + 1
                                 if (ic .eq. 1) then ! That's all we need to know
                                    connected(iev,j) = .true.
                                    go to 100
                                 end if
                              end if
                           end if
                        end do ! End inner loop over phases
                     end do ! End inner loop over events
                  end if
  100             continue
               end do ! End outer loop over phases
            end do ! End outer loop over events
         end if

         
         if (nstep .eq. 0 .or. mt .eq. 0) then
            call mlocout_summary (0, otsh, otsp) ! Summary data for the forward problem
            call mlocout_phase_data (0) ! Phase data for the forward problem
            if (tt5) call tt_local_5 (0, blank30, blank5, dum1, dum2) ! Local data
            if (tt5e) call single_event_tt5_5 (0) ! Single-event tt5 plots
            if (tt5s) call single_station_tt5_5 (0) ! Single-station tt5 plots
            return
         else
            call mlocinv (it) ! Least squares inversion  (solve a * dx = dt)
         end if
         
         ! Convert standard errors in epicenter from km to degrees
         sdxhath(1) = sdxhath(1)*dgkmla
         sdxhath(2) = sdxhath(2)*dgkmlogc(lathgc(it))
         do iev = 1,nev
            sdxhatc(iev,1) = sdxhatc(iev,1)*dgkmla
            sdxhatc(iev,2) = sdxhatc(iev,2)*dgkmlogc(latpgc(iev,it))
         end do
         
         ! Update hypocentroid
         lathgc(it+1) = lathgc(it) + delx0(1,it)*dgkmla
         lonh(it+1) = lonh(it) + delx0(2,it)*dgkmlogc(lathgc(it))
         depthh(it+1) = depthh(it) + delx0(3,it)
         otsh(it+1) = otsh(it) + delx0(4,it)
         call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1))
         
         ! Update cluster event hypocenters
         ! If a cluster vector parameter (e.g., depth) is fixed but the corresponding parameter for the hypocentroid is free,
         ! the cluster vector is updated with the change in the hypocentroid parameter.
         do iev = 1,nev ! Loop over events
         
            if (mindx(iev,1) .ne. 0 .or. latfh) then ! Latitude
               latpgc(iev,it+1) = latpgc(iev,it) + (delx0(1,it) + dxp(iev,1,it))*dgkmla
               call geogra (latpgc(iev,it+1), latp(iev,it+1)) ! Update latitude in geographic coordinates.
            else
               latpgc(iev,it+1) = latpgc(iev,it)
            end if
            
            if (mindx(iev,2) .ne. 0 .or. lonfh) then ! Longitude
               lonp(iev,it+1) = lonp(iev,it) + (delx0(2,it) + dxp(iev,2,it))*dgkmlogc(latpgc(iev,it))
            else
               lonp(iev,it+1) = lonp(iev,it)
            end if
            
            if (mindx(iev,3) .ne. 0 .or. depthfh) then ! Depth
               depthp(iev,it+1) = depthp(iev,it) + delx0(3,it) + dxp(iev,3,it)
               if (depthp(iev,it+1) .le. 0.) then
                  write (msg,'(a,i3,a)') 'mlocset: depth for event ', iev,' went negative - held at previous value'
                  call warnings (trim(msg))
                  depthp(iev,it+1) = depthp(iev,it)
               end if
               ddep = depthp(iev,it+1) - depthp(iev,it)
               if (abs(ddep) .gt. 50.) then
                  depthp(iev,it+1) = depthp(iev,it+1) + ddep/2.
                  write (msg,'(a,i3,a)') 'mlocset: depth change for event ', iev,' was damped'
                  call warnings (trim(msg))
               end if
            else
               depthp(iev,it+1) = depthp(iev,it)
            end if
            
            if (mindx(iev,4) .ne. 0 .or. timefh) then ! OT
               otsp(iev,it+1) = otsp(iev,it) + delx0(4,it) + dxp(iev,4,it)
            else
               otsp(iev,it+1) = otsp(iev,it)
            end if
             
         end do ! End of loop over events
         
         ! Convert OT in seconds to hours-minutes-seconds
         do iev = 1,nev
            call timecr (otsp(iev,it+1), hourp(iev,it+1), minp(iev,it+1), secp(iev,it+1))
         end do
         
         ! Print changes in cluster vectors
         write (io_log,'(a,2x,4(6x,a))') 'iev', 'dtim', 'dlat', 'dlon', 'ddep'
         do iev = 1,nev
            dtim = dxp(iev,4,it)
            dlat = dxp(iev,1,it)*dgkmla
            dlon = dxp(iev,2,it)*dgkmlogc(latpgc(iev,it))
            ddep = dxp(iev,3,it)
            write (io_log,'(i3,2x,4f10.3)') iev, dtim, dlat, dlon, ddep
         end do
         
         ! Recalculate the hypocentroid. This is necessary to keep the hypocentroid synchronized with the cluster vectors
         ! in the case where some parameter (e.g., depth) is fixed for some events and not others.
         x1sum = 0.
         x2sum = 0.
         x3sum = 0.
         x4sum = 0.
         do iev = 1,nev
            x1sum = x1sum + latpgc(iev,it+1)
            x2sum = x2sum + lonp(iev,it+1)
            x3sum = x3sum + depthp(iev,it+1)
            x4sum = x4sum + otsp(iev,it+1)
         end do
         lathgc(it+1) = x1sum/real(nev)
         call geogra (lathgc(it+1),lath(it+1)) ! Update in geographic coordinates
         lonh(it+1) = x2sum/real(nev)
         depthh(it+1) = x3sum/real(nev)
         otsh(it+1) = x4sum/real(nev)
         call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1))
         
         ! Hypocentroid for calculation of empirical path anomalies
         lath_epa = lathgc(it+1)
         lonh_epa = lonh(it+1)
         depthh_epa = depthh(it+1)
         
         ! Print changes in hypocentroid
         dtim = otsh(it+1) - otsh(it)
         dlat = lathgc(it+1) - lathgc(it)
         dlon = lonh(it+1) - lonh(it)
         ddep = depthh(it+1) - depthh(it)
         write (io_log,'(5x,4(6x,a))') 'dtim', 'dlat', 'dlon', 'ddep'
         write (io_log,'(a,4f10.3)') 'hyp: ', dtim, dlat, dlon, ddep
         
         ! Check for convergence
         convrg = .true.
         select case (convergence_test_index)
         
            case (0) ! Convergence tests are done separately on the cluster vectors and hypocentroid
               if (abs(lathgc(it+1) - lathgc(it)) .gt. cl_epi_h) then
                  convrg = .false.
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid latitude change: ', lathgc(it+1) - lathgc(it), ' deg'
               end if
               if (abs(lonh(it+1) - lonh(it)) .gt. cl_epi_h) then
                  convrg = .false.
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid longitude change: ', lonh(it+1) - lonh(it), ' deg'
               end if
               if (abs(depthh(it+1) - depthh(it)) .gt. cl_dep_h) then
                  convrg = .false.
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid depth change: ', depthh(it+1) - depthh(it), ' km'
               end if
               if (abs(otsh(it+1) - otsh(it)) .gt. cl_ot_h) then
                  convrg = .false.
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Origin time change: ', otsh(it+1) - otsh(it), ' sec'
               end if
               do iev = 1,nev
                  if (abs(dxp(iev,1,it)) .gt. cl_epi_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' latitude change: ', dxp(iev,1,it), ' km'
                  end if
                  if (abs(dxp(iev,2,it)) .gt. cl_epi_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' longitude change: ', dxp(iev,2,it), ' km'
                  end if
                  if (abs(dxp(iev,3,it)) .gt. cl_dep_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' depth change: ', dxp(iev,3,it), ' km'
                  end if
                  if (abs(dxp(iev,4,it)) .gt. cl_ot_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' origin time change: ', dxp(iev,4,it), ' sec'
                  end if
               end do
               
            case (1) ! Convergence tests are done on final event hypocenters
               do iev = 1, nev
                  dlat = (latpgc(iev,it+1) - latpgc(iev,it))/dgkmla
                  if (abs(dlat) .gt. cl_epi_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' latitude change: ', dlat, ' km'
                  end if
                  dlon = (lonp(iev,it+1) - lonp(iev,it))/dgkmlogc(latpgc(iev,it))
                  if (abs(dlon) .gt. cl_epi_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' longitude change: ', dlon, ' km'
                  end if
                  ddep = depthp(iev,it+1) - depthp(iev,it)
                  if (abs(ddep) .gt. cl_dep_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' depth change: ', ddep, ' km'
                  end if
                  dtim = otsp(iev,it+1) - otsp(iev,it)
                  if (abs(dtim) .gt. cl_ot_c) then
                     convrg = .false.
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' origin time change: ', dtim, ' sec'
                  end if
               end do
            
         end select
         
         ! Even if convergence criteria are met we require at least one iteration in order to generate values for eci (cluster vectors).
         if (convrg .and. it .eq. 0) convrg = .false.

         if (convrg .or. it .eq. nstep) then
            if (convrg) then
               write (*,'(/a,i1,a/)') 'Converged after ', it, ' iterations'
            else
               write (*,'(/a,i1,a/)') 'Not converged after ', it, ' iterations'
            end if
            
            ! Update geographic latitude for hypocentroid and cluster events
            call geogra (lathgc(it+1), lath(it+1))
            do iev = 1,nev
               call geogra (latpgc(iev,it+1), latp(iev,it+1))
            end do
            
            if (.not. bias_corr) then
               lathgc(it+1) = lathgc(it+1) + bcorr(1)*dgkmla
               call geogra (lathgc(it+1), lath(it+1)) ! Convert to geographic coordinates
               lonh(it+1) = lonh(it+1) + bcorr(2)*dgkmlogc(lathgc(it+1))
               depthh(it+1) = depthh(it+1) + bcorr(3)
               otsh(it+1) = otsh(it+1) + bcorr(4)
               call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1)) 
               do iev = 1,nev
                  latpgc(iev,it+1) = latpgc(iev,it+1) + bcorr(1)*dgkmla
                  call geogra (latpgc(iev,it+1), latp(iev,it+1)) ! Convert to geographic coordinates
                  lonp(iev,it+1) = lonp(iev,it+1) + bcorr(2)*dgkmlogc(latpgc(iev,it+1))
                  depthp(iev,it+1) = depthp(iev,it+1) + bcorr(3)
                  otsp(iev,it+1) = otsp(iev,it+1) + bcorr(4)
               end do
               do iev = 1,nev
                  call timecr (otsp(iev,it+1), hourp(iev,it+1),minp(iev,it+1), secp(iev,it+1))
               end do
            end if
            
            ! Calibration events
            if (calibration) then
            
               ! Calculate shift vector to match calibration locations
               call cal_shift (1, it, otsp)
               do iev = 1,nev ! Add shift vector to all cluster events
                  latp_cal(iev) = latp(iev,it+1) + del_cal_lat
                  lonp_cal(iev) = lonp(iev,it+1) + del_cal_lon
                  depthp_cal(iev) = depthp(iev,it+1) + del_cal_dep
                  ! Don't let focal depth be shifted to a negative value
                  if (depthp_cal(iev) .lt. 0.) then
                     write (msg,'(a,i3,a,f6.2,a)') 'mlocset: calibration shift wants to make depth negative for event',&
                      iev, ' (', depthp_cal(iev), '); reset to 0 km'
                     call warnings (trim(msg))
                     depthp_cal(iev) = 0.0
                  end if
                  otsp_cal(iev) = otsp(iev,it+1) + del_cal_tim
               end do
               ! Update calibrated hypocentroid for calculation of empirical path anomalies.
               lath_epa = lath_epa + del_cal_lat ! Geocentric
               lonh_epa = lonh_epa + del_cal_lon
               depthh_epa = depthh_epa + del_cal_dep
               
               ! Residuals relative to "calibrated" (shifted) locations
               do iev = 1,nev ! Loop over events
                  call geocen (latp_cal(iev), lonp_cal(iev), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                  call depset (depthp_cal(iev), usrc)
                  do i = 1,nst(iev) ! Loop over readings
                     call mlocsteq (iev, i, xlat, xlon, stladg(iev,i), stlndg(iev,i), depthp_cal(iev), ahgts(iev,i),&
                      phase(iev,i), ttcomp, delt_cal(iev,i), azes_cal(iev,i), azse_cal(iev,i), psr, dtpdh,&
                      elcr_cal(iev,i), hgtcr, latp_cal(iev))
                     psd(iev,i) = psr*rpd ! Convert ray parameter to sec/degree
                     
                     ! TT corrections
                     select case (tt_corr)
                     case (0) ! None
                        stacs_cal(iev,i) = 0.0
                        sdstcs_cal(iev,i) = 1.0
                     case (1) ! Station elevation corrections
                        stacs_cal(iev,i) = hgtcr
                        sdstcs_cal(iev,i) = 1.0
                     case (2) ! Patch corrections + station elevation corrections
                        if (delt_cal(iev,i) .ge. 25.) then ! Patch corrections valid for teleseismic only
                           stacs_cal(iev,i) = sca0s(iev,i) + hgtcr
                           sdstcs_cal(iev,i) = sd1s(iev,i)
                        else
                           stacs_cal(iev,i) = hgtcr
                           sdstcs(iev,i) = 1.0
                        end if
                     end select
                     
                     arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                     if (ipah(iev,i) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
                     
                     if (rel_phase(iev,i)) then
                        ttobs = real(ipam(iev,i)*60) + pas(iev,i)
                     else
                        ttobs = arrtime - otsp_cal(iev)
                     end if
                     tto(iev,i) = ttobs
                     s_cal(iev,i) = elcr_cal(iev,i) + stacs_cal(iev,i) ! Path anomalies: ellipticity and station corrections
                     dt_cal(iev,i) = ttobs - ttcomp ! Travel time residuals (uncorrected for path anomalies)
                     
                  end do ! End of loop over readings
               end do ! End of loop over events
               
            end if
            
            exit ! Break out of iteration loop
            
         else ! Not converged
         
            it = it + 1 ! Increment iteration counter
            
            ! Re-identify phases based on new location and OT
            if (phid .and. it .le. 1) then
               nphreid = 0
               do iev = 1,nev
                  do ird = 1,nst(iev)
                     if (phidird(iev,ird)) then
                        call phreid3 (it, iev, ird, nr)
                        nphreid = nphreid + nr
                     end if
                  end do
               end do
               write (msg,'(a,i6,a,i2)') 'mlocset: ', nphreid, ' phases re-identified for iteration ', it
               if (verbose_screen) call fyi (trim(msg))
               write (io_log,'(a)') trim(msg)
            end if
                        
         end if
         
      end do ! End of iteration loop
      
!***********************************************************************

      ! Printed output
      call mlocout_summary (it, otsh, otsp) ! Summary data (.summary), covariance (.cv)
      if (direct_cal) call dircal (it) ! Direct calibration (must be called before mlocout_hdf!)
      call mlocout_hdf (it) ! .hdf and .hdf_cal if there is calibration data
      call mlocout_phase_data (it) ! Phase data (.phase_data), large residuals (.xdat and .lres)
      call mlocout_rderr () ! Empirical reading errors (.rderr)
      if (diffdat) call mlocout_rderr_diff () ! Empirical reading errors for differential time data (.rderr_diff)
      call mlocout_ttsprd (it) ! TT spread for different phases (.ttsprd)
      if (cv_out) call mlocout_cv (it) ! Full covariance matrix
      
      ! GMT scripts
      call fyi ('mlocset: Base plots...')
      if (gmt_version .eq. 5) then
         call base_map (it, 0, .true., .true., .true., '_base') ! Base plot
         if (eplt) call base_map (it, 0, .false., .true., .false., '_ell') ! Confidence ellipses, no event numbers or vectors
         if (splt) call base_map (it, 0, .false., .false., .true., '_seis') ! Seismicity plot, no event numbers or confidence ellipses
         if (direct_cal) call dcal_map (it) ! direct calibration raypath plot
         do ip = 1,nplot_max
            do iev = 1,nev
               if (plot(iev,ip)) then
                  call base_map (it, ip, .true., .true., .true., '_sel') ! Selected events, base plot
                  exit
               end if
            end do
         end do
         if (n_xsec .gt. 0) then
            do i = 1,n_xsec
               call xsec_gmt_5 (it, i)
            end do
         end if      
      end if
      
      call mlocout_kml (it) ! KML file for display in Google Earth
      
      if (comcatout) call mlocout_comcat (it) ! ComCAT format (.comcat)

      if (blocout) call mlocout_bloc (it) ! BAYESLOC format (.bloc)
                  
      if (pukeout) call mlocout_puke (it) ! PUKE format (.puke)
            
      if (md_out) call map_dat (it) ! plot data for GMT (.map_dat)
            
      if (nitomo .gt. 0) then ! Tomography output files (.TOMO_PHASE.ITOMO.tomo)
         do i = 1,nitomo
            call mlocout_tomo (it, tomo_phase(i), itomo(i))
         end do
      end if
      
      if (ttou) then
         if (.not.calibration .and. .not.direct_cal) call warnings ('mlocset: empirical TTs are from an uncalibrated cluster')
         call mlocout_tt (it) ! Output empirical TT data for specific phases
      end if
      
      ! .datf file (event data files with flags and phase IDs as actually used)
      if (datfout) then
         call fyi ('mlocset: .datf file')
         outfil = trim(outfile)//'.datf_mnf'
         if (verbose_screen) then
            write (msg,'(3a,i3)') 'mlocset: opening ', trim(outfil), ' on unit ', io_datf
            call fyi (trim(msg))
         end if
         open (io_datf,file=outfil)
         inquire (unit=io_datf,opened=op)
         if (op) then
            write (io_datf,'(a,3x,a,t121,a)') 'B', 'final disposition of '//trim(outfile), ' '
            write (io_datf,'(a,t121,a)') 'F   MNF v1.3', ' '
            do iev = 1,nev
               call write_mnf_13 (io_datf, iev, it)
            end do
            write (io_datf,'(a,t121,a)') 'EOF', ' '
         else
            msg = 'mlocset: file '//trim(outfil)//' was not opened'
            call warnings (trim(msg))
         end if
         close (io_datf)
      end if
      
      ! Travel time and other plots
      if (verbose_screen) call fyi ('mlocset: travel-time plots')
      if (tt1) call tt_summary_5 (it) ! Summary travel time plot
      if (tt2) call tt_teleseismic_p_5 (it) ! Teleseismic P
      if (tt3) call tt_pkp_caustic_5 (it) ! PKP caustic
      if (tt4) call tt_near_source_5 (it) ! Near source readings
      if (tt5) call tt_local_5 (it, blank30, blank5, dum1, dum2) ! Local distance
      if (tt5e) call single_event_tt5_5 (it) ! Single-event tt5 plots
      if (tt5s) call single_station_tt5_5 (it) ! Single-station tt5 plots
      If (tt6) call tt_local_regional_5 (it) ! Local-regional distance
      if (tt7) call tt_local_regional_s_5 (it) ! Local-regional shear phases
      if (tt8) call tt_rdp_summary (it) ! Relative depth phase summary plot
      if (rdpp) call rdp_single_5 (it) ! Relative depth phase plots for single events
      if (tt9) call tt_s_minus_p_5 (it) ! S-P times
      if (epa_plot) call epa_plot_driver_5 (it) ! Empirical path anomaly plots
      if (fdhp) call focal_depth_histogram (it) ! Histogram of focal depths
      
      ! For ComCat output combine all plots into a single PDF, then delete the postscript files
      if (comcatout) then
         command_line = 'gmt psconvert -TF -Fmloc_plots -A1 '//trim(ccat_folder)//'/*.ps'
         call system (trim(command_line))
         command_line = 'mv mloc_plots.pdf '//trim(ccat_folder)//dirsym//trim(basename)//'_plots.pdf'
         call system (trim(command_line))
         command_line = 'rm '//trim(ccat_folder)//'/*.ps'
         call system (trim(command_line))
      end if
      
      return
      
      end subroutine mlocset


