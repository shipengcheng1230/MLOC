!*****************************************************************************************
subroutine read_mnf (iev, inn)

! Checks MNF format version of a data file and calls the correct subroutine to read it.
! If 'iev' is greater than zero then the subroutine expects to read a standard v1.3 MNF
! arrival time data file. If iev = 0, the subroutine expects to read a v1.5 MNF differential
! time data file.
! The first line of the data file is required to be a format record ('F' in column 1).

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, inn
   character(len=6) :: format_version, current_version_13
   character(len=15) :: format_line
   character(len=132) :: msg
   
   current_version_13 = '1.3.3 '
   
   read (inn,'(a15)') format_line
   if (format_line(1:1) .eq. 'F') then
      if (verbose_screen) then
         write (msg,'(2a)') 'read_mnf: ', trim(format_line)
         call fyi (trim(msg))
      end if
      format_version = format_line(10:15)
      if (format_version(1:3) .eq. '1.3') then ! Standard event file
         if (format_version .eq. current_version_13) then
            call read_mnf_13 (iev, inn, format_version)
         else ! Non-current 1.3-class version
            write (msg,'(3a)') 'read_mnf: MNF version (', trim(format_version), ') is not current'
            call warnings (trim(msg))
            call read_mnf_13 (iev, inn, format_version)
         end if
      else if (format_version(1:3) .eq. '1.5') then
         call read_mnf_15 ()
      else
         if (iev .eq. 0) then
            write (msg,'(4a)') 'read_mnf: unknown format version ', format_version, ' for', diffdatfilnam
         else
            write (msg,'(3a,i3,1x,a)') 'read_mnf: unknown format version ', format_version, ' for event ', iev, evtnam(iev)
         end if
         call oops (trim(msg))
      end if
   else
      if (iev .eq. 0) then
         write (msg,'(2a)') 'read_mnf: no MNF format record found for ', diffdatfilnam
      else
         write (msg,'(a,i3,1x,a)') 'read_mnf: no MNF format record found for event ', iev, evtnam(iev)
      end if
      call oops (trim(msg))
   end if
   
   return
   
end subroutine read_mnf


!*****************************************************************************************
subroutine read_mnf_13 (iev, inn, format_version)

! Read input data files in MNF (mloc native format) format v1.3.
      
   implicit none
   
   include 'mloc.inc'
   
   integer, parameter :: max = 60
   
   integer :: nphase
   real :: tt
   character(len=8) :: phcd
   common /taup/ nphase, tt(max), phcd(max)
   
   integer :: iev, inn, ird, ios, i, j, line_number, nphreid, nr, alpha_smin_in, iptim_base, iptim_spec
!   integer :: scale_length
   real :: p_arrtime, hms2s, smin_in, smaj_in, sece_sigma, res_test
!   real :: t11, t12, t22, 
   real :: usrc(2)
   logical :: h_set, m_set, i_set, skip_station, skip_phase, skip_author, depth_constraint
   character(len=1) :: depth_inp_c_dr
   character(len=5) :: psta, original_code
   character(len=6) :: format_version
   character(len=121) :: linein
   character(len=160) :: msg
               
   nphreid = 0
   mmag(iev) = '  ' ! default, in case there is no magnitude record
   
   ! Loop over all records
   h_set = .false.
   m_set = .false.
   i_set = .false.
   ird = 0
   line_number = 1 ! Format record has already been read
   do
      read (inn,'(a)',iostat=ios) linein
      
      ! Unexpected termination
      if (ios .lt. 0)  call oops ('read_mnf_13: EOF reached unexpectedly')
      
      ! Normal termination
      if (linein(1:3) .eq. 'EOF') exit
      
      ! Write to dat0 file (original data as read in)
      write (io_dat0,'(a)') linein
      
      ! Index for searching
      line_number = line_number + 1
      
      ! Event record
      if (linein(1:1) .eq. 'E')  then
         if (linein(3:3) .eq. '-') then
            write (msg,'(a,i3,1x,a)') 'read_mnf_13: No phase data for event ', iev, trim(evtnam(iev))
            call oops (trim(msg))
         end if
         if (format_version(4:5) .eq. '  ') then ! Event ID may be stored at the end of the line for MNF version '1.3'
            if (len_trim(linein) .eq. 121) then
               evid(iev) = linein(112:121)
               evid_source(iev) = ' '
               i_set = .true.
               write (msg,'(a,i3,1x,3a)') 'read_mnf_13: Event ID for event ', iev, trim(evtnam(iev)),&
                ' taken from event record: ', evid(iev)
               call fyi (trim(msg))
            end if
         end if
         if (len_trim(linein) .ge. 5) then
            event_comment(iev) = linein(5:80)
         else
            event_comment(iev) = ' '
         end if
         ! Take annotation from the first 20 characters of event_comment
         ! This will be over-written by an annotation given previously by the 'anno' command.
         if (event_comment(iev) .ne. ' ') then
            if (annotation(iev) .ne. ' ') then
               write (msg,'(a,i3,1x,2a)') 'read_mnf_13: event comment for event ', iev, trim(evtnam(iev)),&
                ' is being over-written by ANNO command'
               if (verbose_screen) call fyi (trim(msg))
            else
               annotation(iev) = event_comment(iev)(1:20)
               write (msg,'(a,i3,1x,4a)') 'read_mnf_13: annotation for event ', iev, trim(evtnam(iev)),&
                ' read from input file: "', trim(annotation(iev)), '"'
               if (verbose_screen) call fyi (trim(msg))
            end if
         end if
      end if
      
      ! Event ID record
      if (linein(1:1) .eq. 'I') then
         if (linein(3:3) .eq. '=' .or. .not.i_set) then
            evid_source(iev) = linein(5:10)
            j = len_trim(linein)
            if (j .gt. 21) then ! take right-most 10 characters
               i = j - 9
               evid(iev) = linein(i:j)
               i_set = .true.
            else if (j .ge. 12) then
               evid(iev) = adjustr(linein(12:j))
               i_set = .true.
            else
               evid(iev) = ' '
               i_set = .false.
            end if
         end if
      end if
         
      ! Hypocenter record
      if (linein(1:1) .eq. 'H') then
         if (linein(3:3) .eq. '=' .or. .not.h_set) then
            h_set = .true.
            if (verbose_screen) then
               write (msg,'(2a)') 'read_mnf_13: ', trim(linein)
               call fyi (trim(msg))
            end if
            
            ! Minimally-required parameters
            read (linein,'(t5,i4.4)') iyre(iev)
            read (linein,'(t10,i2)') mone(iev)
            read (linein,'(t13,i2)') idye(iev)
            read (linein,'(t16,i2)') hour_inp(iev)
            read (linein,'(t19,i2)') min_inp(iev)
            read (linein,'(t22,f5.2)') sec_inp(iev)
            read (linein,'(t35,f8.4)') lat_inp(iev)
            read (linein,'(t44,f9.4)') lon_inp(iev)
            time_inp(iev) = hms2s(hour_inp(iev), min_inp(iev), sec_inp(iev))
            
            ! Focal depth
            if (linein(73:73) .eq. '.') then
               read (linein,'(t70,f5.1)') depth_inp(iev,1)
            else
               write (msg,'(a,i3,1x,a)') 'read_mnf_13: no focal depth read from preferred hypocenter record for event',&
                iev, trim(evtnam(iev))
               call warnings (trim(msg))
               depth_inp(iev,1) = 10. ! To avoid a crash
            end if
            
            ! Depth code
            ! Turn off cluster default depth for this event if a constrained depth is found in the hypocenter record.
            depth_inp_c(iev) = linein(76:76)
            
            ! Depth uncertainties
            if (linein(81:81) .eq. '.') then
               read (linein(78:82),'(f5.1)') depth_inp(iev,2)
            else
               depth_inp(iev,2) = 99.9
            end if
            if (linein(87:87) .eq. '.') then
               read (linein(84:88),'(f5.1)') depth_inp(iev,3)
            else
               depth_inp(iev,3) = depth_inp(iev,2)
            end if
            
            ! Origin time uncertainty
            if (linein(30:30) .eq. ' ') then
               sece_sigma = 99.9
            else if (linein(30:30) .eq. '.') then
               read (linein,'(t28,f5.2)') sece_sigma
            end if
            
            ! Confidence ellipse for the epicenter
            if (linein(56:56) .ne. ' ' .and. linein(60:60) .eq. '.' .and. linein(66:66) .eq. '.') then
               read (linein(54:56),'(i3)') alpha_smin_in
               read (linein(58:62),'(f5.2)') smin_in
               read (linein(64:68),'(f5.2)') smaj_in
            else
               alpha_smin_in = 0
               smin_in = 99.9
               smaj_in = 99.9
            end if
            
            ! Calibration code
            cal_code_input(iev) = linein(90:93)
            
            ! Check for GTXX classification
            ! This is only informational, user must decide on declaring a calibration event
            if (linein(90:91) .eq. 'GT') then
               write (msg,'(a,i3,1x,4a)') 'read_mnf_13: event ', iev, trim(evtnam(iev)),&
                ' is declared as a GT event (', cal_code_input(iev), ')'
               call fyi (trim(msg))
            end if
            
            ! Author and ORID fields are not read by mloc, but they are written if the 'datf'
            ! command is given. See subroutine write_mnf_13.
            
            ! Attempt to match this event to one from the HDF file.
            ! If found, use the HDF hypocentral parameters for the starting location
            if (read_hdf) call hdf_match (iev)
            
            ! Set starting hypocentral parameters.
            call starting_hypocenter (iev)
            
            ! Initialize tau-p software for calculating travel times
            call depset (depthp(iev,0), usrc)
            
         end if
      end if
      
      ! Depth record, only used if it's been marked as the preferred depth, if it has a "constrained"
      ! depth code, and if the depth from the hypocenter record was not constrained.
      if (linein(1:3) .eq. 'D =') then
         depth_inp_c_dr = linein(11:11)
         if (depth_constraint(depth_inp_c_dr) .and. .not.depth_constraint(depth_inp_c(iev))) then
            depth_inp_c(iev) = depth_inp_c_dr
            read (linein(5:9),'(f5.1)') depth_inp(iev,1)
            if (linein(16:16) .eq. '.') then
               read (linein(13:17),'(f5.1)') depth_inp(iev,2)
            else
               depth_inp(iev,2) = 99.9
            end if
            if (linein(22:22) .eq. '.') then
               read (linein(13:17),'(f5.1)') depth_inp(iev,3)
            else
               depth_inp(iev,3) = depth_inp(iev,2)
            end if
         
            ! Reset starting focal depth (other parameters get set again too).
            call starting_hypocenter (iev)
         
            ! Initialize tau-p software for calculating travel times
            call depset (depthp(iev,0), usrc)
         end if
      end if
      
      ! Magnitude record
      ! If no magnitude record is declared as the preferred one, the first record is used.
      if (linein(1:1) .eq. 'M') then
         if (linein(3:3) .eq. '=' .or. .not.m_set) then
            read (linein(5:8),'(f4.2)') rmag(iev)
            mmag(iev) = linein(10:11)
            m_set = .true.
         end if
         if (debug) write (io_log,'(f4.2,1x,a)') rmag(iev), mmag(iev)
      end if
      
      ! Phase reading record
      if (linein(1:1) .eq. 'P') then
         ird = ird + 1
         mnf_line(iev,ird) = line_number
         diff_line(iev,ird) = 0
         if (ird .eq. ntmax0) call oops ('read_mnf_13: ird exceeds ntmax0')
         nst(iev) = ird
         fcode(iev,ird) = linein(3:3)
         stname(iev,ird) = linein(5:9)
         dist_az(iev,ird) = linein(12:21)
         no_phid(iev,ird) = linein(23:23)
         phase0(iev,ird) = linein(24:31)
         read (linein(33:36),'(i4)') ipay(iev,ird)
         read (linein(38:39),'(i2)') ipamo(iev,ird)
         read (linein(41:42),'(i2)') ipad(iev,ird)
         read (linein(44:45),'(i2)') ipah(iev,ird)
         read (linein(47:48),'(i2)') ipam(iev,ird)
         read (linein(50:55),'(f6.3)') pas(iev,ird)
         if (linein(55:55) .ne. ' ') then
            iptim_base = -3
         else if (linein(54:54) .ne. ' ') then
            iptim_base = -2
         else if (linein(53:53) .ne. ' ') then
            iptim_base = -1
         else
            iptim_base = 0
         end if
         if (linein(58:58) .eq. ' ') then
            iptim(iev,ird) = iptim_base
         else
            read (linein(57:58),'(i2)') iptim_spec
            iptim(iev,ird) = max0(iptim_base,iptim_spec)
         end if
         read (linein(60:64),'(f5.1)') resisc(iev,ird)
         adslcaa(iev,ird)(1:47) = linein(75:121)
         if (read_ad) then
            agency(iev,ird) = linein(75:79)
            deployment(iev,ird) = linein(81:88)
         else
            agency(iev,ird) = ' '
            deployment(iev,ird) = ' '
         end if
         channel(iev,ird) = linein(99:101)
         readsrc(iev,ird) = linein(103:110)
         sad(iev,ird) = stname(iev,ird)//' '//agency(iev,ird)//' '//deployment(iev,ird)
         
         ! Check for changed station code
         original_code = linein(90:94)
         if (original_code .ne. '     ' .and. original_code .ne. stname(iev,ird)) then
            if (verbose_screen) call fyi ('read_mnf_13: '//trim(original_code)//' changed to '//stname(iev,ird))
            if (verbose_log) write (io_stn_log,'(i3,1x,a)') iev, trim(evtnam(iev)), ': ',&
             trim(original_code)//' changed to '//stname(iev,ird)
         end if
         
         ! One-minute errors, only checked for unflagged readings
         if (fcode(iev,ird) .eq. ' ') then
            res_test = abs(resisc(iev,ird))
            if (res_test .ge. 55. .and. res_test .le. 65.) then
               write (msg,'(a,i5)') 'read_mnf_13: possible one-minute error for line ', mnf_line(iev,ird)
               call warnings (trim(msg))
            end if
         end if
                           
         ! Fix a few known issues with certain stations
         call station_check (iev, ird)
         
         if (debug) write (io_log,'(2i5,1x,a,1x,a,i5,4(i3),f7.3,i3,1x,a,1x,f6.1)') ird, mnf_line(iev,ird),&
          fcode(iev,ird), stname(iev,ird), ipay(iev,ird), ipamo(iev,ird), ipad(iev,ird), ipah(iev,ird),&
          ipam(iev,ird), pas(iev,ird), iptim(iev,ird), phase0(iev,ird), resisc(iev,ird)
         
         call stafind (iev, ird) !  Search master list of stations for station coordinates.
         if (debug) write (io_log,'(a,2f10.4)') 'Station coordinates: ', stladg(iev,ird), stlndg(iev,ird)
         
         ! The most recent P arrival is saved for comparison when depth phases are processed, in order to
         ! calculate relative depth phase times.
         if (phase0(iev,ird) .eq. 'P       ') then
            psta = stname(iev,ird)
            p_arrtime = hms2s(ipah(iev,ird), ipam(iev,ird), pas(iev,ird))
            if (ipah(iev,ird) .lt. hourp(iev,0)) p_arrtime = p_arrtime + 86400.
            if (debug) write (io_log,'(a,1x,f10.3)') psta, p_arrtime
         end if
         
         ! Several house-keeping tasks related to phase name
         call phase_utility (iev, ird, psta, p_arrtime)
         
         ! Set data flags   
         call dataflags (dflag, fcode(iev,ird))
         if (debug) write (io_log,'(2a)') 'After setting data flags, fcode = ', fcode(iev,ird)
         
         ! Flag duplicate readings
         call duplicates (iev, ird)
         
         ! Raypath geometric parameters
         call raypath (iev, ird)
         
         ! Skip readings with blank phase codes
         ! Tested against phase0, which is not changed by anything (e.g., phase_utility/pnclean)
         if (skip .and. fcode(iev,ird) .eq. ' ') then
            do i = 1,n_skip
               skip_station = trim(stname(iev,ird)) .eq. trim(skip_params(i,1)) .or. skip_params(i,1) .eq. '*'
               skip_phase = trim(phase0(iev,ird)) .eq. ' ' .and. trim(skip_params(i,2)) .eq. 'blank'
               skip_author = trim(readsrc(iev,ird)) .eq. trim(skip_params(i,3)) .or. skip_params(i,3) .eq. '*'
               if (skip_station .and. skip_phase .and. skip_author) then
                  fcode(iev,ird) = 's'
                  exit
               end if
            end do
         end if
         
         ! Phase re-identification
         if (phid .and. phidird(iev,ird)) then
            call phreid3 (0, iev, ird, nr)
            nphreid = nphreid + nr
            if (phase(iev,ird)(1:7) .eq. 'UNKNOWN') fcode(iev,ird) = 'p'
         end if
         
         ! Explicit command to skip readings on the basis of station code, phase code or author
         if (skip .and. fcode(iev,ird) .eq. ' ') then
            do i = 1,n_skip
               skip_station = trim(stname(iev,ird)) .eq. trim(skip_params(i,1)) .or. skip_params(i,1) .eq. '*'
               skip_phase = trim(phase(iev,ird)) .eq. trim(skip_params(i,2)) .or. skip_params(i,2) .eq. '*'
               skip_author = trim(readsrc(iev,ird)) .eq. trim(skip_params(i,3)) .or. skip_params(i,3) .eq. '*'
               if (skip_station .and. skip_phase .and. skip_author) then
                  fcode(iev,ird) = 's'
                  exit
               end if
            end do
         end if
         
         ! Reading errors
         call readerr (iev, ird)
         if (debug) write (io_log,'(a,f10.3)') 'Reading error: ', sdread(iev,ird)

         ! Initial setting of fltrh. Needed for calculation of open azimuth in model_init.
         if (fcode(iev,ird) .eq. ' ') then
            fltrhflag(iev,ird) = .false.
            fltrh(iev,ird) = .false.
         else
            fltrhflag(iev,ird) = .true.
            fltrh(iev,ird) = .true.
         end if
         if (debug) write (io_log,'(a,1x,l1)') 'fltrhflag = ', fltrhflag(iev,ird)
         if (debug) write (io_log,'(a,1x,l1)') 'fltrh = ', fltrh(iev,ird)        
         
      end if
      
   end do
   
   write (msg,'(3x,a,i6,a)') 'read_mnf_13: ', nphreid, ' phases re-identified for iteration 0'
   if (verbose_screen) call fyi (trim(msg))
   write (io_log, '(a)') trim(msg)
   
   return
   
end subroutine read_mnf_13


!*****************************************************************************************
subroutine write_mnf_13 (io_unit, iev, it)

! Write one event in MNF v1.3 format.
! The old 'header' records (hypocenter, depth and magnitude records) are read from the input file
! and copied forward. If they were in an older format they are translated to v1.3 as far as possible.

   implicit none
   
   include 'mloc.inc'
   
   character(len=3) :: azes_pr
   character(len=4) :: cal_code
   character(len=5) :: dts_pr, fdp, fdm
   character(len=6) :: delt_pr
   character(len=121) :: linein
   character(len=132) :: fmt, fmt0, fmt1, fmt2, fmt3
   character(len=1) :: depth_code, usage
   integer :: io_unit, iev, ird, it, it1, calibration_type, hourpr, minpr
!   integer :: scale_length
   real :: dts, dot, alpha, xl1, xl2, ddep, lat_out, lon_out, lon_test, min_depth_uncertainty, secpr, depth_out
   
   data min_depth_uncertainty/0.1/ ! Required by the format
   
   it1 = it + 1
   
   ! Event line
   write (io_unit,'(a,t5,a,t121,a)') 'E ', event_comment(iev), ' '
      
   ! Depth code
   if (mindx(iev,3) .gt. 0) then ! Free depth solution
      depth_code = 'm' ! Depth set by relocation in mloc with free depth
   else
      depth_code = depset_pr(iev)
   end if
   
   ! Hypocentral uncertainties
   if (calibration) then
      alpha = alphacg(iev)
      xl1 = xl1cg(iev)
      xl2 = xl2cg(iev)
      ddep = sqrt(accv(iev,3,3))
      dot = sqrt(accv(iev,4,4))
   else
      if (direct_cal) then
         alpha = alphadc(iev)
         xl1 = xl1dc(iev)
         xl2 = xl2dc(iev)
         ddep = ddepdc(iev)
         dot = dotdc(iev)
      else
         alpha = alphac(iev)
         xl1 = xl1c(iev)
         xl2 = xl2c(iev)
         ddep = sdxhatc(iev,3)
         dot = sdxhatc(iev,4)
      end if
   end if
   if (alpha .lt. 0.) alpha = alpha + 360.
   
   ! Print variables for uncertainty of focal depth
   fdp = ' '
   fdm = ' '
   if (mindx(iev,3) .gt. 0) then ! Free depth solution
      if (ddep .le. 99.) then
         write (fdp,'(f5.1)') max(min_depth_uncertainty,ddep)
         write (fdm,'(f5.1)') max(min_depth_uncertainty,ddep)
      end if
   else ! Take from default or assigned uncertainties
      if (depthp_plus(iev) .le. 99.) write (fdp,'(f5.1)') max(min_depth_uncertainty,depthp_plus(iev))
      if (depthp_minus(iev) .le. 99.) write (fdm,'(f5.1)') max(min_depth_uncertainty,depthp_minus(iev))
   end if
   
   ! Calibration code.
   ! Indirect calibration takes precedence over direct calibration
   If (calibration) then
      calibration_type = 2
   else if (direct_cal) then
      calibration_type = 1
   else
      calibration_type = 0
   end if
   call calibration_code (iev, calibration_type, cal_code)
   
   ! Hypocenter line
   fmt = '(2a,1x,i4,a,i2.2,a,2(i2.2,1x),i2.2,2f6.2,2f10.4,1x,i3,2f6.2,f6.1,1x,a1,2(1x,a5),1x,a4,1x,a8,1x,a,t121,a)'
   if (calibration) then
      call timecr (otsp_cal(iev), hourpr, minpr, secpr)
      lat_out = latp_cal(iev)
      lon_out = lonp_cal(iev)
      depth_out = depthp_cal(iev)
   else
      hourpr = hourp(iev,it1)
      minpr = minp(iev,it1)
      secpr = secp(iev,it1)
      lat_out = latp(iev,it1)
      lon_out = lonp(iev,it1)
      depth_out = depthp(iev,it1)
   end if
   call set_longitude_range (lon_out, longitude_range)
   usage = '='
!   if (cal_code_input(iev)(1:2) .eq. 'GT') then
!      read (cal_code_input(iev)(3:4),'(i2)') scale_length
!      if (scale_length .le. 1) then ! Use the input GT location as preferred
!         usage = ' '
!      else
!         usage = '='
!      end if
!   else
!      usage = '='
!   end if

   write (io_unit,fmt) 'H ', usage, iyre(iev), '/', mone(iev), '/', idye(iev), hourpr, minpr, secpr,&
    dot, lat_out, lon_out, nint(alpha), xl1, xl2, depth_out, depth_code, fdp, fdm, cal_code, mloc_author,&
    trim(basename), ' '   
     
   ! Retrieve previous header records
   open (io_in,file=infile(iev),status='old')
   read (io_in,'(a)') linein ! Format record
   if (linein(10:12) .eq. '1.3') then
      do
         read (io_in,'(a)') linein
         if (linein(1:1) .eq. 'H') then
            linein(3:3) = ' '
!            if (linein(90:91) .eq. 'GT') then
!               read (linein(92:93),'(i2)') scale_length
!               if (scale_length .le. 1) then ! Use the input GT location as preferred
!                  linein(3:3) = '='
!               else
!                  linein(3:3) = ' '
!               end if
!            else
!               linein(3:3) = ' '
!            end if
            read (linein(44:52),'(f9.4)') lon_test
            call set_longitude_range (lon_test, longitude_range)
            if (linein(49:49) .eq. ' ') then
               write (linein(44:48),'(f5.0)') lon_test
            else if (linein(50:50) .eq. ' ') then
               write (linein(44:49),'(f6.1)') lon_test
            else if (linein(51:51) .eq. ' ') then
               write (linein(44:50),'(f7.2)') lon_test
            else if (linein(52:52) .eq. ' ') then
               write (linein(44:51),'(f8.3)') lon_test
            else
               write (linein(44:52),'(f9.4)') lon_test
            end if
            write (io_unit,'(a121)') linein
         else if (linein(1:1) .eq. 'I') then
            write (io_unit,'(a121)') linein
         else if (linein(1:1) .eq. 'D') then
            linein(3:3) = ' '
            write (io_unit,'(a121)') linein
         else if (linein(1:1) .eq. 'M') then
            write (io_unit,'(a121)') linein
         else if (linein(1:1) .eq. 'P') then ! Bail out when the first phase record is encountered
            exit
         else if (linein(1:1) .eq. 'S') then
            exit
         else if (linein(1:3) .eq. 'EOF') then
            exit
         end if
      end do
   else ! Carry the input hypocenter and magnitude info forward from old formats
      fmt = '(2a,1x,i4,a,i2.2,a,i2.2,1x,i2.2,1x,i2.2,1x,f5.2,1x,5x,2f10.4,6x,6x,4x,f6.1,t95,a8,t121,a)'
      lon_test = lon_inp(iev)
      call set_longitude_range (lon_test, longitude_range)
      write (io_unit,fmt) 'H ', ' ', iyre(iev), '/', mone(iev), '/', idye(iev), hour_inp(iev), min_inp(iev),&
       sec_inp(iev), lat_inp(iev), lon_test, depth_inp(iev,1), hypo_author(iev), ' '
      if (rmag(iev) .gt. 0.)  then
         write (io_unit,'(a4,f3.1,1x,1x,a2,3x,1x,a,t121,a)') 'M = ', rmag(iev), mmag(iev), magnitude_author(iev), ' '
      end if
   end if
   close (io_in)
   
   ! Phase lines
   fmt0 = '(a1,1x,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,4(1x,i2),1x,f3.0,4x,i2,1x,a5,1x,a8,1x,a47)'
   fmt1 = '(a1,1x,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,4(1x,i2),1x,f4.1,3x,i2,1x,a5,1x,a8,1x,a47)'
   fmt2 = '(a1,1x,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,4(1x,i2),1x,f5.2,2x,i2,1x,a5,1x,a8,1x,a47)'
   fmt3 = '(a1,1x,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,4(1x,i2),1x,f6.3,1x,i2,1x,a5,1x,a8,1x,a47)'
   do ird = 1,nst(iev)
      if (calibration) then
         dts = dt_cal(iev,ird) - s_cal(iev,ird)
      else
         dts = dt(iev,ird,it) - s(iev,ird,it)
      end if
      if (abs(dts) .le. 99.9) then
         write (dts_pr,'(f5.1)') dts
      else
         dts_pr = '     '
      end if
      write (delt_pr,'(f6.2)') delt(iev,ird)
      write (azes_pr,'(i3)') nint(azes(iev,ird))
      ! If the station coordinates could not be found, use the delta-azimuth values from the input file
      if (fcode(iev,ird) .eq. 'm') then
         delt_pr = dist_az(iev,ird)(1:6)
         azes_pr = dist_az(iev,ird)(8:10)
      end if
      if (iptim(iev,ird) .eq. -3) then
         fmt = fmt3
      else if (iptim(iev,ird) .eq. -2) then
         fmt = fmt2
      else if (iptim(iev,ird) .eq. -1) then
         fmt = fmt1
      else if (iptim(iev,ird) .eq. 0) then
         fmt = fmt0
      else
         fmt = fmt0
      end if
      write (io_unit,fmt)&
       'P',& ! a1
       fcode(iev,ird),& ! a1
       stname(iev,ird),& ! a5
       delt_pr,& ! f6.2
       azes_pr,& ! i3
       no_phid(iev,ird),& ! a1
       phase(iev,ird),& ! a8
       ipay(iev,ird),& ! i4
       ipamo(iev,ird),& ! i2
       ipad(iev,ird),& ! i2
       ipah(iev,ird),& ! i2
       ipam(iev,ird),& ! i2
       pas(iev,ird),& ! f6.3
       iptim(iev,ird),& ! i2
       dts_pr,& ! a5
       phase0(iev,ird),& ! a8
       adslcaa(iev,ird) ! a47
   end do
   
   ! Stop line
   write (io_unit,'(a,t121,a)') 'STOP', ' '
   
   return
   
end subroutine write_mnf_13


!*****************************************************************************************
subroutine write_mnf_14 (io_unit, iev, it)

! Write one event in MNF v1.4 format. This format is a variant of the standard format that
! is used only for the ".comcat" output file.
! 7/5/2015: v1.41, providing two decimal places of precision for residuals
! 11/16/2017: v1.4.2, adding deployment/network code to station code for each phase reading

   implicit none
   
   include 'mloc.inc'
   
   character(len=2) :: c2
   character(len=3) :: azes_pr, c3
   character(len=4) :: cal_code
   character(len=5) :: fdp, fdm
   character(len=6) :: delt_pr, dts_pr
   character(len=121) :: linein
   character(len=132) :: fmt, fmt0, fmt1, fmt2, fmt3, msg
   character(len=1) :: depth_code, usage, c1, fcode_pr
   integer :: io_unit, iev, ird, it, it1, calibration_type, hourpr, minpr
!   integer :: scale_length
   real :: dts, dot, alpha, xl1, xl2, ddep, lon_test, min_depth_uncertainty, secpr
   
   data min_depth_uncertainty/0.1/ ! Required by the format
   
   it1 = it + 1
   
   ! Event line with event ID
   if (iev .le. 9) then
      write (c1,'(i1)') iev
      write (io_unit,'(a,t5,a,t121,a)') 'E ', 'cec_'//trim(basename)//'_'//c1, ' '
   else if (iev .le. 99) then
      write (c2,'(i2)') iev
      write (io_unit,'(a,t5,a,t121,a)') 'E ', 'cec_'//trim(basename)//'_'//c2, ' '
   else if (iev .le. 999) then
      write (c3,'(i3)') iev
      write (io_unit,'(a,t5,a,t121,a)') 'E ', 'cec_'//trim(basename)//'_'//c3, ' '
   end if
   
   ! Depth code
   if (mindx(iev,3) .gt. 0) then ! Free depth solution
      depth_code = 'm' ! Depth set by relocation in mloc with free depth
   else
      depth_code = depset_pr(iev)
   end if
   
   ! Hypocentral uncertainties
   if (calibration) then
      alpha = alphacg(iev)
      xl1 = xl1cg(iev)
      xl2 = xl2cg(iev)
      ddep = sqrt(accv(iev,3,3))
      dot = sqrt(accv(iev,4,4))
   else
      if (direct_cal) then
         alpha = alphadc(iev)
         xl1 = xl1dc(iev)
         xl2 = xl2dc(iev)
         ddep = ddepdc(iev)
         dot = dotdc(iev)
      else
         alpha = alphac(iev)
         xl1 = xl1c(iev)
         xl2 = xl2c(iev)
         ddep = sdxhatc(iev,3)
         dot = sdxhatc(iev,4)
      end if
   end if
   if (alpha .lt. 0.) alpha = alpha + 360.
   
   ! Print variables for uncertainty of focal depth
   fdp = ' 99.9'
   fdm = ' 99.9'
   if (mindx(iev,3) .gt. 0) then ! Free depth solution
      if (ddep .le. 99.) then
         write (fdp,'(f5.1)') max(min_depth_uncertainty,ddep)
         write (fdm,'(f5.1)') max(min_depth_uncertainty,ddep)
      end if
   else ! Take from default or assigned uncertainties
      if (depthp_plus(iev) .le. 99.) write (fdp,'(f5.1)') max(min_depth_uncertainty,depthp_plus(iev))
      if (depthp_minus(iev) .le. 99.) write (fdm,'(f5.1)') max(min_depth_uncertainty,depthp_minus(iev))
   end if
   
   ! Calibration code.
   ! Indirect calibration takes precedence over direct calibration
   If (calibration) then
      calibration_type = 2
   else if (direct_cal) then
      calibration_type = 1
   else
      calibration_type = 0
   end if
   call calibration_code (iev, calibration_type, cal_code)
   
   ! Hypocenter line
   ! Usage code is not printed in this format because there is only one hypocenter line, the preferred one
   fmt = '(a,t5,i4,a,i2.2,a,2(i2.2,1x),i2.2,2f6.2,2f10.4,1x,i3,2f6.2,f6.1,1x,a1,2(1x,a5),1x,a4,1x,a8,1x,a,t121,a)'
   if (calibration) then
      lon_test = lonp_cal(iev)
      call set_longitude_range (lon_test, longitude_range)
      usage = ' '
      call timecr (otsp_cal(iev), hourpr, minpr, secpr)
      write (io_unit,fmt) 'H', iyre(iev), ' ', mone(iev), ' ', idye(iev), hourpr, minpr, secpr,&
       dot, latp_cal(iev), lon_test, nint(alpha), xl1, xl2, depthp_cal(iev), depth_code, fdp, fdm, cal_code, mloc_author,&
       trim(basename), ' '
   else
      lon_test = lonp(iev,it1)
      call set_longitude_range (lon_test, longitude_range)
      usage = ' '
      write (io_unit,fmt) 'H', iyre(iev), ' ', mone(iev), ' ', idye(iev), hourp(iev,it1), minp(iev,it1), secp(iev,it1),&
       dot, latp(iev,it1), lon_test, nint(alpha), xl1, xl2, depthp(iev,it1), depth_code, fdp, fdm, cal_code, mloc_author,&
       trim(basename), ' '
   end if
       
   ! Retrieve magnitude records, fill in blank fields, remove usage code if present
   open (io_in,file=infile(iev),status='old')
   do
      read (io_in,'(a)') linein
      if (linein(1:1) .eq. 'M') then
         if (linein(3:3) .eq. '=') linein(3:3) = ' '
         if (linein(10:14) .eq. '     ') linein(10:14) = 'UNK  '
         if (linein(16:18) .eq. '   ') linein(16:18) = 'UNK'
         write (io_unit,'(a)') trim(linein(1:110))
      else if (linein(1:1) .eq. 'P') then ! Bail out when the first phase record is encountered
         exit
      else if (linein(1:1) .eq. 'S') then
         exit
      else if (linein(1:3) .eq. 'EOF') then
         exit
      end if
   end do
   close (io_in)

   ! Phase lines
   fmt0 = '(a1,1x,a1,1x,a6,1x,a8,1x,a6,1x,a3,1x,a8,1x,i4,4(1x,i2),1x,f3.0,4x,i2,1x,a6,1x,f5.2)'
   fmt1 = '(a1,1x,a1,1x,a6,1x,a8,1x,a6,1x,a3,1x,a8,1x,i4,4(1x,i2),1x,f4.1,3x,i2,1x,a6,1x,f5.2)'
   fmt2 = '(a1,1x,a1,1x,a6,1x,a8,1x,a6,1x,a3,1x,a8,1x,i4,4(1x,i2),1x,f5.2,2x,i2,1x,a6,1x,f5.2)'
   fmt3 = '(a1,1x,a1,1x,a6,1x,a8,1x,a6,1x,a3,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,i2,1x,a6,1x,f5.2)'
   do ird = 1,nst(iev)
   
      ! Some phase records are skipped for .comcat output
      if (fcode(iev,ird) .eq. 'd') cycle ! skip duplicate readings
      if (fcode(iev,ird) .eq. 'm') cycle ! skip readings from missing stations
      if (fcode(iev,ird) .eq. 'p') cycle ! skip readings for problematic phases
      if (fcode(iev,ird) .eq. 't') cycle ! skip readings with timing problems
        
      if (fcode(iev,ird) .eq. ' ') then
         if (.not.fltrh(iev,ird) .or. .not.fltrc(iev,ird)) then
            fcode_pr = '+'
         else
            fcode_pr = '-'
         end if
      else if (fcode(iev,ird) .eq. 's') then
         fcode_pr = '-'
      else if (fcode(iev,ird) .eq. 'x') then
         fcode_pr = 'x'
      else
         write (msg,'(a,i3,1x,a,1x,a,1x,i5)') 'write_mnf_14: unknown flag "'//fcode(iev,ird)//'" for ',&
          iev, stname(iev,ird), phase(iev,ird), mnf_line(iev,ird)
         call warnings (trim(msg))
         cycle
      end if      
      
      if (calibration) then
         dts = dt_cal(iev,ird) - s_cal(iev,ird)
      else
         dts = dt(iev,ird,it) - s(iev,ird,it)
      end if
      if (abs(dts) .le. 99.9) then
         write (dts_pr,'(f6.2)') dts
      else
         dts_pr = ' 99.99'
      end if
      write (delt_pr,'(f6.2)') delt(iev,ird)
      write (azes_pr,'(i3)') nint(azes(iev,ird))
      if (iptim(iev,ird) .eq. -3) then
         fmt = fmt3
      else if (iptim(iev,ird) .eq. -2) then
         fmt = fmt2
      else if (iptim(iev,ird) .eq. -1) then
         fmt = fmt1
      else if (iptim(iev,ird) .eq. 0) then
         fmt = fmt0
      else
         fmt = fmt0
      end if
      write (io_unit,fmt)&
       'P',& ! a1
       fcode_pr,& ! a1
       stname(iev,ird),& ! a5
       deployment(iev,ird),& ! a8
       delt_pr,& ! f6.2
       azes_pr,& ! i3
       phase(iev,ird),& ! a8
       ipay(iev,ird),& ! i4
       ipamo(iev,ird),& ! i2
       ipad(iev,ird),& ! i2
       ipah(iev,ird),& ! i2
       ipam(iev,ird),& ! i2
       pas(iev,ird),& ! f6.3
       iptim(iev,ird),& ! i2
       dts_pr,& ! a6
       sdread(iev,ird) ! f5.2
   end do
   
   ! Stop line
   write (io_unit,'(a,t121,a)') 'STOP', ' '
   
   return
   
end subroutine write_mnf_14


!*****************************************************************************************
subroutine read_mnf_15 ()

! Read a data file of differential time data in MNF v1.5 format. The input datum is
! "reduced relative arrival time", the time difference between the two arrivals if they are
! treated as occurring on the same day. This is converted to dummy arrival times for the
! associated pair of events. For the template event the theoretical TT is added to the
! origin time of the input data file. For the target event, the reduced relative arrival
! time is added to the template event's dummary arrival time. Dummy arrival times derived
! from differential time data are only used for cluster vectors, never the hypocentroid.

! Phase names for differential time data are not cleaned and phase re-identification is
! disabled. There is no testing for duplicate readings. Differential time data for Lg are
! permitted, but otherwise only phases that are in the tau-p phase list are processed. If
! an event referenced in a differential time datum is not found in the cluster, a warning
! is given.

   implicit none
   
   include 'mloc.inc'

   integer, parameter :: max = 60
   
   integer :: nphase
   real :: tt
   character(len=8) :: phcd
   common /taup/ nphase, tt(max), phcd(max)
   
   character(len=150) :: linein
   character(len=132) :: msg
   character(len=16) :: evtnam_template, evtnam_target
   character(len=10) :: evid_template, evid_target
   character(len=8) :: pha0, deployment0
   character(len=5) :: sta0, agency0
   character(len=1) :: fcode0
   real :: rrat, sig_diff, sig_tt, hms2s
   real :: xlat, sxlat, cxlat, xlon, sxlon, cxlon
   real :: t1, t2, t3, t4, dum1, dum2, dum3, dum4
   real :: dtdd(max), dtdh(max), dddp(max), usrc(2), time_template, time_target, ttime1
   integer :: iev_template, iev_target, ird_template, ird_target
   integer :: nlost, nread, n_line
   integer :: nx(nevmax), nf(nevmax), nfsum, iev, ifrac, i , ios, iptim_base, iptim_spec
   logical :: skipp
   
   nlost = 0 ! Number of data lost because of a failure to associate with cluster events
   ndiff = 0 ! Number of differential time data added
   nread = 0 ! Number of differential times read
   nx = 0 ! Number of differential time data for each event
   nf = 0 ! Number of differential time data that are flagged for each event

   write (io_log,'(a)') 'Dummy times for differential time data'
   
   write (msg,'(2a)') 'read_mnf_15: reading differential times from ', trim(diffdatfilnam)
   call fyi (trim(msg))
      
   n_line = 1 ! The format record has already been read in mnf_read
   
   do
      read (io_diffdat,'(a)',iostat=ios) linein
      
      ! Unexpected termination
      if (ios .lt. 0)  call oops ('read_mnf_15: EOF reached unexpectedly')
      
      ! Normal termination
      if (linein(1:3) .eq. 'EOF') exit
      
      n_line = n_line + 1
      
      if (linein(1:1) .eq. 'D') then ! Differential time record
      
         nread = nread + 1
      
         evtnam_template = linein(5:20)
         evid_template = linein(22:31)
         evtnam_target = linein(33:48)
         evid_target = linein(50:59)
         
         ! Associate with cluster events
         ! Try first with evid then use event designator based on OT
         iev_template = 0
         iev_target = 0
         !Template event
         if (evid_template .ne. '          ') then
            do i = 1,nev
               if (evid_template .eq. evid(i)) then
                  iev_template = i
                  exit
               end if
            end do
         end if
         if (iev_template .eq. 0) then
            do i = 1,nev
               if (evtnam_template .eq. evtnam(i)(1:16)) then
                  iev_template = i
                  exit
               end if
            end do
         end if
         ! Target event
         if (evid_target .ne. '          ') then
            do i = 1,nev
               if (evid_target .eq. evid(i)) then
                  iev_target = i
                  exit
               end if
            end do
         end if
         if (iev_target .eq. 0) then
            do i = 1,nev
               if (evtnam_target .eq. evtnam(i)(1:16)) then
                  iev_target = i
                  exit
               end if
            end do
         end if
         
         if (iev_template .eq. 0 .or. iev_target .eq. 0) then
            if (iev_template .eq. 0) then
               write (msg,'(3a)') 'read_mnf_15: template event', evtnam_template, ' not found in the cluster'
               call warnings (trim(msg))
            end if
            if (iev_target .eq. 0) then
               write (msg,'(3a)') 'read_mnf_15: target event', evtnam_target, ' not found in the cluster'
               call warnings (trim(msg))
            end if
            nlost = nlost + 1
            cycle
         else
            nx(iev_template) = nx(iev_template) + 1
            nx(iev_target) = nx(iev_target) + 1         
         end if
      
         ! Read station, phase, agency, deployment
         sta0 = linein(61:65)
         pha0 = linein(68:75)
         agency0 = linein(114:118)
         deployment0 = linein(120:127)

         ! Increment counter for differential time data.
         ndiff = ndiff + 1
         if (ndiff .gt. ndiffmax) then
            write (msg,'(a,i4,a)') 'read_mnf_15: ndiff exceeds maximum (', ndiffmax, '). Remaining data skipped'
            call warnings (trim(msg))
            ndiff = ndiff - 1
            exit
         end if

         ! Increment counter for number of readings for each event, and set idiff and diff_line
         ! idiff is the index of differential time readings
         ! diff_line is the line number of the data file on which the reading occurs
         ird_template = nst(iev_template) + 1
         nst(iev_template) = ird_template
         idiff(iev_template,ird_template) = ndiff
         diff_line(iev_template,ird_template) = n_line
         ird_target = nst(iev_target) + 1
         nst(iev_target) = ird_target
         idiff(iev_target,ird_target) = ndiff
         diff_line(iev_target,ird_target) = n_line
         
         ! Read station, agency, deployment
         sta0 = linein(61:65)
         agency0 = linein(114:118)
         deployment0 = linein(120:127)
         ! Template event
         stname(iev_template,ird_template) = sta0
         agency(iev_template,ird_template) = agency0
         deployment(iev_template,ird_template) = deployment0
         sad(iev_template,ird_template) = sta0//' '//agency0//' '//deployment0
         ! Target event
         stname(iev_target,ird_target) = sta0
         agency(iev_target,ird_target) = agency0
         deployment(iev_target,ird_target) = deployment0
         sad(iev_target,ird_target) = sta0//' '//agency0//' '//deployment0

         ! Station coordinates
         call stafind (iev_template, ird_template) !  Search master list of stations for station coordinates.
         call stafind (iev_target, ird_target) !  Search master list of stations for station coordinates.
         
         ! Phase name (no cleaning of phase names for differential time readings)
         pha0 = linein(68:75)
         phase0(iev_template,ird_template) = pha0
         phase(iev_template,ird_template) = pha0
         phase0(iev_target,ird_target) = pha0
         phase(iev_target,ird_target) = pha0
         
         ! Phase names for differential times cannot be re-identified
         phidird(iev_template,ird_template) = .false.
         phidird(iev_target,ird_target) = .false.
      
         ! Usage flag
         fcode0 = linein(3:3)
         fcode(iev_template,ird_template) = fcode0
         fcode(iev_target,ird_target) = fcode0
         ! Some phases must be skipped
         if (skipp(phase(iev_template,ird_template))) then
            if (verbose_log) write (io_log,'(3a)') 'read_mnf_15: phase skipped - ',&
             stname(iev_template,ird_template), phase(iev_template,ird_template)
            fcode(iev_template,ird_template) = 'p'
            fcode(iev_target,ird_target) = 'p'
         end if
         
         ! Update the number of flagged readings
         if (fcode0 .ne. ' ') then
            nf(iev_template) = nf(iev_template) + 1
            nf(iev_target) = nf(iev_target) + 1
         end if
      
         ! Set data flags
         call dataflags (dflag, fcode(iev_template,ird_template))
         call dataflags (dflag, fcode(iev_target,ird_target))

         ! Prior residual, not included in v1.5 format
         resisc(iev_template,ird_template) = -99.9
         resisc(iev_target,ird_target) = -99.9
      
         ! Reading source (author)
         if (len(linein) .ge. 149) then
            readsrc(iev_template,ird_template) = linein(142:149)
            readsrc(iev_target,ird_target) = linein(142:149)
         else
            readsrc(iev_template,ird_template) = 'Unknown '
            readsrc(iev_target,ird_target) = 'Unknown '
         end if

         ! These records are not in the main event data file, so they do not get an mnf_line number
         mnf_line(iev_template,ird_template) = 0
         mnf_line(iev_target,ird_target) = 0
      
         adslcaa(iev_template,ird_template)(1:47) = linein(114:149)
         adslcaa(iev_target,ird_target)(1:47) = linein(114:149)
         
         channel(iev_template,ird_template) = linein(138:140)
         channel(iev_target,ird_target) = linein(138:140)
         
         ! Initial setting of fltrh. Needed for calculation of open azimuth in model_init.
         ! Dummy readings for differential times are never used for the hypocentroid.
         if (fcode(iev_template,ird_template) .eq. ' ') then
            fltrhflag(iev_template,ird_template) = .false.
         else
            fltrhflag(iev_template,ird_template) = .true.
         end if
         fltrh(iev_template,ird_template) = .true.
         if (fcode(iev_target,ird_target) .eq. ' ') then
            fltrhflag(iev_target,ird_target) = .false.
         else
            fltrhflag(iev_target,ird_target) = .true.
         end if
         fltrh(iev_target,ird_target) = .true.
         
         ! Reduced relative arrival time
         ! See documentation for MNF v1.5 for definition
         read (linein(77:87),'(f11.4)') rrat

         ! Precision
         if (linein(87:87) .ne. ' ') then
            iptim_base = -4
         else if (linein(86:86) .ne. ' ') then
            iptim_base = -3
         else if (linein(85:85) .ne. ' ') then
            iptim_base = -2
         else if (linein(84:84) .ne. ' ') then
            iptim_base = -1
         else
            iptim_base = 0
         end if
         if (linein(90:90) .eq. ' ') then
            iptim(iev_template,ird_template) = iptim_base
            iptim(iev_target,ird_target) = iptim(iev_template,ird_template)
         else
            read (linein(89:90),'(i2)') iptim_spec
            iptim(iev_template,ird_template) = max0(iptim_base,iptim_spec)
            iptim(iev_target,ird_target) = iptim(iev_template,ird_template)
         end if

         ! Dummy arrival times
      
         call depset (depthp(iev_template,0), usrc)
      
         ! Convert event lat-lon to geocentric coordinates
         call geocen (latp(iev_template,0), lonp(iev_template,0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)      
   
         ! Epicentral distance, azimuth, and back-azimuth
         t1 = xlat*rpd ! Convert to geocentric radians
         t2 = xlon*rpd ! Convert to geocentric radians
         t3 = stladg(iev_template,ird_template)*rpd ! Convert to geocentric radians
         t4 = stlndg(iev_template,ird_template)*rpd ! Convert to geocentric radians
         call delaz (t1, t2, t3, t4, dum1, delt(iev_template,ird_template), dum2, dum3,&
          azes(iev_template,ird_template), dum4, azse(iev_template,ird_template), 1)
         if (delt(iev_template,ird_template) .lt. 0. .or. delt(iev_template,ird_template) .gt. 180.) then
            write (msg,'(a,f10.3)') 'read_mnf_15: illegal value for delta: ', delt(iev_template,ird_template)
            call oops (msg)
         end if

         ! Theoretical travel time
         call trtm (delt(iev_template,ird_template), max, nphase, tt, dtdd, dtdh, dddp, phcd)
         psd(iev_template,ird_template) = 0.
         ttime1 = 0.
         do i = 1,nphase
            if (phcd(i) .eq. phase0(iev_template,ird_template)) then
               ttime1 = tt(i)
               psd(iev_template,ird_template) = dtdd(i)
               exit
            end if
         end do
         
         ! Lg travel time
         if (phase0(iev_template,ird_template)(1:2) .eq. 'Lg') then
            psd(iev_template,ird_template) = lg_b
            ttime1 = lg_a + delt(iev_template,ird_template)*lg_b
         end if
         
         ! Template event dummy arrival time
         time_template = hms2s(hourp(iev_template,0), minp(iev_template,0), secp(iev_template,0)) + ttime1
         call timecr (time_template, ipah(iev_template,ird_template), ipam(iev_template,ird_template),&
          pas(iev_template,ird_template))
      
         ! Target event dummy arrival time
         time_target = time_template + rrat
         call timecr (time_target, ipah(iev_target,ird_target), ipam(iev_target,ird_target), pas(iev_target,ird_target))

         ! Reading error for dummy arrival times
         ! Includes uncertainty of the differetial time measurement plus cluster vector TT error
         if (len(linein) .ge. 97) then
            read (linein(92:97),'(f6.4)') sig_diff
         else
            sig_diff = 0.1 ! default
         end if
         sig_tt = (psd(iev_template,ird_template)*vmr*cscale)/111.2 ! Uncertainty of theoretical TTs for differential time
         sdread(iev_template,ird_template) = sqrt(sig_diff*sig_diff + sig_tt*sig_tt)
         if (sdread(iev_template,ird_template) .lt. 0.05) sdread(iev_template,ird_template) = 0.05 ! Minimum allowed value
         sdread(iev_target,ird_target) = sdread(iev_template,ird_template)

         ! Log
         write (io_log,'(4(a,1x),2(i2,1x,i2,1x,f6.3,1x),f11.4,3f8.3)') evtnam_template, evtnam_target, sta0, pha0,&
          ipah(iev_template,ird_template), ipam(iev_template,ird_template), pas(iev_template,ird_template),&
          ipah(iev_target,ird_target), ipam(iev_target,ird_target), pas(iev_target,ird_target), rrat, sig_diff,&
          sig_tt, sdread(iev_template,ird_template)
         
      end if
      
   end do
   
   nfsum = 0
   do iev = 1,nev
      if (nx(iev) .gt. 0) then
         nfsum = nfsum + nf(iev)
         ifrac = nint(nf(iev)*1.e2/nx(iev))
         write (msg,'(a,i4,a,i4,a,i3,a,i3)') 'read_mnf_15: ', nf(iev), ' of ', nx(iev), ' (',&
          ifrac, '%) differential times flagged for event ', iev
         call fyi (trim(msg))
      end if
   end do

   nfsum = nfsum/2 ! Correction for counting each differential time datum twice, once for each associated event
   ifrac = nint(nfsum*1.e2/nread)
   write (msg,'(a,i4,a,i4,a,i3,a)') 'read_mnf_15: ', nread, ' differential times read; ', nfsum, ' flagged (', ifrac, '%)'
   call fyi (trim(msg))
   write (msg,'(a,i4,a)') 'read_mnf_15: ', ndiff, ' differential times added'
   call fyi (trim(msg))
   if (nlost .gt. 0) then
      write (msg,'(a,i4,a)') 'read_mnf_15: ', nlost, ' unassociated differential time readings skipped'
      call fyi (trim(msg))
   end if
   
   return
   
end subroutine read_mnf_15


!*****************************************************************************************
subroutine hdf_match (iev)

! Match an event to a list read from an HDF file

   implicit none
   
   include 'mloc.inc'
   
   logical :: found_hdf
   real :: time0, hms2s, adt, adt_min
   integer :: iev, j, jf
   character(len=132) :: msg
   
   found_hdf = .false.
   time0 = hms2s(hour_inp(iev), min_inp(iev), sec_inp(iev)) ! OT from input file
   
   ! Try a match based on evids first
   if (evid(iev) .ne. '          ') then
      do j = 1,nrhdf
         if (evid(iev) .eq. hdf_evid(j)) then
            found_hdf = .true.
            jf = j
            if (verbose_screen) then
               write (msg,'(a,i3,a,i3)') 'hdf_match: event ', iev, ' matched in HDF file by evid at line ', jf
               call fyi (trim(msg))
            end if
            exit
         end if
      end do
   end if
   if (.not.found_hdf .and. verbose_screen) then
      write (msg,'(a,i3,a)') 'hdf_match: event ', iev, ' not matched by evid'
      call fyi (trim(msg))
   end if
   
   ! If that didn't work, try a match based on OT
   ! OT in the hypocenter record must be within 10 s of OT in the HDF file.
   ! If more than one event within the 10 s window, the one with smallest time difference is selected.
   if (.not.found_hdf) then
      adt = -999.
      adt_min = 10.
      jf = 0
      do j = 1,nrhdf
         if (iyre(iev) .eq. hdf_yr4(j) .and.&
             mone(iev) .eq. hdf_mon(j) .and.&
             idye(iev) .eq. hdf_day(j)) then
            adt = abs(hdf_time(j)-time0)
            if (adt .lt. adt_min) then
               jf = j
               adt_min = adt
            end if
         end if
      end do
      if (jf .gt. 0) then
         found_hdf = .true.
         if (verbose_screen) then
            write (msg,'(a,i3,a,i3)') 'hdf_match: event ', iev, ' matched in HDF file at line ', jf
            call fyi (trim(msg))
         end if
      end if
   end if
   if (.not.found_hdf) then
      write (msg,'(a,i3,a,i4,4i3,f7.3,1x,a)') 'hdf_match: event ', iev, ' not found in HDF file ', iyre(iev),&
       mone(iev), idye(iev), hour_inp(iev), min_inp(iev), sec_inp(iev), evid(iev)
      call warnings (trim(msg))
   end if
   
   ! Set starting coordinates from values in the HDF file
   if (found_hdf) then
      lat_hdf(iev) = hdf_lat(jf)
      lon_hdf(iev) = hdf_lon(jf)
      depth_hdf(iev,1) = hdf_dep(jf)
      depth_hdf(iev,2) = hdf_dep_plus(jf)
      depth_hdf(iev,3) = hdf_dep_minus(jf)
      if (hdf_dep_code(jf)(2:2) .eq. 'f') then ! free depth solution
         depth_hdf_c(iev) = 'm'
      else
         depth_hdf_c(iev) = hdf_dep_code(jf)(1:1)
      end if
            
      time_hdf(iev) = hdf_time(jf)
      
      if (debug) then
         write (msg,'(a,i3,1x,2f10.4,f5.1,1x,a1,1x,2f5.1,f10.3)') 'hdf_match: ', iev, lat_hdf(iev), lon_hdf(iev),&
          depth_hdf(iev,1), depth_hdf_c(iev), depth_hdf(iev,2), depth_hdf(iev,3), time_hdf(iev)
         call debugger (trim(msg))
      end if
      
   end if
      
   return
   
end subroutine hdf_match


!*****************************************************************************************
subroutine starting_hypocenter (iev)

! Set starting values for the hypocentral parameters.
! For latitude, longitude, and OT, starting locations are taken from the following sources, with
! the order of precedence:
!   1) LAT, LONG or TIME command in the command file
!   2) HDF file, if the event has been found there
!   3) Input file
! For focal depth, there is a complication caused by the frequent use of the DEPC command
! to set a cluster default depth. In cases where the focal depth in the input file (either
! in the preferred hypocenter or depth record) has a flag indicating it has been
! constrained, we do not allow the DEPC command to over-ride that value. It can be over-ridden
! by a DEP_ command in the definition section (after a MEMB command) for that event in the
! command file. 

! Calibration coordinates (calX command) are never taken as starting values.

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev
   real :: time0
!   real :: hms2s
   character(len=132) :: msg
   logical :: depth_constraint
   
   ! Starting latitude
   if (lat_cf(iev) .ge. -90.) then ! Command file
      latp(iev,0) = lat_cf(iev) + hlatshift
   else if (lat_hdf(iev) .ge. -90.) then ! HDF file
      latp(iev,0) = lat_hdf(iev) + hlatshift
   else ! Input file
      latp(iev,0) = lat_inp(iev) + hlatshift
   end if
   
   ! Starting longitude
   if (lon_cf(iev) .ge. -180.) then ! Command file
      lonp(iev,0) = lon_cf(iev) + hlonshift
   else if (lon_hdf(iev) .ge. -180.) then ! HDF file
      lonp(iev,0) = lon_hdf(iev) + hlonshift
   else ! Input file
      lonp(iev,0) = lon_inp(iev) + hlonshift
   end if
   
   ! Starting depth, needed for call to depset and processing phase readings
   if (debug)  write (io_log,'(a)') 'Starting depths'
   if (depth_cf(iev,1) .ge. 0.) then ! Depth settings taken from command file
      depthp(iev,0) = depth_cf(iev,1)
      depset_pr(iev) = depth_cf_c(iev)
      depthp_plus(iev) = depth_cf(iev,2)
      depthp_minus(iev) = depth_cf(iev,3)
      if (debug) then
         write (io_log,'(i3,a,t14,f5.1,1x,a1,2f6.1)') iev, ' cfil:', depthp(iev,0), depset_pr(iev), depthp_plus(iev),&
          depthp_minus(iev)
      end if
   else if (depth_hdf(iev,1) .ge. 0. .and. depth_constraint(depth_hdf_c(iev))) then ! Depth settings taken from a .hdf file if they are constrained
      depthp(iev,0) = depth_hdf(iev,1)
      depset_pr(iev) = depth_hdf_c(iev)
      depthp_plus(iev) = depth_hdf(iev,2)
      depthp_minus(iev) = depth_hdf(iev,3)
      if (debug) then
         write (io_log,'(i3,a,t14,f5.1,1x,a1,2f6.1)') iev, ' hdf:', depthp(iev,0), depset_pr(iev), depthp_plus(iev),&
          depthp_minus(iev)
      end if
   else if (depth_inp(iev,1) .ge. 0. .and. depth_constraint(depth_inp_c(iev))) then ! Depth settings taken from the input data file if depth is constrained
      depthp(iev,0) = depth_inp(iev,1)
      depset_pr(iev) = depth_inp_c(iev)
      depthp_plus(iev) = depth_inp(iev,2)
      depthp_minus(iev) = depth_inp(iev,3)
      if (debug) then
         write (io_log,'(i3,a,t14,f5.1,1x,a1,2f6.1)') iev, ' input:', depthp(iev,0), depset_pr(iev), depthp_plus(iev),&
          depthp_minus(iev)
      end if
   else if (depth_default .ge. 0.) then ! Depth settings taken from the cluster default values
      depthp(iev,0) = depth_default
      depset_pr(iev) = depth_default_c
      depthp_plus(iev) = 99.9
      depthp_minus(iev) = 99.9
      if (debug) then
         write (io_log,'(i3,a,t14,f5.1,1x,a1,2f6.1)') iev, ' cluster:', depthp(iev,0), depset_pr(iev), depthp_plus(iev),&
          depthp_minus(iev)
      end if
   else if (depth_inp(iev,1) .ge. 0.) then ! Depth settings taken from the input data file
      depthp(iev,0) = depth_inp(iev,1)
      depset_pr(iev) = depth_inp_c(iev)
      depthp_plus(iev) = depth_inp(iev,2)
      depthp_minus(iev) = depth_inp(iev,3)
      if (debug) then
         write (io_log,'(i3,a,t14,f5.1,1x,a1,2f6.1)') iev, ' input:', depthp(iev,0), depset_pr(iev), depthp_plus(iev),&
          depthp_minus(iev)
      end if
   else
      write (msg,'(a,i3,a)') 'starting_hypocenter: no depth has been set for event', iev, evtnam(iev)
      call oops (trim(msg))
   end if
   
   ! Starting origin time
   if (time_cf(iev) .ge. 0.) then ! Command file
      time0 = time_cf(iev) + htimeshift
   else if (time_hdf(iev) .ge. 0.) then ! HDF file
      time0 = time_hdf(iev) + htimeshift
   else ! Input file
      time0 = time_inp(iev) + htimeshift
   end if
   call timecr (time0, hourp(iev,0), minp(iev,0), secp(iev,0))
   
   if (debug) write (io_log,'(2f10.4,f6.1,2i3,f6.2)') latp(iev,0), lonp(iev,0), depthp(iev,0), hourp(iev,0),&
    minp(iev,0), secp(iev,0)
    
   return
   
end subroutine starting_hypocenter


!*****************************************************************************************
subroutine phase_utility (iev, ird, psta, p_arrtime)

! Several steps applied on the basis of phase name:
!   Setting a flag to determine if phase re-identification can be done
!   Processing of S-P relative phases
!   Clean-up of phase names
!   Processing for relative depth phases (pP-P and sP-P)
!   Processing PKPdf precursors
!   Flagging certain phases that are always skipped in mloc

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, ird, i
   real :: p_arrtime, dp_arrtime, hms2s
   logical :: skipp
   character(len=5) :: psta
      
   ! Phase re-identification flag
   phidird(iev,ird) = phid
   if (no_phid(iev,ird) .eq. '!') phidird(iev,ird) = .false. ! This phase ID cannot be changed
   if (fcode(iev,ird) .eq. 'm') phidird(iev,ird) = .false. ! Don't re-identify the phase because we don't know where the station is
   if (fcode(iev,ird) .eq. 't') phidird(iev,ird) = .false. ! Don't re-identify the phase because we don't trust the timing
   if (debug) write (io_log,'(a,1x,l1)') 'phase re-identification ', phidird(iev,ird)
   
   ! S-P relative phase
   if (phase0(iev,ird)(1:7) .eq. 'S-P    ') then
      rel_phase(iev,ird) = .true.
      phidird(iev,ird) = .false. ! Don't re-identify this phase
   else
      rel_phase(iev,ird) = .false.
   end if
   if (debug) write (io_log,'(a,1x,l1)') 'relative phase ', rel_phase(iev,ird)
   if (debug) write (io_log,'(a,1x,l1)') 'phase re-identification ', phidird(iev,ird)
   
   ! Clean up phase name (original saved in phase0)
   call pnclean (phase0(iev,ird), phase(iev,ird))
   if (verbose_log) then
      if (phase(iev,ird) .ne. phase0(iev,ird)) then
         write (io_log,'(a,i3,1x,a)') 'pnclean changed '//phase0(iev,ird)//' to '//phase(iev,ird)//' for event ',&
          iev, stname(iev,ird)
      end if
   end if
   
   ! Relative depth phases
   ! pP and sP phases are processed to determine pP-P and sP-P times, using the most recent P arrival from the same station
   if (phase0(iev,ird) .eq. 'pP      ' .or. phase0(iev,ird) .eq. 'sP      ' .or. phase0(iev,ird) .eq. 'pwP     ') then
      if (stname(iev,ird) .eq. psta) then
         dp_arrtime = hms2s(ipah(iev,ird), ipam(iev,ird), pas(iev,ird))
         if (ipah(iev,ird) .lt. hourp(iev,0)) dp_arrtime = dp_arrtime + 86400.
         rel_depth_phase(iev,ird) = dp_arrtime - p_arrtime
         if (debug) write (io_log,'(a,1x,a,f10.3,2i3,f7.3,1x,l1,1x,l1)') stname(iev,ird), phase0(iev,ird),&
          dp_arrtime, ipah(iev,ird), ipam(iev,ird), pas(iev,ird), rel_phase(iev,ird), phidird(iev,ird)
      end if
   end if
   
   ! List of phases that cannot be renamed (command PPRI)
   if (n_no_phreid .gt. 0) then
      do i = 1,n_no_phreid
         if (trim(no_phreid(i)) .eq. trim(phase0(iev,ird))) phidird(iev,ird) = .false. ! This phase ID cannot be changed.
      end do
   end if
   
   ! PKPdf precursors are not used, but their residual will be calculated relative to the PKPdf phase.
   ! They cannot be re-identified.
   if (phase0(iev,ird) .eq. 'PKPdfpre' .or. phase0(iev,ird) .eq. 'PKPpre  ') then
      fcode(iev,ird) = 'p'
      phidird(iev,ird) = .false.
   end if
   
   ! Some phases are always skipped
   if (skipp(phase(iev,ird))) then
      if (verbose_log) write (io_log,'(3a)') 'phase_utility: phase skipped - ', stname(iev,ird), phase(iev,ird)
      fcode(iev,ird) = 'p'
      phidird(iev,ird) = .false.
   end if
      
   return
   
end subroutine phase_utility


!*****************************************************************************************
subroutine duplicates (iev,ird)

! Check for duplicate readings (same station and same time, within 0.1 s).
! This does not apply to S-P readings, for which apparent duplicates are really distinct samples.
! Issue a warning if a dupe is not already flagged as such.

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, ird, i
   real :: t1, t2, tdiff
   
   if (ird .ge. 2 .and. fcode(iev,ird) .eq. ' ') then
      do i = 1,ird-1
         if (stname(iev,ird) .eq. stname(iev,i) .and. phase0(iev,ird) .ne. 'S-P     ') then
            t1 = real(ipah(iev,i))*3600. + float (ipam(iev,i))*60. + pas(iev,i)
            t2 = real(ipah(iev,ird))*3600. + float (ipam(iev,ird))*60. + pas(iev,ird)
            tdiff = t2 - t1
            if (abs(tdiff) .le. 0.1) then
               if (verbose_log) write (io_log,'(i3,i5,1x,a6,1x,a8,2i3,2f5.1,a)') iev, mnf_line(iev,ird),&
                stname(iev,ird), phase0(iev,ird), ipah(iev,ird), ipam(iev,ird), pas(iev,ird), tdiff,&
                ' is a duplicate reading'
               fcode(iev,ird) = 'd'
            end if
         end if
      end do
   end if
   
   if (debug) write (io_log,'(2a)') 'After checking for duplicates, fcode = ', fcode(iev,ird)
   
   return
   
end subroutine duplicates


!*****************************************************************************************
subroutine raypath (iev, ird)

! Calculate various distance and geometric parameters

   implicit none
   
   include 'mloc.inc'
   
   integer, parameter :: max = 60
   integer :: nphase
   real :: tt
   character(len=8) :: phcd
   common /taup/ nphase, tt(max), phcd(max)
   
   integer :: i, iev, ird
   real :: xlat, sxlat, cxlat, xlon, sxlon, cxlon, t1, t2, t3, t4, dum1, dum2, dum3, dum4
   real :: dtdd(max), dtdh(max), dddp(max)
   character(len=132) :: msg
      
   ! Convert event lat-lon to geocentric coordinates
   call geocen (latp(iev,0), lonp(iev,0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
   if (debug) write (io_log,'(a,2f10.4)') 'Geocentric coordinates: ', xlat, xlon      
   
   ! Epicentral distance, azimuth, and back-azimuth
   t1 = xlat*rpd ! Convert to geocentric radians
   t2 = xlon*rpd ! Convert to geocentric radians
   t3 = stladg(iev,ird)*rpd ! Convert to geocentric radians
   t4 = stlndg(iev,ird)*rpd ! Convert to geocentric radians
   call delaz (t1, t2, t3, t4, dum1, delt(iev,ird), dum2, dum3, azes(iev,ird), dum4, azse(iev,ird), 1)
   if (delt(iev,ird) .lt. 0. .or. delt(iev,ird) .gt. 180.) then
      write (msg,'(a,f10.3)') 'raypath: illegal value for delta: ', delt(iev,ird)
      call oops (trim(msg))
   end if
   if (debug) write (io_log,'(a,2f10.4)') 'Distance and azimuth: ', delt(iev,ird), azes(iev,ird)
   
   ! Ray parameter, needed for calculation of surface focus distance
   call trtm (delt(iev,ird), max, nphase, tt, dtdd, dtdh, dddp, phcd)
   psd(iev,ird) = 0.
   do i = 1,nphase
      if (phcd(i) .eq. phase0(iev,ird)) then
         psd(iev,ird) = dtdd(i)
         exit
      end if
   end do
   if (debug) write (io_log,'(a,f10.4)') 'Ray parameter: ', psd(iev,ird)
   
   return
   
end subroutine raypath


!*****************************************************************************************
logical function depth_constraint (c)

! Based on the depth code, returns .true. if the associated focal depth is considered constrained.
! Returns .false. if not.

   implicit none
   
   character(len=1) :: c
      
   if (c .eq. 'c') then ! cluster default depth 
      depth_constraint = .false.
   else if (c .eq. 'd') then ! depth phases
      depth_constraint = .true.
   else if (c .eq. 'e') then ! engineered (man-made explosion)
      depth_constraint = .true.
   else if (c .eq. 'f') then ! fault model (InSAR, GPS, etc.)
      depth_constraint = .true.
   else if (c .eq. 'i') then ! input data file
      depth_constraint = .false.
   else if (c .eq. 'l') then ! local distance readings (more than 2-3 focal depths)
      depth_constraint = .true.
   else if (c .eq. 'm') then ! mloc solution (with free depth)
      depth_constraint = .true.
   else if (c .eq. 'n') then ! near-source station readings
      depth_constraint = .true.
   else if (c .eq. 'r') then ! relocation (outside mloc) with free depth
      depth_constraint = .true.
   else if (c .eq. 'u') then ! unknown
      depth_constraint = .false.
   else if (c .eq. 'w') then ! waveform analysis
      depth_constraint = .true.
   else if (c .eq. ' ') then ! blank
      depth_constraint = .false.
   else
      call warnings ('depth_constraint: unknown depth code "'//trim(c)//'"')
      depth_constraint = .false.
   end if
                      
   return
   
end function depth_constraint

   
   
   
