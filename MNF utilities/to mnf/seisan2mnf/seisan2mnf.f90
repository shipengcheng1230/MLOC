PROGRAM seisan2mnf
      
! Reads a file containing arrival time data for one or more events in SEISAN format as produced by Reza Ghods
! The output is in MNF (mloc native format) format.

! There are two modes for output:
!    Bulletin - All events concatenated in a single file whose name is based on the input file
!    Event - A separate file is created for each event, with the filename based on date and OT (yyyymmdd.hhmm.ss.mnf)
! If bulletin mode is applied to a file containing a single event, the only difference in the output file, other than
! filename, is the presence of a bulletin record at the top of the file.

! Files written by this program are compatible with MNF versions 1.3 through 1.3.3.

   implicit none
   
   include 'seisan2mnf.inc'
   
   character*80 filename, basename
   character*136 linein
   character*120 comment
   character(len=6) :: version
   character*5 c1, c2, code_loc(100), code_int(100) 
   character*1 chinp, m_type
   logical loop, event_on, loop_phase, bulletin, ex, cfil, code_convert
   integer ios, nev, ncode, i
   real distance, residual,diff
 
   version = '1.3.3' ! MNF version
   nev = 0
   event_on = .false.
   cfil = .false.
   comment = ' '

   write (*,'(a/)') 'Release date January 10, 2018, writing MNF v1.3.3 (backward compatible to v1.3)'
   
   write (*,'(a)') 'Enter input filename: '
   read (*,'(a)') filename
   open (1, file=trim(filename), status='old')
   
   ! Bulletin record
   write (*,'(a)') 'Bulletin output?: '
   read (*,'(a)') chinp
   bulletin = chinp .eq. 'y' .or. chinp .eq. 'Y'
   if (bulletin) then ! Output file named after input file
      open (2, file=trim(filename)//'.mnf', status='new')
      write (*,'(a)') 'Bulletin comment: '
      read (*,'(a)') comment
      write (2,'(2a,t121,a)') 'B   ', trim(comment), ' '
      write (2,'(2a,t121,a)') 'F   MNF v', trim(version), ' '
   else
      write (*,'(a)') 'Create mloc command file? '
      read (*,'(a)') chinp
      cfil = chinp .eq. 'y' .or. chinp .eq. 'Y'
      if (cfil) then
         write (*,'(a)') 'Enter command file basename: '
         read (*,'(a)') basename
         open (3, file=trim(basename)//'.cfil', status='new')
      end if
   end if
   
   agency = ' '
   deployment = ' '
   write (*,'(a)') 'Agency: '
   read (*,'(a)') agency
   if (trim(agency) .eq. 'IRSC') then
      deployment = 'ITSN    '
   else if (trim(agency) .eq. 'IIEES') then
      deployment = 'INSN    '
   end if
   
   phase_author = ' '
   write (*,'(a)') 'Phase pick author: '
   read (*,'(a)') phase_author
 
   write (*,'(a)') 'Correct local station codes to internationally registered codes?'
   read (*,'(a)') chinp
   if (chinp .eq. 'y') then
      code_convert = .true.
      write (*,'(a)') 'Enter file name: '
      read (*,'(a)') filename
      open (4,file=trim(filename),status='old')
      read (4,'(a)') linein ! comment line
      write (*,'(a)') trim(linein)
      ncode = 0
      do 
         read (4,'(a5,t10,a5)',iostat=ios) c1, c2
         if (ios .lt. 0) then
            close (4)
            write (*,'(i3,a)') ncode, ' station code translations read'
            exit
         else
            ncode = ncode + 1
            if (ncode .le. 100) then
               code_loc(ncode) = c1
               code_int(ncode) = c2
               write (*,'(i3,2x,a5,2x,a5)') ncode, code_loc(ncode), code_int(ncode)
            else
               write (*,'(//a//)') 'Maximum number of station code translations read!'
            end if
         end if
      end do
   end if
   
   loop = .true.
   do while (loop) ! loop over all lines
   
      read (1,'(a)',iostat=ios) linein
      
      ! EOF
      if (ios .lt. 0) then
         close (1)
         if (bulletin) then
            write (2,'(a,t121,a)') 'EOF', ' '
            close (2)
         end if
         if (cfil) close (3)
         write (*,'(a,i6,a)') 'EOF reached after ', nev, ' events'
         stop
      end if
      
      if (linein(80:80) .eq. 'I') cycle
      if (linein(80:80) .eq. 'E') cycle
      if (linein(80:80) .eq. '2') cycle
      if (linein(80:80) .eq. '3') cycle
      if (linein(80:80) .eq. '6') cycle
      
      if (linein(80:80) .eq. '1' .and. .not. event_on) then ! New event and hypocenter lines
         event_on = .true.
         nev = nev + 1
         evid = ' '
         event_geog = ' '
         event_usage = ' '
         read (linein(1:10),'(1x,i4,1x,2i2)') hypo_year, hypo_month, hypo_day
         read (linein(12:20),'(2i2,1x,f4.1)') hypo_hour, hypo_minute, hypo_seconds
         
         ! Open individual event files if not bulletin output
         if (bulletin) then
            call event_line (2)
         else
            write (filename,'(i4,2i2.2,a,2i2.2,a,i2.2)') hypo_year, hypo_month, hypo_day, '.',&
             hypo_hour, hypo_minute, '.', int(hypo_seconds)
            inquire (file=trim(filename)//'.mnf',exist=ex)
            if (ex) then
               write (*,'(a)') 'file '//trim(filename)//'.mnf already exists'
               close (1)
               stop
            end if
            open (2, file=trim(filename)//'.mnf', status='new')
            write (2,'(2a,t121,a)') 'F   MNF v', trim(version), ' '
            call event_line (2)
            if (cfil) then
               write (3,'(a)') 'memb'
               write (3,'(a)') 'even '//trim(filename)
               write (3,'(a)') 'inpu '//trim(filename)//'.mnf'
            end if
         end if
         
         ! Hypocenter line
         read (linein(24:30),'(f7.3)') latitude
         read (linein(31:38),'(f8.3)') longitude
         read (linein(39:43),'(f5.1)') depth
         hypo_author = linein(46:48)//'     '
         hypo_usage = '='
         hypo_orid = ' '
         call hypocenter_line (2)
                  
         ! Magnitude line
         if (linein(58:58) .eq. '.') then
            read (linein(57:59),'(f3.1)') magnitude
            read (linein(60:60),'(a)') m_type
            if (m_type .eq. 'L') then
               magnitude_scale = 'ML   '
            else if (m_type .eq. 'b') then
               magnitude_scale = 'mb   '
            else if (m_type .eq. 'B') then
               magnitude_scale = 'mB   '
            else if (m_type .eq. 's') then
               magnitude_scale = 'Ms   '
            else if (m_type .eq. 'S') then
               magnitude_scale = 'MS   '
            else if (m_type .eq. 'W') then
               magnitude_scale = 'MW   '
            else if (m_type .eq. 'G') then
               magnitude_scale = 'MbLg '
            else if (m_type .eq. 'c') then
               magnitude_scale = 'Mc   '
            else
               magnitude_scale = '     '
            end if               
            magnitude_usage = '='
            magnitude_author = hypo_author
            magnitude_orid = ' '
            call magnitude_line (2)
         end if
      end if
         
      if (linein(80:80) .eq. '7') then ! Phase lines begin
         loop_phase = .true.
         do while (loop_phase)
            read (1,'(a)') linein
            if (linein(1:4) .eq. '    ') exit
            if (linein(11:11) .eq. ' ') cycle ! Amplitude readings
            station_in = linein(2:6)//' '
            station_out = station_in
            if (code_convert) then
               do i = 1,ncode
                  if (station_out(1:5) .eq. code_loc(i)) then
                     station_out(1:5) = code_int(i)
                     exit
                  end if
               end do
            end if
            phase = linein(11:13)//'     '
            phase_year = hypo_year
            phase_month = hypo_month
            phase_day = hypo_day
            read (linein(19:22),'(2i2)') phase_hour, phase_minute
            if (linein(26:26) .eq. '.') then
               if (linein(28:28) .eq. ' ') then
                  read (linein(24:27),'(f4.1)') phase_seconds
                  iptim_pr = '-1' ! precision of phase readings to nearest tenth of a second
               else
                  read (linein(24:28),'(f5.2)') phase_seconds
                  iptim_pr = '-2' ! precision of phase readings to nearest hundredth of a second
               end if
            else if (linein(25:25) .eq. '.') then
               read (linein(24:26),'(f3.1)') phase_seconds
               iptim_pr = '-1' ! precision of phase readings to nearest tenth of a second
            end if
            ! Phase residual
            if (linein(66:66) .eq. '.') then
               read (linein(64:68),'(f5.2)') residual
               write (residual_pr,'(f5.1)') residual
            else if (linein(67:67) .eq. '.') then
               residual_pr = linein(64:68)
            else
               residual_pr = '     '
            end if
            ! Epicentral distance
            if (linein(73:73) .eq. '.') then
               read (linein(72:75),'(f4.2)') distance
               write (distance_pr,'(f6.2)') distance/111.
            else if (linein(74:74) .eq. '.') then
               read (linein(72:75),'(f4.1)') distance
               write (distance_pr,'(f6.2)') distance/111.
            else if (linein(75:75) .ne. ' ') then
               read (linein(72:75),'(f4.0)') distance
               write (distance_pr,'(f6.2)') distance/111.
            else
               distance_pr = '      '
            end if
            ! Azimuth
            if (linein(79:79) .ne. ' ') then
               azeq_pr = linein(77:79)
            else
               azeq_pr = '   '
            end if
            channel = linein(7:7)//' '//linein(8:8)
            location = '  '
            arrid_pr = ' '
            phase_usage = ' '
            station_flag = ' '
            phase_flag = ' '
            ! S-P phases from Reza Ghods
            if (linein(15:15) .eq. '9') then
               phase_usage = 't' ! Flag to indicate uncalibrated timing
               phase_seconds_p = phase_seconds
               phase_minute_p = phase_minute
               phase_hour_p = phase_hour
               call phase_line (2) ! Write the P reading phase line
               read (1,'(a)') linein  ! read the S phase line. Here I assume that after a 9 phase line an S phase reading will be found
               phase = linein(11:13)//'     '
               read (linein(19:22),'(2i2)') phase_hour, phase_minute
               if (linein(26:26) .eq. '.') then
                  if (linein(28:28) .eq. ' ') then
                     read (linein(24:27),'(f4.1)') phase_seconds
                     iptim_pr = '-1' ! precision of phase readings to nearest tenth of a second
                  else
                     read (linein(24:28),'(f5.2)') phase_seconds
                     iptim_pr = '-2' ! precision of phase readings to nearest hundredth of a second
                  end if
               else if (linein(25:25) .eq. '.') then
                  read (linein(24:26),'(f3.1)') phase_seconds
                  iptim_pr = '-1' ! precision of phase readings to nearest tenth of a second
               end if
               ! Phase residual
               if (linein(66:66) .eq. '.') then
                  read (linein(64:68),'(f5.2)') residual
                  write (residual_pr,'(f5.1)') residual
               else if (linein(67:67) .eq. '.') then
                  residual_pr = linein(64:68)
               else
                  residual_pr = '     '
               end if
               call phase_line (2) ! Write the S reading phase line

               diff = (phase_hour*3600. + phase_minute*60. + phase_seconds) -&
                (phase_hour_p*3600. + phase_minute_p*60. + phase_seconds_p)
               phase_hour = 0
               phase_minute = 0
               phase_seconds = diff
               phase = 'S-P     '
               phase_usage = ' '
            end if

            call phase_line (2) ! Write a phase line
         end do
      end if
      
      if (linein(1:6) .eq. '      ') then ! End of an event
         event_on = .false.
         write (2,'(a,t121,a)') 'STOP', ' '
         if (.not.bulletin) then
            write (2,'(a,t121,a)') 'EOF', ' '
            close (2)
         end if
      end if
      
   end do
   
   stop
      
end program seisan2mnf


!***********************************************************************
subroutine event_line (io_unit)

! Write an MNF event line
      
   implicit none
   
   include 'seisan2mnf.inc'
      
   integer :: io_unit

   write (io_unit,'(a1,1x,a1,1x,a10,1x,a,t121,a)') 'E', event_usage, evid, trim(event_geog), ' '
   
   return
   
end subroutine event_line


!***********************************************************************
subroutine hypocenter_line (io_unit)

! Write an MNF hypocenter line
      
   implicit none
   
   include 'seisan2mnf.inc'

   character(len=120) :: fmt
   character(len=8) :: latitude_pr
   character(len=9) :: longitude_pr
   integer :: io_unit
   
   fmt = '(2a,1x,i4,a,i2.2,a,i2.2,1x,i2.2,1x,i2.2,1x,f5.2,1x,5x,2x,a8,1x,a9,6x,6x,4x,f6.1,t95,a8,t112,a10)'
   
   latitude_pr = ' '
   write (latitude_pr(1:7),'(f7.3)') latitude
   longitude_pr = ' '
   write (longitude_pr(1:8),'(f8.3)') longitude
      
   write (io_unit,fmt)&
    'H ', hypo_usage, hypo_year, '/', hypo_month, '/', hypo_day, hypo_hour, hypo_minute,&
    hypo_seconds, latitude_pr, longitude_pr, depth, hypo_author, hypo_orid
    
   return
   
end subroutine hypocenter_line


!***********************************************************************
subroutine magnitude_line (io_unit)

! Write an MNF magnitude line
      
   implicit none
   
   include 'seisan2mnf.inc'
   
   integer :: io_unit
   
   write (io_unit,'(2a,1x,f3.1,2x,a5,1x,a8,t112,a10)') 'M ', magnitude_usage, magnitude,&
    magnitude_scale, magnitude_author, magnitude_orid
    
   return
   
end subroutine magnitude_line


!***********************************************************************
subroutine phase_line (io_unit)

! Write an MNF phase line
      
   implicit none
   
   include 'seisan2mnf.inc'
         
   integer :: io_unit, iptim
   character(len=120) :: fmt
   character(len=27) :: adslc
   character(len=6) phase_seconds_pr
   
   adslc = agency//'.'//deployment//'.'//station_in(1:5)//'.'//location//'.'//channel
   
   fmt = '(a1,1x,a1,a1,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,1x,4(i2,1x),a6,1x,a2,1x,a5,1x,a8,1x,a27,1x,a8,1x,a10)'
   
   phase_seconds_pr = ' '
   read (iptim_pr,'(i2)') iptim
   select case (iptim)
      case (0)
         write (phase_seconds_pr(1:3),'(f3.0)') phase_seconds
      case (-1)
         write (phase_seconds_pr(1:4),'(f4.1)') phase_seconds
      case (-2)
         write (phase_seconds_pr(1:5),'(f5.2)') phase_seconds
      case (-3)
         write (phase_seconds_pr(1:6),'(f6.3)') phase_seconds
      case default
         write (*,'(2a)') 'Error - illegal value for iptim: ', iptim_pr
         stop
   end select
      
   write (io_unit,fmt)&
    'P',& ! a1
    phase_usage,& ! a1
    station_flag,& ! a1
    station_out,& !a6
    distance_pr,& ! a6, formatted as f6.2
    azeq_pr,& ! a3, formatted as i3
    phase_flag,& ! a1
    phase,& ! a8
    phase_year,&
    phase_month,&
    phase_day,&
    phase_hour,& ! i2
    phase_minute,& ! i2
    phase_seconds_pr,& ! a6 formatted as f6.3
    iptim_pr,& ! a2, formatted as i2
    residual_pr,& ! a5 formatted as f5.1
    phase,& ! a8
    adslc,& ! a27
    phase_author,&  ! a8
    arrid_pr ! a10
    
end subroutine phase_line
