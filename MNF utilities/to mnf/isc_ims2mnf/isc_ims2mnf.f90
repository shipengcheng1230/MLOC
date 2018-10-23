PROGRAM isc_ims2mnf
      
! Reads a file containing arrival time data for one or more events in the IMS1.0 format returned by the ISC Bulletin web search
! when the options for header information and comments are selected. This is the format of the IASPEI Reference Event Bulletin too.
! The output is in MNF (mloc native format) format, v1.3.

! This version supports MNF v1.3.3, putting event ID into a separate record.

! There are two modes for output:
!    Bulletin - All events concatenated in a single file whose name is based on the input file
!    Event - A separate file is created for each event, with the filename based on date and OT (yyyymmdd.hhmm.ss.mnf)
! If bulletin mode is applied to a file containing a single event, the only difference in the output file, other than
! filename, is the presence of a bulletin record at the top of the file.

   implicit none
   
   include 'isc_ims2mnf.inc'
   
   character(len=80) :: filename, basename
   character(len=200) :: linein
   character(len=120) :: comment
   character(len=6) version
   character(len=5) char5
   character(len=1) :: chinp
   logical :: event_on, bulletin, ex, cfil, phase_ok, at, event_written, preferred_hypo
   integer :: ios, nev, phase_seconds_int, iblank, int2
   real :: azeq, smin, smaj
 
   version = '1.3.3 ' ! MNF version
   nev = 0
   event_on = .false.
   cfil = .false.
   preferred_hypo = .false.
   
   write (*,'(a/)') 'Release date April 6, 2018, writing MNF v'//version
   
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
      call bulletin_line (2, comment) ! Write a bulletin line
      call format_line (2, version) ! Write a format line
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
 
   do ! loop over all lines
   
      read (1,'(a)',iostat=ios) linein
      
      ! EOF
      if (linein(1:4) .eq. 'STOP') then
         close (1)
         if (bulletin) then
            write (2,'(a,t121,a)') 'EOF', ' '
            close (2)
         end if
         if (cfil) close (3)
         write (*,'(a,i6,a)') 'EOF reached after ', nev, ' events'
         stop
      end if
            
      if (linein(1:6) .eq. 'Event ') then
         if (.not. event_on) then ! New event
            event_on = .true.
            event_written = .false.
            preferred_hypo = .false.
            nev = nev + 1
            event_usage = ' '
            event_anno = linein(7:len_trim(linein))
            event_anno = adjustl(event_anno)
            ! Extract the evid from the annotation
            iblank = index(event_anno, ' ', .false.)
            evid = ' '
            if (iblank .le. 11) then
               evid = event_anno(1:iblank-1)
               evid = adjustr(evid)
            else
               write (*,'(2a)') 'Warning: evid truncated: ', trim(event_anno)
               evid = event_anno(1:10)
            end if
            ! Remove the evid part of the annotation
            event_anno(1:iblank) = ' '
            event_anno = adjustl(event_anno)
         else ! Finish off an event
            event_on = .false.
            write (2,'(a,t121,a)') 'STOP', ' '
            if (.not.bulletin) then
               write (2,'(a,t121,a)') 'EOF', ' '
               close (2)
            end if
            backspace (1)
         end if
      end if
      
      if (linein(5:5) .eq. '/' .and. linein(8:8) .eq. '/') then ! Hypocenter line  
         
         read (linein(1:10),'(i4,1x,i2,1x,i2)') hypo_year, hypo_month, hypo_day
         read (linein(12:22),'(i2,1x,i2,1x,a5)') hypo_hour, hypo_minute, char5
         if (char5(3:3) .eq. '.') then
            read (char5,'(f5.2)') hypo_seconds
         else
            read (char5(1:2),'(i2)') int2
            hypo_seconds = float(int2)
         end if
         
         ! Open individual event files if not bulletin output
         if (bulletin) then
            if (.not.event_written) then
               call event_line (2)
               call id_line (2)
               event_written = .true.
            end if
         else
            if (.not.event_written) then
               write (filename,'(i4,2i2.2,a,2i2.2,a,i2.2)') hypo_year, hypo_month, hypo_day, '.',&
                hypo_hour, hypo_minute, '.', int(hypo_seconds)
               inquire (file=trim(filename)//'.mnf',exist=ex)
               if (ex) then
                  write (*,'(a)') 'file '//trim(filename)//'.mnf already exists'
                  close (1)
                  stop
               end if
               open (2, file=trim(filename)//'.mnf', status='new')
               call format_line (2, version)
               call event_line (2)
               call id_line (2)
               event_written = .true.
               if (cfil) then
                  write (3,'(a)') 'memb'
                  write (3,'(a)') 'even '//trim(filename)
                  write (3,'(a)') 'inpu '//trim(filename)//'.mnf'
               end if
            end if
         end if
         
         ! Remainder of the hypocenter info
         read (linein(37:44),'(f8.4)') latitude
         read (linein(46:54),'(f9.4)') longitude
         read (linein(72:76),'(f5.1)') depth
         smaj_pr = ' '
         smin_pr = ' '
         az_smin_pr = ' '
         ! The search output formatting is variable. When it's fixed I can bring this section back.
!         if (linein(57:57) .eq. '.' .and. linein(63:63) .eq. '.' .and. linein(70:70) .ne. ' ') then
!            read (linein(55:60),'(f6.1)') smaj
!            write (smaj_pr,'(f5.2)') smaj
!            read (linein(61:66),'(f6.1)') smin
!            write (smin_pr,'(f5.2)') smin
!            az_smin_pr = linein(68:70)
!         end if
         hypo_orid = linein(127:136)
         hypo_author = linein(119:126)
         ! First hypocenter line is preferred
         if (.not.preferred_hypo) then
            hypo_usage = '='
            preferred_hypo = .true.
         else
            hypo_usage = ' '
         end if
         ! Special section for the IASPEI Reference Event Bulletin
!         if (hypo_author(1:6) .eq. 'IASPEI') then
!            hypo_usage = '='
!         else
!            hypo_usage = ' '
!         end if
         call hypocenter_line (2)
      end if
         
      if (linein(1:9) .eq. 'Magnitude') then ! Magnitude lines
         read (1,'(a)') linein
         if (linein(1:1) .ne. ' ') then
            read (linein(8:10),'(f3.1)') magnitude
            magnitude_scale = linein(1:5)
            magnitude_usage = ' '
            magnitude_author = linein(21:28)
            magnitude_orid = linein(29:38)
            call magnitude_line (2)
            do
               read (1,'(a)') linein
               if (linein(1:1) .eq. ' ') then
                  exit
               else
                  read (linein(8:10),'(f3.1)') magnitude
                  magnitude_scale = linein(1:5)
                  magnitude_usage = ' '
                  magnitude_author = linein(21:28)
                  magnitude_orid = linein(29:38)
                  call magnitude_line (2)
               end if
            end do
         end if
      end if
      
      if (linein(1:4) .eq. 'Sta ') then ! Phase readings
         do
            read (1,'(a)') linein
            ! Comment lines within the phase data
            if (linein(2:2) .eq. '(') then
              write (2,'(a1,a)') '#', trim(linein)
              cycle
            end if
            if (linein(1:4) .eq. '    ') then
               event_on = .false.
               write (2,'(a,t121,a)') 'STOP', ' '
               if (.not.bulletin) then
                  write (2,'(a,t121,a)') 'EOF', ' '
                  close (2)
               end if
               exit
            end if
            station = linein(1:5)//' '
            distance_pr = linein(7:12)
            read (linein(14:18),'(f5.1)') azeq
            write (azeq_pr,'(i3)') nint(azeq)
            phase = linein(20:27)
            phase_year = hypo_year
            phase_month = hypo_month
            phase_day = hypo_day
            read (linein(29:40),'(i2,1x,i2)') phase_hour, phase_minute
            if (linein(37:37) .eq. '.') then
               read (linein(35:40),'(f6.3)') phase_seconds
               iptim_pr = '-3'
               if (linein(40:40) .eq. ' ') iptim_pr = '-2'
               if (linein(39:39) .eq. ' ') iptim_pr = '-1'
               if (linein(38:38) .eq. ' ') iptim_pr = ' 0'
            else if (linein(37:37) .eq. ' ') then
               read (linein(35:36),'(i2)') phase_seconds_int
               phase_seconds = float(phase_seconds_int)
               iptim_pr = ' 0'
            else
               write (*,'(2a)') 'Error reading: ', trim(linein)
               close (1)
               stop
            end if
            if (linein(31:31) .eq. ':') then ! Check that arrival time field is not blank
               at = .true.
            else
               at = .false.
            end if
            residual_pr = linein(42:46)
            channel = '   '
            location = '  '
            arrid_pr = ' '//linein(114:122)
            phase_usage = ' '
            station_flag = ' '
            phase_flag = ' '
            phase_author = 'ISC'
            if (phase_ok (phase) .and. at) call phase_line (2) ! Write a phase line
         end do
         
      end if
      
      ! Comment lines
      if (linein(2:2) .eq. '(') then
        write (2,'(a1,a)') '#', trim(linein)
      end if
            
   end do
   
   stop
      
end program isc_ims2mnf

!***********************************************************************
subroutine bulletin_line (io_unit, comment)

! Write an MNF bulletin line
      
   implicit none
      
   character(len=120) :: comment
   integer :: io_unit
   
   write (io_unit,'(2a,t121,a)') 'B   ', trim(comment), ' '
   
   return
   
end subroutine bulletin_line


!***********************************************************************
subroutine format_line (io_unit, version)

! Write an MNF format line
      
   implicit none
      
   integer :: io_unit
   character(len=6) :: version
   
   write (io_unit,'(2a,t121,a)') 'F   MNF v', version, ' '
   
   return
   
end subroutine format_line


!***********************************************************************
subroutine event_line (io_unit)

! Write an MNF event line
      
   implicit none
   
   include 'isc_ims2mnf.inc'
      
   integer :: io_unit

   write (io_unit,'(a1,1x,a1,1x,a,t121,a)') 'E', event_usage, trim(event_anno), ' '
   
   return
   
end subroutine event_line


!***********************************************************************
subroutine id_line (io_unit)

! Write an MNF event ID line
      
   implicit none
   
   include 'isc_ims2mnf.inc'
      
   integer :: io_unit

   write (io_unit,'(a,t5,a,t12,a,t121,a)') 'I', 'ISC', trim(evid), ' '
   
   return
   
end subroutine id_line


!***********************************************************************
subroutine hypocenter_line (io_unit)

! Write an MNF hypocenter line
      
   implicit none
   
   include 'isc_ims2mnf.inc'

   character(len=120) :: fmt
   integer :: io_unit
   
   fmt = '(2a,1x,i4,a,i2.2,a,i2.2,1x,i2.2,1x,i2.2,1x,f5.2,1x,5x,2f10.4,1x,a3,2(1x,a5),f6.1,t95,a8,t112,a10)'
      
   write (io_unit,fmt)&
    'H ', hypo_usage, hypo_year, '/', hypo_month, '/', hypo_day, hypo_hour, hypo_minute,&
    hypo_seconds, latitude, longitude, az_smin_pr, smin_pr, smaj_pr, depth, hypo_author, hypo_orid
    
   return
   
end subroutine hypocenter_line


!***********************************************************************
subroutine magnitude_line (io_unit)

! Write an MNF magnitude line
      
   implicit none
   
   include 'isc_ims2mnf.inc'
   
   integer :: io_unit
   
   write (io_unit,'(2a,1x,f3.1,2x,a5,1x,a8,t112,a10)') 'M ', magnitude_usage, magnitude,&
    magnitude_scale, magnitude_author, magnitude_orid
    
   return
   
end subroutine magnitude_line


!***********************************************************************
subroutine phase_line (io_unit)

! Write an MNF phase line
      
   implicit none
   
   include 'isc_ims2mnf.inc'
         
   integer :: io_unit, iptim
   character(len=120) :: fmt
   character(len=27) :: adslc
   character(len=6) :: phase_seconds_pr
   
   adslc = agency//'.'//deployment//'.'//station(1:5)//'.'//location//'.'//channel
   
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
    station,& !a6
    distance_pr,& ! a6, formatted as f6.2
    azeq_pr,& ! a3, formatted as i3
    phase_flag,& ! a1
    phase,& ! a8
    phase_year,&
    phase_month,&
    phase_day,&
    phase_hour,& ! i2
    phase_minute,& ! i2
    phase_seconds_pr,& ! a6, formatted in several ways
    iptim_pr,& ! a2, formatted as i2
    residual_pr,& ! a5 formatted as f5.1
    phase,& ! a8
    adslc,& ! a27
    phase_author,&  ! a8
    arrid_pr ! a10
    
   return
   
end subroutine phase_line


!*****************************************************************************************
logical function phase_ok (phasein)

   implicit none
   
   character(len=8) :: phasein, bad_phase(48)
   integer :: nbad, i
   
   data bad_phase /'x       ','LZ      ','LN      ','LE      ','LR      ','AML     ','AMB     ','AMS     ',&
                   'pmax    ','MLR     ','IAmb    ','SME     ','SMN     ','LmH     ','LmV     ','PMZ     ',&
                   'PFAKE   ','tx      ','smax    ','LRM     ','MAXIMUM ','LmE     ','LMZ     ','LME     ',&
                   'LMN     ','MLRZ    ','LQ      ','LRZ     ','LRN     ','LRE     ','L       ','FINAL   ',&
                   'UP      ','D       ','MAX     ','QM      ','RM      ','M       ','QSN     ','Trac    ',&
                   'IAML    ','ml      ','XPR1    ','********','********','********','********','********'/
   data nbad /43/
   
   phase_ok = .true.
   do i=1,nbad
      if (phasein .eq. bad_phase(i)) then
         phase_ok = .false.
         exit
      end if
   end do
   
   return
   
end function phase_ok
