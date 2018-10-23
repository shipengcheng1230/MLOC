program any2mnf
      
! Reads input for earthquake hypocenters and arrival time data in several different formats
! and writes MNF event files. The version of MNF that is written is coded in the variable
! "version" that is initialized in the module "mnf_def".

! There are two modes for output:
!    Bulletin - All events concatenated in a single file whose name is based on the input file
!    Event - A separate file is created for each event, with the filename based on date and OT (yyyymmdd.hhmm.ss.mnf)
! If bulletin mode is applied to a file containing a single event, the only difference in the output file, other than
! filename, is the presence of a bulletin record at the top of the file.

! October 15, 2017 by eab
   
   use :: mnf_def
   use :: mnf_proc
   use :: case_proc

   implicit none
      
   character(len=80) :: filename, basename
   character(len=136) :: linein
   character(len=121) :: ev_line, evid_line
   character(len=120) :: comment
   character(len=1) :: chinp
   logical :: loop, loop_hypocenter, event_on, loop_magnitude, loop_phase, bulletin, ex, cfil
   integer :: ios, nev, n_hyp, if_type
   real :: azeq
 
   nev = 0
   event_on = .false.
   cfil = .false.
   
   write (*,'(a/)') 'Release date October 14, 2017, writing MNF v'//trim(version)
   
   write (*,'(4a/)')&
   '1 ISC ISF',&
   '2 NEIC IMS (via getfixed)',&
   '3 NEIC CSV (via getphase)',&
   '4 SEISAN'
   write (*,'(a)') 'Select input file type:'
   read (*,*) if_type
   
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
      call format_line (2) ! Write a format line
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
      
      if (linein(1:6) .eq. 'MSG_ID' .and. .not.event_on) then ! New event
         event_on = .true.
         nev = nev + 1
         evid = change_case(linein(8:17),-1)
      end if
      
      if (linein(1:6) .eq. 'Event ') then ! Event line
         event_annotation = ' '
         event_usage = ' '
         if (bulletin) then
            call event_line (2)
            call id_line (2)
         else
            write (ev_line,'(2a,1x,a,t121,a)') 'E ', event_usage, trim(event_annotation), ' '
            write (evid_line,'(2a,1x,a,t12,a,t121,a)') 'I ', ' ', 'NEIC', evid, ' '
         end if
      end if
      
      if (linein(1:8) .eq. '   Date ') then ! Hypocenter lines
         loop_hypocenter = .true.
         n_hyp = 0
         do while (loop_hypocenter)
            read (1,'(a)') linein
            if (linein(1:10) .eq. '1970/01/01') cycle ! Fix a bug in NEIC's ISF files
            if (linein(1:4) .eq. '    ') exit
            n_hyp = n_hyp + 1
            read (linein(1:10),'(i4,1x,i2,1x,i2)') hypo_year, hypo_month, hypo_day
            read (linein(12:22),'(i2,1x,i2,1x,f5.2)') hypo_hour, hypo_minute, hypo_seconds
            if (.not. bulletin .and. n_hyp .eq. 1) then ! Individual event files are named for the the first encountered hypocenter
               write (filename,'(i4,2i2.2,a,2i2.2,a,i2.2)') hypo_year, hypo_month, hypo_day, '.',&
                hypo_hour, hypo_minute, '.', int(hypo_seconds)
               inquire (file=trim(filename)//'.mnf',exist=ex)
               if (ex) then
                  write (*,'(a)') 'file '//trim(filename)//'.mnf already exists'
                  close (1)
                  stop
               end if
               open (2, file=trim(filename)//'.mnf', status='new')
               call format_line (2)
               write (2,'(a)') ev_line
               write (2,'(a)') evid_line
               if (cfil) then
                  write (3,'(a)') 'memb'
                  write (3,'(a)') 'even '//trim(filename)
                  write (3,'(a)') 'inpu '//trim(filename)//'.mnf'
               end if
            end if
            ot_err = linein(25:29)
            read (linein(37:44),'(f8.4)') latitude
            read (linein(46:54),'(f9.4)') longitude
            smaj = linein(56:60)
            smin = linein(62:66)
            iaz = linein(68:70)
            read (linein(72:76),'(f5.1)') depth
            read (linein(78:82),'(f5.1)') depth_err
            hypo_author = linein(119:126)
            hypo_orid = '          '//linein(129:136)
            if (n_hyp .eq. 1) then ! First encountered hypocenter is preferred
               hypo_usage = '='
            else
               hypo_usage = ' '
            end if
            depth_code = ' '
            call hypocenter_line (2)
         end do
      end if
      
      if (linein(1:9) .eq. 'Magnitude') then ! Magnitude lines
         loop_magnitude = .true.
         do while (loop_magnitude)
            read (1,'(a)') linein
            if (linein(1:4) .eq. '    ') exit
            magnitude_scale = linein(1:5)
            read (linein(8:10),'(f3.1)') magnitude
            magnitude_author = linein(21:28)
            magnitude_orid = '  '//linein(31:38)
            magnitude_usage = ' '
            call magnitude_line (2)
         end do
      end if
      
      if (linein(1:4) .eq. 'Sta ') then ! Phase lines
         loop_phase = .true.
         do while (loop_phase)
            read (1,'(a)') linein
            if (linein(1:4) .eq. '    ') exit
            if (linein(74:74) .ne. 'T') cycle ! Test for arrival time reading, maybe not right
            station = linein(1:5)//' '
            distance_pr = linein(7:12)
            read (linein(14:18),'(f5.1)') azeq
            write (azeq_pr,'(i3)') nint(azeq)
            phase = linein(20:27)
            read (linein(29:40),'(i2,1x,i2,1x,f6.3)') phase_hour, phase_minute, phase_seconds
            if (linein(40:40) .ne. '0') then
               iptim_pr = '-3' ! precision of phase readings to nearest thousandth of a second
            else
               iptim_pr = '-2'
            end if
            residual_pr = linein(42:46)
            agency = ' '
            channel = ' '
            deployment = ' '
            location = ' '
            phase_author = 'NEIC    '
            arrid_pr = ' '
            phase_usage = ' '
            station_flag = ' '
            phase_flag = ' '
            phase_year = hypo_year
            phase_month = hypo_month
            phase_day = hypo_day
            call phase_line (2) ! Write a phase line
         end do
      end if
      
      if (linein(1:4) .eq. 'STOP') then ! End of an event
         event_on = .false.
         write (2,'(a,t121,a)') 'STOP', ' '
         if (.not.bulletin) then
            write (2,'(a,t121,a)') 'EOF', ' '
            close (2)
         end if
      end if
      
   end do
   
   stop
   
end program any2mnf
