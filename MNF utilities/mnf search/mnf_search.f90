program mnf_search

! Search a bulletin in MNF format and select events, optionally create a command file for mloc.
! This version allows output as a new bulletin or as individual event files.
! 2016/12/16: added selection on focal depth (less than some value)
      
   implicit none
   
   integer, parameter :: n_line_max = 7000

   character(len=80) :: filename, basename
   character(len=121) :: linein, buffer(n_line_max)
   character(len=117) :: comment
   character(len=17) :: depth_line
   character(len=15) :: format_line
   character(len=10) :: evid
   character(len=6) :: current_version, evid_source
   character(len=1) :: chinp
   integer :: ios, nev, line, year, month, day, hour, minute, i, n_start, n_end, n_selected, nev_total,&
    nptotal(20), np, i10, npmin, nev_pass_total, output_type
   real :: seconds, latmin, latmax, lonmin, lonmax, lat_in, lon_in, depmax, dep_in, nearest_station, delta
   logical :: cfil, preferred, ex, n_select, epicenter_pass, preferred_depth,&
    bulletin, convert_format, depth_pass, nearest_pass, use_epicenter, use_depth, use_nearest
   
   cfil = .false.
   bulletin = .false.
   n_select = .false.
   buffer = ' '
   n_selected = 0
   current_version = '1.3.3 ' ! Current version of the MNF format
   convert_format = .false.
   evid_source = ' '
   use_epicenter = .false.
   use_depth = .false.
   use_nearest = .false.
   
   write (*,'(a/)') 'Release date October 7, 2018'
         
   write (*,'(a)') 'Enter input filename: '
   read (*,'(a)') filename
   open (1, file=trim(filename), status='old')
   write (*,'(a)') 'Create a new bulletin (1) or individual event files (2)?'
   read (*,*) output_type
   if (output_type .eq. 1) then
      bulletin = .true.
      write (*,'(a)') 'Bulletin comment: '
      read (*,'(a)') comment
   else if (output_type .eq. 2) then
      write (*,'(a)') 'Create mloc command file? '
      read (*,'(a)') chinp
      cfil = chinp .eq. 'y' .or. chinp .eq. 'Y'
      if (cfil) then
         write (*,'(a)') 'Enter command file basename: '
         read (*,'(a)') basename
         open (3, file=trim(basename)//'.cfil', status='new')
      end if
   end if
   
   ! Lat-Lon limits
   latmin = -90.
   latmax = 90.
   lonmin = -180.
   lonmax = 180.
   write (*,'(a)') 'Use lat-lon limits?'
   read (*,'(a)') chinp
   if (chinp .eq. 'y' .or. chinp .eq. 'Y') then
      use_epicenter = .true.
      write (*,'(a)') 'Enter latmin, latmax (-90 to 90):'
      read (*,*) latmin, latmax
      write (*,'(a)') 'Enter lonmin, lonmax (-180 to 180):'
      read (*,*) lonmin, lonmax
   end if
   
   ! Depth limit (events are selected for 0 ≤ depth ≤ depmax)
   depmax = 999.
   write (*,'(a)') 'Use focal depth limit?'
   read (*,'(a)') chinp
   if (chinp .eq. 'y' .or. chinp .eq. 'Y') then
      use_depth = .true.
      write (*,'(a)') 'Enter depmax:'
      read (*,*) depmax
   end if
   
   ! Nearest station
   nearest_station = 180.
   write (*,'(a)') 'Use nearest station distance?'
   read (*,'(a)') chinp
   if (chinp .eq. 'y' .or. chinp .eq. 'Y') then
      use_nearest = .true.
      write (*,'(a)') 'Enter nearest station distance (degrees): '
      read (*,*) nearest_station
   end if
   
   ! Scan the input file
   nev_total = 0
   nev_pass_total = 0
   nptotal = 0
   do
   
      read (1,'(a)',iostat=ios) linein
            
      if (linein(1:1) .eq. 'F') then ! Check if the Bulletin uses the current MNF format
         if (len_trim(linein) .eq. 12) then
            write (*,'(4a)') 'Files will be converted to current version (',&
             trim(current_version), ') from ', linein(10:len_trim(linein))
            convert_format = .true.
         end if
      else if (linein(1:2) .eq. 'E ') then ! New event
         nev_total = nev_total + 1
         epicenter_pass = .false.
         depth_pass = .false.
         nearest_pass = .false.
         np = 0
      else if (linein(1:1) .eq. 'H') then ! Get epicenter to test for search criteria
         if (use_epicenter) then
            read (linein(35:42),'(f8.4)') lat_in
            read (linein(44:52),'(f9.4)') lon_in
            epicenter_pass = lat_in .ge. latmin .and. lat_in .le. latmax .and. lon_in .ge. lonmin .and. lon_in .le. lonmax
         else
            epicenter_pass = .true.
         end if
         if (use_depth) then
            read (linein(70:74),'(f5.1)') dep_in
            depth_pass = dep_in .le. depmax
         else
            depth_pass = .true.
         end if
      else if (linein(1:1) .eq. 'P') then
         np = np + 1
         if (use_nearest) then
            if (.not.nearest_pass) then
               read (linein(12:17),'(f6.2)') delta
               nearest_pass = (delta .le. nearest_station)
            end if
         else
            nearest_pass = .true.
         end if
      else if (linein(1:4) .eq. 'STOP') then
         if (epicenter_pass .and. depth_pass .and. nearest_pass) then
            nev_pass_total = nev_pass_total + 1
            if (np .ge. 10) nptotal(1) = nptotal(1) + 1
            if (np .ge. 20) nptotal(2) = nptotal(2) + 1
            if (np .ge. 30) nptotal(3) = nptotal(3) + 1
            if (np .ge. 40) nptotal(4) = nptotal(4) + 1
            if (np .ge. 50) nptotal(5) = nptotal(5) + 1
            if (np .ge. 60) nptotal(6) = nptotal(6) + 1
            if (np .ge. 70) nptotal(7) = nptotal(7) + 1
            if (np .ge. 80) nptotal(8) = nptotal(8) + 1
            if (np .ge. 90) nptotal(9) = nptotal(9) + 1
            if (np .ge. 100) nptotal(10) = nptotal(10) + 1
            if (np .ge. 110) nptotal(11) = nptotal(11) + 1
            if (np .ge. 120) nptotal(12) = nptotal(12) + 1
            if (np .ge. 130) nptotal(13) = nptotal(13) + 1
            if (np .ge. 140) nptotal(14) = nptotal(14) + 1
            if (np .ge. 150) nptotal(15) = nptotal(15) + 1
            if (np .ge. 160) nptotal(16) = nptotal(16) + 1
            if (np .ge. 170) nptotal(17) = nptotal(17) + 1
            if (np .ge. 180) nptotal(18) = nptotal(18) + 1
            if (np .ge. 190) nptotal(19) = nptotal(19) + 1
            if (np .ge. 200) nptotal(20) = nptotal(20) + 1
        end if
      else if (linein(1:3) .eq. 'EOF') then
         exit
      end if
   end do
   write (*,'(i6,a)') nev_total, ' events read'
   write (*,'(i6,a)') nev_pass_total, ' events pass the search criteria'
   do i = 1,20
      i10 = i * 10
      write (*,'(i4,a,i3,a)') nptotal(i), ' events  that pass the search criteria and have ', i10, ' or more phase readings'
   end do
   rewind (1)
   
   write (*,'(a)') 'Minimum number of phase arrivals: '
   read (*,*) npmin
   
   write (*,'(a)') 'Event number selection: beginning and end numbers: '
   read (*,*) n_start, n_end
   
   nev = 0
   line = 0
   np = 0
   preferred = .false.
   preferred_depth = .false.
   epicenter_pass = .false.
   depth_pass = .false.
   nearest_pass = .false.
   if (bulletin) then
      open (2, file='mnf_search.out', status='new')
      write (2,'(a1,t5,a,t121,a)') 'B', trim(comment), ' '
      write (2,'(a,t121,a)') 'F   MNF v'//current_version, ' '
   end if
   do
      read (1,'(a)',iostat=ios) linein
      
      ! EOF
      if (ios .lt. 0) then
         close (1)
         if (bulletin) write (2,'(a,t121,a)') 'EOF', ' '
         close (2)
         if (cfil) close (3)
         write (*,'(a,i6,a)') 'EOF reached after ', nev, ' events'
         write (*,'(i4,a)') n_selected, ' events selected'
         stop
      end if
      
      if (linein(1:3) .eq. 'EOF') cycle ! Skip embedded EOF's
      
      if (linein(1:1) .eq. 'B') cycle ! Skip bulletin record
      
      if (linein(1:1) .eq. 'F') then
         if (convert_format) then
            format_line = 'F   MNF v'//current_version
         else
            format_line = linein(1:15)
         end if
         cycle
      end if
      
      if (linein(1:4) .eq. 'STOP') then ! Dump buffer
         if (n_select .and. (np .ge. npmin)) then
            n_selected = n_selected + 1
            if (preferred) then
               if (.not.bulletin) open (2, file=trim(filename)//'.mnf', status='new')
               if (cfil) then
                  write (3,'(a)') 'memb'
                  write (3,'(a)') 'even '//trim(filename)
                  write (3,'(a)') 'inpu '//trim(filename)//'.mnf'
                  if (preferred_depth) write (3,'(a)') 'dept '//depth_line
               end if
               if (.not.bulletin) write (2,'(a,t121,a)') format_line, ' '
               do i = 1,line
                  write (2,'(a)') buffer(i)
               end do
               write (2,'(a,t121,a)') 'STOP', ' '
               if (.not.bulletin) then
                  write (2,'(a,t121,a)') 'EOF', ' '
                  close (2)
               end if
            else
               write (*,'(a)') 'Error: no preferred hypocenter found'
               close (1)
               if (cfil) close (3)
               stop
            end if
         end if
         np = 0
         buffer = ' '
         line = 0
         preferred = .false.
         preferred_depth = .false.
         depth_line = ' '
         nearest_pass = .false.
      else ! Load buffer
         if (line .eq. n_line_max) then
            write (*,'(a)') 'Error: maximum number of lines in buffer reached'
            close (1)
            if (cfil) close (3)
            stop
         end if
         line = line + 1
         buffer(line) = linein
         if (linein(1:2) .eq. 'E ') then
            if (convert_format) then
               if (len_trim(linein) .ge. 112) then ! Save the event ID that could be stored at the end of the event record
                  evid = linein(112:121)
                  linein(112:121) = ' '
                  buffer(line) = linein
                  ! Create an ID record
                  line = line + 1
                  write (buffer(line),'(3a,t121,a)') 'I   ', evid_source, evid, ' '
               else
                  evid = ' '
               end if
            end if
         end if
         if (linein(1:3) .eq. 'H =') then ! Get parameters for filename and search criteria from preferred hypocenter record
            preferred = .true.
            if (use_epicenter) then
               read (linein(35:42),'(f8.4)') lat_in
               read (linein(44:52),'(f9.4)') lon_in
               epicenter_pass = lat_in .ge. latmin .and. lat_in .le. latmax .and. lon_in .ge. lonmin .and. lon_in .le. lonmax
            else
               epicenter_pass = .true.
            end if
            if (use_depth) then
               read (linein(70:74),'(f5.1)') dep_in
               depth_pass = dep_in .le. depmax
            else
               depth_pass = .true.
            end if
            read (linein(1:26),'(4x,i4,4(1x,i2),1x,f5.2)') year, month, day, hour, minute, seconds
            write (filename,'(i4,2i2.2,a,2i2.2,a,i2.2)') year, month, day, '.', hour, minute, '.', int(seconds)
            inquire (file=trim(filename)//'.mnf',exist=ex)
            if (ex) then
               write (*,'(a)') 'file '//trim(filename)//'.mnf already exists'
               close (1)
               stop
            end if
            nev = nev + 1
         end if
         if (linein(1:3) .eq. 'D =') then
            preferred_depth = .true.
            depth_line = linein(5:22)
            depth_line(7:7) = ' '
            depth_line(13:13) = ' '
         end if
         if (linein(1:2) .eq. 'P ') then
            np = np + 1
            if (use_nearest) then
               if (.not.nearest_pass) then
                  read (linein(12:17),'(f6.2)') delta
                  nearest_pass = (delta .lt. nearest_station)
               end if
            else
               nearest_pass = .true.
            end if
         end if
         n_select = nev .ge. n_start .and. nev .le. n_end .and. epicenter_pass .and. depth_pass .and. nearest_pass
      end if
            
   end do      
   
   stop
   
end program mnf_search
   