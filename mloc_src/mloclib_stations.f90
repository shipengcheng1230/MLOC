!***********************************************************************************************************************************
subroutine stafind (iev, ird)

! Find coordinates for a given station code. The algorithm now makes use of date codes
! for the start and end of operations of a given station at specified coordinates, allowing
! for changes in station location through time. If the deployment field for a reading is
! non-blank, then the matching algorithm honors it and requires a match for the deployment,
! as well as the station code, in the station list.

   implicit none

   include 'mloc.inc'
   
   integer :: i, iev, ird, k, julday, juldat_test
   logical :: date_range, match, failed, use_deployment
   character(len=132) :: msg

   failed = .false.
   use_deployment = deployment(iev,ird) .ne. '        '

   ! Julian date for current event
   call juldat (iyre(iev), mone(iev), idye(iev), julday, 0)
   juldat_test = iyre(iev)*1000 + julday

   ! Search station data list
   do k = 1,nstat1
      if (use_deployment) then
         if ((stname(iev,ird) .eq. nstr1(k)) .and. (deployment(iev,ird) .eq. sta_deployment(k))) then
            if (date_range(juldat_test,jdate_on(k),jdate_off(k))) then
               stladg(iev,ird) = stalat(k) ! Geocentric coordinates            
               stlndg(iev,ird) = stalon(k) ! Geocentric coordinates            
               ahgts(iev,ird) = ahgtr(k)
               sca0s(iev,ird) = sca0(k)
               sd1s(iev,ird) = sd1(k)
               kcode(iev,ird) = k
               call stations_used (k)
               return
            else
               failed = .true.
            end if
         end if      
      else
         if (stname(iev,ird) .eq. nstr1(k)) then
            if (date_range(juldat_test,jdate_on(k),jdate_off(k))) then
               stladg(iev,ird) = stalat(k) ! Geocentric coordinates            
               stlndg(iev,ird) = stalon(k) ! Geocentric coordinates            
               ahgts(iev,ird) = ahgtr(k)
               sca0s(iev,ird) = sca0(k)
               sd1s(iev,ird) = sd1(k)
               kcode(iev,ird) = k
               call stations_used (k)
               return
            else
               failed = .true.
            end if
         end if
      end if
   end do
   if (failed) then
      n_failed_date_range = n_failed_date_range + 1
      write (msg,'(a,2i4,1x,a,1x,3i8)') 'stafind: failed date_range: ',&
       iev, ird, stname(iev,ird), juldat_test, jdate_on(k), jdate_off(k)
      if (verbose_screen) call warnings (trim(msg))
      write (io_stn_log,'(3x,a)') msg
   end if

   ! station stname(iev,ird) is not in the station set.
   if (verbose_screen) call warnings ('stafind: missing station code: '//stname(iev,ird))
   kcode(iev,ird) = -1
   fcode(iev,ird) = 'm'
   stladg(iev,ird) = 0.
   stlndg(iev,ird) = 0.
   ahgts(iev,ird) = 0.
   sca0s(iev,ird) = 0.
   sd1s(iev,ird) = 0.
   if (nmiss(iev) .lt. n_miss_sta_max) then
      nmiss(iev) = nmiss(iev) + 1
   else
      write (msg,'(a,i4)') 'stafind: maximum number of missing stations reached: ', n_miss_sta_max
      call warnings (trim(msg))
   end if
   
   missta(iev,nmiss(iev)) = sad(iev,ird)
   if (n_miss_sta_total .eq. 0) then
      n_miss_sta(1) = 1
      n_miss_sta_total = 1
      n_miss_sta_list(1) = missta(iev,nmiss(iev))
   else
      match = .false.
      do i = 1,n_miss_sta_total
         if (sad(iev,ird) .eq. n_miss_sta_list(i)) then
            match = .true.
            n_miss_sta(i) = n_miss_sta(i) + 1
         end if
      end do
      if (.not.match) then
         if (n_miss_sta_total .lt. n_miss_sta_max) then
            n_miss_sta_total = n_miss_sta_total + 1
            n_miss_sta(n_miss_sta_total) = 1
            n_miss_sta_list(n_miss_sta_total) = sad(iev,ird)
         else
            write (msg,'(a,i4)') 'stafind: maximum number of missing stations reached: ', n_miss_sta_max
            call warnings (trim(msg))
         end if
      end if
   end if

   return
   
end subroutine stafind
      
      
!***********************************************************************************************************************************
subroutine stations_used (kindex)

! Keeps a record of all stations in the dataset for which station coordinates were found

   implicit none

   include 'mloc.inc'

   integer :: kindex, i
   character(len=132) :: msg

   if (nkstat .gt. 0) then
      if (nkstat .gt. nmax1) then
         msg = 'stations_used: nmax1 exceeded'
         call oops (trim(msg))
      end if
      do i = 1,nkstat
         if (kindex .eq. kstat(i)) return
      end do
      nkstat = nkstat + 1
      kstat(nkstat) = kindex
   else if (nkstat .eq. 0) then
      nkstat = 1
      kstat(1) = kindex
   else
      write (msg,'(a,i8)') 'stations_used: illegal value for kindex: ', kindex
      call oops (trim(msg))
   end if

   return
   
end subroutine stations_used


!***********************************************************************************************************************************
subroutine redsta ()

! Read station locations, elevations, and patch corrections from the main station data file and any
! supplemental station data files that have been declared. Supplemental station data files are read
! first so that they take precedence when searching the list of stations. Several formats are
! defined for supplemental station data files.

! All latitudes and longitudes are returned in geocentric coordinates.

   implicit none

   include 'mloc.inc'

   integer ::i, j, n, k, ios, latdeg, latmin, latsec, londeg, lonmin, lonsec, isstn, station_elevation, lun, lunit, dob
   real :: glat, glon, sxlat, cxlat, sxlon, cxlon, xlat, xlon, xlatmin, xlonmin, xlatsec, xlonsec
   character(len=132) :: msg, line, filename
   character(len=8) :: sta_author_default
   character(len=1) :: char1
     
   n = 0
   n_supp_stn = 0
   isstn = 99
   nstr1 = ' '
   sta_author = ' '
   sta_agency = ' '
   sta_deployment = ' '

   ! Read one or more supplemental station data files, for local (non-registered) stations.
   ! If there are any conflicts, the first instance encountered takes precedence (see stafind).
   ! The first column of the first line contains an index that specifies one of several standard formats.
   ! The rest of the first line can contain a comment about the source of the station data. All station
   ! codes and coordinates defined in supplemental station data files are written to the .log file.
   if (suppstn) then
      do i = 1,n_supp_stn_file
         write (char1,'(i1)') i
         sta_author_default = 'SSTN-'//char1//'  ' ! Author based on index of the supplemental station file
         lun = lunit()
         open (lun,file=suppfilnam(i),status='old')
         if (verbose_screen) write (*,'(/a)') trim(suppfilnam(i))
         write (io_stn_log,'(/a)') trim(suppfilnam(i))
         read (lun,'(a)') line
         if (verbose_screen) write (*,'(a)') trim(line)
         write (io_stn_log,'(a)') trim(line)
         read (line(1:1),'(i1)') isstn
         do k = 1,nmax1
            read (lun,'(a)',iostat=ios) line
            if (ios .lt. 0) exit
            if (line(1:1) .eq. '#') cycle ! Comment line
            n = n + 1
            if (n .gt. nmax1) call oops ('redsta: too many stations')
            select case (isstn)
         
               case (0) ! New master file format from April 28, 2014
                  nstr1(n) = line(1:5)
                  read (line(7:15),'(f9.5)') glat
                  read (line(17:26),'(f10.5)') glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  read (line(28:32),'(i5)') station_elevation
                  read (line(34:37),'(i4)') dob
                  ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                  sta_author(n) = line(39:46)
                  if (read_ad) then
                     sta_agency(n) = line(48:52)
                     sta_deployment(n) = line(54:61)
                  else
                     sta_agency(n) = ' '
                     sta_deployment(n) = ' '
                  end if
                  if (sta_author(n)(1:1) .eq. ' ') sta_author(n) = sta_author_default
                  read (line(66:72),'(i7)') jdate_on(n)
                  read (line(74:80),'(i7)') jdate_off(n)
                  sca0(n) = 0. ! No patch correction
                  sd1(n) = 1.0 ! Don't use patch correction variances.
                  if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,1x,a8)') nstr1(n), glat, stalat(n),&
                   stalon(n), station_elevation, dob, ahgtr(n), sta_author(n)
                  call station_duplicate (n, n-1, glat, duplication(n))

               case (1) ! ISC FFB format (geographic coordinates in deg-min-sec*10, elevation in m)
                  read (line(15:19),'(a5)') nstr1(n)
                  read (line(62:63),'(i2)') latdeg
                  read (line(64:65),'(i2)') latmin
                  read (line(66:68),'(i3)') latsec
                  glat = real(latdeg) + (real(latmin)/60.) + (real(latsec)/36000.)
                  if (line(69:69) .eq. 'S') glat = -glat
                  read (line(70:72),'(i3)') londeg
                  read (line(73:74),'(i2)') lonmin
                  read (line(75:77),'(i3)') lonsec
                  glon = real(londeg) + (real(lonmin)/60.) + (real(lonsec)/36000.)
                  if (line(78:78) .eq. 'W') glon = -glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  read (line(79:82),'(i4)') station_elevation
                  ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  sca0(n) = 0. ! No patch corrections.
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (2) ! SEISAN station list, with 4-character station code starting in column 3
               ! The standard format has been modified by adding optional fields for on- and off-date.
               ! If column 2 is non-blank, reads columns 2:6 for support of 5-character codes.
                  if (line(2:2) .eq. ' ') then
                     read (line(3:6),'(a4)') nstr1(n)(1:4)
                     nstr1(n)(5:5) = ' '
                  else
                     read (line(2:6),'(a5)') nstr1(n)
                  end if
                  read (line(7:8),'(i2)') latdeg
                  read (line(9:13),'(f5.2)') xlatmin
                  glat = real(latdeg) + xlatmin/60.
                  if (line(14:14) .eq. 'S') glat = -glat
                  read (line(15:17),'(i3)') londeg
                  read (line(18:22),'(f5.2)') xlonmin
                  glon = real(londeg) + xlonmin/60.
                  if (line(23:23) .eq. 'W') glon = -glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  read (line(24:27),'(i4)') station_elevation
                  ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  if (len_trim(line) .eq. 40) read (line(34:40),'(i7)') jdate_on(n) ! On-date only
                  if (len_trim(line) .eq. 48) then
                     read (line(34:40),'(i7)') jdate_on(n)
                     read (line(42:48),'(i7)') jdate_off(n)
                  end if
                  sca0(n) = 0. ! No patch corrections.
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4,2i8)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n),&
                   jdate_on(n), jdate_off(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (3) ! Simplified format, decimal degrees, geographic coordinates
               ! Format changed on February 25, 2017 with v10.3.3
               ! Station code in columns 1:5 (5-character station codes are max in v10.4.0)
               ! Decimal geographic latitude in columns 8:15 in f8.4 format
               ! Decimal geographic longitude in columns 17:25 in f9.4 format
               ! Integer elevation in meters in columns 27:31 in i5 format
               ! Integer depth of burial in meters (positive) in columns 33:37 in i5 format
               ! Integer date_on in columns 39:45 in i7 format (4-digit year plus 3-digit day of year)
               ! Integer date_off in columns 47:53 in i7 format (4-digit year plus 3-digit day of year)
                  nstr1(n) = line(1:5)
                  read (line(8:15),'(f8.4)') glat
                  read (line(17:25),'(f9.4)') glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  station_elevation = 0
                  if (len_trim(line) .ge. 31) read (line(27:31),'(i5)') station_elevation ! In meters
                  dob = 0
                  if (len_trim(line) .ge. 37) read (line(33:37),'(i5)') dob ! Depth of burial in meters (positive number)
                  ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  if (len_trim(line) .eq. 45) read (line(39:45),'(i7)') jdate_on(n) ! On-date only
                  if (len_trim(line) .eq. 53) then
                     read (line(39:45),'(i7)') jdate_on(n)
                     read (line(47:53),'(i7)') jdate_off(n)
                  end if
                  sca0(n) = 0. ! No patch corrections
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,2i8)') nstr1(n), glat, stalat(n), stalon(n),&
                   station_elevation, dob, ahgtr(n), jdate_on(n), jdate_off(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (4) ! China Seismic Bureau format
               ! Lower case seismic code in columns 1:3
               ! Elevation in meters in columns 5:8 (i3)
               ! Geographic latitude in degrees-minute-secs in columns 10:21
               ! Longitude in degrees-minute-secs in columns 25:37
                  nstr1(n) = line(1:3)//'  '
                  read (line(5:8),'(i4)') station_elevation
                  ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                  read (line(10:11),'(i2)') latdeg
                  read (line(14:15),'(i2)') latmin
                  read (line(18:21),'(f4.1)') xlatsec
                  glat = real(latdeg) + (real(latmin)/60.) + xlatsec/3600.
                  read (line(25:27),'(i3)') londeg
                  read (line(30:31),'(i2)') lonmin
                  read (line(34:37),'(f4.1)') xlonsec
                  glon = real(londeg) + (real(lonmin)/60.) + xlonsec/3600.
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  sca0(n) = 0. ! No patch corrections
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (5) ! NEIC format
                  nstr1(n) = line(4:8)
                  read (line(40:47),'(f8.4)') glat
                  read (line(48:57),'(f10.4)') glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  read (line(58:62),'(i5)') station_elevation ! In meters
                  dob = 0 ! Depth of burial in meters (positive number)
                  ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                  jdate_on(n) = 0 ! Ignore operational epoch
                  jdate_off(n) = 0 ! Ignore operational epoch
                  sca0(n) = 0. ! No patch corrections
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (read_ad) then
                     sta_agency(n) = 'FDSN '
                     sta_deployment(n) = line(1:2)//'      '
                  else
                     sta_agency(n) = ' '
                     sta_deployment(n) = ' '
                  end if
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (6) ! MSU format
               ! Seismic code in columns 1:5
               ! Elevation in meters in columns 30:33 (i4)
               ! Geographic latitude in degrees-minute-secs in columns 6:15
               ! Longitude in degrees-minute-secs in columns 17:27
                  nstr1(n) = line(1:5)
                  read (line(30:33),'(i4)') station_elevation
                  ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                  read (line(6:7),'(i2)') latdeg
                  read (line(9:10),'(i2)') latmin
                  read (line(12:15),'(f4.1)') xlatsec
                  glat = real(latdeg) + (real(latmin)/60.) + xlatsec/36000.
                  if (line(16:16) .eq. 'S') glat = -glat
                  read (line(17:19),'(i3)') londeg
                  read (line(21:22),'(i2)') lonmin
                  read (line(24:27),'(f4.1)') xlonsec
                  glon = real(londeg) + (real(lonmin)/60.) + xlonsec/36000.
                  if (line(28:28) .eq. 'W') glon = -glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  sca0(n) = 0. ! No patch corrections
                  sd1(n) = 1.0 ! No patch correction variances
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
               
               case (9) ! Same format as the former (up to v10.0.0) master station file
                  nstr1(n) = line(2:6)
                  read (line(8:16),'(f9.3)') stalat(n)
                  read (line(17:26),'(f10.3)') stalon(n)
                  read (line(27:36),'(f10.3)') ahgtr(n)
                  jdate_on(n) = 0
                  jdate_off(n) = 0
                  read (line(37:46),'(f10.3)') sca0(n)
                  sd1(n) = 1.0 ! Don't use patch correction variances.
                  call geogra (stalat(n), glat)
                  sta_author(n) = sta_author_default ! This format should only show up as a supplemental station file now
                  if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  call station_duplicate (n, n-1, glat, duplication(n))
                  
               case (99) ! No index number was supplied
                  call warnings ('redsta: no value supplied for isstn')
                  exit
               
               case default
                  write (msg,'(a,i3)') 'redsta: illegal value for isstn: ', isstn
                  call warnings (trim(msg))
                  exit
               
            end select
         end do
         
         n_supp_stn = n ! Updated only after all entries are read from a supplemental file
         if (verbose_screen) write (*,'(a,i6)') 'Total number of stations = ', n_supp_stn
         close (lun)
      end do
   end if
   if (suppstn) write (io_stn_log,'(/i6,a)') n_supp_stn, ' supplemental stations read'

   ! Read the main station data file.
   ! From April 28, 2014, the file has a completely new format
   ! The former format is still supported as isstn = 9
   filename = trim(station_path)//dirsym//trim(station_master)
   if (verbose_screen) write (*,'(/a)') 'Main station file:'
   write (io_stn_log,'(/a)') 'Main station file:'
   if (verbose_screen) write (*,'(a)') trim(filename)
   write (io_stn_log,'(a)') trim(filename)
   lun = lunit()
   open (lun,file=filename,status='old')
   read (lun,'(a)') line
   read (line(1:1),'(i1)') isstn
   select case (isstn)
   
      case (0) ! New master file format from April 28, 2014
         do k = 1,nmax1
            read (lun,'(a)',iostat=ios) line
            if (ios .lt. 0) exit
            if (line(1:1) .eq. '#') cycle ! Comment line
            n = n + 1
            if (n .gt. nmax1) call oops ('redsta: too many stations')
            nstr1(n) = line(1:5)
            read (line(7:15),'(f9.5)') glat
            read (line(17:26),'(f10.5)') glon
            call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
            stalat(n) = xlat
            stalon(n) = xlon
            read (line(28:32),'(i5)') station_elevation
            read (line(34:37),'(i4)') dob
            ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
            sta_author(n) = line(39:46)
            if (read_ad) then
               sta_agency(n) = 'ISC  '
               sta_deployment(n) = 'IR      '
            else
               sta_agency(n) = ' '
               sta_deployment(n) = ' '
            end if
            read (line(66:72),'(i7)') jdate_on(n)
            read (line(74:80),'(i7)') jdate_off(n)
            sca0(n) = 0. ! No patch correction
            sd1(n) = 1.0 ! Don't use patch correction variances.
            if (suppstn) call station_duplicate (n, n_supp_stn, glat, duplication(n)) ! Only check against stations from supplemental station files
            if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,1x,a8)') nstr1(n), glat, stalat(n), stalon(n),&
             station_elevation, dob, ahgtr(n), sta_author(n)
         end do
         
      case (9) ! Former master station file format, decimal geographic coordinates, elevation in m, operational epochs
         do k = 1,nmax1
            read (lun,'(a)',iostat=ios) line
            if (ios .lt. 0) exit
            if (line(2:2) .ne. ' ') then
               n = n + 1
               if (n .gt. nmax1) call oops ('redsta: too many stations')
               nstr1(n) = line(2:6)
               read (line(60:67),'(f8.4)') glat
               read (line(69:77),'(f9.4)') glon
               call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
               stalat(n) = xlat
               stalon(n) = xlon
               read (line(80:84),'(i5)') station_elevation
               ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
               read (line(90:96),'(i7)') jdate_on(n)
               read (line(99:105),'(i7)') jdate_off(n)
               sca0(n) = 0. ! No patch correction
               sd1(n) = 1.0 ! Don't use patch correction variances.
               sta_author(n) = sta_author_default
               if (suppstn) call station_duplicate (n, n_supp_stn, glat, duplication(n)) ! Only check against stations from supplemental station files
               if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,i6,f8.3)') nstr1(n), glat, stalat(n), stalon(n),&
                station_elevation, ahgtr(n)
            end if
         end do
         
      case default
         write (msg,'(a,i1)') 'redsta: illegal value for isstn for master station file: ', isstn
         call oops (trim(msg))
         
   end select
      
   write (msg,'(a,i6)') 'redsta: total number of stations        ', n
   call fyi (trim(msg))
   write (io_stn_log,'(a,i6)') 'Total number of stations        ', n
   if (suppstn) then
      write (msg,'(a,i6)') 'redsta: supplemental stations           ', n_supp_stn
      call fyi (trim(msg))
      write (msg,'(a,i6)') 'redsta: duplicate stations              ', n_dupe
      call fyi (trim(msg))
      write (msg,'(a,i6)') 'redsta: ...pure duplicates              ', n_dupe_pure
      call fyi (trim(msg))
      write (msg,'(a,i6)') 'redsta: ...with minor differences       ', n_dupe_minor
      call fyi (trim(msg))
      write (msg,'(a,i6)') 'redsta: ...with significant differences ', n_dupe_significant
      call fyi (trim(msg))
      write (msg,'(a,i6)') 'redsta: ...station conflicts            ', n_dupe_conflict
      call fyi (trim(msg))
      write (io_stn_log,'(/a,i6)') 'Duplicate stations              ', n_dupe
      write (io_stn_log,'(a,i6)') '...pure duplicates              ', n_dupe_pure
      write (io_stn_log,'(a,i6)') '...with minor differences       ', n_dupe_minor
      write (io_stn_log,'(a,i6)') '...with significant differences ', n_dupe_significant
      write (io_stn_log,'(a,i6)') '...station conflicts            ', n_dupe_conflict
   end if
   close (lun)
      
   nstat1 = n
      
   return
   
end subroutine redsta
      
      
!***********************************************************************************************************************************
subroutine station_duplicate (n_target, n_test, glat, result)

! Check for duplicate station codes.
! n_target = index of the current station.
! n_test = index of the last station against which to test.

   implicit none

   include 'mloc.inc'

   integer :: j, n_target, n_test
   real :: glat, glatj, diff_lat, diff_lon, diff_hgt, lon1, lon2
   character(len=40) :: result
   
   if (n_test .le. 1) return

   result = ' '
   do j = 1,n_test
      if (nstr1(j) .eq. nstr1(n_target)) then
         n_dupe = n_dupe + 1
         call geogra (stalat(j), glatj)
         diff_lat = glatj - glat
         diff_lon = stalon(j) - stalon(n_target)
         diff_hgt = ahgtr(j) - ahgtr(n_target)
         ! Region of equivalnce is about 100 m
         if (abs(diff_lat) .gt. 1.1 .or. abs(diff_lon) .gt. 1.1) then
            result = ' conflicting station codes'
            n_dupe_conflict = n_dupe_conflict + 1
         else if (abs(diff_lat) .gt. 1.0e-3 .or. abs(diff_lon) .gt. 1.0e-3 .or. abs(diff_hgt) .gt. 1.0e-1) then
            result = ' significant difference in coordinates'
            n_dupe_significant = n_dupe_significant + 1
         else if (abs(diff_lat) .gt. 1.0e-4 .or. abs(diff_lon) .gt. 1.0e-4 .or. abs(diff_hgt) .gt. 1.0e-2) then
            result = ' minor difference in coordinates'
            n_dupe_minor = n_dupe_minor + 1
         else ! Same coordinates
            result = ' '
            n_dupe_pure = n_dupe_pure + 1
         end if
         lon1 = stalon(j)
         call set_longitude_range (lon1, 0)
         lon2 = stalon(n_target)
         call set_longitude_range (lon2, 0)
         write (io_stn_log,'(a,1x,a,1x,a,1x,f10.4,10x,f10.4,f10.3,1x,a)') nstr1(j), sta_agency(j), sta_deployment(j),&
          glatj, lon1, ahgtr(j), sta_author(j) ! First instance
         write (io_stn_log,'(a,1x,a,1x,a,1x,f10.4,10x,f10.4,f10.3,1x,a)') nstr1(n_target), sta_agency(n_target),&
          sta_deployment(n_target), glat, lon2, ahgtr(n_target), sta_author(n_target)
         write (io_stn_log,'(21x,f10.4,10x,f10.4,f10.3,a/)') diff_lat, diff_lon, diff_hgt, result
      end if
   end do

   return
   
end subroutine station_duplicate


!*****************************************************************************************
logical function bdp (station)

! Returns "true" if the station is on a list of stations that are suspected of reporting
! depth phases based on theoretical arrival times from a preliminary location such as the PDE.
! This has been documented for a number of Chinese stations.

   implicit none
   
   include 'mloc.inc'
   
   character(len=5) station
   integer i
   
   bdp = .false.
   
   if (bdp_list) then
      do i = 1,n_bdp
         if (station .eq. bdp_station(i)(2:6)) then
            bdp = .true.
            exit
         end if
      end do
   end if
   
   return

end function bdp


!*****************************************************************************************
subroutine station_check (iev, ird)

! Check for a few known problems with certain stations, mainly derived from Bob Engdahl's
! experiences with the EHB and GEM processing. There are two primary situations:
! 1) Certain stations are known to have timing problems during certain periods, and readings
!    in these periods are flagged (fcode = 't').
! 2) Certain stations were moved significant distances without changing the station code 
!    and a new code for the new location has been assigned retroactively.

   implicit none
   
   include 'mloc.inc'

   integer :: iev, ird
   character(len=132) :: msg
   
   if (stname(iev,ird) .eq. 'DAVO ') then
      if (ipay(iev,ird) .lt. 1999) then
         stname(iev,ird) = 'DAVOS'
      else
         stname(iev,ird) = 'DAVON'
      end if
      call fyi ('station_check: station code DAVO changed to '//stname(iev,ird))
   else if (stname(iev,ird) .eq. 'FCC  ' .and. (ipay(iev,ird) .ge. 1992 .and. ipay(iev,ird) .le. 1993)) then
      fcode(iev,ird) = 't'
      write (msg,'(a,i4)') 'station_check: station FCC had timing problems in ', ipay(iev,ird)
      call fyi (trim(msg))
   else if (stname(iev,ird) .eq. 'FITZ*') then
      stname(iev,ird) = 'FITZ '
      call fyi ('station_check: station code FITZ* changed to FITZ')
   else if (stname(iev,ird) .eq. 'MAG  ' .and. ipay(iev,ird) .gt. 1970) then
      stname(iev,ird) = 'MGD  '
      call fyi ('station_check: station MAG moved without changing codes; subsequent readings changed to MGD')
   else if (stname(iev,ird) .eq. 'NRI  ' .and. ipay(iev,ird) .gt. 1993) then
      stname(iev,ird) = 'NRIS '
      call fyi ('station_check: station NRI moved without changing codes; subsequent readings changed to NRIS')
   else if (stname(iev,ird) .eq. 'TUC  ' .and. (ipay(iev,ird) .ge. 1992 .and. ipay(iev,ird) .le. 1995)) then
      fcode(iev,ird) = 't'
      write (msg,'(a,i4)') 'station_check: station TUC had timing problems in ', ipay(iev,ird)
      call fyi (trim(msg))
   else if (stname(iev,ird) .eq. 'WIT  ' .and. (ipay(iev,ird) .ge. 1964 .and. ipay(iev,ird) .le. 1979)) then
      fcode(iev,ird) = 't'
      write (msg,'(a,i4)') 'station_check: station WIT had timing problems in ', ipay(iev,ird)
      call fyi (trim(msg))
   else if (stname(iev,ird) .eq. 'YKA  ' .and. (ipay(iev,ird) .ge. 1991 .and. ipay(iev,ird) .le. 1995)) then
      fcode(iev,ird) = 't'
      write (msg,'(a,i4)') 'station_check: station YKA had timing problems in ', ipay(iev,ird)
      call fyi (trim(msg))
   end if

   return
end subroutine station_check

      

   
