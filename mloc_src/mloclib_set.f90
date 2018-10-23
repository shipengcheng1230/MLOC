!*****************************************************************************************
subroutine model_init (otsp, otsh, mt)
      
! Initialize location parameters and text strings for output.
! Convert coordinates to geocentric.
! Initial estimate of hypocentroid.
! Number of free parameters for each event and total.
! Initialize bias correction term.
! Initialize ellipticity function subroutine.
      
   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, i, mt
   real :: hms2s, x1sum, x2sum, x3sum, x4sum, sxlat, cxlat, sxlon, cxlon, xlat, xlon, tcor, otsp(nevmax,0:itmax1),&
    otsh(0:itmax1), ellip, openaz, openazlim
!   real :: diflon
   character(len=60) :: depth_set1, depth_set2
!   character(len=132) :: msg
  
   openazlim = 300.
   
   if (verbose_screen) write (*,'(/a/)') 'Starting model and free parameters:'
   
   do iev = 1,nev
   
      call azgap (iev, openaz)
      if (verbose_screen) write (*,'(i3,2(2x,a),a,i3)') iev, trim(evtnam(iev)), trim(infile(iev)),&
       ' / Open azimuth = ', nint(openaz)
      
      latpr(iev) = '   '
      lonpr(iev) = '   '
      depthpr(iev) = '     '
      timepr(iev) = '           '
      
      ! Latitude free or fixed
      if (latf(iev) .and. (openaz .lt. openazlim)) then ! Free
         mindx(iev,1) = 1
         latpr(iev) = 'lat'
      else ! Fixed
         if (verbose_screen) write (*,'(t5,a,f7.3)') 'latitude fixed at ', latp(iev,0)
         mindx(iev,1) = 0
      end if
      
      ! Longitude free or fixed
      if (lonf(iev) .and. (openaz .lt. openazlim)) then ! Free
         mindx(iev,2) = 1
         lonpr(iev) = 'lon'
      else ! Fixed
         if (verbose_screen) write (*,'(t5,a,f8.3)') 'longitude fixed at ', lonp(iev,0)
         mindx(iev,2) = 0
      end if
               
!  Make sure all longitudes are given in the same range.
!      if (iev .gt. 1) then
!         diflon = lonp(iev,0) - lonp(iev-1,0)
!         if (diflon .gt. 100.) then
!            write (msg,'(a,f5.1,a,i3,1x,a,2f7.1)') 'model_init: longitude range violation: difference = ', diflon,&
!             '° for event ', iev, evtnam(iev), lonp(iev,0), lonp(iev-1,0)
!            call oops (trim(msg))
!         end if
!      end if
      
      ! Origin time free or fixed
      if (timef(iev)) then ! Free
         mindx(iev,4) = 1
         timepr(iev) = 'origin time'
      else ! Fixed
         mindx(iev,4) = 0
      end if
      otsp(iev,0) = hms2s(hourp(iev,0), minp(iev,0), secp(iev,0))
      
      ! Depth
      if (depset_pr(iev) .eq. 'c') then
         depth_set1 = ', from cluster default depth'
      else if (depset_pr(iev) .eq. 'd') then
         depth_set1 = ', from depth phases'
      else if (depset_pr(iev) .eq. 'e') then
         depth_set1 = ', from engineering information'
      else if (depset_pr(iev) .eq. 'f') then
         depth_set1 = ', from fault modeling (e.g., InSAR, GPS)'
      else if (depset_pr(iev) .eq. 'l') then
         depth_set1 = ', from local-distance readings'
      else if (depset_pr(iev) .eq. 'm') then
         depth_set1 = ', from mloc solution with free depth'
      else if (depset_pr(iev) .eq. 'n') then
         depth_set1 = ', from near-source readings'
      else if (depset_pr(iev) .eq. 'r') then
         depth_set1 = ', from relocation with free depth outside mloc'
      else if (depset_pr(iev) .eq. 'u') then
         depth_set1 = ', unconstrained'
      else if (depset_pr(iev) .eq. 'w') then
         depth_set1 = ', from waveform analysis'
      else
         depth_set1 = ', from unknown source with code '//depset_pr(iev)
      end if
      if (hdepthshift .lt. 0.01) then
         depth_set2 = ' '
      else
         depth_set2 = ', perturbed in command file'
      end if
      ! Free or fixed
      if (depthf(iev)) then ! Free
         mindx(iev,3) = 1
         depthpr(iev) = 'depth'
      else ! Fixed
         if (verbose_screen) write (*,'(t5,a,f5.1,2a)') 'depth fixed at ', depthp(iev,0), trim(depth_set1), trim(depth_set2)
         mindx(iev,3) = 0
      end if
      
   end do      

   ! Convert event coordinates to geocentric.
   do iev = 1,nev
      call geocen (latp(iev,0), lonp(iev,0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
      latpgc(iev,0) = xlat
      lonp(iev,0) = xlon
   end do

   ! Initial estimate of hypocentroid
   x1sum = 0.
   x2sum = 0.
   x3sum = 0.
   x4sum = 0.
   do iev = 1,nev
      x1sum = x1sum + latp(iev,0)
      x2sum = x2sum + lonp(iev,0)
      x3sum = x3sum + depthp(iev,0)
      x4sum = x4sum + otsp(iev,0)
   end do
   lath(0) = x1sum/real(nev)
   lonh(0) = x2sum/real(nev)
   call geocen (lath(0), lonh(0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
   lathgc(0) = xlat
   depthh(0) = x3sum/real(nev)
   otsh(0) = x4sum/real(nev)
   call timecr (otsh(0), hourh(0), minh(0), sech(0)) 

   ! Number of free parameters for each event and total.
   mt = 0
   do iev = 1,nev
      mtiev(iev) = 0
      do i = 1,4
         if (mindx(iev,i) .ne. 0) then
            mtiev(iev) = mtiev(iev) + 1
            mindx(iev,i) = mtiev(iev)
         end if
      end do
      mt = mt + mtiev(iev)
   end do

   ! Initialize bias correction term
   bcorr = 0.
   
   ! Initialize ellipticity function subroutine
   if (verbose_screen) write (*,'(/a/)') 'Initialize ellip'
   tcor = ellip ('P       ', 30., 30., 30., 30.)
   
   return
   
end subroutine model_init
   

!*****************************************************************************************
subroutine rdp ()
      
! Relative depth phases listed
      
   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, ird
   character(len=8) :: rdp_phase
   
   nrdp = 0
   
   write (io_log,'(/a)') 'Relative depth phases'
   
   do iev = 1,nev
      do ird = 1,nst(iev)
         if (rel_depth_phase(iev,ird) .gt. -100.) then
            rdp_phase = phase(iev,ird)(1:2)//'-P    '
            write (io_log,'(2i4,1x,a,1x,a,f10.2)') iev, ird, stname(iev,ird), rdp_phase, rel_depth_phase(iev,ird)
            nrdp = nrdp + 1
         end if
      end do
   end do
   
   write (io_log,'(/i4,a)') nrdp, ' relative depth phases generated'
   if (verbose_screen) write (*,'(/i4,a/)') nrdp, ' relative depth phases generated'
   
   return
   
end subroutine rdp
      
      
!*****************************************************************************************
subroutine stflt (iev, it)

! Completely rewritten 3/3/06 by eab.
! Logical variables are defined for each reading to determine if it passes criteria on
! duplicate readings, size of residual, and flags for both the cluster vector and the hypocentroid.

! .TRUE. for one of the logical variables means the reading has been filtered for that reason,
! i.e., it should not be used.

! Epicentral distance tests are based on the distance ranges defined by commands CLIM and HLIM.
! Residual size tests are based on the window defined by command WIND. In the range used for
! direct calibration (defined by windloclim) the residual limits for filtering are set explicitly
! and they are rather large because local-distance readings can have large residulas if the
! starting location is poor.
! If a reading is flagged, it is not counted for epicentral distance or residual problems.

   implicit none
   
   include 'mloc.inc'

   integer :: i, ird, iev, it
   real :: rtmp, dtmp, dts
   logical :: ptype, stype
   character(len=160) :: msg

   ndat(iev,it) = 0
   ndatdl(iev,it) = 0
   ndatpr(iev,it) = 0
   ndatfl(iev) = 0
   do ird = 1,nst(iev)
      dtmp = delt(iev,ird)
      rtmp = wind2*ttsprd(iev,ird)
      if (direct_cal .and. dtmp .le. windloclim) rtmp = rtmp * 2.
      
      ! Filter for flags
      fltrcflag(iev,ird) = .false.
      fltrhflag(iev,ird) = .false.
      if (fcode(iev,ird) .eq. 'x') then ! Outlier (large residual), either cluster or absolute.
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. 'p') then ! The phase is problematic.
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. 'd') then ! Duplicate reading.
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. 's') then ! Flagged on the basis of phase name, station and/or author by the SKIP command.
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. 'm') then ! Missing station coordinates
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. 't') then ! Timing is suspect
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      if (fcode(iev,ird) .eq. '*') then ! Engdahl's flag.
         fltrcflag(iev,ird) = .true.
         fltrhflag(iev,ird) = .true.
      end if
      
      ! Filter for epicentral distance
!     if (fltrcflag(iev,ird)) then
!        fltrcdelt(iev,ird) = .false.
!     else
!        fltrcdelt(iev,ird) = .true.
!     end if
      fltrcdelt(iev,ird) = .true.
      do i = 1,nclim
         if (dtmp .ge. clim(i,1) .and. dtmp .le. clim(i,2)) fltrcdelt(iev,ird) = .false.
      end do
!     if (fltrhflag(iev,ird)) then
!        fltrhdelt(iev,ird) = .false.
!     else
!        fltrhdelt(iev,ird) = .true.
!     end if
      fltrhdelt(iev,ird) = .true.
      do i = 1,nhlim
         if (dtmp .ge. hlim(i,1) .and. dtmp .le. hlim(i,2)) fltrhdelt(iev,ird) = .false.
      end do
      
      ! Filter for residual.
      dts = dt(iev,ird,it) - s(iev,ird,it) - ttoff(iev,ird)
      fltrcres(iev,ird) = .true.
      fltrhres(iev,ird) = .true.
      if (abs(dts) .le. rtmp) then
         fltrcres(iev,ird) = .false.
         fltrhres(iev,ird) = .false.
      end if
      
      
      fltrc(iev,ird) = (fltrcdelt(iev,ird) .or. fltrcres(iev,ird) .or. fltrcflag(iev,ird))               
      fltrh(iev,ird) = (fltrhdelt(iev,ird) .or. fltrhres(iev,ird) .or. fltrhflag(iev,ird))
      
      ! Extra filtering when using only P arrivals for hypocentroid.
      if (ponly .and. .not.fltrh(iev,ird)) then
         if (phase(iev,ird) .ne. 'P       ' .and.&
             phase(iev,ird) .ne. 'Pg      ' .and.&
             phase(iev,ird) .ne. 'Pb      ' .and.&
             phase(iev,ird) .ne. 'Pn      ' .and.&
             phase(iev,ird) .ne. 'S-P     ' .and.&
             .not.(depthfh .and. phase(iev,ird) .eq. 'pP      ') .and.&
             .not.(depthfh .and. phase(iev,ird) .eq. 'sP      ')) fltrh(iev,ird) = .true.
      end if
      
      ! Filter for differential time data (never used for hypocentroid).
      if (idiff(iev,ird) .gt. 0) fltrh(iev,ird) = .true.
      
      ! Update counters
      if (.not.fltrh(iev,ird)) ndat(iev,it) = ndat(iev,it) + 1
      if (fltrhdelt(iev,ird)) ndatdl(iev,it) = ndatdl(iev,it) + 1
      if (fltrhres(iev,ird)) then
         ndatpr(iev,it) = ndatpr(iev,it) + 1
         if (verbose_log) write (io_log,'(a,i5,2a,f8.2)') 'stflt:', iev, stname(iev,ird), phase(iev,ird), dts
      end if
      if (fltrhflag(iev,ird)) ndatfl(iev) = ndatfl(iev) + 1
      
      if (debug) then
         write (msg,'(a,i3,i5,1x,a6,1x,a8,1x,f7.2,1x,f7.2,1x,a1,1x,4l1,1x,4l1)') 'stflt: ',&
          iev, ird, stname(iev,ird), phase(iev,ird), dtmp, rtmp, fcode(iev,ird),&
          fltrc(iev,ird), fltrcdelt(iev,ird), fltrcres(iev,ird), fltrcflag(iev,ird),&
          fltrh(iev,ird), fltrhdelt(iev,ird), fltrhres(iev,ird), fltrhflag(iev,ird)
         call debugger (trim(msg))
      end if
      
   end do
   
   return
   
end subroutine stflt


!*****************************************************************************************
subroutine hdf_read ()

! To set the starting locations for a run, the final locations from a previous run can be
! read from an HDF file. This routine reads the HDF file. It will not be compatible with
! very old HDF files because of changes to the format.
   
   implicit none
   
   include 'mloc.inc'
   
   logical :: loop
   integer :: ios, nhdf, hdf_hour, hdf_min
   real :: hdf_sec, hms2s
   character(len=146) :: linein
   character(len=132) :: msg
   
   open (io_rhdf,file=rhdf_filnam,form='formatted',status='old')
   write (io_log,'(/a/a)') 'Starting locations read from:', trim(rhdf_filnam)
   nhdf = 0
   loop = .true.
   do while (loop)
      read (io_rhdf,'(a)',iostat=ios) linein
      if (ios .lt. 0) exit
      write (io_log,'(a)') linein(1:50)
      nhdf = nhdf + 1
      if (nhdf .gt. nhdfmax) then
         call warnings ('hdf_read: nhdfmax exceeded')
         nhdf = nhdf - 1
         exit
      end if
      ! From HDF format v9.8.2 on
      read (linein(1:4),'(i4)') hdf_yr4(nhdf)
      read (linein(6:7),'(i2)') hdf_mon(nhdf)
      read (linein(9:10),'(i2)') hdf_day(nhdf)
      read (linein(12:13),'(i2)') hdf_hour
      read (linein(15:16),'(i2)') hdf_min
      read (linein(18:22),'(f5.2)') hdf_sec
      read (linein(24:32),'(f9.5)') hdf_lat(nhdf)
      read (linein(34:43),'(f10.5)') hdf_lon(nhdf)
      read (linein(45:50),'(f6.2)') hdf_dep(nhdf)
      hdf_dep_code(nhdf) = linein(52:53)
      hdf_evid(nhdf) = linein(67:76)
      read (linein(106:109),'(f4.1)') hdf_dep_plus(nhdf)
      read (linein(111:114),'(f4.1)') hdf_dep_minus(nhdf)
      if (debug) then
         write (msg,'(a,i3,f6.1,1x,a2,2f6.1)') 'hdf_read: ', nhdf, hdf_dep(nhdf), hdf_dep_code(nhdf),&
          hdf_dep_plus(nhdf), hdf_dep_minus(nhdf)
         call debugger (trim(msg))
      end if
      hdf_time(nhdf) = hms2s(hdf_hour, hdf_min, hdf_sec)
   end do
   nrhdf = nhdf
   if (verbose_screen) write (*,'(/i3,a/a)') nrhdf, ' starting locations read from:', trim(rhdf_filnam)
   
   close (io_rhdf)
   
   return
   
end subroutine hdf_read


!*****************************************************************************************
subroutine read_bdp ()

! Read an optional file with a listing of stations suspected of reporting bogus depth phases.
! This information will be reflected in the type 6 ('_tt6') plot of depth phases, in which these
! suspected stations will be plotted with a smaller symbol.

   implicit none
   
   include 'mloc.inc'
   
   character(len=108) :: linein
   character(len=132) :: msg
   integer :: ios
   
   write (io_log,'(/2a)') 'Bogus Depth Phase Stations read from: ', trim(bdp_filnam)
   open (io_bdp,file=bdp_filnam,status='old')
   n_bdp = 0
   do
      read (io_bdp,'(a)',iostat=ios) linein
      if (ios .lt. 0) exit
      if (n_bdp .lt. n_bdp_max) then
         n_bdp = n_bdp + 1
         bdp_station(n_bdp) = linein
         write (io_log,'(a)') linein
      else
         write (msg,'(a,i3)') 'read_bdp: maximum number of entries read: ', n_bdp_max
         call warnings (trim(msg))
         exit
      end if
   end do
   if (verbose_screen) write (*,'(i3,2a)') n_bdp, ' entries read from ', trim(bdp_filnam)
   
   close (io_bdp)
   
   return
   
end subroutine read_bdp


