!***********************************************************************************************************************************
subroutine mlocout_tomo (it, tomo_phase_in, itomo_in)
      
! Output file for use in tomography. "tomo_phase" is a phase name for which data are extracted and "itomo" is
! a flag for the type of readings included:
!  itomo = 1 Extract all readings of the specified phase
!  itomo = 2 Extract only readings which were used for the cluster vectors
!  itomo = 3 Extract empirical path anomalies

   implicit none
   
   include 'mloc.inc'
   
   integer :: hourpr, minpr
   real :: secpr, latout, lonout, depout
   common /hypocenter/ hourpr, minpr, secpr, latout, lonout, depout
   
   integer :: indx(ntmax0), i, ii, it, j, indx2(20), idum(ntmax0), it1, iev, jj, kk, itomo_in
   real :: deltiev(ntmax0), sortime(20)
   logical :: loop
   character(len=1) :: itomo_char
   character(len=5) :: stnam
   character(len=8) :: tomo_phase_in
   character(len=100) :: outfil
   character(len=132) :: msg
         
   write (itomo_char, '(i1)') itomo_in
   outfil = trim(outfile)//'.'//trim(tomo_phase_in)//'.'//itomo_char//'.tomo'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_tomo: opening ', trim(outfil), ' on unit ', io_out
      call fyi (trim(msg))
   end if
   open (io_out,file=outfil,status='new')
   
   it1 = it + 1
   
   if (itomo_in .eq. 1 .or. itomo_in .eq. 2) then
      do iev = 1,nev ! Loop over events
      
         ! Hypocenter info
         ! Uncalibrated cluster or direct calibration
         hourpr = hourp(iev,it1)
         minpr = minp(iev,it1)
         secpr = secp(iev,it1)
         latout = latp(iev,it1)
         lonout = lonp(iev,it1)
         depout = depthp(iev,it1)
         if (calibration) then ! Indirect calibration trumps direct calibration
            call timecr (otsp_cal(iev), hourpr, minpr, secpr)
            latout = latp_cal(iev)
            lonout = lonp_cal(iev)
            depout = depthp_cal(iev)
         end if
         call set_longitude_range (lonout, longitude_range)
         
         ! Index table on delta for sorting phase reading lines
         do i = 1,nst(iev)
            deltiev(i) = delt(iev,i)
         end do
         call indexx (nst(iev), deltiev, indx)
         
         ! Check that the sorting has not disrupted the time order of multiple readings from the same station
         do ii = 1,nst(iev) - 1
            i = indx(ii)
            stnam = stname(iev,i)
            loop = .true.
            j = 1
            do while (loop)
               if (stname(iev,indx(ii+j)) .eq. stnam) then
                  if ((j+ii) .lt. nst(iev)) then
                     loop = .true.
                     j = j + 1
                  else if ((j + ii) .eq. nst(iev)) then ! End of the list
                     loop = .false.
                  else
                     loop = .false.
                  end if
               else ! New station name encountered
                  loop = .false.
               end if
            end do
            if (j .gt. 1) then
               do jj = 1,j
                  kk = indx(ii + jj - 1)
                  sortime(jj) = ipah(iev,kk)*3600. + ipam(iev,kk)*60. + pas(iev,kk)
               end do
               call indexx (j, sortime, indx2)
               do jj = 1,nst(iev)
                  idum(jj) = indx(jj)
               end do
               do jj = 1,j
                  indx(ii + jj - 1) = idum(ii + indx2(jj) - 1)
               end do
            end if
         end do
         
         ! Write phase readings
         
         do ii = 1,nst(iev)
            i = indx(ii)
            if (phase(iev,i) .ne. tomo_phase_in) cycle
            if (itomo_in .eq. 1) then
!            if (fcode (iev,i) .eq. ' ') call tomo2 (it, iev, i)
               if (.not.fltrhres(iev,i) .and. .not.fltrhflag(iev,i) .and. .not.fltrcres(iev,i) .and. &
                .not.fltrcflag(iev,i)) call tomo2 (it, iev, i)
            else if (itomo_in .eq. 2) then
               if (connected(iev,i)) call tomo2 (it, iev, i)
            end if
         end do
         
      end do
   else if (itomo_in .eq. 3) then
      call tomo3 (it, tomo_phase_in)
   end if
   
   close (io_out)
   
   return
   
end subroutine mlocout_tomo
      
      
!***********************************************************************************************************************************
subroutine tomo2 (it, iev, ii)

! Write a tomo format phase reading line.
      
   implicit none
   
   include 'mloc.inc'
   
   integer :: hourpr, minpr
   real :: secpr, latout, lonout, depout, stladg_geog, dts
   common /hypocenter/ hourpr, minpr, secpr, latout, lonout, depout
   
   integer :: it, iev, ii, i, nobsere, igt
   real :: deltout
   character(len=21) :: qtest
   
   deltout = delt(iev,ii)
   if (calibration)  deltout = delt_cal(iev,ii)
   
   ! Travel time residual
   dts = dt(iev,ii,it) - s(iev,ii, it) ! Direct calibration or uncalibrated
   if (calibration) dts = dt_cal(iev,ii) - s_cal(iev,ii) ! Indirect calibration trumps others
   
   call geogra (stladg(iev,ii), stladg_geog) ! Convert from geocentric to geographic latitude
   
   ! Empirical reading error and number of observations
   qtest=stname(iev,ii)//deployment(iev,ii)//phase(iev,ii)
   nobsere = 1
   do i = 1,nqc
      if (qtest .eq. qname1(i)) then
         nobsere = indexq(i)
         exit
      end if
   end do
   
   ! Calibration level (Indirect trumps direct)
   if (calibration) then
      igt = nint(xl2cg(iev))
   else if (direct_cal) then
      igt = nint(xl2dc(iev))
   else
      igt = 99 ! No calibration
   end if
   
   write (io_out,'(f7.3,f7.2,f7.2,1x,i3,1x,i4,2i2.2,1x,2i2.2,f5.2,2x,i4,2i2.2,1x,2i2.2,f5.2,f8.2,'//&
   '2f8.3,f6.1,1x,2f8.3,f6.1,1x,a,1x,i3,1x,i3,1x,a)')&
    deltout,& ! Epicentral distance (deg)
    dts,& ! TT residual (s)
    sdread(iev,ii),& ! Reading error (s)
    nobsere,& ! Number of observations of this station-phase
    iyre(iev),& ! Reading year
    mone(iev),& ! Reading month
    idye(iev),& ! Reading day
    ipah(iev,ii),& ! Reading hour
    ipam(iev,ii),& ! Reading minute
    pas(iev,ii),& ! Reading seconds
    iyre(iev),& ! Origin year
    mone(iev),& ! Origin month
    idye(iev),& ! Origin day
    hourpr,& ! Origin hour
    minpr,& ! Origin minute
    secpr,& ! Origin seconds
    tto(iev,ii),& ! Observed travel time (s)
    latout,& ! Event latitude
    lonout,& ! Event longitude
    depout,& ! Event depth (km)
    stladg_geog,& ! Station latitude
    stlndg(iev,ii),& ! Station longitude
    ahgts(iev,ii),& ! Station elevation
    qtest,& ! Station-phase name
    iev,& ! Event number
    igt,& ! Calibration level (99 for no uncalibrated)
    trim(basename) ! Title
    
   return
   
end subroutine tomo2


!***********************************************************************************************************************************
subroutine tomo3 (it, tomo_phase_in)

! Tomography output file of empirical path anomalies. Based on subroutine mlocout_rderr.
! The full travel time associated with an empirical path anomaly is calculated by adding the
! empirical path anomaly to the theoretical travel time between the hypocentroid and the
! appropriate station. In this case, the hypocentroid is not of the entire cluster, but the
! hypocentroid of those events which contribute to this station-phase. The entries for origin
! time and reading time are ignored.
      
   implicit none
   
   include 'mloc.inc'
   
   integer, parameter :: max = 60
   
   integer :: i, j, nphase, igt, it, iev
   real :: tt(max), dtdd(max), dtdh(max), dddp(max), ttepa, usrc(2), x1sum, x2sum, x3sum
   real :: qlatdg, qlondg, delt1, t1, t2, t3, t4, dum1, dum2, dum3, dum4, dum5, dum6, lat_epa_gc, lat_epa, lon_epa, depth_epa
   character(len=8) :: tomo_phase_in, phcd(max)
   
   write (io_log,'(/a)') 'tomo3:'
   
   do i = 1,nqc ! Loop over station-phases with data
      if (indexq(i) .ge. 2 .and. idiff0(i) .eq. 0) then
         if (qname1(i)(14:21) .ne. tomo_phase_in) cycle
         
         ! Calculate hypocentroid of events contributing to this empirical path anomally
         x1sum = 0.
         x2sum = 0.
         x3sum = 0.
         do j = 1,indexq(i)
            iev = iqiev(i,j)
            x1sum = x1sum + latpgc(iev,it+1)
            x2sum = x2sum + lonp(iev,it+1)
            x3sum = x3sum + depthp(iev,it+1)
         end do
         lat_epa_gc = x1sum/real(indexq(i))
         call geogra (lat_epa_gc,lat_epa) ! Update in geographic coordinates
         lon_epa = x2sum/real(indexq(i))
         depth_epa = x3sum/real(indexq(i))
         call set_longitude_range (lon_epa, longitude_range)
         write (io_log,'(2i4,1x,a,1x,2f10.3,f6.1)') i, indexq(i), qname1(i), lat_epa, lon_epa, depth_epa
         write (io_log,'(4x,10i4)') (iqiev(i,j),j=1,indexq(i))
         
         ! Epicentral distance
         t1 = lat_epa_gc*rpd  ! Convert to geocentric radians
         t2 = lon_epa*rpd  ! Convert to geocentric radians
         t3 = qlat(i)*rpd  ! Convert to geocentric radians
         t4 = qlon(i)*rpd  ! Convert to geocentric radians
         call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, dum4, dum5, dum6, 1)
         call geogra (qlat(i), qlatdg)
         qlondg = qlon(i)
         call set_longitude_range (qlondg, 0)
         
         !  Theoretical travel-time
         call depset (depth_epa, usrc)
         if (.not.locmod) then
            call trtm (delt1, max, nphase, tt, dtdd, dtdh, dddp, phcd)
         else
            if (verbose_log) write (io_log,'(a)') 'Calling tt_mixed_model from tomo3'
            call tt_mixed_model (depth_epa, delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
!               if (delt1 .le. dlimlocmod .and. depth_epa .le. zlimlocmod) then
!                  call ttloc2 (depth_epa, delt1, 'D', nphase, tt, dtdd, dtdh, dddp, phcd, ierr, 30000)
!               else
!                  call trtm (delt1, max, nphase, tt, dtdd, dtdh, dddp, phcd)
!               end if
         end if
         if (nphase .eq. 0) then
            if (verbose_log) write (io_log,'(a,f10.3)') ' mlocsteq: no phases returned for delt1 = ', delt1
            nphase = 1
            phcd(1) = 'CRAP    '
         end if
         ttepa = 999.99
         do j = 1,nphase
            if (tomo_phase_in .eq. phcd(j)) then
               ttepa = tt(j) + epa(i)
               exit
            end if
         end do
         
         ! Calibration level (Indirect trumps direct)
         if (calibration) then
            igt = nint(s2gt)
         else if (direct_cal) then
            igt = nint(xl2h)
         else
            igt = 99 ! No calibration
         end if
         
         write (io_out,'(f7.3,f7.2,f7.2,1x,i3,1x,i4,2i2.2,1x,2i2.2,f5.2,2x,i4,2i2.2,1x,2i2.2,f5.2,f8.2,'//&
         '2f8.3,f6.1,1x,2f8.3,f6.1,1x,a,1x,i3,1x,i3,1x,a)')&
          delt1,& ! Epicentral distance (deg)
          epa(i),& ! Empirical path anomaly (s)
          ere(i),& ! Empirical reading error (s)
          indexq(i),& ! Number of observations of this station-phase
          9999,& ! Reading year (N/A for empirical path anomalies)
          99,& ! Reading month (N/A for empirical path anomalies)
          99,& ! Reading day (N/A for empirical path anomalies)
          99,& ! Reading hour (N/A for empirical path anomalies)
          99,& ! Reading minute (N/A for empirical path anomalies)
          99.99,& ! Reading seconds (N/A for empirical path anomalies)
          9999,& ! Origin year (N/A for empirical path anomalies)
          99,& ! Origin month (N/A for empirical path anomalies)
          99,& ! Origin day (N/A for empirical path anomalies)
          99,& ! Origin hour (N/A for empirical path anomalies)
          99,& ! Origin minute (N/A for empirical path anomalies)
          99.99,& ! Origin seconds (N/A for empirical path anomalies)
          ttepa,& ! Empirical path anomaly travel time (s)
          lat_epa,& ! Hypocentroid latitude
          lon_epa,& ! Hypocentroid longitude
          depth_epa,& ! Hypocentroid depth (km)
          qlatdg,& ! Station latitude
          qlondg,& ! Station longitude
          qelev(i),& ! Station elevation
          qname1(i),& ! Station-deployment-phase name
          0,& ! Event number = 0 for empirical path anomalies
          igt,& ! Hypocentroid calibration level (99 = no calibration)
          trim(basename) ! Title
         
      end if
   end do ! End of loop over station-phases with data
   
   return
   
end subroutine tomo3
