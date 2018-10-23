!***********************************************************************************************************************************
      subroutine mlocout_puke (it)
      
! Output file in PUKE format (based on the PICK format). The PUKE format is extended from the PICK format because some
! fields, for both the hypocenter line (e.g., open azimuth) and the phase reading line (e.g., defining phase), are relevant
! to both the hypocentroid and the cluster vector. Output of some parameters is based on whether calibration has been done,
! and which method of calibration was used.

! The PUKE format still does not carry information about the HD analysis as a whole. For example,
! it does not carry the open azimuth for hypocentroid or cluster vectors over all events.

      implicit none
      
      include 'mloc.inc'
      
      real deltiev(ntmax0)
      integer indx(ntmax0), i, ii, it, j, indx2(20), idum(ntmax0), iev, jj, kk
      character*100 outfil
      logical loop
      character*5 stnam
      real sortime(20)
      character(len=132) :: msg

      outfil = trim(outfile)//'.puke'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_puke: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      open (io_out,file=outfil,status='new')
      
      do iev = 1,nev ! Loop over events

         call puke1 (it, iev) ! Write hypocenter line

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
         
         do ii = 1,nst(iev)
            i = indx(ii)
            call puke2 (it, iev, i) ! Write phase reading line         
         end do
         
         write (io_out,'()') ! Blank line separates events
      
      end do
            
      close (io_out)
!      if (lg_out) close (io_outlg)
      
      return
      end


!***********************************************************************************************************************************
      subroutine puke1 (it, iev)
      
      ! Write a puke format hypocenter line.
      
      include 'mloc.inc'
      
      real azesiev(ntmax0), alpha, xl1, xl2, ddep, dot, secpr, latout, lonout, depout
      integer indx(ntmax0), i, ii, nsdfh, nsdfc, hourpr, minpr
      character*4 cal_code
      character*5 stnam
      integer it1, ntoth, ntotc, it, iev, ilast, jj, j, calibration_type
      real saaz1, saaz2, azdif, openazh, openazc, rstadelh, cstadelh, rstadelc, cstadelc, fdp, fdm
      
      it1 = it + 1
      
      ! Uncalibrated cluster, uses cluster vector confidence ellipse
      xl1 = xl1c(iev)
      xl2 = xl2c(iev)
      alpha = alphac(iev)
      ddep = sdxhatc(iev,3)
      dot = sdxhatc(iev,4)
      hourpr = hourp(iev,it1)
      minpr = minp(iev,it1)
      secpr = secp(iev,it1)
      latout = latp(iev,it1)
      lonout = lonp(iev,it1)
      depout = depthp(iev,it1)
      calibration_type = 0
      if (direct_cal) then ! Direct calibration
         xl1 = xl1dc(iev)
         xl2 = xl2dc(iev)
         alpha = alphadc(iev)
         ddep = ddepdc(iev)
         dot = dotdc(iev)
         hourpr = hourp(iev,it1)
         minpr = minp(iev,it1)
         secpr = secp(iev,it1)
         latout = latp(iev,it1)
         lonout = lonp(iev,it1)
         depout = depthp(iev,it1)
         calibration_type = 1
      end if
      if (calibration) then ! Indirect calibration trumps direct calibration
         xl1 = xl1cg(iev)
         xl2 = xl2cg(iev)
         alpha = alphacg(iev)
         ddep = sqrt(accv(iev,3,3))
         dot = sqrt(accv(iev,4,4))
         call timecr (otsp_cal(iev), hourpr, minpr, secpr)
         latout = latp_cal(iev)
         lonout = lonp_cal(iev)
         depout = depthp_cal(iev)
         calibration_type = 2
      end if
      call set_longitude_range (lonout, longitude_range)
      
      ! Calibration code
      call calibration_code (iev, calibration_type, cal_code)
      
      ! Make sure the semi-axis azimuths are between 0 and 360 degrees
      saaz1 = alpha
      if (saaz1 .lt. 0.) saaz1 = saaz1 + 360.
      saaz2 = saaz1 + 90.
      if (saaz2 .gt. 360.) saaz2 = saaz2 - 360.
      
      ! Defining phases, based on use for hypocentroid
      ntoth = 0 ! Number of defining phases
      nsdfh = 0 ! Number of stations with defining phases
      stnam = '     '
      do i = 1,nst(iev)
         if (.not.fltrh(iev,i)) then
            ntoth = ntoth + 1
            if (stname(iev,i) .ne. stnam) then
               nsdfh = nsdfh + 1
               stnam = stname(iev,i)
            end if
          end if
      end do
      
      ! Defining phases, based on use for cluster vector
      ntotc = 0 ! Number of defining phases
      nsdfc = 0 ! Number of stations with defining phases
      stnam = '     '
      do i = 1,nst(iev)
         if (connected(iev,i)) then
            ntotc = ntotc + 1
            if (stname(iev,i) .ne. stnam) then
               nsdfc = nsdfc + 1
               stnam = stname(iev,i)
            end if
          end if
      end do
      
      ! Open azimuth based on defining phases for hypocentroid
      do i = 1,nst(iev) ! Index table on azimuth for sorting
         azesiev(i) = azes(iev,i)
      end do
      call indexx (nst(iev), azesiev, indx)
      openazh = 0.
      do ii = 1,nst(iev)-1
         i = indx(ii)
         ilast = i
         if (.not.fltrh(iev,i)) then
            do jj = ii+1,nst(iev)
               j = indx(jj)
               if (.not.fltrh(iev,j)) then
                  azdif = azesiev(j) - azesiev(i)
                  openazh = max(openazh,azdif)
                  ilast = j
                  exit
               end if
            end do
         end if
      end do
      ! Gap between the largest and smallest azimuth
      do jj = 1,nst(iev)-1
         j = indx(jj)
         if (.not.fltrh(iev,j)) then
            azdif = azesiev(j) + 360. - azesiev(ilast)
            openazh = max(openazh,azdif)
            exit
         end if
      end do
      if (nsdfh .le. 1) openazh = 360. ! Catch events with 0 or 1 defining phases for hypocentroid
      
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
      
      ! Closest and furthest station, based on defining stations for hypocentroid
      rstadelh = 180.
      cstadelh = 0.
      do i = 1,nst(iev)
         if (.not.fltrh(iev,i)) then
            rstadelh = min(rstadelh,delt(iev,i))
            cstadelh = max(cstadelh,delt(iev,i))
         end if
      end do
      
      ! Closest and furthest station, based on defining stations for cluster vector
      rstadelc = 180.
      cstadelc = 0.
      do i = 1,nst(iev)
         if (connected(iev,i)) then
            rstadelc = min(rstadelc,delt(iev,i))
            cstadelc = max(cstadelc,delt(iev,i))
         end if
      end do

      ! Uncertainty of focal depth
      if (mindx(iev,3) .gt. 0) then ! Free depth solution
         fdp = ddep
         fdm = ddep
      else ! Take from default or assigned uncertainties
         fdp = depthp_plus(iev)
         fdm = depthp_minus(iev)
      end if
      
      ! Write
      write (io_out,'(a4,1x,i4,2i2.2,1x,2i2.2,f5.2,1x,f4.2,f8.3,f9.3,f6.1,2f5.1,f6.2,2(f6.1,f5.1),2(2i4,3f6.1),f4.1,a2)')& 
       cal_code,& ! GTCNU code
       iyre(iev),& ! Origin year.
       mone(iev),& ! Origin month.
       idye(iev),& ! Origin day.
       hourpr,& ! Origin hour.
       minpr,& ! Origin minute.
       secpr,& ! Origin seconds.
       dot,& ! Standard error in origin time (s).
       latout,& ! Geographic latitude.
       lonout,& ! Geographic longitude.
       depout,& ! Depth (km).
       fdp,& ! Error in depth, positive (km).
       fdm,& ! Error in depth, minus (km).
       eciev(iev,it)/ndatc(iev,it),& ! Standard error of observations in cluster analysis (s).
       saaz1,& ! Semi-axis azimuth (deg).
       xl1,& ! Semi-axis length (km).
       saaz2,& ! Semi-axis azimuth (deg).
       xl2,& ! Semi-axis length (km).
       ntoth,& ! Number of defining phases for hypocentroid.
       nsdfh,& ! Number of stations with defining phases for hypocentroid.
       openazh,& ! Open azimuth based on defining phases for hypocentroid.
       rstadelh,& ! Closest (defining) station distance (deg) for hypocentroid.
       cstadelh,& ! Furthest (defining) station distance (deg) for hypocentroid.
       ntotc,& ! Number of defining phases for cluster vector.
       nsdfc,& ! Number of stations with defining phases for cluster vector.
       openazc,& ! Open azimuth based on defining phases for cluster vector.
       rstadelc,& ! Closest (defining) station distance (deg) for cluster vector.
       cstadelc,& ! Furthest (defining) station distance (deg) for cluster vector.
       rmag(iev),& ! preferred magnitude
       mmag(iev) ! preferred magnitude scale
      
      return
      end
      

!***********************************************************************************************************************************
subroutine puke2 (it, iev, ii)

! Write a puke format phase reading line.

   implicit none

   include 'mloc.inc'

   character(len=1) :: defh, defc
   integer :: it, iev, ii
   real :: deltout, azesout, dts, stladg_geog, stlndg180

   ! Epicentral distance and azimuth from event to station
   deltout = delt(iev,ii)
   azesout = azes(iev,ii)
   if (calibration) then
      deltout = delt_cal(iev,ii)
      azesout = azes_cal(iev,ii)
   end if

   ! Defining phase for hypocentroid
   defh = 'n'
   if (.not.fltrh(iev,ii)) defh = 'y'

   ! Defining phase for cluster vector
   defc = 'n'
   if (connected(iev,ii)) defc = 'y'

   ! Travel time residual
   dts = dt(iev,ii,it) - s(iev,ii,it) ! Direct calibration or uncalibrated
   if (calibration) dts = dt_cal(iev,ii) - s_cal(iev,ii) ! Indirect calibration trumps others
   if (phase(iev,ii)(1:7) .eq. 'UNKNOWN') dts = 999.
   
   ! Convert from geocentric to geographic latitude
   call geogra (stladg(iev,ii), stladg_geog) 
   
   ! Station longitude runs from -180 to 180
   stlndg180 = stlndg(iev,ii)
   call set_longitude_range (stlndg180, 0)

   write (io_out,'(a5,1x,f8.4,1x,f9.4,1x,i5,1x,f6.2,1x,i3,1x,a8,1x,i4,2i2.2,1x,2i2.2,f6.3,1x,f6.2,1x,f8.2,1x,f8.2,1x,a8,1x,2a1)')& 
    stname(iev,ii),& ! Station code
    stladg_geog,& ! Station latitude
    stlndg180,& ! Station longitude
    nint(ahgts(iev,ii)*1.0e3),& ! Station elevation, m
    deltout,& ! Epicentral distance (deg)
    nint(azesout),& ! Azimuth, event to station
    phase(iev,ii),& ! Phase name
    iyre(iev),& ! Year
    mone(iev),& ! Month
    idye(iev),& ! Day
    ipah(iev,ii),& ! Hour
    ipam(iev,ii),& ! Minute
    pas(iev,ii),& ! Seconds
    sdread(iev,ii),& ! Reading error (s)
    tto(iev,ii),& ! Observed travel time (s)
    dts,& ! TT residual (s) or TT in the case where a travel time cannot be calculated
    readsrc(iev,ii),& ! Reading source
    defh,& ! Defining phase for hypocentroid (y/n)
    defc ! Defining phase for cluster vector (y/n)
 
   return
   
end subroutine puke2

