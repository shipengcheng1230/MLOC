!***********************************************************************************************************************************
subroutine mlocout_rderr ()

! Summary file of updated reading errors.
! The measure of spread is 'sn', a robust, efficient scale estimator (Croux & Rousseeuw, 1992)
! The average is also estimated, for empirical path anomalies (but these are uncalibrated).
! The last column is the reading error actually used in this run.
! This is only done for cluster mode runs.
! 11/17/2017: format changed to add deployment code
      
   implicit none
   
   include 'mloc.inc'
   
   real :: sumdata, qlatdg, qlondg, delt1, azeqst, t1, t2, t3, t4, dum1, dum2, dum3, dum5, dum6
   real :: lath_epa_g, lonh_epa_g
   
   character*100 outfil
   integer i, j, k
   real res_data(n_qres_max)
   character(len=132) :: msg
   
   outfil=trim(outfile)//'.rderr'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_rderr: opening ', trim(outfil), ' on unit ', io_out
      call fyi (trim(msg))
   end if
   open (io_out,file=outfil,status='new')
   
   ! First line carries the hypocentroid (calibrated if available) for plotting empirical path anomalies
   call geogra (lath_epa,lath_epa_g) ! hypocentroid in geographic coordinates
   lonh_epa_g = lonh_epa
   call set_longitude_range (lonh_epa_g, 0)
   write (io_out,'(a,3f10.3)') 'Hypocentroid: ', lath_epa_g, lonh_epa_g, depthh_epa

   do i = 1,nqc ! Loop over station-phases with data
      if (indexq(i) .ge. 2 .and. idiff0(i) .eq. 0) then
         sumdata = 0.
         k = 0
         do j = 1,indexq(i)
            ! On rare occasions, the calibration shift during indirect calibration causes the epicentral distance
            ! of a reading to fall off the phase branch and the residual becomes very large, upsetting the calculation of 
            ! mean residual (epa) and, if the sample size is small, the empirical reading error (ere). The test on
            ! absolute value of qres is meant to remove such readings from the calculation.
            if (abs(qres(i,j)) .lt. 20.) then
               k = k + 1
               res_data(k) = qres(i,j)
               sumdata = sumdata + res_data(k)
               if (debug) write (io_log,'(a,3x,i3,2f10.3)') qname1(i), j, qres(i,j), sumdata
            end if
         end do
         if (k .gt. 1) then
            epa(i) = sumdata/real(k)
            call croux (res_data, min(k,1000), ere(i)) ! Robust scale estimator Sn
         else
            epa(i) = 0.0
            ere(i) = 1.0
            write (msg,'(a,i1,2a)') 'mlocout_rderr: k = ', k, ' for ', qname1(i)
            call warnings (trim(msg))
         end if
         if (debug) write (io_log,'(a,3x,2f10.3)') qname1(i), epa(i), ere(i)
         t1 = lath_epa*rpd  ! Convert to geocentric radians
         t2 = lonh_epa*rpd  ! Convert to geocentric radians
         t3 = qlat(i)*rpd  ! Convert to geocentric radians
         t4 = qlon(i)*rpd  ! Convert to geocentric radians
         call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, azeqst, dum5, dum6, 1)
         call geogra (qlat(i), qlatdg)
         qlondg = qlon(i)
         call set_longitude_range (qlondg, 0)
         write (io_out,'(a,1x,i3,7f10.3)') qname1(i), k, epa(i), ere(i), rderr0(i), qlatdg, qlondg, delt1, azeqst
         if (abs(epa(i)) .gt. 10.0 .or. ere(i) .gt. 3.0) then
            write (msg,'(2a,2f10.3)') 'mlocout_rderr: possible outliers ', qname1(i), epa(i), ere(i)
            call fyi (trim(msg))
         end if
      end if
   end do ! End of loop over station-phases with data
   close (io_out)
   
   return
   
end subroutine mlocout_rderr


!***********************************************************************************************************************************
subroutine mlocout_rderr_diff ()

! Summary file of updated reading errors for differential time data. The estimate is made from
! all differential times for the same station and phase. Empirical path anomalies are not calculated.
! The measure of spread is 'sn', a robust, efficient scale estimator (Croux & Rousseeuw, 1992)
      
   implicit none
   
   include 'mloc.inc'
      
   character(len=100) :: outfil
   character(len=132) :: msg
   integer :: i, j, k
   real :: res_data(n_qres_max)
   
   outfil=trim(outfile)//'.rderr_diff'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_rderr_diff: opening ', trim(outfil), ' on unit ', io_out
      call fyi (trim(msg))
   end if
   open (io_out,file=outfil,status='new')
   
   do i = 1,nqc ! Loop over station-phases with data
      if (indexq(i) .ge. 2 .and. idiff0(i) .ne. 0) then
         k = 0
         do j = 1,indexq(i)
            ! On rare occasions, the calibration shift during indirect calibration causes the epicentral distance
            ! of a reading to fall off the phase branch and the residual becomes very large, upsetting the calculation of 
            ! the empirical reading error (ere). The test on
            ! absolute value of qres is meant to remove such readings from the calculation.
            if (abs(qres(i,j)) .lt. 20.) then
               k = k + 1
               res_data(k) = qres(i,j)
               if (debug) write (io_log,'(a,3x,i3,f10.3)') qname1(i), j, qres(i,j)
            end if
         end do
         if (k .gt. 1) then
            call croux (res_data, min(k,1000), ere(i)) ! Robust scale estimator Sn
         else
            ere(i) = 1.0
            write (msg,'(a,i1,2a)') 'mlocout_rderr_diff: k = ', k, ' for ', qname1(i)
            call warnings (trim(msg))
         end if
         if (debug) write (io_log,'(a,3x,f10.3)') qname1(i), ere(i)
         write (io_out,'(a,1x,i3,2f10.3)') qname1(i), k, ere(i), rderr0(i)
         if (ere(i) .gt. 0.5) then ! Tolerance for outliers is smaller for differential time data
            write (msg,'(2a,f10.3)') 'mlocout_rderr_diff: possible outliers ', qname1(i), ere(i)
            call fyi (trim(msg))
         end if
      end if
   end do ! End of loop over station-phases with data
   close (io_out)
   
   return
   
end subroutine mlocout_rderr_diff
