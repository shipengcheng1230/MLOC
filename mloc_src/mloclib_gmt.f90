!*****************************************************************************************
subroutine ttdat (phtest, i1, i2, dincrement, vred, hp)
   
! vred = inverse reduction velocity, in seconds per degree.
! i1 and i2 are indices for epicentral distance, which is discretized at "dincrement" deg.
! hp is the depth for which travel times should be calculated, only needed for custom crust model

   implicit none
   
   include 'mloc.inc'

   integer max
   parameter (max=60)
   real tt(max), dtdd(max), dtdh(max), dddp(max), hp, tt_lg
   integer nphase
   character*8 phcd(max), phtest
   logical stype
   
   integer i, j, i1, i2
   real d, ttred, vred, dincrement
   
   do i = i1,i2
      d = real(i) * dincrement
      if (.not.locmod) then
         call trtm (d, max, nphase, tt, dtdd, dtdh, dddp, phcd)
      else
         if (verbose_log) write (io_log,'(a)') 'Calling tt_mixed_model from ttdat'
         call tt_mixed_model (hp, d, nphase, tt, dtdd, dtdh, dddp, phcd)
      end if
      do j = 1,nphase
         if (phcd(j) .eq. phtest) then
            ttred = tt(j) - vred*d
            write (io_gmt,'(2f10.3)') ttred, d
            exit
         end if
      end do
      
      ! S-P relative phase
      if (phtest(1:3) .eq. 'S-P') then
         ! Find the first S phase
         do j = 1,nphase
            if (stype(phcd(j))) exit
         end do
         ttred = tt(j) - tt(1) - vred*d
         write (io_gmt,'(2f10.3)') ttred, d
      end if
      
      ! Lg travel time
      if (phtest .eq. 'Lg      ') then
         tt_lg = lg_a + d*lg_b
         ttred = tt_lg - vred*d
         write (io_gmt,'(2f10.3)') ttred, d
      end if
      
   end do
   
   return
   
end subroutine ttdat


!*****************************************************************************************
subroutine map_boundaries (it1, iplot, xlonmin, xlonmax, ylatmin, ylatmax)
   
! Creates a string for defining the map region in the GMT plot script
   
   implicit none
   
   include 'mloc.inc'

   integer :: it1, iplot, iev
   real :: ylatmin, ylatmax, xlonmin, xlonmax, lon_test
   
   xlonmin = 360.
   xlonmax = -360.
   ylatmin = 90.
   ylatmax = -90.
   do iev = 1,nev
      if (plot(iev,iplot)) then
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         if (latp(iev,it1) .lt. ylatmin) ylatmin = latp(iev,it1)
         if (latp(iev,it1) .gt. ylatmax) ylatmax = latp(iev,it1)
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         if (lon_test .lt. xlonmin) xlonmin = lon_test
         if (lon_test .gt. xlonmax) xlonmax = lon_test
         if (latp(iev,0) .lt. ylatmin) ylatmin = latp(iev,0)
         if (latp(iev,0) .gt. ylatmax) ylatmax = latp(iev,0)
         lon_test = lonp(iev,0)
         call set_longitude_range (lon_test, longitude_range)
         if (lon_test .lt. xlonmin) xlonmin = lon_test
         if (lon_test .gt. xlonmax) xlonmax = lon_test
         if (calibration) then
            if (latp_cal(iev) .lt. ylatmin) ylatmin = latp_cal(iev)
            if (latp_cal(iev) .gt. ylatmax) ylatmax = latp_cal(iev)
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            if (lon_test .lt. xlonmin) xlonmin = lon_test
            if (lon_test .gt. xlonmax) xlonmax = lon_test
         end if
      end if
   end do

   ylatmin = ylatmin - 0.1
   ylatmax = ylatmax + 0.1
   xlonmin = xlonmin - 0.1
   xlonmax = xlonmax + 0.1
   
   return
   
end subroutine map_boundaries


!*****************************************************************************************
subroutine run_gmt_script (gmt_script, psfile)

! Standard processing of GMT scripts for plotting:
   
! Execute the GMT script (in mloc_working directory)
! Convert the postscript file to PDF
! Move the PDF file back into the cluster directory
! Delete the postscript file, or move certain ones into the _comcat folder if ComCat output is being created.
   
   implicit none
   
   include 'mloc.inc'
   
   character(len=160) :: command_line
   character(len=132) :: pdf_file
   character*(*) gmt_script, psfile
   integer :: lps
   
   lps = len_trim(psfile)
   pdf_file = psfile(1:lps-2)//'pdf'
   
   command_line = 'chmod +x '//gmt_script
   call system (trim(command_line))
   command_line = './'//gmt_script
   call system (trim(command_line))
   command_line = 'gmt psconvert -Tf -A1 '//trim(psfile)
   call system (trim(command_line))
   command_line = 'mv '//trim(pdf_file)//' '//trim(datadir)//dirsym//trim(pdf_file)
   call system (trim(command_line))
   
   ! If ComCat output is requested some postscript files are moved to the _comcat folder
   ! after conversion to PDF files. Otherwise, they are deleted.  
   if (comcatout) then
      if (index(psfile, '_epa_') .ne. 0) then
         call system ('rm '//trim(psfile))
      else if (index(psfile, '_rdp') .ne. 0) then ! Catches summary plot and individual event plots
         call system ('rm '//trim(psfile))   
      else if (index(psfile, '_tt5e_local') .ne. 0) then
         call system ('rm '//trim(psfile))   
      else if (index(psfile, '_tt5s_local') .ne. 0) then
         call system ('rm '//trim(psfile))   
      else 
         command_line = 'mv '//trim(psfile)//' '//trim(ccat_folder)//dirsym//trim(psfile)
         call system (trim(command_line))
      end if
   else
      call system ('rm '//trim(psfile))
   end if
   
   return
   
end subroutine run_gmt_script


!*****************************************************************************************
logical function pass (fcode)
   
   character*1 fcode
   
   pass = (fcode .ne. 'x' .and. &
           fcode .ne. 'd' .and. &
           fcode .ne. 'p' .and. &
           fcode .ne. 's' .and. &
           fcode .ne. 'a' .and. &
           fcode .ne. 'm' .and. &
           fcode .ne. 't')
   
   return
   
end function pass


!*****************************************************************************************
subroutine map_dat (it)

! Creates a file with filename suffix '.map_dat' that can be read by a GMT script to plot the epicenters. For example:

!    awk '{print $2, $1, 0.10}' $datfile | psxy $projection $region -Sc -W1 -G255/0/0 -O -K >> $psfile

! where the variables "datfile", "projection", "region" and "psfile" are defined previously. The file created by this subroutine
! constitutes the "datfile". It can be concatenated with similar files from other clusters to make a regional map.
! The coordinates are given as lat-lon, so the awk/print command is used to reverse the order for GMT. The size of the symbol (0.10) is
! given in the awk command also. "-Sc -W1 -G255/0/0" defines a red filled circle. The data file also contains event number, date,
! origin time and cluster sequence number to aid in interpretation, especially when multiple clusters are combined for plotting.
! The awk command only pulls out the first two columns.

   implicit none
   
   include 'mloc.inc'
   
   character(len=132) :: outfil, msg
   integer :: iev, it, it1
   
   it1 = it + 1
   
   outfil = trim(outfile)//'.map_dat'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'map_dat: opening ', trim(outfil), ' on unit ', io_gmt
      call fyi (trim(msg))
   end if
   open (io_gmt,file=outfil,status='new')
   
   do iev = 1,nev
      if (calibration) then
         write (io_gmt,'(2f10.4,1x,i3,1x,a,1x,a)') latp_cal(iev), lonp_cal(iev), iev, evtnam(iev)(1:16), trim(basename)
      else
         write (io_gmt,'(2f10.4,1x,i3,1x,a,1x,a)') latp(iev,it1), lonp(iev,it1), iev, evtnam(iev)(1:16), trim(basename)
      end if
   end do
   
   close (io_gmt)
   
   return
   
end subroutine map_dat
