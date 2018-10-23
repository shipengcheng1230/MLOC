!*****************************************************************************************
subroutine mlocout_tt (it)

! Outout of empirical TT data for specific phases (command "ttou")

   implicit none
   
   include 'mloc.inc'

   integer :: it, i, j, ii, jj, iev, indx(ntmax0)
   real :: deltiev(ntmax0), dts, delt_pr, azes_pr
   character(len=132) :: file_folder, tt_file, msg
   character(len=160) :: command_line, tt_file_path
   
   ! Make a folder for the files
   file_folder = trim(datadir)//dirsym//trim(basename)//'_tt'
   command_line = 'mkdir '//trim(file_folder)
   call system (trim(command_line))
      
   ! Loop over requested phases
   do j = 1, n_ttou
      if (j .gt. 1) then
         do jj = 1, j-1
            if (ttou_phase(jj) .eq. ttou_phase(j)) then
               msg = 'mlocout_tt: repeated phase (,'//trim(ttou_phase(j))//')'
               call warnings (trim(msg))
               cycle
            end if
         end do
      end if
      tt_file = trim(basename)//'_tt_'//trim(ttou_phase(j))//'.dat'
      tt_file_path = trim(file_folder)//dirsym//trim(tt_file)
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_tt: opening ', trim(tt_file_path), ' on unit ', io_tt
         call fyi (trim(msg))
      end if
      open (io_tt,file=tt_file_path,status='new')
      
      if (ttou_phase(j) .eq. 'S-P     ') then ! S-P is handled separately
         call mlocout_sp (it)
      else
         ! Loop over events
         do iev = 1, nev
            ! Index table on delta for sorting
            do i = 1,nst(iev)
               if (calibration) then
                  deltiev(i) = delt_cal(iev,i)
               else
                  deltiev(i) = delt(iev,i)
               end if
            end do
            call indexx (nst(iev), deltiev, indx)
            ! Loop over readings, sorted by delta
            do ii = 1,nst(iev)
               i = indx(ii)
               if (calibration) then
                  dts = dt_cal(iev,i) - s_cal(iev,i)
                  delt_pr = delt_cal(iev,i)
                  azes_pr = azes_cal(iev,i)
               else
                  dts = dt(iev,i,it) - s(iev,i,it)
                  delt_pr = delt(iev,i)
                  azes_pr = azes(iev,i)
               end if
               if (phase(iev,i) .eq. ttou_phase(j))&
                write (io_tt,'(i3,1x,a16,1x,f5.1,1x,i5,1x,a20,1x,f8.3,1x,i3,1x,3f10.2,3x,a)')&
                iev, evtnam(iev)(1:16), depthp(iev,it), mnf_line(iev,i), sad(iev,i), delt_pr,&
                nint(azes_pr), tto(iev,i), dts, sdread(iev,i), fcode(iev,i)
            end do
         end do
      end if
      
      close (io_tt)
   end do
   
   return

end subroutine mlocout_tt


!*****************************************************************************************
subroutine mlocout_sp (it)

! Output file with S-P data

   implicit none
   
   include 'mloc.inc'

   character(len=132) :: msg
   integer :: iev, it, ird, nevsp, nsp, nspx, lastiev
      
   nevsp = 0
   nsp = 0
   nspx = 0
   lastiev = 0
      
   do iev = 1,nev ! Loop over events
      do ird = 1,nst(iev) ! Loop over readings
         if (phase(iev,ird) .eq. 'S-P     ') then
            call sp_write (iev, ird, it)
            if (iev .gt. lastiev) then
               nevsp = nevsp + 1
               lastiev = iev
            end if
            nsp = nsp + 1
            if (fcode (iev,ird) .ne. ' ') nspx = nspx + 1
         end if
      end do
   end do
      
   write (msg,'(a,i3,a)') 'mlocout_sp: ', nevsp, ' events with S-P readings'
   call fyi (trim(msg))
   write (io_log,'(/a)') trim(msg)
   write (msg,'(a,i3,a,i3,a)') 'mlocout_sp: ', nsp, ' S-P readings (', nspx, ' flagged)'
   call fyi (trim(msg))
   write (io_log,'(a)') trim(msg)
   
   return
   
end subroutine mlocout_sp


!*****************************************************************************************
subroutine sp_write (iev, ird, it)

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, ird, it, it1, hourpr, minpr
   real lonout, latout, depout, scale_length, xl1, xl2, dot, secpr, deltout, azesout, azseout,&
    dts, stladg_geog, stlndg180, fdp, fdm, ddep
   character(len=164) :: fmt
   
   it1 = it + 1
   
   ! Uncalibrated cluster, uses cluster vector confidence ellipse
   xl1 = xl1c(iev)
   xl2 = xl2c(iev)
   ddep = sdxhatc(iev,3)
   dot = sdxhatc(iev,4)
   hourpr = hourp(iev,it1)
   minpr = minp(iev,it1)
   secpr = secp(iev,it1)
   latout = latp(iev,it1)
   lonout = lonp(iev,it1)
   depout = depthp(iev,it1)
   if (direct_cal) then ! Direct calibration
      xl1 = xl1dc(iev)
      xl2 = xl2dc(iev)
      ddep = ddepdc(iev)
      dot = dotdc(iev)
      hourpr = hourp(iev,it1)
      minpr = minp(iev,it1)
      secpr = secp(iev,it1)
      latout = latp(iev,it1)
      lonout = lonp(iev,it1)
      depout = depthp(iev,it1)
   end if
   if (calibration) then ! Indirect calibration trumps direct calibration
      xl1 = xl1cg(iev)
      xl2 = xl2cg(iev)
      ddep = sqrt(accv(iev,3,3))
      dot = sqrt(accv(iev,4,4))
      call timecr (otsp_cal(iev), hourpr, minpr, secpr)
      latout = latp_cal(iev)
      lonout = lonp_cal(iev)
      depout = depthp_cal(iev)
   end if
   call set_longitude_range (lonout, longitude_range)
      
   ! Epicentral distance and azimuth from event to station
   deltout = delt(iev,ird)
   azesout = azes(iev,ird)
   azseout = azse(iev,ird)
   if (calibration) then
      deltout = delt_cal(iev,ird)
      azesout = azes_cal(iev,ird)
      azseout = azse_cal(iev,ird)
   end if   

   ! Scale length ! Semi-major axis of the confidence ellipse
   if (direct_cal) then ! Direct calibration
      scale_length = xl2dc(iev)
   else if (calibration) then ! Indirect calibration
      scale_length = xl2cg(iev)
   else
      scale_length = xl2c(iev) ! Uncalibrated
   end if
   scale_length = scale_length / 111.2 ! Convert km to degrees

   ! Travel time residual
   dts = dt(iev,ird,it) - s(iev,ird,it) ! Direct calibration or uncalibrated
   if (calibration) dts = dt_cal(iev,ird) - s_cal(iev,ird) ! Indirect calibration trumps others

   ! Convert from geocentric to geographic latitude
   call geogra (stladg(iev,ird), stladg_geog) 
   
   ! Station longitude runs from -180 to 180
   stlndg180 = stlndg(iev,ird)
   call set_longitude_range (stlndg180, 0)

   ! Uncertainty of focal depth
   if (mindx(iev,3) .gt. 0) then ! Free depth solution
      fdp = ddep
      fdm = ddep
   else ! Take from default or assigned uncertainties
      fdp = depthp_plus(iev)
      fdm = depthp_minus(iev)
   end if
   
   fmt = '(i3,1x,a5,1x,f6.2,1x,f5.2,1x,f6.2,f6.3,2f6.1,1x,f5.3,1x,i4,2i2.2,1x,2(i2,1x),f4.1,1x,&
    &f4.2,f9.4,f10.4,f6.1,2f5.1,f9.4,f10.4,f6.3,1x,a1)'
   write (io_tt,fmt)&
    iev,& ! Event number 
    stname(iev,ird),& ! Station code
    tto(iev,ird),& ! S-P time, s
    sdread(iev,ird),& ! Reading error (s)
    dts,& ! TT residual (s)
    deltout,& ! Epicentral distance (deg)
    azseout,& ! Back-azimuth, station to event
    azesout,& ! Azimuth, event to station
    scale_length,& ! Uncertainty in epicentral distance, degrees
    iyre(iev),& ! Origin year.
    mone(iev),& ! Origin month.
    idye(iev),& ! Origin day.
    hourpr,& ! Origin hour.
    minpr,& ! Origin minute.
    secpr,& ! Origin seconds.
    dot,& ! Standard error in origin time (s).
    latout,& ! Event latitude
    lonout,& ! Event longitude
    depout,& ! Depth
    fdp,& ! Error in depth, positive (km).
    fdm,& ! Error in depth, minus (km).
    stladg_geog,& ! Station latitude
    stlndg180,& ! Station longitude
    ahgts(iev,ird),& ! Station elevation, km
    fcode(iev,ird) ! flag

   return
   
end subroutine sp_write
