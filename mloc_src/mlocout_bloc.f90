!*****************************************************************************************
subroutine mlocout_bloc (it)
      
! Output file in BAYESLOC format. 

   implicit none
  
   include 'mloc.inc'
  
   integer, parameter :: real_kind = selected_real_kind(16,30)
   real(kind=real_kind) :: xot, xot_idc

   integer :: it, iev, it1, calibration_type, ret, hourpr, minpr
   real :: lon_test, ddep, dot, ddep_idc, dot_idc
   real :: alpha, xl1, xl2, secpr, latout, lonout, depout
   real :: alpha_idc, xl1_idc, xl2_idc, secpr_idc, latout_idc, lonout_idc, depout_idc
   character(len=1) :: cal2
   character(len=4) :: cal_code, fdu_pr
   character(len=100) :: outfil
   character(len=102) :: header
   character(len=132) :: msg, fmt

   outfil = trim(outfile)//'.bloc'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_bloc: opening ', trim(outfil), ' on unit ', io_bayes
      call fyi (trim(msg))
   end if
   open (io_bayes,file=outfil,status='new')

   it1 = it + 1
   fmt = '(i4,2i2.2,a1,2i2.2,a1,i2.2,f10.5,f11.5,3x,f10.3,f7.1,3x,a4,x,f20.3,f6.2)' 
   header = 'ev_id lat_mean lon_mean dist_sd depth_mean depth_sd time_mean time_sd'
   
   if (calibration) then ! Indirect calibration
      write (io_bayes,'(a)') 'Indirect calibration'
      write (io_bayes,'(a)') header
      do iev = 1,nev
      
         ! Hypocentral uncertainties
         xl1_idc = xl1cg(iev)
         xl2_idc = xl2cg(iev)
         alpha_idc = alphacg(iev)
         ddep_idc = sqrt(accv(iev,3,3))
         dot_idc = sqrt(accv(iev,4,4))
         call timecr (otsp_cal(iev), hourpr, minpr, secpr_idc)
         latout_idc = latp_cal(iev)
         lonout_idc = lonp_cal(iev)
         depout_idc = depthp_cal(iev)
         calibration_type = 2


         ! Calibration code
         call calibration_code (iev, calibration_type, cal_code)
         cal2 = cal_code(2:2)
         call uctolc (cal2,-1)  
         
         call focal_depth_uncertainty (iev, mindx(iev,3), ddep_idc, fdu_pr)    
            
         ! Unix time
         call unix_time (iyre(iev), mone(iev), idye(iev), hourp(iev,it1), minp(iev,it1), secpr_idc, xot_idc)
          
         write (io_bayes,fmt)&
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          '.',&
          hourpr,& ! Origin hour.
          minpr,& ! Origin minute.
          '.',&
          int(secpr_idc),& ! Origin seconds.
          latout_idc, & ! Geographic latitude.
          lonout_idc,& ! Geographic longitude.
          eqr(iev),& !Uncertanity in lat & long (km)
          depout_idc,& ! Final depth (km).
          fdu_pr,& ! Uncertainty in depth (km)
          xot_idc,& !Origin time in epoch (sec)
          dot_idc  ! Uncertainty in origin time (sec).
      end do
      close (io_bayes)
      return
   end if
   
   ! Output for direct calibration or uncalibrated
   if (direct_cal) then
      write (io_bayes,'(a)') 'Direct calibration'
   else
      write (io_bayes,'(a)') 'Uncalibrated'
   end if
   write (io_bayes,'(a)') header

   do iev = 1,nev
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
      else ! Uncalibrated
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
      end if
      

      ! Calibration code
      call calibration_code (iev, calibration_type, cal_code)
      cal2 = cal_code(2:2)
      call uctolc (cal2,-1)
         
      call focal_depth_uncertainty (iev, mindx(iev,3), ddep, fdu_pr)
            
      ! Unix time
      call unix_time (iyre(iev), mone(iev), idye(iev), hourp(iev,it1), minp(iev,it1), secp(iev,it1), xot)
       
      write (io_bayes,fmt)&
       iyre(iev),& ! Origin year.
       mone(iev),& ! Origin month.
       idye(iev),& ! Origin day.
       '.',&
       hourpr,& ! Origin hour.
       minpr,& ! Origin minute.
       '.',&
       int(secpr),& ! Origin seconds.
       latout,& ! Geographic latitude.
       lonout,& ! Geographic longitude.
       eqr(iev),& !Uncertanity in lat & long (km)
       depout,& ! Final depth (km).
       fdu_pr,& ! Uncertainty in depth (km)
       xot,& !Origin time in epoch (sec)
       dot ! Uncertainty in origin time (sec).

   end do

   close (io_bayes)      
   return

end subroutine mlocout_bloc


!*****************************************************************************************
subroutine focal_depth_uncertainty (iev, mindx3, ddep, fdu_pr)

   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, mindx3
   real :: fdu, ddep, min_depth_uncertainty, default_depth_uncertainty
   character(len=4) :: fdu_pr

   data min_depth_uncertainty/0.1/ ! To avoid rounding to zero on output
   data default_depth_uncertainty/5./ ! For events set to cluster default depth

   if (mindx3 .gt. 0) then ! Free depth solution
         fdu = ddep
   else ! Take from default or assigned uncertainties
      if (depthp_plus(iev) .le. 99.) then
         fdu = max(depthp_plus(iev),depthp_minus(iev))
      else
         fdu = default_depth_uncertainty ! default value for events set to cluster default depth
      end if
   end if
   write (fdu_pr,'(f4.1)') max(fdu,min_depth_uncertainty)
   
   return
   
end subroutine focal_depth_uncertainty

