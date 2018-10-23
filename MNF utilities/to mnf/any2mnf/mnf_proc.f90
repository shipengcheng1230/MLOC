module mnf_proc

! Procedures to write MNF files

   use :: mnf_def

contains

   subroutine bulletin_line (io_unit, comment)

   ! Write an MNF bulletin line
   
      implicit none
   
      character(len=120) :: comment
      integer :: io_unit

      write (io_unit,'(2a,t121,a)') 'B   ', trim(comment), ' '

      return

   end subroutine bulletin_line


   subroutine format_line (io_unit)

   ! Write an MNF format line

      implicit none
   
      integer :: io_unit

      write (io_unit,'(2a,t121,a)') 'F   MNF v', version, ' '

      return

   end subroutine format_line


   subroutine event_line (io_unit)

   ! Write an MNF event line
   
      implicit none
      
      integer :: io_unit

      write (io_unit,'(a1,1x,a1,1x,a,t121,a)') 'E', event_usage, trim(event_annotation), ' '

      return

   end subroutine event_line


   subroutine id_line (io_unit)

   ! Write an MNF event ID line
   
      implicit none
      
      integer :: io_unit

      write (io_unit,'(a,t5,a,t12,a,t121,a)') 'I', 'NEIC', evid, ' '

      return

   end subroutine id_line


   subroutine hypocenter_line (io_unit)

   ! Write an MNF hypocenter line
   
      implicit none

      character(len=120) :: fmt
      character(len=5) :: depth_err_pr
      character(len=4) :: gtcnu
      integer :: io_unit

      fmt = '(2a,1x,i4,2(a,i2.2),1x,2(i2.2,1x),f5.2,1x,a5,1x,f9.4,f10.4,1x,a3,1x,2(a5,1x),f5.1,1x,a1,1x,2(a5,1x),a4,1x,a8,1x,a18)'

      if (depth_err .gt. 0.1) then
         write (depth_err_pr,'(f5.1)') depth_err
      else
         depth_err_pr = '     '
      end if
      gtcnu = '    '
   
      write (io_unit,fmt)&
       'H ', hypo_usage, hypo_year, '/', hypo_month, '/', hypo_day, hypo_hour, hypo_minute, hypo_seconds, ot_err,&
        latitude, longitude, iaz, smin, smaj, depth, depth_code, depth_err_pr, depth_err_pr, gtcnu, hypo_author, hypo_orid
 
      return

   end subroutine hypocenter_line


   subroutine magnitude_line (io_unit)

   ! Write an MNF magnitude line
   
      implicit none
   
      integer :: io_unit

      write (io_unit,'(2a,1x,f3.1,2x,a5,1x,a8,1x,a10,t121,a)') 'M ', magnitude_usage, magnitude, magnitude_scale,&
       magnitude_author, magnitude_orid, ' '
 
      return

   end subroutine magnitude_line


   subroutine phase_line (io_unit)

   ! Write an MNF phase line
   
      implicit none
         
      integer :: io_unit
      character(len=120) :: fmt
      character(len=27) :: adslc

      adslc = agency//'.'//deployment//'.'//station(1:5)//'.'//location//'.'//channel

      if (iptim_pr .eq. '-3') then
         fmt = '(a2,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,1x,4(i2,1x),f6.3,1x,a2,1x,a5,1x,a8,1x,a27,1x,a8,1x,a10)'
      else if (iptim_pr .eq. '-2') then
         fmt = '(a2,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,1x,4(i2,1x),f5.2,2x,a2,1x,a5,1x,a8,1x,a27,1x,a8,1x,a10)'
      else if (iptim_pr .eq. '-1') then
         fmt = '(a2,a1,1x,a6,1x,a6,1x,a3,1x,a1,a8,1x,i4,1x,4(i2,1x),f4.1,3x,a2,1x,a5,1x,a8,1x,a27,1x,a8,1x,a10)'
      end if
   
      write (io_unit,fmt)&
       'P ',& ! a1
       phase_usage,& ! a1
       station,& !a6
       distance_pr,& ! a6, formatted as f6.2
       azeq_pr,& ! a3, formatted as i3
       phase_flag,& ! a1
       phase,& ! a8
       phase_year,&
       phase_month,&
       phase_day,&
       phase_hour,& ! i2
       phase_minute,& ! i2
       phase_seconds,& ! Format depends on iptim_pr
       iptim_pr,& ! a2, formatted as i2
       residual_pr,& ! a5 formatted as f5.1
       phase,& ! a8
       adslc,& ! a27
       phase_author,&  ! a8
       arrid_pr ! a10
 
      return

   end subroutine phase_line


   logical function phase_ok (phasein)

      implicit none

      character(len=8) :: phasein, bad_phase(32)
      integer :: nbad, i

      data bad_phase /'x       ','LZ      ','LN      ','LE      ','LR      ','AML     ','AMB     ','AMS     ',&
                      'pmax    ','MLR     ','IAmb    ','SME     ','SMN     ','LmH     ','LmV     ','PMZ     ',&
                      'PFAKE   ','tx      ','smax    ','LRM     ','MAXIMUM ','LmE     ','LMZ     ','LME     ',&
                      'LMN     ','MLRZ    ','LQ      ','LRZ     ','LRN     ','LRE     ','********','********'/
      data nbad /30/

      phase_ok = .true.
      do i=1,nbad
         if (phasein .eq. bad_phase(i)) then
            phase_ok = .false.
            exit
         end if
      end do

      return

   end function phase_ok

end module mnf_proc