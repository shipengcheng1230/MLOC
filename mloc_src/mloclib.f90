!*****************************************************************************************
subroutine calibration_code (iev, calibration_type, cal_code)

! Constructs the GTCU code for calibration level
!  calibration_type = 0 for uncalibrated
!                     1 for direct calibration
!                     2 for indirect calibration

   implicit none
   
   include 'mloc.inc'
   
   character(len=4) :: cal_code
   character(len=2) :: scale_length_pr
   logical :: epicenter_free, depth_calibrated
   integer :: iev, igt, calibration_type
   
   ! Calibrated events must have both latitude and longitude as free parameters
   if (mindx(iev,1) .gt. 0 .and. mindx(iev,2) .gt. 0) then
      epicenter_free = .true.
   else
      epicenter_free = .false.
   end if
   
   ! Depth calibration.
   ! If depth is a free parameter it is assumed to be calibrated, but certain values of
   ! depset_pr are also taken to mean a calibrated depth. Depth is also considered to be
   ! calibrated if an event is a calibration event with cal_par = 'h'.
   depth_calibrated = .false.
   if (cal_event(iev,3) .and. (cal_par(iev) .eq. 'f' .or. cal_par(iev) .eq. 'h')) depth_calibrated = .true.
   if (mindx(iev,3) .gt. 0) then ! Free depth
      depth_calibrated = .true.
   else ! Depth fixed but calibrated in some way
      if (depset_pr(iev) .eq. 'd') then ! Depth phases
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'e') then ! engineering information
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'f') then ! from fault modeling (InSAR, GPS, etc.)
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'l') then ! from local-distance readings
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'm') then ! from an mloc run with free depth
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'n') then ! from near-source readings
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'r') then ! from relocation outside of mloc with free depth
         depth_calibrated = .true.
      else if (depset_pr(iev) .eq. 'w') then ! from waveform analysis
         depth_calibrated = .true.
      end if
   end if
   
   ! Scale length ! Nearest integer to the semi-major axis length, capped at 99.
   ! Always written as a 2-digit integer.
   select case (calibration_type)
   
      case (0) ! Uncalibrated
         igt = nint(xl2c(iev))
         
      case (1) ! Direct Calibration
         igt = nint(xl2dc(iev))
         
      case (2) ! Indirect Calibration
         igt = nint(xl2cg(iev))
         
   end select
   
   if (igt .le. 99) then
      write (scale_length_pr,'(i2.2)') igt
   else
      scale_length_pr = '99'
   end if
   
   select case (calibration_type)
   
      case (0) ! Uncalibrated
         if (depth_calibrated) then
            cal_code = 'UF'//scale_length_pr
         else
            cal_code = 'UE'//scale_length_pr
         end if
      
      case (1) ! Direct calibration
         if (epicenter_free) then
            if (depth_calibrated) then
               cal_code = 'CH'//scale_length_pr
            else
               cal_code = 'CT'//scale_length_pr
            end if
         else
            if (depth_calibrated) then
               cal_code = 'UF'//scale_length_pr
            else
               cal_code = 'UE'//scale_length_pr
            end if
         end if
      
      case (2) ! Indirect calibration
         if (epicenter_free) then
            if (depth_calibrated) then
               if (ot_cal) then
                  cal_code = 'CH'//scale_length_pr
               else
                  cal_code = 'CF'//scale_length_pr
               end if
            else
               if (ot_cal .or. cal_par(iev) .eq. 't') then
                  cal_code = 'CT'//scale_length_pr
               else
                  cal_code = 'CE'//scale_length_pr
               end if
            end if
         else
            if (depth_calibrated) then
               cal_code = 'UF'//scale_length_pr
            else
               cal_code = 'UE'//scale_length_pr
            end if
         end if
      
   end select
      
   return
   
end subroutine calibration_code


!*****************************************************************************************
subroutine readerr (iev, ird)

! Sets default reading errors for phase arrivals. There is a hierarchy of methods:
! 1) A set of phase-specific default reading errors is automatically read during start-up.
! 2) If a reading is not found in that phase list, the old hard-wired algorithm is used.
! 3) Whatever value results from these first two methods is over-ridden by the value of
!    empirical reading error from a previous run (.rderr file), if it has been specified.
!    A minimum allowed value is still enforced.
! 4) Reading errors for local-distance stations may be further reset by the 'rels' command,
!    with no minimum value enforced.
! 5) Finally, a check is made to ensure that the reading error is not less than the standard
!    deviation of the uniform rectangular continuous distribution over the range defined
!    by the precision of the reading.
      
   implicit none
   
   include 'mloc.inc'
   
   integer :: iev, ird, i
   real :: readerror
   logical :: stype, ptype
   character(len=21) :: qtest
   character(len=132) :: msg
   
   ! Reading errors for differential time data are not set here.
   if (idiff(iev,ird) .gt. 0) return

   if (.not.data_weight) then ! No weighting by reading precision.
      sdread(iev,ird) = 1.0
      return
   end if
   
   ! Initialize
   sdread(iev,ird) = 0.
   
   ! See if a value can be found in the phase-specific default reading error list.
   do i = 1,n_psdre
      if (phase(iev,ird) .eq. psdre_phase(i)) then
         sdread(iev,ird) = psdre(i)
      else ! Use the old hard-wired values
   
         ! First cut
         if (phase(iev,ird) .eq. 'P       ' .or. phase(iev,ird)(1:3) .eq. 'PKP') then
            sdread(iev,ird) = 0.5
         else
            sdread(iev,ird) = 0.8
         end if
   
         ! Depth phases
         if (phase(iev,ird)(1:1) .eq. 'p') then
            sdread(iev,ird) = 1.0
         else if (phase(iev,ird)(1:1) .eq. 's') then
            sdread(iev,ird) = 1.5
         end if
   
         ! S-P
         if (phase(iev,ird)(1:3) .eq. 'S-P') sdread(iev,ird) = 0.20
     
         ! Increase reading error for phases with S-type.
         if (stype(phase(iev,ird))) sdread(iev,ird) = sdread(iev,ird)*2.0
         
         ! Lg
         if (phase(iev,ird)(1:2) .eq. 'Lg') sdread(iev,ird) = 4.0
   
         ! T-phase
         if (phase(iev,ird)(1:2) .eq. 'T ') sdread(iev,ird) = 4.0
      
      end if
   end do
   
   ! Set reading errors according to observed distribution of
   ! residuals from a previous run. This is done by searching in a summary file ("rderr" file)
   ! of station/phase residuals. A minimum allowed reading error is still enforced.
   if (read_rderr) then
      qtest = stname(iev,ird)//deployment(iev,ird)//phase(iev,ird)
      readerror = -1.
      do i = 1,nqname
         if (qtest .eq. qname(i)) then
            readerror = rdsigma(i)
            if (phase(iev,ird)(2:2) .eq. 'g' .or. phase(iev,ird)(2:2) .eq. 'b') then
               sdread(iev,ird) = amax1(readerror,rderr_min_loc)
            else if (phase(iev,ird)(1:1) .eq. 'p' .or. phase(iev,ird)(1:1) .eq. 's') then
               sdread(iev,ird) = amax1(readerror,rderr_min_depth)
            else
               sdread(iev,ird) = amax1(readerror,rderr_min)
            end if
            exit
         end if
      end do
   end if

   ! Set reading errors for crustal phases to fixed values.
   if (rels_set) then
      if ((phase(iev,ird)(2:2) .eq. 'g' .or. phase(iev,ird)(2:2) .eq. 'b') .and. delt(iev,ird) .le. rderr_loc_delt) then
         if (ptype(phase(iev,ird))) sdread(iev,ird) = rderr_loc_p
         if (stype(phase(iev,ird))) sdread(iev,ird) = rderr_loc_s
      end if
   end if
   
   ! Reading errors cannot be smaller than the standard deviation of the uniform rectangular
   ! continuous distribution over the range defined by the precision of the reading. This is
   ! irrelevant for iptim=-1 or iptim=-2, since the theoretical standard deviation is smaller
   ! than any plausible minimum reading error.
   if (iptim(iev,ird) .eq. 0) then ! Precision to nearest second.
      sdread(iev,ird) = max(sdread(iev,ird),0.29)
   else if (iptim(iev,ird) .eq. 1) then ! Precision to nearest 10 seconds.
      sdread(iev,ird) = max(sdread(iev,ird),2.9)
   else if (iptim(iev,ird) .eq. 2) then ! Precision to nearest minute.
      sdread(iev,ird) = max(sdread(iev,ird),17.0)
   else if (iptim(iev,ird) .eq. 3) then ! Precision to nearest tenth of a minute.
      sdread(iev,ird) = max(sdread(iev,ird),1.7)
   end if
   
   ! Catch anything that fell through the cracks with zero reading error
   if (sdread(iev,ird) .lt. 0.01) then
      write (msg,'(a,i3,a,1x,a,a,e12.4,a)') 'readerr: ', iev, stname(iev,ird), phase(iev,ird),&
       ': reading error of ', sdread(iev,ird), ' not allowed, set to 1.0'
      call warnings (trim(msg))
      sdread(iev,ird) = 1.0
   end if
         
   return
   
end subroutine readerr
      

!***********************************************************************************************************************************
subroutine getrderr ()

! Get empirical reading errors from a previous run (.rderr file), or use default values.
! As of v8.3, the first line of a .rderr file contains the hypocentroid used for plotting empirical path anomalies.
      
   implicit none
   
   include 'mloc.inc'
   
   integer :: i, ios
   logical :: op
   character(len=12) :: line
   
   if (read_rderr) then
      
      ! Terminal version
      open (io_rderr, file=rderrfname, status='old', form='formatted')
      
      ! MRWE version
      !open (io_rderr,file="",status='old',form='formatted')
      
      inquire (unit=io_rderr,opened=op,name=rderrfname)
      if (op) then
         ! First line is the hypocentroid used for plotting empirical path anomalies. Skip it.
         read (io_rderr,'(a)') line
         if (line(1:12) .ne. 'Hypocentroid') then
            call warnings ('getrderr: attempt to read an invalid .rderr file! Using default values.')
            close (io_rderr)
            read_rderr = .false.
            rderrfname = 'default'
            return
         end if
         do i = 1,nqmax
            read (io_rderr,'(a,t37,f10.3)',iostat=ios) qname(i), rdsigma(i)
            if (ios .lt. 0) exit
         end do
         close (io_rderr)
         nqname = i - 1
         if (verbose_screen) write (*,'(/i4,2a)') nqname, ' station-phase reading errors read from:', trim(rderrfname)
      else
         call warnings ('getrderr: failed to open file '//trim(rderrfname)//'...using default values')
         read_rderr = .false.
         rderrfname = 'default'
      end if
   else
      read_rderr = .false.
      rderrfname = 'default'
   end if
   
   return
   
end subroutine getrderr


!***********************************************************************************************************************************
integer function lenb (string)

!  finds index of last nonblank character in string

   implicit none

   integer :: i, n
   character(len=*) :: string

   n = len(string)
   do i = n,1,-1
      if (string(i:i) .ne. ' ') then
         lenb = i
         return
      end if
   end do
   lenb = 0

   return
   
end function lenb


!***********************************************************************************************************************************
integer function lennb (string)

!  finds index of first nonblank character in string

   implicit none

   integer :: i
   character(len=*) :: string

   do i = 1,len(string)
      if (string(i:i) .ne. ' ') then
         lennb = i
         return
      end if
   end do
   lennb = 1

   return
   
end function lennb


!***********************************************************************************************************************************
integer function lunit ()

! find an unopened fortran unit number between 10 and 100.
      
   implicit none
   
   logical :: lopen
   
   do lunit = 10,99
      inquire (unit=lunit,opened=lopen)
      if (.not. lopen) return
   end do
   lunit = 0
   call warnings ('lunit: failed to find an unopened unit between 10 and 99')
   
   return
   
end function lunit
      
      
!***********************************************************************************************************************************
subroutine dataflags (dflag, fcode)

! Set data flags
      
   implicit none
   
   logical :: dflag
   character(len=1) :: fcode, blank
   
   data blank/' '/

   ! My flags
   if (.not. dflag) then
      if (fcode .eq. 'x') fcode = blank ! Don't use because of large residual
      if (fcode .eq. 'p') fcode = blank ! Don't use because of the phase or lack of station information
      if (fcode .eq. 'd') fcode = blank ! Don't use because of duplicate reading
      if (fcode .eq. 's') fcode = blank ! Don't use because this station is to skipped
      if (fcode .eq. 'a') fcode = blank ! Don't use because readings by this author are to be skipped
      if (fcode .eq. 'm') fcode = blank ! Don't use because this station is missing in the station list
      if (fcode .eq. 't') fcode = blank ! Don't use because timing is suspect
   end if

   return
   
end subroutine dataflags


!***********************************************************************************************************************************
!      subroutine du
!      
!      ! Implements Istvan Bondar's dU metric for network geometry
!      
!      implicit none
!      
!      include 'mloc.inc'
!      
!      unifsum = 0
!      eazsum = 0
!      do i = 0,nsta-1
!         j = 
!         unifsum = unifsum + real(i)*360.0/real(nsta)
!         eazsum = eaz_sum + eaz(i)
!      end do
!      unif_avg = unifsum/real(nsta)
!      eaz_avg = eazsum/real(nsta)
!      
!      return
!      end


!*****************************************************************************************
Pure Function to_upper (str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lower case
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

End Function to_upper


!*****************************************************************************************
Pure Function to_lower (str) Result (string)

!   ==============================
!   Changes a string to lower case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Change each letter to lower case if it is uppercase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

End Function to_lower
