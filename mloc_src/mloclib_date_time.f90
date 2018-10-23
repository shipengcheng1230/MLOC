!***********************************************************************************************************************************
real function hms2s (i, j, atmp)

! Convert hours-minutes-seconds to seconds

   implicit none

   integer :: i, j
   real :: atmp

   hms2s = real(3600*i + 60*j) + atmp

   return
   
end function hms2s
      

!***********************************************************************************************************************************
subroutine timecr (ots0, hour, min, sec)

!  convert number of seconds since beginning of the day (ots0) to
!  hour-minute-seconds time. if ots0 < 0, hour < 0. if ots0 > 86400,
!  hour > 23.

   implicit none

   integer :: hour, min
   real :: sec, ots0, ots, dsec
   
   ots = ots0

   dsec = ots
   hour = 0
   min = 0
   if (dsec .ge. 60.) then
      hour = 0
      min = 0
      do while (dsec .ge. 60.)
         dsec = dsec - 60.
         min = min + 1
         if (min .eq. 60) then
            min = 0
            hour = hour + 1
         end if
      end do
   else if (dsec .lt. 0.) then
      hour = -1
      min = 59
      do while (dsec .lt. -60.)
         dsec = dsec + 60.
         min = min - 1
         if (min .eq. -1) then
            min = 59
            hour = hour - 1
         end if
      end do
      dsec = 60. + dsec
   end if
 
   sec = dsec

   return
   
end subroutine timecr

!***********************************************************************************************************************************
subroutine juldat (year, month, day, julday, iflg)

! Conversion between Gregorian and day-of-year (Julian day). The conversion direction
! is controlled by input parameter "iflg". Based on code by R. Buland, re-structured
! by eab.

! iflg = 0:
!   Given the usual year, month, and day (as integer numbers) subroutine
!   juldat returns integer Julian day JULDAY.  It properly accounts for
!   the leap years through 2099.

! iflg = 1:
!   Given the integer year and julian day, subroutine gredat returns the
!   correct integer month and day of the month.  Correctly treats
!   leap years through 2099.

   implicit none

   integer :: year, day, dpm(12), feb, month, iflg, im, julday, i
   character(len=132) :: message
   
   data dpm/31,28,31,30,31,30,31,31,30,31,30,31/,feb/28/

   select case (iflg)
      case (0)
         dpm(2) = feb
         if ((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or. (mod(year,400) .eq. 0)) dpm(2) = feb + 1
         if (month .ge. 1 .and. month .le. 12 .and. day .ge. 1 .and. day .le. dpm(month)) then
            julday = day
            if (month .ge. 2) then
               im = month - 1
               do i = 1,im
                  julday = julday + dpm(i)
               end do
            end if
         else
            write (message,'(a,2i8)') 'juldat: illegal value for MONTH or DAY: ', month, day
            call oops (message)
         end if
      
      case (1)
         if (julday .ge. 1 .and. julday .le. (366 - min0(mod(year,4),1))) then
            dpm(2) = feb
            if (mod(year,4) .eq. 0) dpm(2) = feb + 1
            day = julday
            month = 1
            do i = 1,12
               if (day .le. dpm(i)) exit
               day = day - dpm(i)
               month=month + 1
            end do
         else
            write (message,'(a,i8)') 'juldat: illegal value for JULDAY: ', julday
            call oops (message)
         end if
      
      case default
         write (message,'(a,i8)') 'juldat: illegal value for IFLG: ', iflg
         call oops (message)
      
   end select

   return
   
end subroutine juldat
     
!***********************************************************************************************************************************
logical function date_range (julday_test, julday_start, julday_end)

! Test if a Julian date is within a given date range.
! Date format is YYYYDDD

   implicit none

   integer :: julday_test, julday_start, julday_end
   character(len=132) :: message

   date_range = .false.

   if (julday_start .eq. 0 .and. julday_end .eq. 0) then ! No info on date range
      date_range = .true.
   else if (julday_start .eq. 0 .and. julday_end .gt. 0) then ! No info on start of date range
      if (julday_test .le. julday_end) date_range = .true.
   else if (julday_start .gt. 0 .and. julday_end .eq. 0) then ! No info on end of date range
      if (julday_test .ge. julday_start) date_range = .true.
   else if (julday_start .gt. 0 .and. julday_end .gt. 0) then ! Both ends of date range are defined
      if (julday_test .ge. julday_start .and. julday_test .le. julday_end) date_range = .true.
   else
      write (message,'(a,3i8)') 'date_range: illegal dates ', julday_test, julday_start, julday_end
      call oops (message)
   end if

   return
   
end function date_range
      
      
!*****************************************************************************************
subroutine unix_time (year, month, day, hr, minute, sec, eot)

! Unix time converted from 'calendar' or 'human' time
! Referenced to the Unix 'epoch': 00:00:00 UTC on 1 January 1970
! Times for dates before the epoch are negative.
! Fractional seconds are supported.
! This algorithm is valid for dates from 1/1/1902.
! Although the Unix epoch is January 1, 1970, the algorithm calculates number of days to
! the target date from January 1, 1900, and then subtracts the number of days between
! January 1, 1900 and January 1, 1970.

! Written by EAB (2016/03/10)

   implicit none
   
   integer, parameter :: real_kind = selected_real_kind(16,30)
   real(kind=real_kind) :: eot
   
   integer :: dpm(12) ! days per month (non-leap-year)
   integer :: spd ! seconds per day
   integer :: offset_days ! Number of days between 1/1/1900 and 1/1/1970
   integer :: year, month, day, hr, minute
   integer :: i, year_days, month_days, day_days, whole_days, whole_days_corrected
   real :: sec
   
   data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
   data spd/86400/
   data offset_days/25568/
   
   if (year .lt. 1902) then
      write (*,'(a)') 'unix_time: invalid date'
      eot = 0.0_real_kind
      return
   end if
   
   ! Number of days in whole years up to the year of the date of interest, starting from 1900.
   ! 1900 was not a leap-year.
   year_days = 0
   if (year .ge. 1901) then
      do i = 1900,year-1
         year_days = year_days + 365
         if (mod(i,4) .eq. 0 .and. year .ne. 1900) year_days = year_days + 1 ! leap-year
      end do
   end if
   
   ! Number of days in whole months of the target year, up to the date of interest
   month_days = 0
   if (month .ge. 2) then
      do i = 1,month-1
         if (mod(year,4) .eq. 0 .and. year .ne. 1900 .and. i .eq. 2) then ! February of a leap year
            month_days = month_days + dpm(i) + 1
         else
            month_days = month_days + dpm(i)
         end if
      end do
   end if
      
   day_days = day - 1 ! Number of days up to the current day of interest  
   whole_days = year_days + month_days + day_days ! Number of whole days since 1/1/1900
   whole_days_corrected = whole_days - offset_days ! Correct for the offset from 1/1/1900 to 1/1/1970
   
   ! Convert to seconds and add the time elapsed on the current day
   eot = real(whole_days_corrected*spd,real_kind) + real(hr*3600) + real(minute*60) + sec
   
   return
   
end subroutine unix_time

