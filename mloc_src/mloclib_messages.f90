!***********************************************************************************************************************************
subroutine oops (msg)

! Error report and stop

   implicit none

   character(len=*) :: msg

   write (*,'(/a)') 'Oops! Fatal error in '//trim(msg)

   stop
   
end subroutine oops

      
!***********************************************************************************************************************************
subroutine warnings (msg)

! Warnings

   implicit none

   include 'mloc.inc'

   character(len=*) :: msg
      
   write (*,'(t3,a)') 'Warning from '//trim(msg)
      
   return
   
end subroutine warnings

      
!***********************************************************************************************************************************
subroutine fyi (msg)

! Informational messages

   implicit none
      
   character(len=*) :: msg
      
   write (*,'(t3,a)') 'FYI from '//trim(msg)
      
   return
   
end subroutine fyi
      
!***********************************************************************************************************************************
subroutine debugger (msg)

! Debugging information

   implicit none
      
   include 'mloc.inc'
   
   character(len=*) :: msg
      
   write (io_log,'(t3,a)') 'Debugging '//trim(msg)
      
   return
   
end subroutine debugger


!*****************************************************************************************
subroutine logit (msg)

! Write to the log file. This is needed for routines that do not have access to the io_log
! unit number through mloc.inc

   implicit none

   include 'mloc.inc'
   
   character(len=*) :: msg
   
   write (io_log,'(a)') trim(msg)

   return
   
end subroutine logit


