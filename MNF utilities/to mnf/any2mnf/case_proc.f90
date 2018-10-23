module case_proc

! Procedure to change case

contains

   pure function change_case (string_in, mode) result (string_out)
   
   ! Change the case of a character string, based on 'mode':
   !  mode < 0 : to lower case
   !  mode = 0 : no change
   !  mode > 0 : to upper case

       implicit none
       
       character(*), intent(in) :: string_in
       integer, intent(in) :: mode
       character(len(string_in)) :: string_out
       integer :: ic, i

       character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

       string_out = string_in
       do i = 1, len_trim(string_in)
           if (mode < 0) then ! to lower case
              ic = index(cap, string_in(i:i))
              if (ic > 0) string_out(i:i) = low(ic:ic)
           else if (mode > 0) then ! to upper case
              ic = index(low, string_in(i:i))
              if (ic > 0) string_out(i:i) = cap(ic:ic)
           end if
       end do

   end function change_case
   
end module case_proc
