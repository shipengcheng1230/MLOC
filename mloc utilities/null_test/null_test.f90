program null_test

   character(len=4) :: prefix
   
   prefix = char(0)
   write (*,'(a)') char(0)//'basemap'
   prefix = 'gmt '
   write (*,'(a)') 'gmt '//'basemap'
   
   stop
end program null_test
