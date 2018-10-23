program xdat

! Version of xdat for MNF v1.1 files.

   implicit none

   integer, parameter :: line_max = 20000 ! Maximum number fo lines in an MNF file
   integer :: iev, imnf0, nline, i, ios
   character(len=20) :: infile0
   character(len=6) :: stname0, stname
   character(len=8) :: phase0, phase, readsrc0, readsrc
   character(len=110) :: filein, command_line
   character(len=121) :: buffer(line_max), line

   write (*,'(a/)') 'Release date November 7, 2017'

   write (*,'(a)',advance='no') ' XDAT file: '
   read (*,'(a)') filein
   open (8,file=filein,status='old')

   do
      read (8,'(i3,1x,a6,2(1x,a8),t31,i4,t38,a20)',iostat=ios) iev, stname0, readsrc0, phase0, imnf0, infile0
      if (ios .lt. 0) exit
      print *, iev, stname0, readsrc0, phase0, imnf0, infile0
   
      ! Do not process differential time data (imnf0 = 0)
      if (imnf0 .eq. 0) cycle
   
      ! Read all lines into buffer
      open (9,file=trim(infile0),status='old')
      nline = 0
      do
         read (9,'(a)',iostat=ios) line
         if (ios .lt. 0) exit
         nline = nline + 1
         buffer(nline) = line
         if (nline .eq. line_max) then
            write (*,'(//a,i5,a)') 'Maximum number of lines (', line_max, ') reached!'
            close (8)
            close (9)
            stop
         end if
      end do
      write (*,'(i4,2a)') nline, ' lines read from ', infile0
      close (9)
   
      ! Edit the appropriate line of the buffer, delete the original data file, write a new one from the buffer
      if (imnf0 .lt. nline) then
         buffer(imnf0)(3:3) = 'x'
         write (*,'(a)') buffer(imnf0)
         command_line = 'rm '//trim(infile0)
         call system (trim(command_line))
         open (9,file=trim(infile0),status='new')
         do i = 1,nline
            write (9,'(a)') buffer(i)
         end do
         close (9)
      else
         write (*,'(a,i3,1x,a6,2(1x,a8),1x,i4)') 'Not found: ', iev, stname0, readsrc0, phase0, imnf0
      end if
              
   end do

   close (8)
        
   stop
   
end program xdat
        