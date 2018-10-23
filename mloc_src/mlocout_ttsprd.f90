!***********************************************************************************************************************************
subroutine mlocout_ttsprd (it)

! Runs through all possible phases in the tau-p software, collects the readings,
! calculates the spread (as robust scale estimator Sn), and writes an output file.
      
   implicit none
   
   include 'mloc.inc'
   
   integer, parameter :: nphases_max = 200
   integer, parameter :: nphdat_max = 30000
   
   integer :: nph(nphases_max), nphases, it, iev, ird, j, i, jout, lunit
   real :: phdat(nphases_max,nphdat_max), d(nphdat_max), sn, dts, total, average
   logical :: newphase
   character(len=8) :: ph(nphases_max)
   character(len=100) :: outfil
   character(len=132) :: msg
   
   ph(1) = 'P       '
   nphases = 1
   nph(1) = 0
   
   do iev = 1,nev
      do ird = 1,nst(iev)
         if (idiff(iev,ird) .gt. 0) cycle ! Skip differential time data for this purpose, since absolute times are dummy values
         if (.not.fltrh(iev,ird) .or. .not.fltrc(iev,ird)) then ! Only consider data used for the cluster vectors or hypocentroid
            dts = dt(iev,ird,it) - s(iev,ird,it)
            do j = 1,nphases
               newphase = .true.
               if (phase(iev,ird)(1:7) .eq. ph(j)(1:7)) then
                  if (nph(j) .lt. nphdat_max) then
                     nph(j) = nph(j) + 1
                     phdat(j,nph(j)) = dts
                     newphase = .false.
                  else
                     write (msg,'(a,i5)') 'mlocout_ttsprd: maximum number of data reached: ', nphdat_max
                     call warnings (trim(msg))
                     if (verbose_log) write (io_log,'(a)') msg
                  end if
                  exit
               end if
            end do
            if (newphase) then
               if (nphases .lt. nphases_max) then
                  nphases = nphases + 1
                  ph(nphases) = phase(iev,ird)
                  nph(nphases) = 1
                  phdat(nphases,nph(nphases)) = dts
               else
                  write (msg,'(a,i3)') 'mlocout_ttsprd: maximum number of phases reached: ', nphases_max
                  call warnings (trim(msg))
                  if (verbose_log) write (io_log,'(a)') msg
               end if
            end if
         end if
      end do
   end do

   ! Output file
   outfil = trim(outfile)//'.ttsprd'
   jout = lunit()
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_ttsprd: opening ', trim(outfil), ' on unit ', jout
      call fyi (trim(msg))
   end if
   open (jout,file=outfil,status='new')
   do i = 1,nphases
      total = 0.
      if (nph(i) .ge. 5) then
         do j = 1,nph(i)
            d(j) = phdat(i,j)
            total = total + d(j)
         end do
         average = total/real(nph(i))
         call croux (d, min(nph(i),1000), sn)
         write (jout,'(a8,i8,2f10.3)') ph(i), nph(i), sn, average
      end if
   end do
   close (jout)
   
   return
   
end subroutine mlocout_ttsprd
      