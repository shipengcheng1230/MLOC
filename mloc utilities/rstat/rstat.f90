!************************************************************************
      program rstat

! This program analyzes the cluster vector residuals from program
! 'mloc'.  The individual residuals are read from the 'phase_data'
! file.  Every instance of a specified station and phase is listed.

! This version works with mloc6, which uses an index (idiff) to keep track
! of differential time data.
! 11/18/2017: updated to work with phase_data files from v10.4.0
! 6/29/2018: Fixed a bug reading the phase_data file when differential time data are used

   implicit none
 
   logical loop, good, diff, case_sens, use_deployment, match1, match2
   integer, parameter :: ndiffmax = 1500
   integer iter, i, j, n, iev, idiff0, idiff, ios, nres(ndiffmax), dot
   real data(1000), sum, mean, sn, dts, res(ndiffmax,4), iev_idiff(ndiffmax,2), resout, xmean
   character(len=1) :: char1
   character(len=5) :: sta0
   character(len=8) :: pha0, deployment0
   character(len=14) :: sd
   character(len=60) :: datafile
   character(len=80) :: file2
   character(len=165) :: charin
   
   write (*,'(/a/)') 'Release date June 29, 2018'

   write (*,'(a)',advance='no') 'Enter file name for phase_data: '
   read (*,'(a)') file2
   open (2, file=file2, status='old', form='formatted')
   write (*,'(a)',advance='no') 'Enter number of iterations: '
   read (*,*) iter
   write (*,'(a)',advance='no') 'Is there differential time data? y or n: '
   read (*,'(a)') char1
   diff = (char1 .eq. 'y')
   write (*,'(a)',advance='no') 'Case sensitive station names? y or n: '
   read (*,'(a)') char1
   case_sens = (char1 .eq. 'y')
   write (*,'(a)',advance='no') 'Use deployment codes? y or n: '
   read (*,'(a)') char1
   use_deployment = (char1 .eq. 'y')
   
   ! Pre-calculate the differential time residuals
   if (diff) then
      nres = 0
      res = 0.
      do
         read (2,'(a)',iostat=ios) charin
         if (ios .lt. 0) exit
         if (charin(1:8) .eq. ' CLUSTER') read (charin(16:18),'(i3)') iev
         if (charin(162:165) .ne. '****' .and. charin(162:165) .ne. '    ' .and. charin(162:165) .ne. 'IDIF') then
            read (charin(162:165),'(i4)') idiff
            if (idiff .gt. 0 .and. idiff .le. ndiffmax) then
               nres(idiff) = nres(idiff) + 1
               if (nres(idiff) .eq. 3) then
                  write (*,'(a)') 'Error: more than two samples for a given idiff:'
                  write (*,'(a)') charin
                  stop
               end if
               if (iter .eq. 0) then
                  read (charin(69:75),'(f7.2)') res(idiff,nres(idiff))
               else if (iter .eq. 1) then
                  read (charin(76:82),'(f7.2)') res(idiff,nres(idiff))
               else if (iter .eq. 2) then
                  read (charin(83:89),'(f7.2)') res(idiff,nres(idiff))
               else if (iter .eq. 3) then
                  read (charin(90:96),'(f7.2)') res(idiff,nres(idiff))
               else if (iter .eq. 4) then
                  read (charin(97:103),'(f7.2)') res(idiff,nres(idiff))
               end if
               iev_idiff(idiff,nres(idiff)) = iev
               if (nres(idiff) .eq. 2) then
                  xmean = (res(idiff,1) + res(idiff,2))*0.5
                  res(idiff,3) = res(idiff,1) - xmean
                  res(idiff,4) = res(idiff,2) - xmean
                  print *, idiff, iev_idiff(idiff,1), iev_idiff(idiff,2), res(idiff,1), res(idiff,2), res(idiff,3), res(idiff,4)
               end if
            end if
         end if
      end do
   end if

   do
      if (use_deployment) then
         write (*,'(/a)',advance='no') 'Enter station code "dot" deployment code (no blanks, "q.q" to quit): '
         read (*,'(a)') sd
         if (sd .eq. 'q.q') then
            close (2)
            stop
         end if
         dot = scan (sd, '.')
         i = dot -1
         j = dot + 1
         if (dot .ge. 4 .and. dot .le. 6) then
            sta0 = sd(1:i)
            deployment0 = sd(j:len(trim(sd)))
         else
            write (*,'(/a,i2/)') 'Error, invalid dot position: ', dot
            cycle
         end if
      else
         write (*,'(/a)',advance='no') 'Enter station name (q to quit): '
         read (*,'(a)') sta0
         deployment0 = ' '
         if (sta0 .eq. 'q') then
            close (2)
            stop
         end if
      end if
      
      ! Convert to upper case station name
      if (.not.case_sens .and. ichar(sta0(1:1)) .ge. 97) then
         do i=1,len(trim(sta0))
            sta0(i:i)=char(ichar(sta0(i:i))-32)
         end do
      end if
            
      write (*,'(a)',advance='no') 'Enter phase name (* for all phases): '
      read (*,'(a)') pha0
      idiff0 = 0
      if (diff) then
         write (*,'(a)',advance='no') 'Enter idiff0 (-1 for all differential times): '
         read (*,*) idiff0
      end if
   
      rewind 2
      write (*,'(/t30,a,t46,a,t52,a,t59,a)') 'rderr  delta', 'dts', 'wgt', 'eci'
      good = .true.
      n = 0
      do
         read (2,'(a)',iostat=ios) charin
         if (ios .lt. 0) exit
         if (charin(65:67) .eq. 'BAD') good=.false.
         if (charin(1:8) .eq. ' CLUSTER') then
            read (charin(16:18),'(i3)') iev
            read (2,'(1x,a)') datafile
         end if
         if (charin(162:165) .eq. '****' .or. charin(162:165) .eq. '    ' .or. charin(162:165) .eq. 'IDIF') then
            cycle
         else
            read (charin(162:165),'(i4)') idiff
         end if
         match1 = (charin(2:6) .eq. sta0 .and. (charin(19:26) .eq. pha0 .or. pha0(1:1) .eq. '*')) ! Match station and phase
         match2 = (charin(8:15) .eq. deployment0 .or. .not.use_deployment) ! Match deployment
         if (diff .and. (idiff .gt. 0)) then
            if (match1 .and. match2) then
               if (idiff .eq. idiff0 .or. idiff0 .lt. 0) then
                  if (iev_idiff(idiff,1) .eq. iev) resout = res(idiff,3)
                  if (iev_idiff(idiff,2) .eq. iev) resout = res(idiff,4)
                  if (good) then
                     n = n + 1
                     data(n) = resout
                  end if
                  write (*,'(i3,a,f7.2,2x,3a,2x,a)') iev, charin(1:38), resout, charin(53:56), charin(126:160),&
                   charin(162:165), trim(datafile)
               end if
            end if
         else if (.not.diff .or. idiff0 .eq. 0) then
            if (match1 .and. match2) then
               read (charin(162:165),'(i4)') idiff
               if (idiff .gt. 0) cycle
               if (iter .eq. 0) then
                  read (charin(69:75),'(f7.2)') dts
               else if (iter .eq. 1) then
                  read (charin(76:82),'(f7.2)') dts
               else if (iter .eq. 2) then
                  read (charin(83:89),'(f7.2)') dts
               else if (iter .eq. 3) then
                  read (charin(90:96),'(f7.2)') dts
               else if (iter .eq. 4) then
                  read (charin(97:103),'(f7.2)') dts
               end if 
               if (good) then
                  n = n + 1
                  data(n) = dts
               end if
               if (charin(17:17) .ne. 'd') then ! Don't print duplicates
                  write (*,'(i3,a,f7.2,2x,2a,2x,a)') iev, charin(1:38), dts, charin(53:56), charin(126:160), trim(datafile)
               end if
            end if
         end if
      end do
   
      if (n .ge. 2) then
         sum = 0.
         do i = 1,n
            sum = sum + data(i)
         end do
         mean = sum/float(n)
         write (*,'(a,f6.3)') 'Mean = ', mean
         call croux (data, n, sn)
         write (*,'(a,f6.3)') 'Sn = ', sn
         write (*,'(a,i4,a)') 'On ', n, ' readings.'
      else
         write (*,'(a,i1)') 'n = ', n
      end if
   
   end do

end

!***********************************************************************
      subroutine croux (x, nin, sn)
      
      ! Implementation of the "naive" algorithm for Sn in:
      ! "Time Efficient algorithms for two highly robust estimators of scale" by Croux & Rousseuw
      ! Computational Statistics, V1 (1992), Dodge and Whittaker, ed., Physica-Verlag, Heidleberg, pp. 411-428.
      ! Sn is a robust estimator for the spread of a sample distribution that does not need an estimate
      ! of central location. It is well-behaved even with small sample size.
      
      ! The original formulation behaves badly with n=3 and two values close to each other. 
      ! the trivial difference (i = j) always yeilds a zero value and if one of the other
      ! differences is also small, the estimate of Sn implodes. After consulting
      ! Christophe Croux, I altered the algorithm for n=3 so that the inner loop takes the
      ! average instead of the lomed of the three differences.
      
      ! It follows that the constant cn(3) should be recalculated. I did the same experiment as
      ! reported in Croux & Rousseuw with the new algorithm and found cn(3) = 1.172.
      
      implicit none
      
      integer nmax
      parameter (nmax=1000) ! Maximum size of input array "x"
      
      integer i, j, n, nin, ihimed, ilomed, indx(nmax)
      dimension x(nin)
      real a1(nmax), a2(nmax), sn, cn, x
      external cn
      
      if (nin .le. nmax .and. nin .ge. 2) then
         n = nin
      else if (nin .lt. 2) then
         write (*,'(a,i6)') 'croux: illegal value for n: ', nin
         sn = 1.0
         return   
      else if (nin .gt. nmax) then
         write (*,'(a,i6)') 'croux: nin set to ', nmax
         n = nmax
      end if
      
      ! Equation 1
      do i = 1,n
         do j = 1,n
            a1(j) = abs(x(i) - x(j))
         end do
         call indexx(n, a1, indx)
         ihimed = indx((n/2)+1)
         a2(i) = a1(ihimed)
         if (n .eq. 3) a2(i) = (a1(1)+a1(2)+a1(3))/3.
      end do
      
      call indexx(n, a2, indx)
      ilomed = indx((n+1)/2)
      sn = cn(n)*1.1926*a2(ilomed)
      
      return
      end
      
      
!***********************************************************************
      real function cn (n)
      
      ! Small sample correction terms for subroutine croux
      
      integer n
      
      if (n .eq. 2) then
         cn = 0.743
      else if (n .eq. 3) then
         ! cn = 1.851   
         cn = 1.172 ! Special correction for modified algorithm using average of differences
      else if (n .eq. 4) then
         cn = 0.954
      else if (n .eq. 5) then
         cn = 1.351
      else if (n .eq. 6) then
         cn = 0.993
      else if (n .eq. 7) then
         cn = 1.198
      else if (n .eq. 8) then
         cn = 1.005
      else if (n .eq. 9) then
         cn = 1.131
      else if (n .ge. 10) then
         if (mod(n,2) .eq. 0.) then ! n even
            cn = 1.
         else ! n odd
            cn = float(n)/(float(n)-0.9)
         end if
      else
         cn = 1.
      end if
      
      return
      end


!***********************************************************************
      subroutine indexx (n, arrin, indx)

! indexes array 'arrin' of length 'n', i.e., outputs the array 'indx'
! such that arrin(indx(j)) is in ascending order for j=1,2,...,n.  the
! input quantities n and arrin are not changed.

      dimension arrin(n), indx(n)
 
      do j=1,n
         indx(j)=j
      end do
      if (n .eq. 1) return
      l=n/2+1
      ir=n
   10 continue
      if (l .gt. 1) then
         l=l-1
         indxt=indx(l)
         q=arrin(indxt)
      else
         indxt=indx(ir)
         q=arrin(indxt)
         indx(ir)=indx(1)
         ir=ir-1
         if (ir .eq. 1) then
            indx(1)=indxt
            return
         end if
      end if
      i=l
      j=l+l
   20 if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i)=indx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         go to 20 
      end if 
      indx(i)=indxt
      go to 10

      end
