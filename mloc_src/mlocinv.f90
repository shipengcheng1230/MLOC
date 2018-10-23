!***********************************************************************************************************************************
subroutine mlocinv (it)

! generalized inverse for hypocentroid and cluster vectors.

! January 26, 1989 by eric bergman
! Converted to use double precision 6/29/95 by eab.
! Major change to the calculation of confidence ellipses for cluster vectors, 2/3/05 by eab.

   implicit none
   save
   
   include 'mloc.inc'

   double precision  ahat(ntmax1,mtmax)
   real              al, bl, kcrit
   real              alphah1
   real              alphah12
   real              alphah2
   double precision  aqbahat(mtmax,mtmax)
   double precision  awa(4,4)
   real              az
   double precision  a0(nqmax,4)
   double precision  dtemp
   double precision  dthat(ntmax1)
   double precision  dx(mtmax)
   double precision  dx01(4)
   double precision  dx02(4)
   double precision  ehath(ntmax1)
   double precision  ehati
   double precision  e2
   real              f
   integer           i
   integer           ibin
   integer           ibinq(nqmax)
   integer           icond
   integer           iev     ! Event counter
   integer           indx(mtmax)
   integer           irank(mtmax)
   integer           it
   integer           j
   integer           k
   integer           kbayes
   integer           kst
   integer           m
   integer           mf
   integer           mindxh(4)
   real              minimum_weight
   integer           mt
   integer           mtc
   integer           mth
   integer           mttotal
   integer           mt0
   integer           nf
   integer           nq ! number of distinct station-phases used
   integer           nt ! number of TT residuals used
   integer           ntiev(ntmax1) ! event number associated with each TT residual
   real              pc
   real              ph
   double precision  pqahat(ntmax1)
   double precision  pwa0(nqmax)
   double precision  qbadx(ntmax1)
   double precision  qbahat(ntmax1,mtmax)
   double precision  qbdthat(ntmax1)
   double precision  qm(mtmax)
   double precision  qmtr(mtmax)
   character*13      qhname(nqmax)
   character*21      qcname(nqmax)
   double precision  q4(4)
   double precision  shatsqh
   real              sn(nqmax)
   real              temp
   double precision  tm(mtmax)
   double precision  tmn(mtmax,ntmax1)
   double precision  tn(ntmax1)
   double precision  tn4(ntmax1,4)
   real              totimp
   real              totimpiev
   double precision  tq(nqmax)
   double precision  t22a(2,2)
   double precision  t4m(4,mtmax)
   double precision  t4m2(4,mtmax)
   double precision  t4n1(4,ntmax1)
   double precision  t4n2(4,ntmax1)
   double precision  t4q(4,nqmax)
   double precision  umm(mtmax,mtmax)
   double precision  unm(ntmax1,mtmax)
   double precision  uq4(nqmax,4)
   logical           used
   double precision  vc(mtmax)
   double precision  vh(4)
   double precision  vhath1(4,4)
   double precision  vhath2(4,4)
   double precision  vmm(mtmax,mtmax)
   double precision  v44(4,4)
   double precision  wa0(nqmax,4)
   double precision  wa0d(4,nqmax)
   real              xl1h1
   real              xl1h12
   real              xl1h2
   real              xl2h1
   real              xl2h12
   real              xl2h2
   real              t11, t12, t22
   real              dxp_sum(4)
   real              dxp_mean(4)
   double precision  tf      
   character(132):: msg
          
   ! pc and ph are the probabilities used for calculating confidence
   ! ellipses for cluster vectors and the hypocentroid, respectively.

   data pc/0.90/, ph/0.90/
   
   tf = dble(tikhonov_factor)
   
   ! cluster vectors

   write (*,'(a)') '  cluster vectors:'
   if (nev .gt. 1) then
   
      ! Initialization

      ntiev = 0
      dthat = 0.0d0
      ahat = 0.0d0
      qcname = ' '
      ntqi = 0
      ntq = 0
      do iev = 1,nevmax
         eciev(iev,it) = 0.
         ndatc(iev,it) = 0
      end do
      mt = 0
      nq = 0
      nt = 0
      
      do iev = 1,nev
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               nt = nt + 1
               if (nt .gt. ntmax1) then
                  write (msg,'(a,i5,a,i3)') 'mlocinv: NT exceeds NTMAX1 (', ntmax1, ') for event ', iev
                  call oops (trim(msg))
               end if
               ntiev(nt) = iev
               ndatc(iev,it) = ndatc(iev,it) + 1 ! Number of observations from each event used for cluster vector
               if (nt .gt. 1) then
                  used = .false.
                  kst = 0
                  do while (.not.used)
                     kst = kst + 1
                     if (kst .le. nq) then
                        if (stname(iev,j)//deployment(iev,j)//phase(iev,j) .eq. qcname(kst)&
                         .and. idiff(iev,j) .eq. idiff0(kst)) then
                           jb(nt) = kst
                           used = .true.
                        end if
                     else
                        nq = nq + 1
                        if (nq .gt. nqmax) then
                           write (msg,'(2(a,i3))') 'mlocinv: NQ exceeds NQMAX (', nqmax, ') for event ', iev
                           call oops (trim(msg))
                        end if 
                        qcname(nq) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                        qlat(nq) = stladg(iev,j) ! Geocentric coordinates
                        qlon(nq) = stlndg(iev,j) ! Geocentric coordinates
                        qelev(nq) = ahgts(iev,j)
                        idiff0(nq) = idiff(iev,j)
                        if (data_weight) then
                           rderr0(nq) = sdread(iev,j)
                        else
                           rderr0(nq) = 1.
                        end if
                        jb(nt) = nq
                        used = .true.
                     end if 
                  end do
               else
                  nq = 1
                  qcname(1) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                  qlat(1) = stladg(iev,j) ! Geocentric coordinates
                  qlon(1) = stlndg(iev,j) ! Geocentric coordinates
                  qelev(1) = ahgts(iev,j)
                  idiff0(1) = idiff(iev,j)
                  if (data_weight) then
                     rderr0(1) = sdread(iev,j)
                  else
                     rderr0(1) = 1.
                  end if
                  jb(1) = 1
               end if
               if (data_weight) then
                  vnhat(nt) = sdread(iev,j)*sdread(iev,j)
               else
                  vnhat(nt) = 1.0
               end if
               sighatj(nt) = sqrt(vnhat(nt))
               temp = 1.0/sighatj(nt)
               dthat(nt) = dble(dt(iev,j,it)*temp)
               do k = 1,4
                  if (mindx(iev,k) .ne. 0) then
                     m = mt + mindx(iev,k)
                     ahat(nt,m) = dble(a(iev,j,k)*temp)
                  end if
               end do
            end if
         end do
         mt = mt + mtiev(iev)
         if (mt .gt. mtmax) then
            write (msg,'(2(a,i3))') 'mlocinv: MT exceeds MTMAX (', mtmax, ') for event ', iev
            call oops (trim(msg))
         end if
         if (ndatc(iev,it) .lt. mtiev(iev)) then
            write (msg,'(a,i3)') 'mlocinv: Fewer data than free parameters for event ', iev
            call oops (trim(msg))         
         end if
      end do
      do j = 1,nq
         i = 0
         do k = 1,nt
            if (jb(k) .eq. j) then
               i = i + 1
               if (i .gt. ntqmax) then
                  write (msg,'(3(a,i4,3a))') 'mlocinv: max number (', ntqmax, ') of instances of station-phase ',&
                   qcname(j), ' exceeded'
                  call oops (trim(msg))
               end if
               ntq(i,j) = k
            end if
         end do
         ntqi(j) = i
      end do
      
      ntc = nt
      nqc = nq
      write (*,'(t7,a,i5,a)') ' NT = ', nt, ' (# data)'
      write (*,'(t7,a,i5,a)') ' MT = ', mt, ' (# free parameters)'
      write (*,'(t7,a,i5,a)') ' NQ = ', nq, ' (# indep. station-phases)'
      do i = 1,nqc
         qname1(i) = qcname(i)
      end do

      ! solve delxc = (qbhat * ahat)**dagger * qbhat * dthat (eq 72)

      ! diagonal elements of wq2(eq 67)

      do i = 1,nq
         wq2(i) = 0.
         do j = 1,ntqi(i)
            wq2(i) = wq2(i) + (1./vnhat(ntq(j,i)))
         end do
      end do
      
      call postq (ahat, ntmax1, mtmax, nt, mt, qbahat)
      call postqv (dthat, ntmax1, nt, qbdthat)
      write (*,'(a)') '   ...calling dsvd...'
      call dsvd (qbahat, nt, mt, ntmax1, mtmax, unm, vmm, qm)
      call dindexx (mt, qm, indx) 
      call rank (mt, indx, irank)
      do j = 1,mt
         if (irank(j) .gt. mtiev(1)) then
            qmtr(j) = ((qm(j)*qm(j)) + (tf*tf))/qm(j) ! Tikhonov regularization
            tm(j) = 1.0d0/qmtr(j)
         else
            tm(j)=0.0d0
         end if
      end do
      icond = int(qmtr(indx(mt))/qmtr(indx(mtiev(1)+1)))
      write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
      if (debug) then ! write out the singular values
         write (io_log,'(/a,e12.6)') 'Tikhonov factor :', tf
         write (io_log,'(a)') 'Singular values for the cluster vectors:'
         write (io_log,'(a)') '     j irank  indx     qm           qmtr'
         do j = 1,mt
            write (io_log,'(3i6,1x,e12.6,1x,e12.6)') j, irank(j), indx(j), qm(j), qmtr(j)
         end do
         write (io_log,'(a,i6)') 'indx(mt) = ', indx(mt)
         write (io_log,'(a,i6)') 'indx(mtiev(1)+1) = ', indx(mtiev(1)+1)
         write (io_log,'(a,i6)') 'indx(mtiev(1)) = ', indx(mtiev(1))
         write (io_log,'(2(e12.6,a),i12)') qmtr(indx(mt)), ' / ', qmtr(indx(mtiev(1)+1)), ' = ', icond
         write (io_log,'(/a)') 'tm:'
         write (io_log,'(12(e12.6,1x))') (tm(j),j=1,mt)
      end if
      do j = 1,mt
         dtemp = tm(j)
         do i = 1,mt
            vmm(i,j) = vmm(i,j)*dtemp
         end do
      end do 
      call ddot3 (vmm, unm, tmn, mtmax, mtmax, ntmax1, mtmax, mtmax, ntmax1, mt, mt, nt)
      call ddot1 (tmn, qbdthat, dx, mtmax, ntmax1, ntmax1, 1, mtmax, 1, mt, nt, 1)

      ! data importances

      do i = 1,nt
         pqahat(i) = 0.0d0
         do k = 1,mt
            pqahat(i) = pqahat(i) + (qbahat(i,k)*tmn(k,i))
         end do
      end do
      dimpc = 0.
      dimpciev = 0.
!         do ibin = 1,12
!            dimpc(ibin) = 0.
!            do iev = 1,nevmax
!               dimpciev(iev,ibin) = 0.
!            end do
!         end do
      totimp = 0.
      k = 0
      do iev = 1,nev
         totimpiev = 0.
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               k = k + 1
               dtmpc(iev,j) = sngl(pqahat(k))
               totimp = totimp + sngl(pqahat(k))
               totimpiev = totimpiev + sngl(pqahat(k))
            else
               dtmpc(iev,j) = 0.
            end if
            az = azes(iev,j)
            if (az .lt. 0.) az = az + 360.
            ibin = int(az/30.) + 1
            if (nint(az) .eq. 360 .and. ibin .eq. 13) ibin = 1
            if (ibin .lt. 1 .or. ibin .gt. 12) then
               write (msg,'(a,i3,f10.3,2i5)') 'mlocinv: illegal value for ibin = ', ibin, az, iev, j
               call warnings (trim(msg))
               write (io_log,'(a)') msg
               ibin = 1
            end if
            dimpc(ibin) = dimpc(ibin) + dtmpc(iev,j)
            dimpciev(iev,ibin) = dimpciev(iev,ibin) + dtmpc(iev,j)
         end do
         do ibin = 1,12
            dimpciev(iev,ibin) = dimpciev(iev,ibin)/totimpiev
         end do
      end do
      do ibin = 1,12
         dimpc(ibin) = dimpc(ibin)/totimp
      end do

      ! variance matrix of the estimation process:
      ! vhatc = (ahat**tr * qbhat * ahat)**dagger (eq 78)

      print  *,'    variance matrix'
      call ddot2 (ahat, qbahat, aqbahat, ntmax1, mtmax, ntmax1, mtmax, mtmax, mtmax, mt, nt, mt) 
      call dsvd (aqbahat, mt, mt, mtmax, mtmax, umm, vmm, qm)
      call dindexx (mt, qm, indx) 
      call rank (mt, indx, irank)
      icond = int(qm(indx(mt))/qm(indx(mtiev(1)+1)))
      write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
      if (debug) then ! write out the singular values
         write (io_log,'(/a)') 'Singular values for the variance matrix:'
         write (io_log,'(12(e12.6,1x))') (qm(j),j=1,mt)
         write (io_log,'(a,i6)') 'indx(mt) = ', indx(mt)
         write (io_log,'(a,i6)') 'indx(mtiev(1)+1) = ', indx(mtiev(1)+1)
         write (io_log,'(a,i6)') 'indx(mtiev(1)) = ', indx(mtiev(1))
         write (io_log,'(2(e12.6,a),i12)') qm(indx(mt)), ' / ', qm(indx(mtiev(1)+1)), ' = ', icond
      end if
      do j = 1,mt
         do i = 1,mt
            dtemp = 0.
            do k = 1,mt 
               if (irank(k) .gt. mtiev(1)) dtemp = dtemp + (umm(i,k)*umm(j,k)/qm(k))
            end do
            vhatc(i,j) = dtemp
         end do
      end do 

      ! error.
      ! ehatc = qbhat * dthat - qbhat * ahat * dx            (eq 65)
      ! shatc**2 = (ehatc*ehatc)/(nt-nq-(nev-1)mtiev)        (eq 81)

      call ddot1 (qbahat, dx, qbadx, ntmax1, mtmax, mtmax, 1, ntmax1, 1, nt, mt, 1)
      dtemp = 0.0d0
      do i = 1,nt
         ehati = qbdthat(i) - qbadx(i)
         e2 = ehati*ehati
         dtemp = dtemp + e2
         eciev(ntiev(i),it) = eciev(ntiev(i),it) + sngl(e2)
      end do
      ehatsqc(it) = sngl(dtemp)
      shatsqc = dtemp/(nt-nq-mt+mtiev(1))
      shatc(it) = sngl(dsqrt(shatsqc))

      ! cluster vector error for output

      k = 0
      do iev = 1,nev
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               k = k + 1
               eci(iev,j) = sngl(qbdthat(k) - qbadx(k))
            else
               eci(iev,j) = 0.
            end if
         end do
      end do

      ! confidence ellipses for lat-long

      ! This section is tricky! Naively taking the result in EQ 80 suggests
      !         mf=2*(nev-1)
      !         nf=nt-nq-(nev-1)*2
      ! which leads to impossibly large confidence ellipses because they scale
      ! up as the number of events in the cluster. However, we need to treat this
      ! problem as one of "marginal" confidence ellipses, as in the discussion leading
      ! to EQ 38. Therefore, it is valid to use mf=2 instead of mf=2*(nev-1). Also,
      ! the relevant number of degrees of freedom for the uncertainty of a cluster
      ! vector is the number of readings used to estimate that cluster vector, not the
      ! total number in the inversion.

      ! By the same logic, the squared error term used to scale "kcrit" is not based
      ! on the entire set of cluster vectors, but on the squared error term for the readings
      ! actually used, normalized by the number of readings.

      ! Until 2/3/05 I used mf=2 and nf=nt-nq-mf, and used shatsqc as the error term:
      !         nf=nt-nq-mf
      !         call fstat1 (mf, nf, pc, f)
      !         kcrit=mf*sngl(shatsqc)*f
      ! This is not too bad, but gives essentially the same kcrit (typically about 8.2) for each
      ! cluster vector's confidence ellipse scaling. I now think there's a better way to handle this
      ! calculation. I now use the normalized cluster error (eciev(iev,it)/ndatc(iev,it) for
      ! each event and the degrees of freedom are based on ndatc. mf still is 2 of course, as it
      ! always would be for a confidence ellipse on the epicenter. The result of this change is
      ! confidence ellipses will be a little larger than before for the poorest events (kcritc > 8.26),
      ! and smaller than before for the better located events (kcritc < 8.26).

      ! Since ndatc can be a small number in this algorithm, the argument for using a Bayesian approach
      ! re-emerges. I have chosen a conservative value of K=3, giving a s.d. of 0.4 on the estimated
      ! "normalized cluster error". 

      mttotal = 0
      mf = 2
      kbayes = 3
      do iev = 1,nev
         nf = kbayes + ndatc(iev,it) - mf
         call fstat1 (mf, nf, pc, f)
         shatsqci(iev) = (real(kbayes) + eciev(iev,it))/real(nf)
         kcritc(iev) = mf*shatsqci(iev)*f 
         if (mindx(iev,1) .eq. 0 .or. mindx(iev,2) .eq. 0) then
            write (msg,'(a,i3)') 'mlocinv: no confidence ellipse calculated for event ', iev
            call warnings (trim(msg))
            alphac(iev)=0.
            xl1c(iev)=0.
            xl2c(iev)=0.
            alic(iev) = 0.
            blic(iev) = 0.
            ccv(iev,1,1) = 0. ! used for output file .cv
            ccv(iev,1,2) = 0.
            ccv(iev,1,3) = 0.
            ccv(iev,1,4) = 0.
            ccv(iev,2,2) = 0.
            ccv(iev,2,3) = 0.
            ccv(iev,2,4) = 0.
            ccv(iev,3,3) = 0.
            ccv(iev,3,4) = 0.
            ccv(iev,4,4) = 0.
         else
            t22a(1,1) = vhatc(mttotal+1,mttotal+1)
            t22a(2,1) = -vhatc(mttotal+2,mttotal+1) ! geocentric latitude reverses sign.
            t22a(1,2) = -vhatc(mttotal+1,mttotal+2) ! geocentric latitude reverses sign.
            t22a(2,2) = vhatc(mttotal+2,mttotal+2)
            ccv(iev,1,1) = sngl(t22a(1,1)) ! used for output file .cv
            ccv(iev,1,2) = sngl(t22a(1,2))
            ccv(iev,2,1) = sngl(t22a(2,1))
            ccv(iev,2,2) = sngl(t22a(2,2))
            if (mindx(iev,3) .ne. 0) then
               ccv(iev,1,3) = sngl(-vhatc(mttotal+1,mttotal+3))
               ccv(iev,2,3) = sngl(vhatc(mttotal+2,mttotal+3))
               ccv(iev,3,3) = sngl(vhatc(mttotal+3,mttotal+3))
               ccv(iev,3,4) = sngl(vhatc(mttotal+3,mttotal+4))
            else
               ccv(iev,1,3) = 0.
               ccv(iev,2,3) = 0.
               ccv(iev,3,3) = 0.
               ccv(iev,3,4) = 0.
            end if
            if (mindx(iev,4) .ne. 0 .and. mindx(iev,3) .ne. 0) then
               ccv(iev,1,4) = sngl(-vhatc(mttotal+1,mttotal+4))
               ccv(iev,2,4) = sngl(vhatc(mttotal+2,mttotal+4))
               ccv(iev,4,4) = sngl(vhatc(mttotal+4,mttotal+4))
            else if (mindx(iev,4) .ne. 0 .and. mindx(iev,3) .eq. 0) then
               ccv(iev,1,4) = sngl(-vhatc(mttotal+1,mttotal+3))
               ccv(iev,2,4) = sngl(vhatc(mttotal+2,mttotal+3))
               ccv(iev,4,4) = sngl(vhatc(mttotal+3,mttotal+3))
            else
               ccv(iev,1,4) = 0.
               ccv(iev,2,4) = 0.
               ccv(iev,4,4) = 0.
            end if
            call delips (t22a, alphac(iev), al, bl)
            alic(iev) = 1./al
            blic(iev) = 1./bl
            xl1c(iev) = sqrt(kcritc(iev)*alic(iev))
            xl2c(iev) = sqrt(kcritc(iev)*blic(iev))

            ! Convert 90% confidence ellipse back to matrix form under the asumption of kcrit = 1 and
            ! add the uncertainty from radius_cvff, then get the new confidence ellipse.
            call ell2cv (1.0, alphac(iev), xl1c(iev), xl2c(iev), t11, t12, t22)
            t22a(1,1) = dble(t11) + dble(radius_cvff*radius_cvff)
            t22a(1,2) = dble(t12)
            t22a(2,1) = dble(t12)
            t22a(2,2) = dble(t22) + dble(radius_cvff*radius_cvff)
            call delips (t22a, alphac(iev), al, bl)
            alic(iev) = 1./al
            blic(iev) = 1./bl
            xl1c(iev) = sqrt(alic(iev))
            xl2c(iev) = sqrt(blic(iev))
            ! Get the modified raw covariance matrix back
            call ell2cv (kcritc(iev), alphac(iev), xl1c(iev), xl2c(iev), t11, t12, t22)
            ccv(iev,1,1) = t11
            ccv(iev,1,2) = t12
            ccv(iev,2,1) = t12
            ccv(iev,2,2) = t22
         end if
         mt0 = 0
         if (mindx(iev,1) .ne. 0) mt0 = mt0 + 1
         if (mindx(iev,2) .ne. 0) mt0 = mt0 + 1
         if (mindx(iev,3) .ne. 0) then
            mt0 = mt0 + 1
            ccv(iev,3,3) = sngl(vhatc(mttotal+mt0,mttotal+mt0))
         end if
         if (mindx(iev,4) .ne. 0) then
            mt0 = mt0 + 1
            ccv(iev,4,4) = sngl(vhatc(mttotal+mt0,mttotal+mt0))
         end if
         mttotal = mttotal + mtiev(iev)
      end do

      ! scale solution variance by shatsq

      do i = 1,mt
         vc(i) = shatsqc*vhatc(i,i)
      end do

      ! convert to 4-vectors.
      ! Check for non-zero mean displacement of clsuter vectors.
      ! In theory it should be zero, but in practice some datasets produce a non-zero mean
      ! in one or more parameters that is large enough to prevent convergence.

      mt0 = 0
      dxp_sum = 0.
      dxp_mean = 0.
      if (debug) write (io_log,'(/a)') 'dxp:'
      do iev = 1,nev
         do k = 1,4
            if (mindx(iev,k) .ne. 0) then
               m = mt0 + mindx(iev,k)
               dxp(iev,k,it) = sngl(dx(m))
               sdxhatc(iev,k) = sngl(dsqrt(vc(m)))
            else
               dxp(iev,k,it) = 0.0
               sdxhatc(iev,k) = 0.0
            end if
            dxp_sum(k) = dxp_sum(k) + dxp(iev,k,it)
         end do
         mt0 = mt0 + mtiev(iev)
         if (debug) then
            write (io_log,'(i4,4f10.3)') iev, (dxp(iev,k,it),k=1,4)
         end if
      end do
      do i = 1,4
         dxp_mean(i) = dxp_sum(i)/real(nev)
      end do
      write (io_log,'(a)') 'Average:'
      write (io_log,'(4x,4f10.3)') (dxp_mean(k),k=1,4)
      write (msg,'(a,4f10.3)') 'mlocinv: removing mean cluster vector changes = ', (dxp_mean(k),k=1,4)
      call fyi (trim(msg))
      
      ! Remove the mean, to aid convergence
      do iev = 1,nev
         do k = 1,4
            if (mindx(iev,k) .ne. 0) dxp(iev,k,it) = dxp(iev,k,it) - dxp_mean(k)
         end do
      end do
              
   end if

   ! hypocentroid

   write (*,'(a)') '  hypocentroid:'
   
   ! Initialization.

   ntiev = 0
   dthat = 0.0d0
   ahat = 0.0d0
   qhname = ' ' 
   ntqi = 0
   ntq = 0
   a0 = 0.0d0
   do j = 1,nevmax
      ehiev(j,it) = 0.
   end do
   
   ! Free parameters
   mth = 0
   if (latfh) then
      mth = mth + 1
      mindxh(1) = mth
   else
      mindxh(1) = 0
   end if
   if (lonfh) then
      mth = mth + 1
      mindxh(2) = mth
   else
      mindxh(2) = 0
   end if
   if (depthfh) then
      mth = mth + 1
      mindxh(3) = mth
   else
      mindxh(3) = 0
   end if
   if (timefh) then
      mth = mth + 1
      mindxh(4) = mth
   else
      mindxh(4) = 0
   end if
   
   mtc = 0
   nq = 0
   nt = 0
   if (verbose_log) write (io_log,'(t32,a,t41,a,t51,a,t60,a,t70,a,t80,a)') 'dt', 'vnhat', 'rderr', 'ttsprd', 'weight',&
    'derivatives'
   do iev = 1,nev
      if (verbose_log) write (io_log,'(a,i3)') 'Event ', iev
      do j = 1,nst(iev)
         if (.not.fltrh(iev,j)) then
            nt = nt + 1
            if (nt .gt. ntmax1) then
               write (msg,'(a,i5,a,i3)') 'mlocinv: NT exceeds NTMAX1 (', ntmax1, ') for event ', iev
               call oops (trim(msg))
            end if 
            ntiev(nt) = iev
            if (nt .gt. 1) then
               used = .false.
               kst = 0
               do while (.not.used)
                  kst = kst + 1
                  if (kst .le. nq) then
                     if (stname(iev,j)//phase(iev,j) .eq. qhname(kst)) then
                        jb(nt) = kst
                        used = .true.
                     end if
                  else
                     nq = nq + 1
                     if (nq .gt. nqmax) then
                        write (msg,'(2(a,i4))') 'mlocinv: NQ exceeds NQMAX (', nqmax, ') for event ', iev
                        call oops (trim(msg))
                     end if
                     qhname(nq) = stname(iev,j)//phase(iev,j)
                     sn(nq) = s(iev,j,it)
                     jb(nt) = nq
                     do k = 1,4
                        if (mindxh(k) .ne. 0) then
                           m = mindxh(k)
                           a0(nq,m) = dble(a(iev,j,k))
                        end if
                     end do
                     used = .true.
                  end if 
               end do
            else
               nq = 1
               qhname(1) = stname(iev,j)//phase(iev,j)
               sn(1) = s(iev,j,it)
               jb(1) = 1
               do k = 1,4
                  if (mindxh(k) .ne. 0) then
                     m = mindxh(k)
                     a0(1,m) = dble(a(iev,j,k))
                  end if
               end do
            end if
            if (data_weight) then
               ! Check for zero weight. No reading that gets to this point should have zero weight,
               ! because it will throw a divide-by-zero exception here.
               minimum_weight = 0.01
               if (weight(iev,j) .lt. minimum_weight) then
                  write (msg,'(a,e12.4,a,i3,1x,a,1x,a,1x,i5,a,f4.2)') 'mlocinv: ~zero weight ', weight(iev,j),&
                   ' for ', iev, stname(iev,j), phase(iev,j), mnf_line(iev, j), ' set to ', minimum_weight 
                  call warnings (trim(msg))
                  weight(iev,j) = minimum_weight
               end if
               if (.not.pttt) then
                  vnhat(nt) = (sdread(iev,j)**2+ttsprd(iev,j)**2)/weight(iev,j)
               else
                  vnhat(nt) = (sdread(iev,j)**2)/weight(iev,j)
               end if
            else
               vnhat(nt) = 1.0
            end if
            if (verbose_log) write (io_log,'(2a,1x,a,1x,9f10.3)') 'mlocinv: ', stname(iev,j), phase(iev,j), dt(iev,j,it),&
             vnhat(nt), sdread(iev,j), ttsprd(iev,j), weight(iev,j), a(iev,j,1), a(iev,j,2), a(iev,j,3), a(iev,j,4)
            sighatj(nt) = sqrt(vnhat(nt))
            temp = 1.0/sighatj(nt)
            dthat(nt) = dble(dt(iev,j,it)*temp)
            do k = 1,4
               if (mindx(iev,k) .ne. 0) then
                  m = mtc+mindx(iev,k)
                  ahat(nt,m) = dble(a(iev,j,k)*temp)
               end if
            end do
         else
            if (verbose_log) write (io_log,'(2a,1x,a,1x,f10.3,10x,7f10.3)') 'mlocinv: ', stname(iev,j), phase(iev,j),&
             dt(iev,j,it), sdread(iev,j), ttsprd(iev,j), weight(iev,j), a(iev,j,1), a(iev,j,2), a(iev,j,3), a(iev,j,4)
         end if
      end do
      mtc = mtc + mtiev(iev)
   end do
   do j = 1,nq
      i = 0
      do k = 1,nt
         if (jb(k) .eq. j) then
            i = i + 1
            ntq(i,j) = k
         end if
      end do
      ntqi(j) = i
   end do

   write (*,'(t7,a,i5,a)') ' NT  = ', nt, ' (# data)'
   write (*,'(t7,a,i5,a)') ' MTH = ', mth, ' (# free parameters)'
   write (*,'(t7,a,i5,a)') ' NQ  = ', nq, ' (# indep. station-phases)'
   
   if (nq .lt. mth) then
      write (msg,'(a,3i6)') 'mlocinv: number of data is less than the number of free parameters ', nt, mth, nq
      call warnings (trim(msg))
      delx0(1,it) = 0.
      delx0(2,it) = 0.
      delx0(3,it) = 0.
      delx0(4,it) = 0.
      dxp(iev,1,it) = 0.
      dxp(iev,2,it) = 0.
      dxp(iev,3,it) = 0.
      dxp(iev,4,it) = 0.
      sdxhath(1) = 0.
      sdxhath(2) = 0.
      sdxhath(3) = 0.
      sdxhath(4) = 0.
      return
   end if

   ! diagonal elements of wq and wq2(eq 67)

   do i = 1,nq
      temp = 0.
      do j = 1,ntqi(i)
         temp = temp + (1.0/vnhat(ntq(j,i)))
      end do
      wq2(i) = temp
      wq(i) = sqrt(temp)
   end do

   ! form w*a0 and (w*a0)**dagger = v * q**-1 * u**tr

   do j = 1,mth
      do i = 1,nq
         wa0(i,j) = dble(wq(i))*a0(i,j)
      end do
   end do 
   write (*,'(a)') '   ...calling dsvd...'
   call dsvd (wa0, nq, mth, nqmax, 4, uq4, v44, q4)
   call dindexx (mth, q4, indx) 
   icond = int(q4(indx(mth))/q4(indx(1)))
   write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
   do j = 1,mth
      dtemp = 1.0d0/q4(j)
      do i = 1,mth
         v44(i,j) = v44(i,j)*dtemp
      end do
   end do
   call ddot3 (v44, uq4, wa0d, 4, 4, nqmax, 4, 4, nqmax, mth, mth, nq)  

   ! data importances 
   ! pwa0 = wa0 * wa0d    (eq 102)

   do i = 1,nq
      dtemp = 0.0d0
      do j = 1,mth
         dtemp = dtemp + (wa0(i,j)*wa0d(j,i))
      end do
      pwa0(i) = dtemp
   end do
   do i = 1,nq
      ibinq(i) = 0
   end do
   do ibin = 1,12
      dimph(ibin) = 0.
   end do
   k = 0
   do iev = 1,nev
      do j = 1,nst(iev)
         if (.not.fltrh(iev,j)) then
            k = k + 1
            dtmph(iev,j) = sngl(pwa0(jb(k)))
            az = azes(iev,j)
            if (az .lt. 0.) az = az + 360.
            ibinq(jb(k)) = int(az/30.) + 1
         else
            dtmph(iev,j) = 0.
         end if
      end do
   end do
   totimp = 0.
   do i = 1,nq
      if (ibinq(i) .gt. 12) ibinq(i) = 12
      if (ibinq(i) .lt. 1) ibinq(i) = 1
      dimph(ibinq(i)) = dimph(ibinq(i)) + sngl(pwa0(i))
      totimp = totimp + sngl(pwa0(i))
   end do
   do ibin = 1,12
      dimph(ibin) = dimph(ibin)/totimp
   end do 

   ! dx01 = (w*a0)**dagger * {w**-1 * bhat**tr * dthat - w*sn} (eqns 83, 101)

   do i = 1,nq
      dtemp = 0.
      do j = 1,ntqi(i)
         k = ntq(j,i)
         dtemp = dtemp + dthat(k)/dble((sighatj(k)*wq(i)))
      end do 
      tq(i) = dtemp - dble((wq(i)*sn(i)))
   end do
   call ddot1 (wa0d, tq, dx01, 4, nqmax, nqmax, 1, 4, 1, mth, nq, 1)
   write (io_log,'(a,4f10.3)') ' weighted epicentroid: ', (dx01(i),i=1,mth)

   ! dx02 = (w*a0)**dagger *{w**-1 * bhat**tr * ahat * dx} (eqns 84, 101)

   call ddot1 (ahat, dx, tn, ntmax1, mtmax, mtmax, 1, ntmax1, 1, nt, mtc, 1)
   do i = 1,nq
      dtemp = 0.0d0
      do j = 1,ntqi(i)
         k = ntq(j,i)
         dtemp = dtemp + tn(k)/dble((sighatj(k)*wq(i)))
      end do 
      tq(i) = dtemp
   end do
   call ddot1 (wa0d, tq, dx02, 4, nqmax, nqmax, 1, 4, 1, mth, nq, 1)
   write (io_log,'(a,4f10.3)') '      bias correction: ', (dx02(i),i=1,mth)

   ! vhath1 = (a0**tr * wq2 * a0)**-1 (eq 90)

   do j = 1,mth
      do i = 1,mth
         dtemp = 0.0d0
         do k = 1,nq
            dtemp = dtemp + (a0(k,i)*dble(wq2(k))*a0(k,j))
         end do
         awa(i,j) = dtemp
      end do
   end do
   call dmatinv2 (awa, vhath1, 4, 4, mth, mth)

   ! vhath2 = (vhath1 * a0**tr * bhat**tr) * ahat * vhatc * ahat**tr *
   !          (vhath1 * a0**tr * bhat**tr)**tr (eq 90, 101)

   call ddot3 (vhath1, a0, t4q, 4, 4, nqmax, 4, 4, nqmax, mth, mth, nq)
   do j = 1,nt
      dtemp = 1.0d0/dble(sighatj(j))
      do i = 1,mth
         t4n1(i,j) = t4q(i,jb(j))*dtemp
      end do
   end do
   call ddot1 (t4n1, ahat, t4m, 4, ntmax1, ntmax1, mtmax, 4, mtmax, mth, nt, mtc)
   call ddot1 (t4m, vhatc, t4m2, 4, mtmax, mtmax, mtmax, 4, mtmax, mth, mtc, mtc)
   call ddot3 (t4m2, ahat, t4n2, 4, mtmax, ntmax1, mtmax, 4, ntmax1, mth, mtc, nt)
   call ddot3 (t4n2, t4n1, vhath2, 4, ntmax1, 4, ntmax1, 4, 4, mth, nt, mth)
   
   ! error
   ! ehath = (pbhat * dthat) - (bhat * sn) - (bhat * a0 * dx01)    (eq 92)
   ! shath = sqrt(ehath**2/(nq-mth))               (eq 96)

   call postpv (dthat, ntmax1, nt, ehath)
   do i = 1,nt
      ehath(i) = ehath(i) - dble(sn(jb(i))/sighatj(i))
   end do
   do j = 1,mth
      do i = 1,nt
         tn4(i,j) = a0(jb(i),j)/dble(sighatj(i))
      end do
   end do
   call ddot1 (tn4, dx01, tn, ntmax1, 4, 4, 1, ntmax1, 1, nt, mth, 1)
   dtemp = 0.0d0
   do i = 1,nt
      ehath(i) = ehath(i) - tn(i)
      e2 = ehath(i)*ehath(i)
      dtemp = dtemp + e2
      ehiev(ntiev(i),it) = ehiev(ntiev(i),it) + sngl(e2)
   end do
   ehatsqh(it) = sngl(dtemp)
   nqmth(it) = nq-mth
   shatsqh = dtemp/(nq-mth)
   shath(it) = sngl(dsqrt(shatsqh))

   ! Confidence ellipses for bias-corrected epicentroid (eq. 100),
   ! weighted epicentroid (eq. 97), and bias correction term (eq. 98).

   write (io_log,'(i3,a)') int(ph*100), '% confidence ellipses:'
   if (mindxh(1) .eq. 0 .or. mindxh(2) .eq. 0) then
      call warnings ('mlocinv: no confidence ellipse calculated for epicentroid')
      alphah = 0.
      xl1h = 0.
      xl2h = 0.
   else
      mf = 2
      nf = max(nt-(nev*mf), 2) ! To prevent negative nf when using a small amount of local readings
      call fstat1 (mf, nf, ph, f)
      t22a(1,1) = shatsqh*vhath1(1,1) + shatsqc*vhath2(1,1)
      t22a(2,1) = -(shatsqh*vhath1(2,1) + shatsqc*vhath2(2,1)) ! geocentric latitude reverses sign.
      t22a(1,2) = -(shatsqh*vhath1(1,2) + shatsqc*vhath2(1,2)) ! geocentric latitude reverses sign.
      t22a(2,2) = shatsqh*vhath1(2,2) + shatsqc*vhath2(2,2)
      call delips (t22a, alphah12, al, bl)
      kcrit  =mf*f
      xl1h12 = sqrt(kcrit/al)
      xl2h12 = sqrt(kcrit/bl)
      write (io_log,'(a,f6.1,a,2f10.3)') ' Bias-corrected: alpha = ', alphah12, ' axes = ', xl1h12, xl2h12
      if (bias_corr) then
         hcv(1,1) = sngl(t22a(1,1)) ! used for output file .cv
         hcv(1,2) = sngl(t22a(1,2))
         hcv(2,2) = sngl(t22a(2,2))
         if (mindxh(3) .ne. 0) then
            hcv(1,3) = sngl(-(shatsqh*vhath1(1,3) + shatsqc*vhath2(1,3)))
            hcv(2,3) = sngl(shatsqh*vhath1(2,3) + shatsqc*vhath2(2,3))
            hcv(3,3) = sngl(shatsqh*vhath1(3,3) + shatsqc*vhath2(3,3))
            hcv(3,4) = sngl(shatsqh*vhath1(3,4) + shatsqc*vhath2(3,4))
         else
            hcv(1,3) = 0.
            hcv(2,3) = 0.
            hcv(3,3) = 0.
            hcv(3,4) = 0.
         end if
         if (mindxh(4) .ne. 0 .and. mindxh(3) .ne. 0) then
            hcv(1,4) = sngl(-(shatsqh*vhath1(1,4) + shatsqc*vhath2(1,4)))
            hcv(2,4) = sngl(shatsqh*vhath1(2,4) + shatsqc*vhath2(2,4))
            hcv(4,4) = sngl(shatsqh*vhath1(4,4) + shatsqc*vhath2(4,4))
         else if (mindxh(4) .ne. 0 .and. mindxh(3) .eq. 0) then
            hcv(1,4) = sngl(-(shatsqh*vhath1(1,3) + shatsqc*vhath2(1,3)))
            hcv(2,4) = sngl(shatsqh*vhath1(2,3) + shatsqc*vhath2(2,3))
            hcv(4,4) = sngl(shatsqh*vhath1(3,3) + shatsqc*vhath2(3,3))
         else
            hcv(1,4) = 0.
            hcv(2,4) = 0.
            hcv(4,4) = 0.
         end if
      end if
      nf = nq - mf
      call fstat1 (mf, nf, ph, f)
      t22a(1,1) = vhath1(1,1)
      t22a(2,1) = -vhath1(2,1) ! geocentric latitude reverses sign.
      t22a(1,2) = -vhath1(1,2) ! geocentric latitude reverses sign.
      t22a(2,2) = vhath1(2,2)
      call delips (t22a, alphah1, al, bl)
      kcrit = mf*sngl(shatsqh)*f
      xl1h1 = sqrt(kcrit/al)
      xl2h1 = sqrt(kcrit/bl)
      write (io_log,'(a,f6.1,a,2f10.3)') '       Weighted: alpha = ', alphah1, ' axes = ', xl1h1, xl2h1
      if (.not.bias_corr) then
         hcv(1,1) = sngl(t22a(1,1)) ! used for output file .cv
         hcv(1,2) = sngl(t22a(1,2))
         hcv(2,2) = sngl(t22a(2,2))
         if (mindxh(3) .ne. 0) then
            hcv(1,3) = sngl(-vhath1(1,3))
            hcv(2,3) = sngl(vhath1(2,3))
            hcv(3,3) = sngl(vhath1(3,3))
            hcv(3,4) = sngl(vhath1(3,4))
         else
            hcv(1,3) = 0.
            hcv(2,3) = 0.
            hcv(3,3) = 0.
            hcv(3,4) = 0.
         end if
         if (mindxh(4) .ne. 0 .and. mindxh(3) .ne. 0) then
            hcv(1,4) = sngl(-vhath1(1,4))
            hcv(2,4) = sngl(vhath1(2,4))
            hcv(4,4) = sngl(vhath1(4,4))
         else if (mindxh(4) .ne. 0 .and. mindxh(3) .eq. 0) then
            hcv(1,4) = sngl(-vhath1(1,3))
            hcv(2,4) = sngl(vhath1(2,3))
            hcv(4,4) = sngl(vhath1(3,3))
         else
            hcv(1,4) = 0.
            hcv(2,4) = 0.
            hcv(4,4) = 0.
         end if
      end if         
      if (nev .gt. 1) then
         nf = max(nt-nq-(nev-1)*mf, 2) ! To prevent negative nf when using a small amount of local readings
         call fstat1 (mf, nf, ph, f)
         t22a(1,1) = vhath2(1,1)
         t22a(2,1) = -vhath2(2,1) ! geocentric latitude reverses sign.
         t22a(1,2) = -vhath2(1,2) ! geocentric latitude reverses sign.
         t22a(2,2) = vhath2(2,2)
         call delips (t22a, alphah2, al, bl)
         kcrit = mf*sngl(shatsqc)*f
         xl1h2 = sqrt(kcrit/al)
         xl2h2 = sqrt(kcrit/bl)
         write (io_log,'(a,f6.1,a,2f10.3)') '     Correction: alpha = ', alphah2, ' axes = ', xl1h2, xl2h2
      end if
      if (bias_corr) then
         alphah = alphah12
         xl1h = xl1h12
         xl2h = xl2h12
      else
         alphah = alphah1
         xl1h = xl1h1
         xl2h = xl2h1
      end if
   end if

   ! Convert to 4-vector.  
   ! Keep total of bias-correction terms (dx02).
   do k = 1,4
      if (mindxh(k) .ne. 0) then
         m = mindxh(k)
         delx0(k,it) = sngl(dx01(m) - dx02(m))
         bcorr(k) = bcorr(k) + sngl(dx02(m))
         if (bias_corr) then
            vh(m) = shatsqh*vhath1(m,m) + shatsqc*vhath2(m,m)
         else
            vh(m) = shatsqh*vhath1(m,m)
         end if 
         sdxhath(k) = sngl(dsqrt(vh(m)))
      else
         delx0(k,it) = 0.
         sdxhath(k) = 0.
      end if
   end do

   return
   
end subroutine mlocinv

