! Subroutines to calculate travel times at local distances through a user-specified crustal model.
! Based on Johannes Schweitzer's hyposat code.

!***********************************************************************************************************************************
subroutine rd_loc_mod (iunit, rmaxout, zmaxout, ierr)

! Reads crustal velocity model
! First line starts with epicentral distance limit to which this model applies (f4.1). The remainder of the first line
! can be used for a comment about the source of the model.

   implicit none
   
   include 'hyposat_loc.inc'

   character(len=34) :: line
   character(len=132) :: msg
   integer :: i, iunit, ios, ierr
!   integer :: j
   real :: zmaxout, rmaxout

   rmax = 0.
   zmax = 0.
   ierr = 0
   
   read (iunit,'(a)') line
   read (line(1:4),'(f4.1)') rmax
   rmaxout = rmax
   
   ! Read model layers
   do i = 1,max_layers
      if (i .eq. max_layers) then
         write (msg,'(a,i2)') 'rd_loc_mod: maximum number of model layers reached: ', max_layers
         call warnings (trim(msg))
      end if
      read (iunit,'(a)',iostat=ios) line
      if (ios .eq. 0) then
         read (line,'(3f10.3,a4)') h(i), v0(1,i), v0(2,i), azo(i)
      else if (ios .lt. 0) then ! EOF
         jmod = i - 1
!         write (*,'(i2,a)') jmod, ' lines read'
!         do j = 1,jmod
!            write (*,'(3f10.3,a4)') h(j), v0(1,j), v0(2,j), azo(j)
!         end do
         imoh = 600
         icon = 600
         ipd = 0
         isd = 0
         zoold = -999.
         zmax = h(jmod)
         zmaxout = zmax
!         write (*,'(a,2f10.2)') 'rmaxout, zmaxout = ', rmaxout, zmaxout
         return
      end if
   end do
   
   return
   
end subroutine rd_loc_mod


!***********************************************************************************************************************************
subroutine ttloc2 (ho, dis, czo1, nphas2, ttc, dtdd, dtdh, dddp, phcd, ierr, indph)

!     A subroutine to calculate travel-time tables for hyposat.f
!              johannes schweitzer
!     (developed from laufps.f (version 28. april 1997)) 11. june  1997
!     nov 14, 1997  indph corrected for multiples and surface focus.
!     march 1999    changes due to crust 5.1 model input
!     october 2000  including converted surface reflections and converted reflections from the discontinuities conrad and moho
!     last changes/corrections 27. july 2001
! Modified for use in mloc by eric bergman, january 2006.
!     input:
!             ho      = source depth in km (fixed or not, see fixho)
!             dis     = receiver distance in deg
!             czo1    = character to inform about fixed or not fixed 
!                       depth (D == no-fixed, other == fixed)
!                       only P-onsets at receiver
!             indph   = 10000 for only direct P-onsets
!                     = 11000 for P, pP
!                     = 11100 for P, pP, PP 
!                     = 11110 for P, pP, PP, PbP, PmP 
!                     = 11111 for P, pP, PP, PbP, PmP, sP, SmP
!                       only S-onsets at receiver
!             indph   = 20000 for only direct S-onsets
!                     = 22000 for S, sS
!                     = 22200 for S, sS, SS 
!                     = 22220 for S, sS, SS, SbS, SmS 
!                     = 22222 for S, sS, SS, SbS, SmS, pS, PmS 
!                       P- and S-onsets at receiver
!             indph   = 30000 for P and S (direct only)
!                     = 33000 for P, S, pP, sS
!                     = 33300 for P, S, pP, sS, PP, SS
!                     = 33330 for P, S, pP, sS, PP, SS, PbP, PmP, SbS, SmS
!                     = 33333 for P, S, pP, sS, PP, SS, PbP, PmP, SbS, SmS, pS, sP, PmS, SmP
!    in common model :
!             v0(1,i) =  p velocity in layer i
!             v0(2,i) =  s velocity in layer i
!             h(i)    =  depth of layer i
!             jmod    =  number of layers
!             azo(i)  =  conrad/moho indicator
!             zmax    =  maximum depth for which this model can be used
!     output 
!             nphas2  = number of found onsets at receiver distance
!             phcd    = array with names of found onsets
!             ttc     = travel times of onsets in [sec]
!             dtdd    = ray parameters of onsets in [sec/deg]
!             dtdh    = partial derivative d(travel time)/d(ho) in [sec/km]
!             dddp    = partial derivative d(ray parameter)/d(dis) in [(sec/deg**2]
!             dpdh    = partial derivative d(ray parameter)/d(ho) in [(sec/deg)/km]
!             ierr    = 0     everything o.k.
!                       else  some error occurred.

   implicit none

!     mphas  = max number of phases for tt-calculations (see in
!              tauget_mod and hyposat_loc )
!
!              parameter maxp in ttlim.h must have the same value!
!
   integer, parameter ::  mphas = 60

   character(len=8) :: phcd(mphas)
   real :: ttc(mphas), dtdd(mphas), dtdh(mphas)

   include 'hyposat_loc.inc'
   
   character(len=1) :: czo1
   character(len=4) :: az(max_layers)
   character(len=8) :: phas1
   character(len=132) :: msg

   integer :: phnum, ibn, i, j, k, nzo, jh2, ij, i2, izo, m, indph1, indphx, indph2, indph3, indph4, indph5, imul, iq4, mult,&
    iq5, iq6, ib2, ib1, ibm, kk, kc, iphase, iph, klr, nphas, nphas1, ni, j1, k1, ndisc(max_layers), indx(np), j2, nphas2, ierr,&
    indph, lenb
  
   real :: re, d, b, zdiff, ffa, vmax, c, g1, g2, vv, r, t, e, p, o, f, q, g3, dt, del1, fct1, fct2, fct3, tt1, pa2, dtdh1,&
    dtdh2, dpdh1, dpdh2, dpdd1, dpdd2, dh1, dh2, dd1, dd2, z(max_layers), zo(3), pm(3,3), tm(3,3), dpdh(mphas), dddp(mphas),&
    dh3, dd3, ho, dis, tti(np), eps
   
   logical :: fixho, kp, ks, zof, surf, mul, surfc

   data phlist /'Pg      ','Pb      ','Pn      ','P       ','pPg     ','pPb     ','pPn     ','pP      ','PbP     ','PmP     ',&
                'PgPg    ','PbPb    ','PnPn    ','PP      ','pSg     ','pSb     ','pSn     ','pS      ','PbS     ','PmS     ',&
                'Sg      ','Sb      ','Sn      ','S       ','sSg     ','sSb     ','sSn     ','sS      ','SbS     ','SmS     ',&
                'SgSg    ','SbSb    ','SnSn    ','SS      ','sPg     ','sPb     ','sPn     ','sP      ','SbP     ','SmP     '/

   eps = 2.*epsilon(dh3)
   
   piref = 4.0*atan(1.0)
   pim = piref/180.0
   re = 6371.0
   aa = pim*re
   ib = 20
   ibn = ib*10
   ierr = 0

   dtdh = 0.
   dtdd = 0.
   dddp = 0.
   dpdh = 0.

   ! Reset onset table.

   if (ho .gt. zmax) then
      ierr = 99
      write (msg,'(a,2f6.1)') 'ttloc2: depth greater than maximum model depth!', ho, zmax
      call oops (trim(msg))
   end if
   ion = 0
   zof = .false.
   if (czo1 .eq. 'D') then
      fixho = .false.
      nzo = 3
      jh1 = 2
      jh2 = 3
      zo(2) = ho
      zo(1) = zo(2) - 1.0
      if (zo(1) .lt. 0.0) zo(1) = 0.0
      zo(3) = zo(2) + 1.0
   else
      fixho = .true.
      nzo = 1
      jh1 = 1
      jh2 = 1
      zo(1) = ho
      if (abs(zo(1)-zoold) .lt. eps) zof = .true.
   end if
   del(2) = dis
   del(1) = del(2) - 0.010
   if (del(1) .lt. 0.0) del(1) = 0.0
   del(3) = del(2) + 0.010

   do iql = 1,nzo

      ! First work to establish the model. We will calculate at first P then S phases. (k-loop)

      do k = 1,2
         if (k .eq. 1) then
            kp = .true.
            ks = .false.
         else
            kp = .false.
            ks = .true.
         end if
         if (zof) cycle
         ij = 0
         izo = 0
         do i = 1,jmod
            i2  = i + 1
            ij = ij + 1
            d = v0(k,i)
            b = re/(re-h(i))
            if (abs(zo(iql)-h(i)) .lt. eps .and. izo .eq. 0) then
               iqq = ij
               izo = 1
            end if
            az(ij) = azo(i)
            v(k,ij) = d*b
            if (v(1,ij) .ge. 10.0 .and. ipd .ne. 0) ipd = i
            if (v(2,ij) .ge. 4.70 .and. isd .ne. 0) isd = i
            z(ij) = re*log(b)
            if (h(i2) .gt. zo(iql) .and. h(i) .lt. zo(iql) .and. izo .eq. 0) then
               d = (zo(iql)-h(i))*(v0(k,i2)-v0(k,i))/(h(i2)-h(i)) + v0(k,i)
               b = re/(re-zo(iql))
               ij = ij + 1
               iqq = ij
               az(ij) = ''
               izo = 1
               v(k,ij) = d*b
               if (v(1,ij) .ge. 10.0 .and. ipd .ne. 0) ipd = i
               if (v(2,ij) .ge. 4.70 .and. isd .ne. 0) isd = i
               z(ij) = re*log(b)
            end if
         end do
         j = ij
         m = j - 1
         do i = 1,m
            i2 = i + 1
            v2(k,i) = v(k,i2)
            if (az(i) .eq. 'CONR')  icon = i
            if (az(i) .eq. 'MOHO')  imoh = i
            if (abs(v2(k,i)-v(k,i)) .lt. eps) then
               v2(k,i) = 1.00010*v(k,i)
               v(k,i2) = v2(k,i)
            end if
            zdiff = z(i2) - z(i)
            ndisc(i) = 0
            if (abs(zdiff) .lt. eps)  then
               zdiff = 0.00010
               z(i2) = z(i2) + zdiff
               ndisc(i) = 1
            endif
            g(k,i) = (v2(k,i)-v(k,i))/zdiff
         end do
      end do

      do k = 1,2
         if (k .eq. 1) then
           kp = .true.
           ks = .false.
         else
           kp = .false.
           ks = .true.
         end if
         indph1 = indph  / 10000
         indphx = indph  - indph1*10000
         indph2 = indph  / 1000
         indphx = indph  - indph2*1000
         indph3 = indphx / 100
         indphx = indphx - indph3*100
         indph4 = indphx / 10
         indph5 = indphx - indph4*10
         if (indph1 .ne. 3) then
           if (kp .and. indph1 .ne. 1) cycle
           if (ks .and. indph1 .ne. 2) cycle
         end if
         vhq(k) = v(k,iqq)*v(k,iqq)
         pa1(k) = (re-z(iqq))*pim/v(k,iqq)

         conv = .false.
         
         ! Direct waves (if source deeper than 0.) plus defining requested direct phases
         if (iqq .ne. 1)  then
            if (kp) phaseref(1:1) = 'P'
            if (ks) phaseref(1:1) = 'S'
            if (iqq .le. imoh) then
               if (iqq .le. icon) phaseref(2:) = 'g     '
               if (iqq .gt. icon) phaseref(2:) = 'b     '
            else
               phaseref(2:) = 'n     '
            end if
            if ((kp .and. iqq .ge. ipd .and. ipd .ne. 0) .or. (ks .and. iqq .ge. isd .and. isd .ne. 0) ) phaseref(2:) = '      '
            call reflex (iqq, k)
         end if

         ! Body waves 

         imul = 1111
         mul = .false.
         surf = .false.
         surfc = .false.
         iq4 = 0
         mult = 0
         ffa = 0.0

1100     iq5 = 1
         iq6 = 0
         if (iqq .gt. 1) iq5 = iqq
         if (iq4 .gt. 1) iq5 = iq4
         if (imul .lt. m .and. imul .gt. iq6) iq6 = imul
         iq5 = max0(iq6,iq5)
         vmax = v(k,iq5)

         do i = 1,iq5
            fa(1,i) = 2.0
            fa(2,i) = 0.0
            if (i .ge. imul) fa(1,i) = ffa
            if (i .lt. iqq)  then
               fa(1,i)=fa(1,i)-1.0
               cycle
            end if
            if (i .lt. iq4 .and. surf) fa(1,i) = fa(1,i) + 1.0
            if (i .lt. iq4 .and. surfc) then
               fa(2,i) = fa(2,i) + 1.0
               if (k .eq. 1 .and. vmax .lt. v(2,i)) vmax = v(2,i)
               if (k .eq. 2 .and. vmax .lt. v(1,i)) vmax = v(1,i)
            end if
            if (vmax .lt. v(k,i)) vmax = v(k,i)
         end do

         do i = iq5,m
            ib2 = ib
            ib1 = 1
            ibm = 0
            fa(1,i) = 2.0
            if (i .ge. imul) fa(1,i) = ffa
            if (i .eq. ndisc(i)) cycle
            if (kp) phaseref(1:1) = 'P'
            if (ks) phaseref(1:1) = 'S'
            if (i .lt. imoh) then
               if (i .le. icon) phaseref(2:) = 'g     '
               if (i .gt. icon) phaseref(2:) = 'b     '
            else
               phaseref(2:) = 'n     '
            end if
            if ((kp .and. i .ge. ipd .and. ipd .ne. 0) .or. (ks .and. i .ge. isd .and. isd .ne. 0) ) phaseref(2:) = '      '
            d = v2(k,i)
            if (d .le. vmax) cycle
            
            if (imul .ge. m .or. i .ge. imul) then
               c = v(k,i)
               if (c .lt. vmax) c = vmax
1350           g1 = ib2-1
               b = (d-c)/g1
               do i2 = ib1,ib2
                  g2 = i2-1
                  vv = c + g2*b
                  r = 0.0
                  t = 0.0
                  do kk = 1,i
                     e = v(k,kk)
                     g1 = e/vv
                     p = sqrt(abs(1.0-g1*g1))
                     o = 1.0/g(k,kk)
                     if (kk .ge. i)  then
                        f = vv
                        q = 0.0
                     else
                        f = v2(k,kk)
                        g3 = f/vv
                        q = sqrt(abs(1.0-g3*g3))
                     end if
                     r = r + fa(1,kk)*(p-q)*o
                     dt = fa(1,kk)*log(f*(1.0+p)/(e*(1.0+q)))*o
                     t = t + dt
                  end do

                  ! Extension for converted surface reflections (i.e. pS or sP)

                  if (surfc) then
                     if (k .eq. 1) kc = 2
                     if (k .eq. 2) kc = 1
                     do kk = 1,iq4
                        if (fa(2,kk) .lt. 0.90) cycle
                        e = v(kc,kk)
                        g1 = e/vv
                        p = sqrt(abs(1.0-g1*g1))
                        o = 1.0/g(kc,kk)
                        f = v2(kc,kk)
                        g3 = f/vv
                        q = sqrt(abs(1.0-g3*g3))
                        r = r + fa(2,kk)*(p-q)*o
                        dt = fa(2,kk)*log(f*(1.0+p)/(e*(1.0+q)))*o
                        t = t + dt
                     end do
                  end if
                  rr(2) = r*vv/aa
                  tt(2) = t
                  pa(2) = aa/vv
                  phas1 = phaseref
                  iphase = lenb(phaseref)
                  if (mul) phas1 = phaseref(1:iphase)//phaseref(1:iphase)
                  if (surf) then
                     if (kp) phas1 = 'p' // phaseref(1:iphase)
                     if (ks) phas1 = 's' // phaseref(1:iphase)
                  else if (surfc) then
                     if (kp) phas1 = 's' // phaseref(1:iphase)
                     if (ks) phas1 = 'p' // phaseref(1:iphase)
                  end if
                  iph = phnum(phas1)
                  do klr = 1,3
                     if (iql .ne. jh1 .and. klr .ne. 2) cycle
                     del1 = del(klr)
                     if (abs(rr(2)-del1) .lt. eps) then
                        if (ibm .eq. 0 ) then
                           ib2 = ibn
                           ib1 = (i2-1)*10
                           if (ib1 .le. 0) ib1 = 1
                           ibm = 1
                           go to 1350
                        end if
                        ion(iph,iql,klr) = ion(iph,iql,klr) + 1
                        if (ion(iph,iql,klr) .eq. 1) then
                           ttp(iph,iql,klr) = tt(2)
                           ppp(iph,iql,klr) = pa(2)
                        else
                           if (tt(2) .lt. ttp(iph,iql,klr)) then
                              ttp(iph,iql,klr) = tt(2)
                              ppp(iph,iql,klr) = pa(2)
                           end if
                        end if
                        cycle
                     end if
                     if (i2 .le. 1) cycle
                     fct1 = del1 - rr(1)
                     fct2 = del1 - rr(2)
                     if (fct1*fct2 .lt. 0.0) then
                        if (ibm .eq. 0 ) then
                           ib2 = ibn
                           ib1 = (i2-1)*10
                           ibm = 1
                           go to 1350
                        end if
                        fct3 = fct1/(rr(2) - rr(1))
                        tt1 = fct3*(tt(2) - tt(1)) + tt(1)
                        pa2 = fct3*(pa(2) - pa(1)) + pa(1)
                        ion(iph,iql,klr) = ion(iph,iql,klr) + 1
                        if (ion(iph,iql,klr) .eq. 1) then
                           ttp(iph,iql,klr) = tt1
                           ppp(iph,iql,klr) = pa2
                        else
                           if (tt1 .lt. ttp(iph,iql,klr)) then
                              ttp(iph,iql,klr) = tt1
                              ppp(iph,iql,klr) = pa2
                           end if
                        end if
                     end if
                  end do
                  rr(1) = rr(2)
                  tt(1) = tt(2)
                  pa(1) = pa(2)
               end do
            end if
            
            vmax = d

         end do

         if (mul) go to 3500

         ! Now the surface reflections of the body waves will be calculated (i.e., pP, sS...).

         if (surf .or. iqq .eq. 1) then
            if (surf) then
               iqq = iq4
               iq4 = 1
               surf = .false.
            end if
            go to 3100
         end if
         if (indph2 .ne. 3) then
           if (kp .and. indph2 .ne. 1) go to 3100
           if (ks .and. indph2 .ne. 2) go to 3100
         endif
         surf = .true.
         iq4 = iqq
         iqq = 1
         go to 1100

         ! Now the converted surface reflections of the body waves will be calculated (i.e., pS, sP...).

3100        continue

         if (surfc .or. iqq .eq. 1) then
            if (surfc) then
               iqq = iq4
               iq4 = 1
               surfc = .false.
            end if
            go to 3200
         end if

         if (indph5 .ne. 3) then
            if (kp .and. indph5 .ne. 1) go to 3200
            if (ks .and. indph5 .ne. 2) go to 3200
         end if

         surfc = .true.
         iq4 = iqq
         iqq = 1

         go to 1100

         ! Now the multiple phases will be done (e.g., PgPg, PnPn or SnSn...)
         ! At the moment only single multiples are calculated (e.g., no PPP or SnSnSn ...).

3200        continue

         if (indph3 .ne. 3) then
            if (kp .and. indph3 .ne. 1) go to 3500
            if (ks .and. indph3 .ne. 2) go to 3500
         end if

         imul = 1
         mult = 1
         mul = .true.

         ffa = mult*2 + 2
         go to 1100

         ! End of the body-phase and direct-wave loop

3500     continue

         ! Reflections for P and S from the two possible layers: the Conrad and the Moho discontinuities.

         if (indph4 .ne. 3) then
            if (kp .and. indph4 .ne. 1) go to 6500
            if (ks .and. indph4 .ne. 2) go to 6500
         end if
         do i = iqq,j
            if (i .ne. icon .and. i .ne. imoh) cycle
            if (i .eq. icon) then
               if (kp) phaseref = 'PbP'
               if (ks) phaseref = 'SbS'
            else if (i .eq. imoh) then
               if (kp) phaseref = 'PmP'
               if (ks) phaseref = 'SmS'
            end if
            i2=i
            call reflex (i2, k)
         end do
6500     continue

         ! Converted reflections for p and s from the two possible layers: the Conrad and the Moho discontinuities.
         ! PbS, SbP, PmS, SmP

         if (indph5 .ne. 3) then
            if (kp .and. indph5 .ne. 1) cycle
            if (ks .and. indph5 .ne. 2) cycle
         end if

         conv = .true.

         do i = iqq,j
            if (i .ne. icon .and. i .ne. imoh) cycle
            if (i .eq. icon) then
               if (kp) phaseref = 'SbP'
               if (ks) phaseref = 'PbS'
            else if (i.eq.imoh) then
               if (kp) phaseref = 'SmP'
               if (ks) phaseref = 'PmS'
            end if
            i2=i
            call reflex (i2, k)
         end do

         conv = .false.

      end do

   end do
   
   ! Finally we have to do some interpolations

   nphas = mphas + 1
   nphas1 = 0
   do i = 1,np
      indx(i) = 0
      tti(i) = 0.
      dtdh1 = 0.0
      dtdh2 = 0.0
      dpdh1 = 0.0
      dpdh2 = 0.0
      dpdd1 = 0.0
      dpdd2 = 0.0
      dh1 = 0.0
      dh2 = 0.0
      dd1 = 0.0
      dd2 = 0.0
      ni = 0
      do j = 1,nzo
         do k = 1,3
            tm(j,k) = 0.0
            pm(j,k) = 0.0
            if (ion(i,j,k) .eq. 0) go to 8455
            tm(j,k) = ttp(i,j,k)
            pm(j,k) = ppp(i,j,k)
            if (fixho) go to 8450
            j1 = j - 1
            if (k .eq. 2) then
               if (j .eq. 1) then
                  dtdh1 = tm(j,k)
                  dpdh1 = pm(j,k)
                  dh1   = zo(j) 
               else if (j .eq. 2 .and. ion(i,j1,k) .eq. 0) then
                  dtdh1 = tm(j,k)
                  dpdh1 = pm(j,k)
                  dh1   = zo(j) 
               else if (j .eq. 3 .and. ion(i,j,k) .ne. 0) then
                  dtdh2 = tm(j,k)
                  dpdh2 = pm(j,k)
                  dh2   = zo(j)
               else if (j .eq. 3 .and. ion(i,j,k) .eq. 0) then
                  dtdh2 = tm(j1,k)
                  dpdh2 = pm(j1,k)
                  dh2   = zo(j1)
               end if
            end if

8450        k1 = k - 1

            if (j .eq. jh1) then
               if (k .eq. 1) then
                  dpdd1 = pm(j,k)
                  dd1 = del(k)
               else if (k .eq. 2 .and. ion(i,j,k1) .eq. 0) then
                  dpdd1 = pm(j,k)
                  dd1 = del(k)
               else if (k .eq. 3 .and. ion(i,j,k) .ne. 0) then
                  dpdd2 = pm(j,k)
                  dd2  = del(k)
               else if (k .eq. 3 .and. ion(i,j,k) .eq. 0) then
                  dpdd2 = pm(j,k1)
                  dd2   = del(k1)
               end if
            end if

            if (j .eq. jh1 .and. k .eq. 2) then
               ni = 1
               nphas1= nphas1 + 1
               nphas = nphas - 1
               tti(nphas1) = tm(j,k)
               dtdd(nphas) = pm(j,k)
               phcd(nphas) = phlist(i)
               if (lenb(phcd(nphas)) .ge. 3) tti(nphas1) = tti(nphas1) + 0.00001
            end if

8455           if (j .eq. jh2 .and. k .eq. 3 .and. ni .ne. 0) then
               dtdh(nphas) = 0.
               dpdh(nphas) = 0.0
               dddp(nphas) = 0.0
               if (.not.fixho) then
                  dh3 = dh2 - dh1
                  if (abs(dh3) .gt. eps) then
                     dtdh(nphas) = (dtdh2-dtdh1)/dh3
                     dpdh(nphas) = (dpdh2-dpdh1)/dh3
!                     if (phcd(nphas)(1:2) .eq. 'Pg') print *, nphas, dtdh1, dtdh2, dh3, dtdh(nphas)
                  end if
               end if
               dd3 = dd2 - dd1
               ! If the next statement is uncommented, there is a bug that overwrites dtdh
!                  if (dd3 .ne. 0.0) dddp(nphas) = (dpdd2-dpdd1)/dd3
!                  if (phcd(nphas)(1:2) .eq. 'Pg') print *, nphas, dpdd1, dpdd2, dd3, dddp(nphas)
!                  print *, i, j, k, dtdh
            end if
            
         end do
      end do
!         if (nphas1 .gt. 0) then
!            print *, i, nphas, nphas1, phcd(nphas), tti(nphas1), dtdd(nphas), dtdh(nphas), dpdh(nphas), dddp(nphas),
!     #       dh3, dtdh1, dtdh2
!         end if
      
   end do
   
   call indexx (nphas1, tti, indx)

   do i = 1,nphas1
      j = indx(i)
      j2 = mphas + 1 - j
      ttc(i) = tti(j)
      phcd(i) = phcd(j2)
      dtdd(i) = dtdd(j2)
      dtdh(i) = dtdh(j2)
      dddp(i) = dddp(j2)
      dpdh(i) = dpdh(j2)
!      print *, i, j, j2, dis, phcd(i), ttc(i), dtdd(i), dtdh(i), dpdh(i), dddp(i)
   end do
   nphas2 = nphas1
   zoold = ho
   
   return

end subroutine ttloc2


!***********************************************************************************************************************************
subroutine reflex (ii, k)
      
   implicit none

   include 'hyposat_loc.inc'

   real :: vmax, b, rvv, t, r, e, g1, p, o, f, q, dt, fi, del1, fct1, fct2, fct3, tt1, pa2, eps
   integer :: phnum, ii, k, l, kc, i, ibc, kk, klr, iph

   eps = 2.*epsilon(del1)

   iph = phnum(phaseref)
   l = ii - 1
   
   if (conv) then
      vmax = max(v(1,iqq),v(2,iqq))
      if (k .eq. 1) kc = 2
      if (k .eq. 2) kc = 1
      do i = 1,l
         fa(1,i) = 1.0
         fa(2,i) = 0.0
         if (v(k,i) .gt. vmax) vmax = v(k,i)
         if (v2(k,i) .gt. vmax) vmax = v2(k,i)
         if (i .ge. iqq) then
            fa(2,i) = 1.0
            if (v(kc,i) .gt. vmax) vmax = v(kc,i)
            if (v2(kc,i) .gt. vmax) vmax = v2(kc,i)
         end if
      end do
   else
      vmax = v(k,iqq)
      do i = 1,l
         fa(1,i) = 2.0
         if (i .lt. iqq) fa(1,i) = 1.0
         if (v(k,i) .gt. vmax) vmax = v(k,i)
         if (v2(k,i) .gt. vmax) vmax = v2(k,i)
      end do
   end if

   ibc = ib*30
   if (ibc .gt. 600) ibc = 600
   rr(1) = 0.0
   rr(2) = 0.0
   b = piref/(2.0*(ibc-1))

   do i = 1,ibc
      rvv = sin(b*(i-1))/vmax
      t = 0.0
      r = 0.0

      do kk = 1,l
         e = v(k,kk)
         g1 = e*rvv
         p = sqrt(abs(1.0-g1*g1))
         o = 1.0/g(k,kk)
         f = v2(k,kk)
         g1 = f*rvv
         q = sqrt(abs(1.0-g1*g1))
         r = r + fa(1,kk)*(p-q)*o
         dt = fa(1,kk)*log(f*(1.0+p)/(e*(1.0+q)))*o
         t = t + dt
      end do

      if (conv) then
         do kk = 1,l
            if (fa(2,kk) .lt. 0.90) cycle
            e = v(kc,kk)
            g1 = e*rvv
            p = sqrt(abs(1.0-g1*g1))
            o = 1.0/g(kc,kk)
            f = v2(kc,kk)
            g1 = f*rvv
            q = sqrt(abs(1.0-g1*g1))
            r = r + fa(2,kk)*(p-q)*o
            dt = fa(2,kk)*log(f*(1.0+p)/(e*(1.0+q)))*o
            t = t + dt
         end do
      end if

      tt(2) = t

      if (i .gt. 1) then
         rr(2) = r/rvv
         p = sqrt(abs(1.0/(rvv*rvv*vhq(k))-1.0))
         if (abs(p) .lt. eps) then
            fi = 0.50*piref
         else
            fi = atan(1.0/p)
         end if
      else
         fi = 0.0
      end if

      if (ii .eq. iqq) fi = piref - fi

      rr(2) = rr(2)/aa
      pa(2) = sin(fi)*pa1(k)

      do 1400 klr = 1,3

         if (iql .ne. jh1 .and. klr .ne. 2) go to 1400

         del1 = del(klr)

         if (abs(rr(2)-del1) .lt. eps) then
            ion(iph,iql,klr) = ion(iph,iql,klr) + 1
            if (ion(iph,iql,klr) .eq. 1) then
               ttp(iph,iql,klr) = tt(2)
               ppp(iph,iql,klr) = pa(2)
            else
               if (tt(2) .lt. ttp(iph,iql,klr)) then
                  ttp(iph,iql,klr) = tt(2)
                  ppp(iph,iql,klr) = pa(2)
               end if
            end if
            go to 1400
         end if

         if (i .eq. 1) go to 1400

         fct1 = del1 - rr(1)
         fct2 = del1 - rr(2)
   
         if (fct1*fct2 .lt. 0.0) then
            fct3 = fct1/(rr(2)-rr(1))
            tt1 = fct3*(tt(2)-tt(1)) + tt(1)
            pa2 = fct3*(pa(2)-pa(1)) + pa(1)
            ion(iph,iql,klr) = ion(iph,iql,klr) + 1
            if (ion(iph,iql,klr) .eq. 1) then
               ttp(iph,iql,klr) = tt1
               ppp(iph,iql,klr) = pa2
            else
               if (tt1 .lt. ttp(iph,iql,klr)) then
                  ttp(iph,iql,klr) = tt1
                  ppp(iph,iql,klr) = pa2
               end if
            end if
         end if

1400  continue

      tt(1) = tt(2)
      rr(1) = rr(2)
      pa(1) = pa(2)

   end do

   return
   
end subroutine reflex
      
      
!***********************************************************************************************************************************
function phnum (phasein)
      
   implicit none
   
   include 'hyposat_loc.inc'
   
   integer :: phnum, i
   character(len=8) phasein
   
   phnum = 999
   do i = 1,np
      if (phasein .eq. phlist(i)) then
         phnum = i
         return
      end if
   end do

   return
      
end function phnum
