!***********************************************************************************************************************************      
subroutine phreid3 (it, iev, ird, nr)
      
! Phase identification for single readings. No concept of groups of readings.
! This version uses a traditional "best fit" approach, but modified to take advantage of
! information on the probability distribution functions (PDFs) of different phases.
! The target arrival time is tested against all phases in the theoretical TT model for the
! corresponding focal depth and epicentral distance, regardless arrival time order. Probability
! is calculated for each possible phase ID, based on the candidate phases's PDF
! and the choice is made on the basis of highest probability. There are a couple of rules:
!  1) Depth phases are handled separately
!  2) Phase type (P or S) is honored if a phase name has been provided
! A phase name is not changed unless the probability of the new phase ID is 0.05 or greater.
! The necessary data of coefficients to describe the PDFs of all the needed
! phases is not yet available. Therefore, at this time all PDS are the same and the choice boils
! down to the classical "smallest residual" criterion.

   implicit none
   
   include 'mloc.inc'
         
   integer :: max
   parameter (max=60)

   common /taup/ nphase, tt(max), phcd(max)
   integer :: nphase
   real :: tt
   character(len=8) :: phcd

   common /taupnd/ nphasend, ttnd(max), phcdnd(max)
   integer :: nphasend
   real :: ttnd
   character(len=8) :: phcdnd

   character(len=8) :: phasein, phaseout
   real :: otsp, arrtime, ttobs, hms2s, usrc(2), dtdd(max), dtdh(max), dddp(max), probi, probmax, prob
   integer :: i, iev, it, nr, ird
   logical :: ptype, stype
   
   nr = 0
         
   ! Theoretical travel-time
   call depset (depthp(iev,it), usrc)
   if (.not.locmod) then
      call trtm (delt(iev,ird), max, nphase, tt, dtdd, dtdh, dddp, phcd)
   else
      if (verbose_log) write (io_log,'(a)') 'phreid3: calling tt_mixed_model'
      call tt_mixed_model (depthp(iev,it), delt(iev,ird), nphase, tt, dtdd, dtdh, dddp, phcd)
   end if
         
   call nodupe () ! Get the equivalent phase list without duplicates
   if (nphasend .eq. 0) then
      write (io_plog,'(a,f10.3)') 'phreid3: no phases returned delta = ', delt(iev,ird)
      nphasend = 1
      phcdnd(1) = 'CRAP    '
   end if
   
   ! Some cases in which phase identification is not done. 
   if (phase(iev,ird) .eq. phcd(1)) return
   if (phase(iev,ird)(1:2) .eq. 'T ') return
   if (phase(iev,ird)(1:3) .eq. 'pwP') return
   
   phasein = phase(iev,ird)
   if (phase(iev,ird) .eq. 'Lg      ' .and. delt(iev,ird) .le. lg_min) phasein = 'Sg      '
   
   ! Observed travel time in seconds, relative to current origin time
   otsp = hms2s(hourp(iev,it),minp(iev,it),secp(iev,it))
   arrtime = hms2s(ipah(iev,ird),ipam(iev,ird),pas(iev,ird))
   if (ipah(iev,ird) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
   ttobs = arrtime - otsp
   
   if (nphase .eq. 0) then ! Phase ID not changed.
      phaseout = phasein
   else if (nphase .ge. 1) then
      ! Probability of candidate phases
      ! Depth phases can only be identified as depth phases.
      probmax = 0.
      if (phasein(1:1) .eq. 'p' .or. phasein(1:1) .eq. 's') then
         do i = 1,nphase
            if (phcd(i)(1:1) .ne. 'p' .and. phcd(i)(1:1) .ne. 's') cycle
            probi = prob(phcd(i), 0, ttobs, delt(iev,ird), depthp(iev,it))
            if (probi .gt. probmax) then
               probmax = probi
               phaseout = phcd(i)
            end if
         end do
      else
         do i = 1,nphase
            if (phcd(i)(1:1) .eq. 'p' .or. phcd(i)(1:1) .eq. 's') cycle
            if (ptype(phasein) .and. .not.ptype(phcd(i))) cycle
            if (stype(phasein) .and. .not.stype(phcd(i))) cycle
            probi = prob(phcd(i), 0, ttobs, delt(iev,ird), depthp(iev,it))
            if (probi .gt. probmax) then
               probmax = probi
               phaseout = phcd(i)
            end if
         end do
         ! Lg is not in the standard phase list so check it separately. Also check if an
         ! Sg might be re-identified as Lg
         if (phasein .eq. 'Lg      ' .or. (phasein .eq. 'Sg      ' .and. delt(iev,ird) .gt. lg_min)) then
            probi = prob('Lg      ', 0, ttobs, delt(iev,ird), depthp(iev,it))
            if (probi .gt. probmax) then
               probmax = probi
               phaseout = 'Lg      '
            end if
         end if
      end if
      if (probmax .lt. 0.05) phaseout = phasein ! Leave it unchanged
      if (verbose_log) write (io_plog,'(2a,2f8.2,2f8.3,1x,a)') 'phreid3: ', phasein, delt(iev,ird), ttobs, probi, probmax, phaseout
   end if
      
   phase(iev,ird) = phaseout
   
   if (phaseout .ne. phasein) then
      call readerr (iev,ird) ! Assign the correct reading error for the new phase ID
      nr = nr + 1
      write (io_plog,'(a,i3,1x,a,2(1x,a),f7.2,1x,a)') 'phreid3: ', iev, stname(iev,ird),&
       phasein, phaseout, delt(iev,ird), infile20(iev)
   end if
                            
   return
   
end subroutine phreid3


!***********************************************************************************************************************************      
real function prob (phasein, itype, t, d, h)
      
! For a given phase arrival, return the probability of being observed.
! In general, the calculation is done on the basis of:
!   1) relative observability of that phase
!   2) median (absolute or relative to AK135) as a function of delta
!   3) spread as a function of delta
!   4) a composit PDF which defines the shape of the distribution 
! The details of the algorithm can be different for each phase

! Input: itype - 0 for travel time, 1 for residual
!        t - time, either travel time or a residual
!        d - epicentral distance
!        h - focal depth

   implicit none
   
   include 'mloc.inc'

   integer :: max ! number of possible phases in return list and tau-p list
   parameter (max=60)
   
   common /taup/ nphase, tt(max), phcd(max)
   integer :: nphase
   real :: tt
   character(len=8) :: phcd
   
   character(len=8) :: phasein
   character(len=1) :: c
   integer :: mb, itype, i
   real :: b(8), xmed, psprd, t, r, d, robs, x, f, h, ttlg, tttp
   logical :: ptype, stype
         
   data mb/8/, b/-0.07,  0.88, -0.40,  0.58,  0.13,  2.00,  1.10,  0.10/

   if (itype .eq. 0) then ! Convert travel time to residual
      do i = 1,nphase
         if (phcd(i) .eq. phasein) then
            r = t - tt(i)
            go to 10 ! to avoid later arrivals of the same phase
         end if
      end do
      if (phasein(1:2) .eq. 'Lg') then
         ttlg = lg_a + d*lg_b
         r = t - ttlg
      end if
      if (phasein(1:2) .eq. 'T ') then
         tttp = tphase_a + d*tphase_b
         r = t - tttp
      end if
10    continue
   else
      r = t
   end if
   
   ! Default values
   xmed = 0.
   if (ptype(phasein)) then
      psprd = 2.5
   else if (stype(phasein)) then
      psprd = 4.
   else
      psprd = 3.
   end if
   c = phasein(1:1)
   robs = 1.
   if (c .eq. 'p' .or. c .eq. 's') then
      if (h .le. 15) then
         robs = robs*0.4
      else
         robs = robs*0.6
      end if
   end if
   
   ! Use offset from a .ttsprd file if it exists. This helps remove bias due to baseline offset
   ! for different phases.
   ! It is problematic to use the spread until I have distance-dependent parameterization. Use
   ! of a single value for a phase such as Pn, which covers a large epicentral distance range
   ! and has greatly increasing spread with distance, creates problems at both ends. Even
   ! the offset creates problems in certain situations (direct calibration) if it is not
   ! distance-dependent
   
!      if (read_ttsprd) then
!         do i = 1,nsprd
!            if (phasein .eq. sprdph(i)) then
!               xmed = ttoffset(i)
!               psprd = sprd(i)
!               exit
!            end if
!         end do
!      end if
      
   x = (r - xmed)/psprd
   call gswt (x, mb, b, f)
   prob = f*robs
   
   if (verbose_log) write (io_plog,'(2a,6f8.3)') 'prob: ', phasein, psprd, r, x, f, robs, prob
   
   return
   
end function prob

      
!***********************************************************************************************************************************      
subroutine gswt (x, m, b, f0)

! Derived from GSWT2, by removing calculation for derivatives, and entrypoint.
! gswt evaluates a function at a given location. The function is composed
! of one or more PDFs (a Gaussian plus one or two Cauchy PDFs).

! Input:      
! x  = the location at which the function is evaluated.
! m  = number of coefficients needed to describe the desired combination of PDFs which is evaluated:
!    = 2 for fitting only a Gaussian
!    = 5 for fitting a Gaussian plus a Cauchy
!    = 8 for fitting a Gaussian plus 2 Cauchy distributions
! b  = coefficients for the PDFs

! Output:
! f0 = the predicted value at x for the desired combination of PDFs

   implicit none
   
   real :: b(8), pii, rpi2, arg0, x, t0, g0, f0, arg1, t1, g1, arg2, t2, g2
   integer :: m
   
   data pii/0.3183098/ ! 1/pi
   data rpi2/0.3989422/ ! 1/sqrt(2pi)

   ! Guassian PDF (see Evans et al., p. 145)
   ! b(1) = location parameter, the mean
   ! b(2) = scale parameter, >0, the standard deviation

   arg0 = (x - b(1))/b(2)
   if (abs(arg0) .le. 4.) then
      t0 = rpi2*exp(-.5*arg0*arg0)/b(2) ! Gaussian PDF
   else ! To avoid FPU underflow
      t0 = 0.
   end if
   if (m .gt. 2) go to 1
   g0 = t0
   f0 = g0
   return
      
   ! Cauchy PDF (see Evans et al., p. 48)
   ! b(3) = location parameter a, the median
   ! b(4) = scale parameter b, b>0
   ! b(5) is the relative percentage of this Cauchy
      
 1 arg1 = (x - b(3))/b(4)
   t1 = pii/(b(4)*(1. + arg1*arg1)) ! Cauchy PDF
   g1 = b(5)*t1
   if (m .gt. 5) go to 2
   g0 = (1. - b(5))*t0 ! Scale the Gaussian back
   f0 = g0+g1
   return
      
   ! A second Cauchy PDF
   ! b(6) = location parameter a, the median
   ! b(7) = scale parameter b, b>0
   ! b(8) is the relative percentage of this Cauchy
      
 2 arg2 = (x - b(6))/b(7)
   t2 = pii/(b(7)*(1. + arg2*arg2)) ! Cauchy PDF
   g2 = b(8)*t2
   g0 = (1. - b(5) - b(8))*t0 ! Scale the Gaussian back
   f0 = g0 + g1 + g2
            
   return
   
end subroutine gswt

!***********************************************************************************************************************************      
logical function ptype (phasein)

! returns .true. if phasein is a P phase
! The first leg of depth phases is ignored.

   implicit none

   integer :: lenb, k
   logical :: dphase
   character(len=2) :: c2
   character(len=3) :: c3
   character(len=4) :: c4
   character(len=5) :: c5
   character(len=8) :: phasein

   if (phasein .eq. '        ') then
      ptype = .false.
      return
   end if
      
   if (phasein .eq. 'UNKNOWNP') then ! Unknown P-type phase
      ptype = .true.
      return
   end if

   if (phasein(1:2) .eq. 'T ') then ! T is considered an P phase.
      ptype = .true.
      return
   end if

   k = lenb(phasein)
   c2 = '  '
   c3 = '   '
   c4 = '    '
   c5 = '     '

   dphase = .false.
   if (phasein(1:1) .eq. 'p' .or. phasein(1:1) .eq. 's') dphase = .true. ! Depth phase
   if (dphase) then
      if (k .ge. 3) c2 = phasein(2:3)
      if (k .ge. 4) c3 = phasein(2:4)
      if (k .ge. 5) c4 = phasein(2:5)
      if (k .ge. 6) c5 = phasein(2:6)
   else
      if (k .ge. 2) c2 = phasein(1:2)
      if (k .ge. 3) c3 = phasein(1:3)
      if (k .ge. 4) c4 = phasein(1:4)
      if (k .ge. 5) c5 = phasein(1:5)
   end if

   if (phasein(k:k) .eq. 'P') then ! Check last leg
      ptype = .true.
   else if (c2 .eq. 'Pn' .or.&
            c2 .eq. 'Pb' .or.&
            c2 .eq. 'P*' .or.&
            c2 .eq. 'Pg') then
      ptype = .true.
   else if (c3 .eq. 'PKP' .or.&
            c3 .eq. 'SKP') then
      ptype = .true.
   else if (c4 .eq. "P'P'" .or.&
            c4 .eq. 'PKKP' .or.&
            c4 .eq. 'SKKP') then
      ptype = .true.
   else if (c5 .eq. 'Pdiff') then
      ptype = .true.
   else
      ptype = .false.
   end if

   return
   
end function ptype

            
!***********************************************************************************************************************************      
logical function stype (phasein)

! Returns .true. if phasein is an S phase.
! The first leg of depth phases is ignored.

   implicit none

   integer :: lenb, k
   logical :: dphase
   character(len=2) :: c2
   character(len=3) :: c3
   character(len=4) :: c4
   character(len=5) :: c5
   character(len=8) :: phasein

   if (phasein .eq. '        ') then
      stype = .false.
      return
   end if
      
   if (phasein .eq. 'UNKNOWNS') then ! Unknown S-type phase
      stype = .true.
      return
   end if

   if (phasein(1:2) .eq. 'Lg') then ! Lg is considered an S phase.
      stype = .true.
      return
   end if

   if (phasein(1:2) .eq. 'T ') then ! T is considered an P phase.
      stype = .false.
      return
   end if

   k = lenb(phasein)
   c2 = '  '
   c3 = '   '
   c4 = '    '
   c5 = '     '

   dphase = .false.
   if (phasein(1:1) .eq. 'p' .or.&
       phasein(1:1) .eq. 's') dphase = .true. ! Depth phase
   if (dphase) then
      if (k .ge. 3) c2 = phasein(2:3)
      if (k .ge. 4) c3 = phasein(2:4)
      if (k .ge. 5) c4 = phasein(2:5)
      if (k .ge. 6) c5 = phasein(2:6)
   else
      if (k .ge. 2) c2 = phasein(1:2)
      if (k .ge. 3) c3 = phasein(1:3)
      if (k .ge. 4) c4 = phasein(1:4)
      if (k .ge. 5) c5 = phasein(1:5)
   end if

   if (phasein(k:k) .eq. 'S') then ! Check last leg
      stype = .true.
   else if (c2 .eq. 'Sn' .or.&
            c2 .eq. 'Sb' .or.&
            c2 .eq. 'S*' .or.&
            c2 .eq. 'Sg') then
      stype = .true.
   else if (c3 .eq. 'SKS' .or.&
            c3 .eq. 'PKS') then
      stype = .true.
   else if (c4 .eq. "S'S'" .or.&
            c4 .eq. 'SKKS' .or.&
            c4 .eq. 'PKKS') then
      stype = .true.
   else if (c5 .eq. 'Sdiff') then
      stype = .true.
   else
      stype = .false.
   end if

   return
   
end function stype

      
!***********************************************************************************************************************************      
subroutine nodupe ()

! Removes duplicate phases from the tau-p phase list (as given in the
! common block /taup/, and returns
! the equivalent phase list in the common block /taupnd/.

   implicit none

   integer :: max
   parameter (max=60)

   common /taup/ nphase, tt(max), phcd(max)
   integer :: nphase
   real :: tt
   character(len=8) :: phcd

   common /taupnd/ nphasend, ttnd(max), phcdnd(max)
   integer :: nphasend
   real :: ttnd
   character(len=8) :: phcdnd

   integer :: i, j, nd
   logical :: dupe

   ! Initialize
   do i = 1,max
      phcdnd(i) = '        '
      nphasend = 0
      ttnd(i) = 0.
   end do

   if (nphase .eq. 1) then
      phcdnd(1) = phcd(1)
      nphasend = 1
      ttnd(1) = tt(1)
   else if (nphase .gt. 1) then
      phcdnd(1) = phcd(1)
      ttnd(1) = tt(1)
      nd = 1
      do i = 2,nphase
         dupe = .false.
         do j = 1,nd
            if (phcdnd(j) .eq. phcd(i)) dupe = .true.
         end do
         if (.not.dupe) then
            nd = nd + 1
            phcdnd(nd) = phcd(i)
            ttnd(nd) = tt(i)
         end if
      end do
      nphasend = nd
   end if

   return
   
end subroutine nodupe
      
      
!***********************************************************************************************************************************
character(len=8) function iscphase (n)
      
!  Given a phase code number, returns the ISC phase identification.

   implicit none

   integer :: n
   character(len=8) :: phase(0:100)
   
   data phase/'P       ','PP      ','PPP     ','PCP     ','PKP     ',&
              'PKP2    ','PKPPKP  ','PCPPKP  ','PS      ','PPS     ',&
              'PCS     ','PKS     ','PKKS    ','PCSPKP  ','PKPPKS  ',&
              'PKPSKS  ','PKKP    ','3PKP    ','PKIKP   ','PP2     ',&
              'PPP2    ','PKS2    ','PSS     ','PSS2    ','SSP2    ',&
              'PCPPKP2 ','PCSPKP2 ','SS2     ','PKKP2   ','PKKS2   ',&
              'SCSPKP3 ','SCSPKP2 ','SCSP2   ','SKSP2   ','SSS2    ',&
              'S       ','SS      ','SSS     ','SCS     ','SKS     ',&
              'SKKS    ','SKKKS   ','SCSPKP  ','SKSSKS  ','SCSP    ',&
              'SKSP    ','SCP     ','SP      ','SKP     ','SKKP    ',&
              'SKPPKP  ','SSP     ','SKP2    ','SKS2    ','SKKS2   ',&
              'SKKS3   ','SKKKS2  ','*SPKP2  ','*PPCP   ','*PPKP   ',&
              '*PP     ','*PPP    ','*SP     ','*SPKP   ','*SS     ',&
              '*SSS    ','*SPP    ','*SPCP   ','*SSCS   ','*PPKP2  ',&
              'P*      ','S*      ','PG      ','SG      ','PN      ',&
              'SN      ','PGPG    ','SGSG    ','LR      ','LQ      ',&
              'L       ','PKKP3   ','PKKS3   ','SPP     ','PHASE84 ',&
              'P DIFF  ','QM      ','RM      ','T       ','T(MAX)  ',&
              'NORTH   ','SOUTH   ','EAST    ','WEST    ','UP      ',&
              'DOWN    ','E       ','I       ','MAXIMUM ','FINAL   ',&
              '        '/
  
   iscphase = phase(n)
   
   return
   
end function iscphase
     

!***********************************************************************************************************************************
logical function crustal_phase (phasein)
      
! Returns "true" if the candidate phase matches an entry in the list of crustal phases.
! Basic list of crustal phases from Storchak et al., The IASPEI Standard Seismic Phase List,
! Seismological Research Letters, v. 74, No. 6, 2003. Crustal surface waves (Lg, Rg) 
! are not considered, nor are multiple Moho reflections, e.g., PmPN. Some additional
! phases have been added (PbPb, SbSb, depth phases).
      
   implicit none
   
   integer :: nphases
   parameter (nphases=28)
   
   integer :: i
   character(len=8) :: phasein, phases(nphases)
   
   data phases /'Pg      ','Pb      ','Pn      ','Sg      ','Sb      ','Sn      ',&
                'pPg     ','pPb     ','pPn     ','pSg     ','pSb     ','pSn     ',&
                'sPg     ','sPb     ','sPn     ','sSg     ','sSb     ','sSn     ',&
                'PgPg    ','PbPb    ','PnPn    ','SgSg    ','SbSb    ','SnSn    ',&
                'PmP     ','PmS     ','SmS     ','SmP     '/
   
   crustal_phase = .false.
   do i = 1,nphases
      if (phasein .eq. phases(i)) then
         crustal_phase = .true.
         exit
      end if
   end do
   
   return
   
end function crustal_phase


!***********************************************************************************************************************************
subroutine phwgt (phcd, delta, wgt)
   
   implicit none
   
   integer :: mxbr, nb, mxo
   parameter (mxbr=45,nb=180,mxo=8)
   
   integer :: k
   real :: delta, wgt, spread, x, evlply
   character(len=8) :: phcd
   character(len=20) :: modnam

   common /resprm/ nbr, phnm(mxbr), ms(mxbr), me(mxbr), xm(mxbr), mo(3,mxbr), prm(mxo,3,mxbr)
   character(len=8) :: phnm
   real :: xm, prm
   integer :: nbr, ms, me, mo

   if(phcd.ne.'START') go to 1
   call getprm (modnam)
   return
    
1  k=0
   if(phcd(1:3).eq.'Pb '.or.phcd(1:3).eq.'Pg '.or. phcd(1:3).eq.'Pn '.or.phcd(1:3).eq.'P  ') k=1
   if((phcd(1:3).eq.'Pn '.or.phcd(1:2).eq.'P ').and. delta.gt.16.5.and.delta.le.30.5) k=2
   if(phcd(1:2).eq.'P '.and. delta.gt.30.5.and.delta.le.91.5) k=3
   if(phcd(1:2).eq.'P '.and.delta.gt.91.5) k=4
   if(phcd(1:6).eq.'Pdiff ') k=5
   if(phcd(1:6).eq.'PKiKP ') k=6
   if(phcd(1:6).eq.'PKPdf ') k=7
   if(phcd(1:6).eq.'PKPab ') k=8
!  if(phcd(1:6).eq.'PKPbc ') k=9
!  if(phcd(1:7).eq.'PKPdiff') k=10
!  if(phcd(1:6).eq.'PKPpre') k=11
   if(phcd(1:3).eq.'pP '.or.phcd(1:4).eq.'pPn ') k=12
   if(phcd(1:7).eq.'pPdiff ') k=13
   if(phcd(1:3).eq.'sP '.or.phcd(1:4).eq.'sPn '.or. phcd(1:7).eq.'sPdiff ') k=14
   if(phcd(1:4).eq.'pwP '.or.phcd(1:5).eq.'pwPn '.or. phcd(1:8).eq.'pwPdiff') k=15
   if(phcd(1:3).eq.'Sg ') k=20
   if(phcd(1:3).eq.'Sb ') k=21
   if((phcd(1:3).eq.'Sn '.or.phcd(1:3).eq.'S  ').and. delta.le.18.5) k=22
   if((phcd(1:3).eq.'Sn '.or.phcd(1:2).eq.'S ').and. delta.gt.18.5.and.delta.le.27.5) k=23
   if(phcd(1:2).eq.'S '.and.delta.gt.27.5) k=24
   if(phcd(1:6).eq.'Sdiff ') k=25

   if(k.ne.0) go to 2
   write(6,*) phcd,delta
   wgt=0.0
   return

2  if(delta.lt.ms(k)) delta=ms(k)
   if(delta.gt.me(k)) delta=me(k)
   x=delta-xm(k)
!      med_dt=evlply(x,mo(2,k),prm(1,2,k))
   spread=evlply(x,mo(3,k),prm(1,3,k))
   wgt=wgt/(spread**2)
!      hitcnt=evlply(x,mo(1,k),prm(1,1,k))
!
!      print *,'Evlply:  phcd,phnm,k,delta,x,med_dt,spread,hitcnt =',phcd,phnm(k),k,delta,x,med_dt,spread,hitcnt
!      go to 10
!

   return
   
end subroutine phwgt


!***********************************************************************************************************************************
subroutine phase_alias (phase, delta, ip)
      
! Check for alternative phase names
! Input: phase, delta
! Output: ip (index of phcod)

   implicit none
   
   integer :: ip
   real :: delta
   character(len=*) :: phase

   if (phase(1:3).eq.'Pn ') then ! phase='P       '
     ip=2
   else if(phase(1:3).eq.'Sn ') then ! phase='S       '
     ip=33
   else if(phase(1:4).eq.'pPn ') then ! phase='pP      '
     ip=8
   else if(phase(1:4).eq.'pwP ') then ! phase='pP      '
     ip=8
   else if(phase(1:5).eq.'pwPn ') then ! phase='pP      '
     ip=8
   else if(phase(1:4).eq.'sPn ') then ! phase='sP      '
     ip=13
   else if(phase(1:4).eq.'pSn ') then ! phase='pS      '
     ip=37
   else if(phase(1:4).eq.'sSn ') then ! phase='sS      '
     ip=40
   else if(phase(1:4).eq.'SPn ') then ! phase='SP      '
     ip=55
   else if(phase(1:4).eq.'SnP ') then ! phase='SP      '
     ip=55
   else if(phase(1:4).eq.'PSn ') then ! phase='PS      '
     ip=56
   else if(phase(1:5).eq.'PnPn ') then ! phase='PP      '
     ip=30
   else if(phase(1:5).eq.'SnSn ') then ! phase='SS      '
     ip=53
   else if(phase(1:2).eq.'p ') then ! phase='Pup     '
     ip=1  
   else if(phase(1:2).eq.'s ') then ! phase='Sup     '
     ip=32 
   else if(phase(1:6).eq."P'P'ab") then ! phase="P'P'    '
     ip=31
   else if(phase(1:6).eq."P'P'bc") then ! phase="P'P'    '
     ip=31
   else if(phase(1:6).eq."P'P'df") then ! phase="P'P'    '
     ip=31
   else if(phase(1:6).eq."S'S'ac") then ! phase="S'S'    '
     ip=54
   else if(phase(1:6).eq."S'S'df") then ! phase="S'S'    '
     ip=54
   else if(delta.le.100.0.and.phase.eq.'pPdiff  ') then ! phase='pP      '
     ip=8
   else if(delta.le.100.0.and.phase.eq.'pwPdiff ') then ! phase='pP      '
     ip=8
   else if(delta.le.100.0.and.phase.eq.'sPdiff  ') then ! phase='sP      '
     ip=13
   else if(delta.le.100.0.and.phase.eq.'pSdiff  ') then ! phase='pS      '
     ip=37
   else if(delta.le.100.0.and.phase.eq.'sSdiff  ') then ! phase='sS      '
     ip=40
   else if(delta.le.165.0.and.phase.eq.'PKPdiff ') then ! phase='PKPbc   '
     ip=5
   else if(delta.le.165.0.and.phase.eq.'pPKPdiff') then ! phase='pPKPbc  '
     ip=10
   else if(delta.le.165.0.and.phase.eq.'sPKPdiff') then ! phase='sPKPbc  '
     ip=15
   else
     ip=-1
   end if
   
   return
   
end subroutine phase_alias


!***********************************************************************************************************************************
subroutine pnclean (pin, pout)

! Clean up phase names.
! Remove initial 'X', 'e', 'i', and 'q' and trailing 'c', ,'d' or 'r' from phase name.
! Remove the !-flag in character position 8 that fixes phase ID.
! Also removes parentheses

   implicit none

   integer :: lenb, j, k, i
   logical :: trim
   character(len=1) :: c
   character(len=8) :: pin, pout, b8
   
   data b8/'        '/

   pout = pin
   if (pout .eq. b8 .or. pout .eq. 'UNK     ') then
      pout = 'UNKNOWN '
      return
   end if
   if (pout(8:8) .eq. '!') pout(8:8) = ' '
   if (pout .eq. 'PKPpre  ') pout = 'PKPdfpre'
   if (pout .eq. 'PKPdfpre') return
   if (pout .eq. 'S-P     ') return

   j = lenb(pout)
   if (j .ge. 1) then

      ! Trim junk off the beginning
      trim = .true.
      do while (trim)
         if (pout(1:1) .eq. 'X' .or.& 
             pout(1:1) .eq. ' ' .or.& 
             pout(1:1) .eq. 'e' .or.& 
             pout(1:1) .eq. 'E' .or.& 
             pout(1:1) .eq. 'i' .or.& 
             pout(1:1) .eq. 'I' .or.& 
             pout(1:1) .eq. 'q' .or.& 
             pout(1:1) .eq. '(') then
            if (j .ge. 2) then
               do i = 2,j
                  k = i - 1
                  pout(k:k) = pout(i:i)
               end do
               pout(j:j) = ' '
               j = j - 1
            else
               pout = 'UNKNOWN '
               return
            end if
         else
            trim = .false.
         end if
      end do
         
      ! Trim junk off the end
      if (j .ge. 2) then
         k = j - 1
         if (pout(k:j) .ne. 'bc' .and. pout(k:j) .ne. 'ac') then ! Watch out for ac and bc branches
            c = pout(j:j)
            if (c .eq. 'c' .or. c .eq. 'd' .or. c .eq. 'r' .or. c .eq. ')') pout(j:j) = ' '
         end if
      else if (j .eq. 1) then
         c = pout(1:1)
         if (c .ne. 'P') return
         if (c .ne. 'S') return
         if (c .ne. 'L') return
         if (c .ne. 'M') return
         if (c .ne. 'Q') return
         if (c .ne. 'R') return
         pout = 'UNKNOWN '
         return
      end if
   
      j = lenb(pout)
   
      ! Fix second character of some crustal phases
      if (j .eq. 2) then
         if (pout(2:2) .eq. 'N') pout(2:2) = 'n'
         if (pout(2:2) .eq. 'G') pout(2:2) = 'g'
         if (pout(2:2) .eq. 'B') pout(2:2) = 'b'
         if (pout(2:2) .eq. '*') pout(2:2) = 'g'
      end if
   
      ! Outer core reflections
      if (j .eq. 3 .and. pout(2:2) .eq. 'C') pout(2:2) = 'c'

      ! Inner-core reflections
      if (j .eq. 5 .and. pout(3:3) .eq. 'I') pout(3:3) = 'i'
   
      ! Fix any remaining blanks inside the phase name (trailing blanks are OK)
      do i = 1,j
         if (pout(i:i) .eq. ' ') pout(i:i) = '_'
      end do
   
   end if

   return
   
end subroutine pnclean


!***********************************************************************************************************************************
logical function skipp (phasein)

! Returns a logical value controlling whether a phase should be kept
! in the dataset or not. Mostly used to get rid of various forms of
! surface waves, and amplitude measurements.

! Some exotic phases are skipped because they are not in the ak135 phase set.

   implicit none

   character(len=1) :: c
   character(len=8) :: phasein

   skipp = .false.
   if (phasein(1:2) .eq. 'Lg')  return
   if (phasein(1:2) .eq. 'T ') return

   if (phasein(1:3) .eq. 'amp') then
      skipp = .true.
      return
   end if

   if (phasein(1:3) .eq. 'PPP' .or. phasein(1:3) .eq. 'SSS') then ! Not in ak135 phase set
      skipp = .true.
      return
   end if

   if (phasein(1:4) .eq. 'P3KP') then ! Not in ak135 phase set
      skipp = .true.
      return
   end if

   c = phasein(1:1)
   if (c .eq. 'L') then ! Catches, e.g., L, LR. LQ
      skipp = .true.
   else if (c .eq. 'M') then ! Catches, e.g.,  M, MLR, MAXIMUM
      skipp = .true.
   else if (c .eq. 'R') then ! Catches, e.g., R, RM
      skipp = .true.
   else if (c .eq. 'Q') then ! Catches, e.g., Q, QM
      skipp = .true.
   else if (c .eq. 'A') then ! Catches AMP
      skipp = .true.
   else if (c .eq. 'F') then ! Catches FINAL
      skipp = .true.
   end if

   if (trim(phasein) .eq. 'Pmax') skipp = .true.
   if (trim(phasein) .eq. 'Smax') skipp = .true.
   if (trim(phasein) .eq. 'max') skipp = .true.

   return
   
end function skipp


!***********************************************************************
real function evlply (x, m, b)

! Evlply evaluates the polynomial.

   implicit none

!  dimension b(1)
   real :: b(10), x, xx
   integer :: m, j

   xx = 1.
   evlply = b(1)
   do j = 2,m
      xx = xx*x
      evlply = evlply + xx*b(j)
   end do

   return
   
end function evlply


!***********************************************************************
subroutine getprm (modnam)

! Read the parameter file.

   implicit none

   integer :: mxcm, mxbr, mxo
   parameter (mxcm=10, mxbr=45, mxo=8)

   character(len=*) :: modnam
   character(len=20) :: name
   integer :: ios, lunit, i, j, k, nprm
   !logical log

   common /resprm/ nbr, phnm(mxbr), ms(mxbr), me(mxbr), xm(mxbr), mo(3,mxbr), prm(mxo,3,mxbr)
   character(len=8) :: phnm
   real :: xm, prm
   integer :: nbr, ms, me, mo

   !data nprm/2/

   nprm = lunit()
   open (nprm,file='res_ac.prm',access='sequential',form='formatted', status='old')

   read (nprm,'(1x,a)') name
   !modnam = modnam(1:ntrbl(modnam))//name
   modnam = name
   !print *,'Getprm:  modnam = ',modnam

   do i = 1,mxbr
      read (nprm,'(1x,a8,2i8,f8.1)',iostat=ios) phnm(i), ms(i), me(i), xm(i)
      if (ios .lt. 0) go to 1
      !print *,'Getprm:  i phnm ms me xm =', i, ' ', phnm(i), ' ', ms(i), me(i), xm(i)
      do k = 1,3
         read (nprm,'(6x,i5,1p,4e15.7:/(11x,1p,4e15.7))') mo(k,i), (prm(j,k,i),j=1,mo(k,i))
         !print *,'Getprm:  k mo prm =', k, mo(k,i), (prm(j,k,i), j=1,mo(k,i))
      end do
   end do
   i = mxbr + 1

 1 nbr = i - 1
   !print *,'Getprm:  nbr =', nbr
   close (nprm)
   
   return
   
end subroutine getprm
      
