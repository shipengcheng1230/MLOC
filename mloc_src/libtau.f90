!***********************************************************************************************************************************      
subroutine tabin (in, taup_path, dirsym, modname)
      
! Read .hed and .tbl files for tau-p calculations.
      
   implicit none
   
   include 'ttlim.inc'

   character(len=*) :: modname
   character(len=8) :: phdif(6)
   character(len=40) :: basename
   character(len=80) :: taup_path
   character(len=1) :: dirsym
   character(len=132) :: msg
   integer :: in, nasgr, nl, len2, nph, i, j, k, l, ind

   double precision :: pm, zm
   integer :: ndex, mt
   common /umdc/ pm(jsrc,2), zm(jsrc,2), ndex(jsrc,2), mt(2)
   
   double precision :: us, pt, tau, xlim, xbrn, dbrn 
   real :: xn, pn, tn, dn, hn
   integer :: jndx, idel, mbr1, mbr2
   common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
    idel(jbrn,3), mbr1, mbr2
       
   double precision :: zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
   real :: odep, fcs
   integer :: nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
   common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
    coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
    isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
   
   character(len=8) :: phcd
   common /pcdc/ phcd(jbrn)      

   phdif(1) = 'P'
   phdif(2) = 'S'
   phdif(3) = 'pP'
   phdif(4) = 'sP'
   phdif(5) = 'pS'
   phdif(6) = 'sS'

   ! Diagnostic output from tau-p software. This is normally not printed.
!   open (iottim,file='ttim1.lis',status='unknown')

   nin = in
   
   basename = trim(taup_path)//dirsym//trim(modname)
   if (verbose_taup) call fyi ('Opening '//trim(basename)//'.hed')
   open (nin,file=trim(basename)//'.hed',status='old',form='unformatted')
   rewind (nin)
   
   if (verbose_taup) then
      write (msg,'(a,i4)') 'tabin: jseg = ', jseg
      call fyi (trim(msg))
   end if
   
   read (nin) nasgr, nl, len2, xn, pn, tn, mt, nseg, nbrn, ku, km, fcs, nafl, indx, kndx
   read (nin) pm, zm, ndex
   read (nin) pu, pux
   read (nin) phcd, px, xt, jndx
   read (nin) pt, taut
   read (nin) coef
   
   close (nin)
   
   if (verbose_taup) call fyi ('tabin: opening '//trim(basename)//'.tbl')
   open (nin,file=trim(basename)//'.tbl',status='old',form='unformatted', access='direct', recl=nasgr)
   ! rewind (nin)

   do nph = 1,2
      pu(ku(nph)+1,nph) = pm(1,nph)
   end do
   tn = 1./tn
   dn = 3.1415927/(180.*pn*xn)
   odep = -1.
   ki = 0
   msrc(1) = 0
   msrc(2) = 0
   k = 1
   do 3 i = 1,nbrn
      jidx(i) = jndx(i,2)
      do j = 1,2
         dbrn(i,j) = -1d0
      end do
8     if (jndx(i,2) .le. indx(k,2)) go to 7
      k = k+1
      go to 8
7     if (nafl(k,2) .gt. 0) go to 9
      ind = nafl(k,1)
      l = 0
      do j = jndx(i,1),jndx(i,2)
         l = l + 1
         tp(l,ind) = pt(j)
      end do
9     if (nafl(k,1) .gt. 0 .and. (phcd(i)(1:1) .eq. 'P'.or. phcd(i)(1:1) .eq. 'S')) go to 3
      do j = 1,6
         if (phcd(i) .eq. phdif(j)) go to 6
      end do
      go to 3
6     dbrn(i,1) = 1d0
      phdif(j) = ' '
3  continue

   return
   
end subroutine tabin
      
      
!***********************************************************************************************************************************      
block data tabin_init

   implicit none
   
   include 'ttlim.inc'

   double precision :: zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
   real :: odep, fcs
   integer :: nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
   common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
    coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
    isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
   
   data tauc,xc/jtsm*0d0,jxsm*0d0/
      
end

      
!***********************************************************************************************************************************      
subroutine depset(dep,usrc)
      
   implicit none
   save 
   include 'ttlim.inc'
   
   logical :: dop, dos
   real :: usrc(2), dep, rdep, eps
   integer :: i, ind, k, j, int, nph
   
   double precision :: pm, zm
   integer :: ndex, mt
   common /umdc/ pm(jsrc,2), zm(jsrc,2), ndex(jsrc,2), mt(2)
   
   double precision :: us, pt, tau, xlim, xbrn, dbrn
   real :: xn, pn, tn, dn, hn
   integer :: jndx, idel, mbr1, mbr2
   common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn,&
    jndx(jbrn,2), idel(jbrn,3), mbr1, mbr2
    
   double precision :: zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
   real :: odep, fcs
   integer :: nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
   common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
    coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
    isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
   
   character(len=8) :: phcd
   common /pcdc/ phcd(jbrn)
   
   logical :: segmsk, prnt
   common /prtflc/ segmsk(jseg), prnt(2)
   
   eps = 2.*epsilon(dep)
   
   if (abs(max(dep,.011)-odep) .gt. eps) go to 1
   dop = .false.
   dos = .false.
   do 2 i = 1,nseg
      if (.not.segmsk(i) .or. iidx(i) .gt. 0) go to 2
      if (iabs(nafl(i,1)) .le. 1) dop = .true.
      if (iabs(nafl(i,1)) .ge. 2) dos = .true.
2  continue
   if (.not.dop .and. .not.dos) return
   go to 3

1  nph0 = 0
   int0(1) = 0
   int0(2) = 0
   mbr1 = nbrn + 1
   mbr2 = 0
   dop = .false.
   dos = .false.
   do 4 i = 1,nseg
      if (.not.segmsk(i)) go to 4
      if (iabs(nafl(i,1)) .le. 1) dop = .true.
      if (iabs(nafl(i,1)) .ge. 2) dos = .true.
4  continue
   do 5 i = 1,nseg
      if (nafl(i,2) .gt. 0 .or. odep .lt. 0.) go to 5
      ind = nafl(i,1)
      k = 0
      do 15 j = indx(i,1),indx(i,2)
         k = k + 1
15       pt(j) = tp(k,ind)
5  iidx(i) = -1
   do 6 i = 1,nbrn
6     jndx(i,2) = -1
   if (ki .le. 0) go to 7
   do 8 i = 1,ki
      j = kk(i)
8     pt(j)=pk(i)
   ki = 0
   
   ! Sample the model at the source depth.
7  odep = amax1(dep,.011)
   rdep = dep
   if (rdep .lt. .011) rdep = 0.
   zs = amin1(alog(amax1(1.-rdep*xn,1e-30)),0.)
   hn = 1./(pn*(1.-rdep*xn))
   if (prnt(1) .or. prnt(2)) write (10,'(/1x,a,f7.2/)') 'Depth =', dep

3  if (nph0 .gt. 1) go to 12
   if (dop) call depcor(1)
   if (dos) call depcor(2)
   go to 14
12 if (dos) call depcor(2)
   if (dop) call depcor(1)

   ! Interpolate all tau branches.

14 j = 1
   do 9 i = 1,nseg
      if (.not.segmsk(i)) go to 9
      nph = iabs(nafl(i,1))
!      print *, 'i iidx nph msrc nafl =', i, iidx(i), nph,msrc(nph), nafl(i,1)
      if (iidx(i) .gt. 0 .or. (msrc(nph) .le. 0 .and. nafl(i,1) .gt. 0)) go to 9
      iidx(i) = 1
      if (nafl(i,2) .le. 0) int = nafl(i,1)
      if (nafl(i,2) .gt. 0 .and. nafl(i,2) .eq. iabs(nafl(i,1))) int = nafl(i,2) + 2
      if (nafl(i,2) .gt. 0 .and. nafl(i,2) .ne. iabs(nafl(i,1))) int = iabs(nafl(i,1)) + 4
      if (nafl(i,2) .gt. 0 .and. nafl(i,2) .ne. nafl(i,3)) int = nafl(i,2) + 6
11    if (jndx(j,1) .ge. indx(i,1)) go to 10
      j = j + 1
      go to 11
10    idel(j,3) = nafl(i,1)
!      print *,'spfit:  j int =', j, int 
      call spfit (j,int)
      mbr1 = min0(mbr1,j)
      mbr2 = max0(mbr2,j)
      if (j .ge. nbrn) go to 9
      j = j + 1
!      print *, 'j jidx indx jndx =', j, jidx(j), indx(i,2), jndx(j,2)
      if (jidx(j) .le. indx(i,2) .and. jndx(j,2) .gt. 0) go to 10
9  continue

!   write (10,*) 'mbr1 mbr2', mbr1, mbr2
!   write (10,*) 'msrc isrc odep zs us', msrc, isrc, odep, sngl(zs), sngl(us(1)), sngl(us(2))
!   write (10,'(/10x,i5/(1x,3i5,f12.6))') ki, (i, iidx(i), kk(i), pk(i), i = 1,nseg)

   usrc(1) = sngl(us(1)/pn)
   usrc(2) = sngl(us(2)/pn)
   
   return
   
end


!***********************************************************************************************************************************      
block data depset_init
      
   implicit none
   
   include 'ttlim.inc'

   logical :: segmsk, prnt
   common /prtflc/ segmsk(jseg), prnt(2)
   
   data segmsk,prnt/jseg*.true.,2*.false./
      
end
      
      
!***********************************************************************************************************************************      
      subroutine depcor(nph)
      
      implicit none
      save
      include 'ttlim.inc'
      
      logical noend, noext
      double precision tup(jrec), umod, zmod, tauus1(2), tauus2(2), xus1(2), xus2(2), ttau, tx, sgn, umin, dtol, u0, u1,&
       z0, z1, fac, du, eps
      real tol, ztol
      integer lpower, nph, ks, i, n1, k2, k, k1, ms, mu, is, j, mod, iph, l, lp, kph, i1, i2, m
      
      common /umdc/ pm(jsrc,2), zm(jsrc,2), ndex(jsrc,2), mt(2)
      double precision pm, zm
      integer ndex, mt
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
       coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
       isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
      double precision zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
      real odep, fcs
      integer nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
      
      common/pcdc/phcd(jbrn)
      character*8 phcd      
      
      common /pdec/ ua(5,2), taua(5,2), deplim, ka
      double precision ua, taua
      real deplim
      integer ka
      
      common /prtflc/ segmsk(jseg), prnt(2)
      logical segmsk, prnt
      
      equivalence (tauc,tup)
      
      data tol,dtol,lpower/.01,1d-6,7/
      
      eps = 2.*epsilon(umin)
 
!     write(10,*)'depcor:  nph nph0',nph,nph0
      if(nph.eq.nph0) go to 1
      nph0=nph
      us(nph)=umod(zs,isrc,nph)
!   If we are in a high slowness zone, find the slowness of the lid.
      umin=us(nph)
      ks=isrc(nph)
!     write(10,*)'ks us',ks,sngl(umin)
      do 2 i=1,ks
      if(pm(i,nph).gt.umin) go to 2
      umin=pm(i,nph)
 2    continue
!   Find where the source slowness falls in the ray parameter array.
      n1=ku(nph)+1
      do 3 i=2,n1
      if(pu(i,nph).gt.umin) go to 4
 3    continue
      k2=n1
      if(abs(pu(n1,nph)-umin) .lt. eps) go to 50
      call abort1('Source slowness too large.')
 4    k2=i
!50   write(10,*)'k2 umin',k2,sngl(umin)
!
!   Read in the appropriate depth correction values.
!
 50   noext=.false.
      sgn=1d0
      if(msrc(nph).eq.0) msrc(nph)=1
!   See if the source depth coincides with a model sample
      ztol=xn*tol/(1.-xn*odep)
      if(dabs(zs-zm(ks+1,nph)).gt.ztol) go to 5
      ks=ks+1
      go to 6
 5    if(dabs(zs-zm(ks,nph)).gt.ztol) go to 7
!   If so flag the fact and make sure that the right integrals are
!   available.
 6    noext=.true.
      if(msrc(nph).eq.ks) go to 8
      call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
      go to 11
!   If it is necessary to interpolate, see if appropriate integrals
!   have already been read in.
 7    if(msrc(nph).ne.ks+1) go to 9
      ks=ks+1
      sgn=-1d0
      go to 8
 9    if(msrc(nph).eq.ks) go to 8
!   If not, read in integrals for the model depth nearest the source
!   depth.
      if(dabs(zm(ks,nph)-zs).le.dabs(zm(ks+1,nph)-zs)) go to 10
      ks=ks+1
      sgn=-1d0
 10   call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
!   Move the depth correction values to a less temporary area.
 11   do 31 i=1,ku(nph)
 31   tauu(i,nph)=tup(i)
      k=ku(nph)
      do 12 i=1,km(nph)
      k=k+1
      xc(i)=tup(k)
 12   xu(i,nph)=tup(k)
!     write(10,*)'bkin',ks,sngl(sgn),sngl(tauu(1,nph)),sngl(xu(1,nph))
!
!   Fiddle pointers.
!
 8    msrc(nph)=ks
!     write(10,*)'msrc sgn',msrc(nph),sngl(sgn)
      noend=.false.
      if(dabs(umin-pu(k2-1,nph)).le.dtol*umin) k2=k2-1
      if(dabs(umin-pu(k2,nph)).le.dtol*umin) noend=.true.
      if(msrc(nph).le.1.and.noext) msrc(nph)=0
      k1=k2-1
      if(noend) k1=k2
!     write(10,*)'noend noext k2 k1',noend,noext,k2,k1
      if(noext) go to 14
!
!   Correct the integrals for the depth interval [zm(msrc),zs].
!
      ms=msrc(nph)
      if(sgn)15,16,16
 16   u0=pm(ms,nph)
      z0=zm(ms,nph)
      u1=us(nph)
      z1=zs
      go to 17
 15   u0=us(nph)
      z0=zs
      u1=pm(ms,nph)
      z1=zm(ms,nph)
 17   mu=1
!     write(10,*)'u0 z0',sngl(u0),sngl(z0)
!     write(10,*)'u1 z1',sngl(u1),sngl(z1)
      do 18 k=1,k1
      call tauint(pu(k,nph),u0,u1,z0,z1,ttau,tx)
      tauc(k)=tauu(k,nph)+sgn*ttau
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 18
      xc(mu)=xu(mu,nph)+sgn*tx
!     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 18   continue
      go to 39
!   If there is no correction, copy the depth corrections to working
!   storage.
 14   mu=1
      do 40 k=1,k1
      tauc(k)=tauu(k,nph)
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 40
      xc(mu)=xu(mu,nph)
!     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 40   continue
!
!   Calculate integrals for the ray bottoming at the source depth.
!
 39   xus1(nph)=0d0
      xus2(nph)=0d0
      mu=mu-1
      if(dabs(umin-us(nph)).gt.dtol.and.dabs(umin-pux(mu,nph)).le.dtol) mu=mu-1
!   This loop may be skipped only for surface focus as range is not
!   available for all ray parameters.
      if(msrc(nph).le.0) go to 1
      is=isrc(nph)
      tauus2(nph)=0d0
      if(dabs(pux(mu,nph)-umin).gt.dtol.or.dabs(us(nph)-umin).gt.dtol) go to 48
!   If we happen to be right at a discontinuity, range is available.
      tauus1(nph)=tauc(k1)
      xus1(nph)=xc(mu)
!     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)), sngl(xus1(nph)),'  *'
      go to 33
!   Integrate from the surface to the source.
 48   tauus1(nph)=0d0
      j=1
      if(is.lt.2) go to 42
      do 19 i=2,is
      call tauint(umin,pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
 19   j=i
!     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)), sngl(xus1(nph))
 42   if(dabs(zm(is,nph)-zs).le.dtol) go to 33
!   Unless the source is right on a sample slowness, one more partial
!   integral is needed.
      call tauint(umin,pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
!     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)), sngl(xus1(nph))
 33   if(pm(is+1,nph).lt.umin) go to 41
!   If we are in a high slowness zone, we will also need to integrate
!   down to the turning point of the shallowest down-going ray.
      u1=us(nph)
      z1=zs
      do 35 i=is+1,mt(nph)
      u0=u1
      z0=z1
      u1=pm(i,nph)
      z1=zm(i,nph)
      if(u1.lt.umin) go to 36
      call tauint(umin,u0,u1,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
 35   xus2(nph)=xus2(nph)+tx
!36   write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)), sngl(xus2(nph))
 36   z1=zmod(umin,i-1,nph)
      if(dabs(z0-z1).le.dtol) go to 41
!   Unless the turning point is right on a sample slowness, one more
!   partial integral is needed.
      call tauint(umin,u0,umin,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
      xus2(nph)=xus2(nph)+tx
!     write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)), sngl(xus2(nph))
!
!   Take care of converted phases.
!
 41   iph=mod(nph,2)+1
      xus1(iph)=0d0
      xus2(iph)=0d0
      tauus1(iph)=0d0
      tauus2(iph)=0d0
      go to (59,61),nph
 61   if(umin.gt.pu(ku(1)+1,1)) go to 53
!
!   If we are doing an S-wave depth correction, we may need range and
!   tau for the P-wave which turns at the S-wave source slowness.  This
!   would bd needed for sPg and SPg when the source is in the deep mantle.
!
      do 44 j=1,nbrn
      if((phcd(j)(1:2).ne.'sP'.and.phcd(j)(1:2).ne.'SP').or.px(j,2).le.0d0) go to 44
!     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1), px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 44   continue
      go to 53
!
!   If we are doing an P-wave depth correction, we may need range and
!   tau for the S-wave which turns at the P-wave source slowness.  This
!   would be needed for pS and PS.
!
 59   do 60 j=1,nbrn
      if((phcd(j)(1:2).ne.'pS'.and.phcd(j)(1:2).ne.'PS').or. px(j,2).le.0d0) go to 60
!     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1), px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 60   continue
      go to 53
!
!   Do the integral.
 45   j=1
!     write(10,*)'Depcor:  do pS or sP integral - iph =',iph
      do 46 i=2,mt(iph)
      if(umin.ge.pm(i,iph)) go to 47
      call tauint(umin,pm(j,iph),pm(i,iph),zm(j,iph),zm(i,iph),ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
 46   j=i
 47   z1=zmod(umin,j,iph)
      if(dabs(zm(j,iph)-z1).le.dtol) go to 53
!   Unless the turning point is right on a sample slowness, one more
!   partial integral is needed.
      call tauint(umin,pm(j,iph),umin,zm(j,iph),z1,ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
!     write(10,*)'is ks tauusp xusp',j,ks,sngl(tauus1(iph)), sngl(xus1(iph))
!
 53   ua(1,nph)=-1d0
!     if(odep.ge.deplim.or.odep.le..1) go to 43
      if(odep.ge.deplim) go to 43
      do 57 i=1,nseg
      if(.not.segmsk(i)) go to 57
      if(nafl(i,1).eq.nph.and.nafl(i,2).eq.0.and.iidx(i).le.0) go to 58
 57   continue
      go to 43
!
!   If the source is very shallow, we will need to insert some extra
!   ray parameter samples into the up-going branches.
!
 58   du=amin1(1e-5+(odep-.4)*2e-5,1e-5)
!     write(10,*)'Add:  nph is ka odep du us =',nph,is,ka,odep, sngl(du),sngl(us(nph))
      lp=lpower
      k=0
      do 56 l=ka,1,-1
      k=k+1
      ua(k,nph)=us(nph)-(l**lp)*du
      lp=lp-1
      taua(k,nph)=0d0
      j=1
      if(is.lt.2) go to 54
      do 55 i=2,is
      call tauint(ua(k,nph),pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph), ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
 55   j=i
!     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 54   if(dabs(zm(is,nph)-zs).le.dtol) go to 56
!   Unless the source is right on a sample slowness, one more partial
!   integral is needed.
      call tauint(ua(k,nph),pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
!     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 56   continue
      go to 43
!
!   Construct tau for all branches.
!
 1    mu=mu+1
 43   j=1
!     write(10,*)'mu',mu
!     write(10,*)'killer loop:'
      do 20 i=1,nseg
      if(.not.segmsk(i)) go to 20
!     write(10,*)'i iidx nafl nph',i,iidx(i),nafl(i,1),nph
      if(iidx(i).gt.0.or.iabs(nafl(i,1)).ne.nph.or.(msrc(nph).le.0.and.nafl(i,1).gt.0)) go to 20
!
      iph=nafl(i,2)
      kph=nafl(i,3)
!   Handle up-going P and S.
      if(iph.le.0) iph=nph
      if(kph.le.0) kph=nph
      sgn=isign(1,nafl(i,1))
      i1=indx(i,1)
      i2=indx(i,2)
!     write(10,*)'i1 i2 sgn iph',i1,i2,sngl(sgn),iph
      m=1
      do 21 k=i1,i2
      if(pt(k).gt.umin) go to 22
 23   if(dabs(pt(k)-pu(m,nph)).le.dtol) go to 21
      m=m+1
      go to 23
 21   tau(1,k)=taut(k)+sgn*tauc(m)
      k=i2
!     write(10,*)'k m',k,m
      go to 24
!22   write(10,*)'k m',k,m
 22   if(dabs(pt(k-1)-umin).le.dtol) k=k-1
      ki=ki+1
      kk(ki)=k
      pk(ki)=pt(k)
      pt(k)=umin
      fac=fcs(i,1)
!     write(10,*)'ki fac',ki,sngl(fac)
      tau(1,k)=fac*(tauus1(iph)+tauus2(iph)+tauus1(kph)+tauus2(kph))+sgn*tauus1(nph)
!     write(10,*)'&&&&& nph iph kph tauus1 tauus2 tau =',nph,iph,kph,sngl(tauus1(1)),sngl(tauus1(2)),sngl(tauus2(1)),&
!      sngl(tauus2(2)),sngl(tau(1,k))
 24   m=1
 26   if(jndx(j,1).ge.indx(i,1)) go to 25
      j=j+1
      go to 26
 25   jndx(j,2)=min0(jidx(j),k)
      if(jndx(j,1).lt.jndx(j,2)) go to 37
      jndx(j,2)=-1
      go to 20
!37   write(10,*)'j jndx jidx',j,jndx(j,1),jndx(j,2),jidx(j),' ', phcd(j)
 37   do 30 l=1,2
 28   if(dabs(pux(m,nph)-px(j,l)).le.dtol) go to 27
      if(m.ge.mu) go to 29
      m=m+1
      go to 28
 27   xbrn(j,l)=xt(j,l)+sgn*xc(m)
!     write(10,*)'x up:  j l m  ',j,l,m
      go to 30
 29   xbrn(j,l)=fac*(xus1(iph)+xus2(iph)+xus1(kph)+xus2(kph))+sgn*xus1(nph)
!     write(10,*)'x up:  j l end',j,l
!     write(10,*)'&&&&& nph iph kph xus1 xus2 xbrn =',nph,iph,kph,sngl(xus1(1)),sngl(xus1(2)),sngl(xus2(1)),&
!      sngl(xus2(2)),sngl(xbrn(j,l))
 30   continue
      if(j.ge.nbrn) go to 20
      j=j+1
      if(jndx(j,1).le.k) go to 25
 20   continue
      return
      end


!***********************************************************************************************************************************      
      block data depcor_init
      
      implicit none

      common /pdec/ ua(5,2), taua(5,2), deplim, ka
      double precision ua, taua
      real deplim
      integer ka
      
      data deplim,ka/1.1,4/
      
      end
      
      
!***********************************************************************************************************************************      
      double precision function umod(zs,isrc,nph)
      
      implicit none
      save 
      include 'ttlim.inc'
      
      character*31 msg
      double precision zs, uend, dtol, zmod
      real dep
      integer isrc(2), m1, nph, i, j, js
      
      common /umdc/ pm(jsrc,2), zm(jsrc,2), ndex(jsrc,2), mt(2)
      double precision pm, zm
      integer ndex, mt
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      data dtol/1d-6/
!
      m1=mt(nph)
      do 1 i=2,m1
      if(zm(i,nph).le.zs) go to 2
 1    continue
      dep=sngl(1d0-dexp(zs))/xn
      write(msg,100)dep
      write(6,100)dep
 100  format('Source depth (',f6.1,') too deep.')
      call abort1(msg)
 2    if(dabs(zs-zm(i,nph)).le.dtol.and.dabs(zm(i,nph)-zm(i+1,nph)).le. dtol) go to 3
      j=i-1
      isrc(nph)=j
      umod=pm(j,nph)+(pm(i,nph)-pm(j,nph))*(dexp(zs-zm(j,nph))-1d0)/ (dexp(zm(i,nph)-zm(j,nph))-1d0)
      return
 3    isrc(nph)=i
      umod=pm(i+1,nph)
      return
!
      entry zmod(uend,js,nph)
      i=js+1
      zmod=zm(js,nph)+dlog(dmax1((uend-pm(js,nph))*(dexp(zm(i,nph)- zm(js,nph))-1d0)/(pm(i,nph)-pm(js,nph))+1d0,1d-30))
      return
      end


!***********************************************************************************************************************************      
      subroutine bkin (lu, nrec, len, buf)
      
      ! Reads a block of len double precision words into array buf(len)
      ! from record nrec of the direct access unformatted file connected to
      ! logical unit lu.
      
      implicit none
      save
      
      double precision buf(len), tmp
      integer lu, nrec, len, i
      
      if (nrec .le. 0) then ! If the record doesn't exist, zero fill the buffer
         do i = 1,len
            buf(i) = 0d0
         end do
      else
         read (lu,rec=nrec) buf
         tmp = buf(1)
      end if
      
      return
      end
      
      
!***********************************************************************************************************************************      
      subroutine tauint(ptk,ptj,pti,zj,zi,tau,x)
      
      implicit none
      save
!
! $$$$$ calls warning $$$$$
!
!   Tauint evaluates the intercept (tau) and distance (x) integrals  for
!   the spherical earth assuming that slowness is linear between radii
!   for which the model is known.  The partial integrals are performed
!   for ray slowness ptk between model radii with slownesses ptj and pti
!   with equivalent flat earth depths zj and zi respectively.  The partial
!   integrals are returned in tau and x.  Note that ptk, ptj, pti, zj, zi,
!   tau, and x are all double precision.
!
      character*71 msg
      double precision ptk, ptj, pti, zj, zi, tau, x, xx, b, sqk, sqi, sqj, sqb, eps
      
      eps = 2.*epsilon(ptk)
!
      if(dabs(zj-zi).le.1d-9) go to 13
      if(dabs(ptj-pti).gt.1d-9) go to 10
      if(dabs(ptk-pti).le.1d-9) go to 13
      b=dabs(zj-zi)
      sqj=dsqrt(dabs(ptj*ptj-ptk*ptk))
      tau=b*sqj
      x=b*ptk/sqj
      go to 4
 10   if(ptk.gt.1d-9.or.pti.gt.1d-9) go to 1
!   Handle the straight through ray.
      tau=ptj
      x=1.5707963267948966d0
      go to 4
 1    b=ptj-(pti-ptj)/(dexp(zi-zj)-1d0)
      if(ptk.gt.1d-9) go to 2
      tau=-(pti-ptj+b*dlog(pti/ptj)-b*dlog(dmax1((ptj-b)*pti/((pti-b)*ptj),1d-30)))
      x=0d0
      go to 4
 2    if(abs(ptk-pti) .lt. eps) go to 3
      if(abs(ptk-ptj) .lt. eps) go to 11
      sqk=ptk*ptk
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(sqb.gt.1d-30) go to 5
      xx=0d0
      x=ptk*(dsqrt(dabs((pti+b)/(pti-b)))-dsqrt(dabs((ptj+b)/(ptj-b))))/b
      go to 6
 5    if(b*b.lt.sqk) go to 7
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*(sqb*sqj+b*ptj-sqk)),1d-30))
      x=ptk*xx/sqb
      go to 6
 7    xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptk*dabs(pti-b)),1d0),-1d0))-dasin(dmax1(dmin1((b*ptj-sqk)/(ptk*dabs(ptj-b)),1d0),-1d0))
      x=-ptk*xx/sqb
 6    tau=-(sqi-sqj+b*dlog((pti+sqi)/(ptj+sqj))-sqb*xx)
      go to 4
 3    sqk=pti*pti
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 8
      xx=dlog(dmax1((ptj-b)*(b*pti-sqk)/((pti-b)*(sqb*sqj+b*ptj-sqk)),1d-30))
      x=pti*xx/sqb
      go to 9
 8    xx=dsign(1.5707963267948966d0,b-pti)-dasin(dmax1(dmin1((b*ptj-sqk)/(pti*dabs(ptj-b)),1d0),-1d0))
      x=-pti*xx/sqb
 9    tau=-(b*dlog(pti/(ptj+sqj))-sqj-sqb*xx)
      go to 4
 11   sqk=ptj*ptj
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 12
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*(b*ptj-sqk)),1d-30))
      x=ptj*xx/sqb
      go to 14
 12   xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptj*dabs(pti-b)),1d0),-1d0))-dsign(1.5707963267948966d0,b-ptj)
      x=-ptj*xx/sqb
 14   tau=-(b*dlog((pti+sqi)/ptj)+sqi-sqb*xx)
!
!   Handle various error conditions.
!
 4    if(x.ge.-1d-10) go to 15
      write(msg,100)ptk,ptj,pti,tau,x
 100  format('Bad range: ',1p,5d12.4)
      call warning (msg)
 15   if(tau.ge.-1d-10) go to 16
      write(msg,101)ptk,ptj,pti,tau,x
 101  format('Bad tau: ',1p,5d12.4)
      call warning (msg(1:69))
 16   return
!   Trap null integrals and handle them properly.
 13   tau=0d0
      x=0d0
      return
      end


!***********************************************************************************************************************************      
      subroutine spfit(jb,int)
      
      implicit none
      save 
      include 'ttlim.inc'
      
      character*3 disc
      logical newgrd,makgrd
!     logical log
      double precision pmn,dmn,dmx,hm,shm,thm,p0,p1,tau0,tau1,x0,x1,pe,pe0,spe0,scpe0,pe1,spe1,scpe1,dpe,dtau,dbrnch,cn,x180,x360,&
       dtol,ptol,xmin
      integer i, i1, i2, k, j, jb, is, int, nn, mxcnt, mncnt, ios, mod
       
      common /umdc/ pm(jsrc,2), zm(jsrc,2), ndex(jsrc,2), mt(2)
      double precision pm, zm
      integer ndex, mt
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
       coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
       isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
      double precision zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
      real odep, fcs
      integer nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
      
      common/pcdc/phcd(jbrn)
      character*8 phcd
      
      common/prtflc/segmsk(jseg),prnt(2)
      logical segmsk, prnt
      
      data dbrnch,cn,x180,x360,xmin,dtol,ptol/2.5307274d0,57.295779d0,3.1415927d0,6.2831853d0,3.92403d-3,1d-6,2d-6/
!
      if(prnt(1)) write(10,102)
      i1=jndx(jb,1)
      i2=jndx(jb,2)
!     write(10,*)'Spfit:  jb i1 i2 pt =',jb,i1,i2,sngl(pt(i1)), sngl(pt(i2))
      if(i2-i1.gt.1.or.dabs(pt(i2)-pt(i1)).gt.ptol) go to 14
      jndx(jb,2)=-1
      return
 14   newgrd=.false.
      makgrd=.false.
      if(dabs(px(jb,2)-pt(i2)).gt.dtol) newgrd=.true.
!     write(10,*)'Spfit:  px newgrd =',sngl(px(jb,2)),newgrd
      if(.not.newgrd) go to 10
      k=mod(int-1,2)+1
      if(int.ne.int0(k)) makgrd=.true.
!     write(10,*)'Spfit:  int k int0 makgrd =',int,k,int0(k),makgrd
      if(int.gt.2) go to 12
!     call query('Enter xmin:',log)
!     read *,xmin
!     xmin=xmin*xn
      xmin=xn*amin1(amax1(2.*odep,2.),25.)
!     write(10,*)'Spfit:  xmin =',xmin,xmin/xn
      call pdecu(i1,i2,xbrn(jb,1),xbrn(jb,2),xmin,int,i2)
      jndx(jb,2)=i2
 12   nn=i2-i1+1
      if(makgrd) call tauspl(1,nn,pt(i1),tcoef(1,1,k))
!     write(10,301,iostat=ios)jb,k,nn,int,newgrd,makgrd, xbrn(jb,1),xbrn(jb,2),(i,pt(i-1+i1),tau(1,i-1+i1),&
!      (tcoef(j,i,k),j=1,5),i=1,nn)
!301  format(/1x,4i3,2l3,2f12.8/(1x,i5,0p,2f12.8,1p,5d10.2))
      call fitspl(1,nn,tau(1,i1),xbrn(jb,1),xbrn(jb,2),tcoef(1,1,k))
      int0(k)=int
      go to 11
 10   call fitspl(i1,i2,tau,xbrn(jb,1),xbrn(jb,2),coef)
 11   pmn=pt(i1)
      dmn=xbrn(jb,1)
      dmx=dmn
      mxcnt=0
      mncnt=0
!     call appx(i1,i2,xbrn(jb,1),xbrn(jb,2))
!     write(10,300)(i,pt(i),(tau(j,i),j=1,3),i=i1,i2)
!300  format(/(1x,i5,4f12.6))
      pe=pt(i2)
      p1=pt(i1)
      tau1=tau(1,i1)
      x1=tau(2,i1)
      pe1=pe-p1
      spe1=dsqrt(dabs(pe1))
      scpe1=pe1*spe1
      j=i1
      is=i1+1
      do 2 i=is,i2
      p0=p1
      p1=pt(i)
      tau0=tau1
      tau1=tau(1,i)
      x0=x1
      x1=tau(2,i)
      dpe=p0-p1
      dtau=tau1-tau0
      pe0=pe1
      pe1=pe-p1
      spe0=spe1
      spe1=dsqrt(dabs(pe1))
      scpe0=scpe1
      scpe1=pe1*spe1
      tau(4,j)=(2d0*dtau-dpe*(x1+x0))/(.5d0*(scpe1-scpe0)-1.5d0*spe1*spe0*(spe1-spe0))
      tau(3,j)=(dtau-dpe*x0-(scpe1+.5d0*scpe0-1.5d0*pe1*spe0)*tau(4,j))/(dpe*dpe)
      tau(2,j)=(dtau-(pe1*pe1-pe0*pe0)*tau(3,j)-(scpe1-scpe0)*tau(4,j))/dpe
      tau(1,j)=tau0-scpe0*tau(4,j)-pe0*(pe0*tau(3,j)+tau(2,j))
      xlim(1,j)=dmin1(x0,x1)
      xlim(2,j)=dmax1(x0,x1)
      if(xlim(1,j).ge.dmn) go to 5
      dmn=xlim(1,j)
      pmn=pt(j)
      if(x1.lt.x0) pmn=pt(i)
 5    disc=' '
      if(dabs(tau(3,j)).le.1d-30) go to 4
      shm=-.375d0*tau(4,j)/tau(3,j)
      hm=shm*shm
      if(shm.le.0d0.or.(hm.le.pe1.or.hm.ge.pe0)) go to 4
 7    thm=tau(2,j)+shm*(2d0*shm*tau(3,j)+1.5d0*tau(4,j))
      xlim(1,j)=dmin1(xlim(1,j),thm)
      xlim(2,j)=dmax1(xlim(2,j),thm)
      if(thm.ge.dmn) go to 6
      dmn=thm
      pmn=pe-hm
 6    disc='max'
      if(tau(4,j).lt.0d0) disc='min'
      if(disc.eq.'max') mxcnt=mxcnt+1
      if(disc.eq.'min') mncnt=mncnt+1
 4    if(prnt(1)) write(10,100,iostat=ios)disc,j,pt(j),(tau(k,j),k=1,4),(cn*xlim(k,j),k=1,2)
 100  format(1x,a,i5,f10.6,1p,4e10.2,0p,2f7.2)
      dmx=dmax1(dmx,xlim(2,j))
 2    j=i
!     if(prnt(1)) write(10,100,iostat=ios)'   ',j,pt(j)
      xbrn(jb,1)=dmn
      xbrn(jb,2)=dmx
      xbrn(jb,3)=pmn
      idel(jb,1)=1
      idel(jb,2)=1
      if(xbrn(jb,1).gt.x180) idel(jb,1)=2
      if(xbrn(jb,2).gt.x180) idel(jb,2)=2
      if(xbrn(jb,1).gt.x360) idel(jb,1)=3
      if(xbrn(jb,2).gt.x360) idel(jb,2)=3
      if(int.gt.2) go to 1
      phcd(jb)=phcd(jb)(1:1)
      i=jb
      do 8 j=1,nbrn
      i=mod(i,nbrn)+1
      if(phcd(i)(1:1).eq.phcd(jb).and.phcd(i)(2:2).ne.'P'.and.(pe.ge.px(i,1).and.pe.le.px(i,2))) go to 9
 8    continue
      go to 1
 9    phcd(jb)=phcd(i)
      if(dabs(pt(i2)-pt(jndx(i,1))).le.dtol) phcd(jb)=phcd(i-1)
 1    if(prnt(1).and.prnt(2)) write(10,102)
 102  format()
      if(dbrn(jb,1).le.0d0) go to 3
      dbrn(jb,1)=dmx
      dbrn(jb,2)=dbrnch
      if(prnt(2)) write(10,101,iostat=ios)phcd(jb),(jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),&
       (cn*dbrn(jb,k),k=1,2),(idel(jb,k),k=1,3),int,newgrd,makgrd
 101  format(1x,a,2i5,2f8.2,f8.4,2f8.2,4i3,2l2)
      go to 15
 3    if(prnt(2)) write(10,103,iostat=ios)phcd(jb),(jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),&
       (idel(jb,k),k=1,3),int,newgrd,makgrd
 103  format(1x,a,2i5,2f8.2,f8.4,16x,4i3,2l2)
 15   if(mxcnt.gt.mncnt.or.mncnt.gt.mxcnt+1) call warning ('Bad interpolation on '//phcd(jb))
      return
      end


!***********************************************************************************************************************************      
      subroutine pdecu(i1,i2,x0,x1,xmin,int,len)
      
      implicit none
      save 
      include 'ttlim.inc'
      
      double precision x0, x1, xmin, dx, dx2, sgn, rnd, xm, axm, x, h1, h2, hh, xs
      integer i, j, k, m, n, int, i1, i2, is, ie, len
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      common/pdec/ua(5,2),taua(5,2),deplim,ka
      double precision ua, taua
      real deplim
      integer ka
!
!     write(10,*)'Pdecu:  ua =',sngl(ua(1,int))
      if(ua(1,int).le.0d0) go to 17
!     write(10,*)'Pdecu:  fill in new grid'
      k=i1+1
      do 18 i=1,ka
      pt(k)=ua(i,int)
      tau(1,k)=taua(i,int)
 18   k=k+1
      pt(k)=pt(i2)
      tau(1,k)=tau(1,i2)
      go to 19
!
 17   is=i1+1
      ie=i2-1
      xs=x1
      do 11 i=ie,i1,-1
      x=xs
      if(i.ne.i1) go to 12
      xs=x0
      go to 14
 12   h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      xs=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 14   if(dabs(x-xs).le.xmin) go to 15
 11   continue
      len=i2
      return
 15   ie=i
      if(dabs(x-xs).gt..75d0*xmin.or.ie.eq.i2) go to 16
      xs=x
      ie=ie+1
 16   n=max0(idint(dabs(xs-x0)/xmin+.8d0),1)
      dx=(xs-x0)/n
      dx2=dabs(.5d0*dx)
      sgn=dsign(1d0,dx)
      rnd=0d0
      if(sgn.gt.0d0) rnd=1d0
      xm=x0+dx
      k=i1
      m=is
      axm=1d10
      do 1 i=is,ie
      if(i.lt.ie) go to 8
      x=xs
      go to 5
 8    h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      x=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 5    if(sgn*(x-xm).le.dx2) go to 2
      if(k.lt.m) go to 3
      do 4 j=m,k
 4    pt(j)=-1d0
 3    m=k+2
      k=i-1
      axm=1d10
 7    xm=xm+dx*idint((x-xm-dx2)/dx+rnd)
 2    if(dabs(x-xm).ge.axm) go to 1
      axm=dabs(x-xm)
      k=i-1
 1    continue
      if(k.lt.m) go to 9
      do 6 j=m,k
 6    pt(j)=-1d0
 9    k=i1
      do 10 i=is,i2
      if(pt(i).lt.0d0) go to 10
      k=k+1
      pt(k)=pt(i)
      tau(1,k)=tau(1,i)
 10   continue
 19   len=k
!     write(10,300)(i,pt(i),tau(1,i),i=i1,len)
!300  format(/(1x,i5,0p,f12.6,1p,d15.4))
      return
      end


!***********************************************************************************************************************************      
      subroutine tauspl(i1,i2,pt,coef)
!
! $$$$$ calls only library routines $$$$$
!
!   Given ray parameter grid pt;i (pt sub i), i=i1,i1+1,...,i2, tauspl
!   determines the i2-i1+3 basis functions for interpolation I such
!   that:
!
!      tau(p) = a;1,i + Dp * a;2,i + Dp**2 * a;3,i + Dp**(3/2) * a;4,i
!
!   where Dp = pt;n - p, pt;i <= p < pt;i+1, and the a;j,i's are
!   interpolation coefficients.  Rather than returning the coefficients,
!   a;j,i, which necessarily depend on tau(pt;i), i=i1,i1+1,...,i2 and
!   x(pt;i) (= -d tau(p)/d p | pt;i), i=i1,i2, tauspl returns the
!   contribution of each basis function and its derivitive at each
!   sample.  Each basis function is non-zero at three grid points,
!   therefore, each grid point will have contributions (function values
!   and derivitives) from three basis functions.  Due to the basis
!   function normalization, one of the function values will always be
!   one and is not returned in array coef with the other values.
!   Rewritten on 23 December 1983 by R. Buland.
!
      implicit none
      save
      
      double precision pt(i2), coef(5,i2)
      double precision del(5), sdel(5), deli(5), d3h(4), d1h(4), dih(4), d(4), ali, alr, b3h, b1h, bih, th0p, th2p, th3p, th2m
      integer i, j, k, l, m, i1, i2, is, n2
!
      n2=i2-i1-1
      if(n2.le.-1) return
      is=i1+1
!
!   To achieve the requisite stability, proceed by constructing basis
!   functions G;i, i=0,1,...,n+1.  G;i will be non-zero only on the
!   interval [p;i-2,p;i+2] and will be continuous with continuous first
!   and second derivitives.  G;i(p;i-2) and G;i(p;i+2) are constrained
!   to be zero with zero first and second derivitives.  G;i(p;i) is
!   normalized to unity.
!
!   Set up temporary variables appropriate for G;-1.  Note that to get
!   started, the ray parameter grid is extrapolated to yeild p;i, i=-2,
!   -1,0,1,...,n.
      del(2)=pt(i2)-pt(i1)+3d0*(pt(is)-pt(i1))
      sdel(2)=dsqrt(dabs(del(2)))
      deli(2)=1d0/sdel(2)
      m=2
      do 1 k=3,5
      del(k)=pt(i2)-pt(i1)+(5-k)*(pt(is)-pt(i1))
      sdel(k)=dsqrt(dabs(del(k)))
      deli(k)=1d0/sdel(k)
      d3h(m)=del(k)*sdel(k)-del(m)*sdel(m)
      d1h(m)=sdel(k)-sdel(m)
      dih(m)=deli(k)-deli(m)
 1    m=k
      l=i1-1
      if(n2.le.0) go to 10
!   Loop over G;i, i=0,1,...,n-3.
      do 2 i=1,n2
      m=1
!   Update temporary variables for G;i-1.
      do 3 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 3
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 3    m=k
      l=l+1
      del(5)=pt(i2)-pt(l+1)
      sdel(5)=dsqrt(dabs(del(5)))
      deli(5)=1d0/sdel(5)
      d3h(4)=del(5)*sdel(5)-del(4)*sdel(4)
      d1h(4)=sdel(5)-sdel(4)
      dih(4)=deli(5)-deli(4)
!   Construct G;i-1.
      ali=1d0/(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*dih(1)*del(3))*del(3))
      alr=ali*(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*del(3)*deli(2)-sdel(3))*del(3))
      b3h=d3h(2)+alr*d3h(1)
      b1h=d1h(2)+alr*d1h(1)
      bih=dih(2)+alr*dih(1)
      th0p=d1h(1)*b3h-d3h(1)*b1h
      th2p=d1h(3)*b3h-d3h(3)*b1h
      th3p=d1h(4)*b3h-d3h(4)*b1h
      th2m=dih(3)*b3h-d3h(3)*bih
!   The d;i's completely define G;i-1.
      d(4)=ali*((dih(1)*b3h-d3h(1)*bih)*th2p-th2m*th0p)/((dih(4)*b3h-d3h(4)*bih)*th2p-th2m*th3p)
      d(3)=(th0p*ali-th3p*d(4))/th2p
      d(2)=(d3h(1)*ali-d3h(3)*d(3)-d3h(4)*d(4))/b3h
      d(1)=alr*d(2)-ali
!   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
!   G;i-1(p;i-1) need not be constructed as it is normalized to unity.
      coef(1,l)=(.125d0*del(5)*sdel(5)-(.75d0*sdel(5)+.375d0*deli(5)*del(4)-sdel(4))*del(4))*d(4)
      if(i.ge.3) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+.375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
!   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      coef(3,l)=-.75d0*(sdel(5)+deli(5)*del(4)-2d0*sdel(4))*d(4)
      if(i.ge.2) coef(4,l-1)=-.75d0*((sdel(2)+deli(2)*del(3)-2d0*sdel(3))*d(2)-(d1h(1)+dih(1)*del(3))*d(1))
      if(i.ge.3) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-2d0*sdel(2))*d(1)
 2    continue
!   Loop over G;i, i=n-2,n-1,n,n+1.  These cases must be handled
!   seperately because of the singularities in the second derivitive
!   at p;n.
 10   do 4 j=1,4
      m=1
!   Update temporary variables for G;i-1.
      do 5 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 5
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 5    m=k
      l=l+1
      del(5)=0d0
      sdel(5)=0d0
      deli(5)=0d0
!   Construction of the d;i's is different for each case.  In cases
!   G;i, i=n-1,n,n+1, G;i is truncated at p;n to avoid patching across
!   the singularity in the second derivitive.
      if(j.lt.4) go to 6
!   For G;n+1 constrain G;n+1(p;n) to be .25.
      d(1)=2d0/(del(1)*sdel(1))
      go to 9
!   For G;i, i=n-2,n-1,n, the condition dG;i(p)/dp|p;i = 0 has been
!   substituted for the second derivitive continuity condition that
!   can no longer be satisfied.
 6    alr=(sdel(2)+deli(2)*del(3)-2d0*sdel(3))/(d1h(1)+dih(1)*del(3))
      d(2)=1d0/(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*deli(2)*del(3)-sdel(3))*del(3)-(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*&
       dih(1)*del(3))*del(3))*alr)
      d(1)=alr*d(2)
      if(j-2)8,7,9
!   For G;n-1 constrain G;n-1(p;n) to be .25.
 7    d(3)=(2d0+d3h(2)*d(2)+d3h(1)*d(1))/(del(3)*sdel(3))
      go to 9
!   No additional constraints are required for G;n-2.
 8    d(3)=-((d3h(2)-d1h(2)*del(4))*d(2)+(d3h(1)-d1h(1)*del(4))*d(1))/(d3h(3)-d1h(3)*del(4))
      d(4)=(d3h(3)*d(3)+d3h(2)*d(2)+d3h(1)*d(1))/(del(4)*sdel(4))
!   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
 9    if(j.le.2) coef(1,l)=(.125d0*del(3)*sdel(3)-(.75d0*sdel(3)+.375d0*deli(3)*del(4)-sdel(4))*del(4))*d(3)-&
       (.125d0*d3h(2)-(.75d0*d1h(2)+.375d0*dih(2)*del(4))*del(4))*d(2)-(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*dih(1)*del(4))*del(4))*&
       d(1)
      if(l-i1.gt.1) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+.375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
!   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      if(j.le.2) coef(3,l)=-.75d0*((sdel(3)+deli(3)*del(4)-2d0*sdel(4))*d(3)-(d1h(2)+dih(2)*del(4))*d(2)-(d1h(1)+&
       dih(1)*del(4))*d(1))
      if(j.le.3.and.l-i1.gt.0) coef(4,l-1)=0d0
      if(l-i1.gt.1) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-2d0*sdel(2))*d(1)
 4    continue
      return
      end


!***********************************************************************************************************************************      
      subroutine fitspl(i1,i2,tau,x1,xn,coef)
!
! $$$$$ calls only library routines $$$$$
!
!   Given ray parameter grid p;i (p sub i), i=1,2,...,n, corresponding
!   tau;i values, and x;1 and x;n (x;i = -dtau/dp|p;i); tauspl finds
!   interpolation I such that:  tau(p) = a;1,i + Dp * a;2,i + Dp**2 *
!   a;3,i + Dp**(3/2) * a;4,i where Dp = p;n - p and p;i <= p < p;i+1.
!   Interpolation I has the following properties:  1) x;1, x;n, and
!   tau;i, i=1,2,...,n are fit exactly, 2) the first and second
!   derivitives with respect to p are continuous everywhere, and
!   3) because of the paramaterization d**2 tau/dp**2|p;n is infinite.
!   Thus, interpolation I models the asymptotic behavior of tau(p)
!   when tau(p;n) is a branch end due to a discontinuity in the
!   velocity model.  Note that array a must be dimensioned at least
!   a(4,n) though the interpolation coefficients will be returned in
!   the first n-1 columns.  The remaining column is used as scratch
!   space and returned as all zeros.  Programmed on 16 August 1982 by
!   R. Buland.
!
      implicit none
      save 
      
      double precision tau(4,i2), x1, xn, coef(5,i2), a(2,100), ap(3), b(100), alr, g1, gn
      integer i, j, n, i1, i2, ie, is, n1
!
      if(i2-i1)13,1,2
 1    tau(2,i1)=x1
 13   return
 2    n=0
      do 3 i=i1,i2
      n=n+1
      b(n)=tau(1,i)
      do 3 j=1,2
 3    a(j,n)=coef(j,i)
      do 4 j=1,3
 4    ap(j)=coef(j+2,i2)
      n1=n-1
!
!   Arrays ap(*,1), a, and ap(*,2) comprise n+2 x n+2 penta-diagonal
!   matrix A.  Let x1, tau, and xn comprise corresponding n+2 vector b.
!   Then, A * g = b, may be solved for n+2 vector g such that
!   interpolation I is given by I(p) = sum(i=0,n+1) g;i * G;i(p).
!
!   Eliminate the lower triangular portion of A to form A'.  A
!   corresponding transformation applied to vector b is stored in
!   a(4,*).
      alr=a(1,1)/coef(3,i1)
      a(1,1)=1d0-coef(4,i1)*alr
      a(2,1)=a(2,1)-coef(5,i1)*alr
      b(1)=b(1)-x1*alr
      j=1
      do 5 i=2,n
      alr=a(1,i)/a(1,j)
      a(1,i)=1d0-a(2,j)*alr
      b(i)=b(i)-b(j)*alr
 5    j=i
      alr=ap(1)/a(1,n1)
      ap(2)=ap(2)-a(2,n1)*alr
      gn=xn-b(n1)*alr
      alr=ap(2)/a(1,n)
!   Back solve the upper triangular portion of A' for coefficients g;i.
!   When finished, storage g(2), a(4,*), g(5) will comprise vector g.
      gn=(gn-b(n)*alr)/(ap(3)-a(2,n)*alr)
      b(n)=(b(n)-gn*a(2,n))/a(1,n)
      j=n
      do 6 i=n1,1,-1
      b(i)=(b(i)-b(j)*a(2,i))/a(1,i)
 6    j=i
      g1=(x1-coef(4,i1)*b(1)-coef(5,i1)*b(2))/coef(3,i1)
!
      tau(2,i1)=x1
      is=i1+1
      ie=i2-1
      j=1
      do 7 i=is,ie
      j=j+1
 7    tau(2,i)=coef(3,i)*b(j-1)+coef(4,i)*b(j)+coef(5,i)*b(j+1)
      tau(2,i2)=xn
      return
      end


!***********************************************************************************************************************************      
      subroutine trtm(delta,max,n,tt,dtdd,dtdh,dddp,phnm)
      
      implicit none
      save 
      include 'ttlim.inc'
      
      character*(*) phnm(max)
      character*8 ctmp(60)
      double precision x(3), cn, dtol, pi, pi2
      real tt(max), dtdd(max), dtdh(max), dddp(max), tmp(60,4), atol, delta
      integer iptr(60), i, j, k, n, max
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      data cn,dtol,atol,pi,pi2/.017453292519943296d0,1d-6,.005,3.1415926535897932d0,6.2831853071795865d0/
!
      n=0
      if(mbr2.le.0) return
      x(1)=dmod(dabs(cn*delta),pi2)
      if(x(1).gt.pi) x(1)=pi2-x(1)
      x(2)=pi2-x(1)
      x(3)=x(1)+pi2
      if(dabs(x(1)).gt.dtol) go to 9
      x(1)=dtol
      x(3)=-10d0
 9    if(dabs(x(1)-pi).gt.dtol) go to 7
      x(1)=pi-dtol
      x(2)=-10d0
 7    do 1 j=mbr1,mbr2
 1    if(jndx(j,2).gt.0) call findtt(j,x,max,n,tmp,tmp(1,2),tmp(1,3),tmp(1,4),ctmp)
      if(n-1)3,4,5
 4    iptr(1)=1
      go to 6
 5    call r4sort(n,tmp,iptr)
 6    k=0
      do 2 i=1,n
      j=iptr(i)
      if(k.le.0) go to 8
      if(phnm(k).eq.ctmp(j).and.abs(tt(k)-tmp(j,1)).le.atol) go to 2
 8    k=k+1
      tt(k)=tmp(j,1)
      dtdd(k)=tmp(j,2)
      dtdh(k)=tmp(j,3)
      dddp(k)=tmp(j,4)
      phnm(k)=ctmp(j)
 2    continue
      n=k
 3    return
      end


!***********************************************************************************************************************************      
      subroutine findtt(jb,x0,max,n,tt,dtdd,dtdh,dddp,phnm)
      
      implicit none
      save 
      include 'ttlim.inc'
      
      character*(*) phnm(max)
      character*67 msg
      double precision x, x0(3), p0, p1, arg, dp, dps, delp, tol, ps, deps
      real tt(max), dtdd(max), dtdh(max), dddp(max), hsgn, dsgn, dpn, dp0
      integer nph, jb, i, ij, ie, is, le, n, in, j, jj, ln, max
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      common/pcdc/phcd(jbrn)
      character*8 phcd
      
      data tol/3d-6/,deps/1d-10/
!
      nph=iabs(idel(jb,3))
      hsgn=isign(1,idel(jb,3))*hn
      dsgn=(-1.)**idel(jb,1)*dn
      dpn=-1./tn
      do 10 ij=idel(jb,1),idel(jb,2)
      x=x0(ij)
      dsgn=-dsgn
      if(x.lt.xbrn(jb,1).or.x.gt.xbrn(jb,2)) go to 12
      j=jndx(jb,1)
      is=j+1
      ie=jndx(jb,2)
      do 1 i=is,ie
      if(x.le.xlim(1,j).or.x.gt.xlim(2,j)) go to 8
      le=n
      p0=pt(ie)-pt(j)
      p1=pt(ie)-pt(i)
      delp=dmax1(tol*(pt(i)-pt(j)),1d-3)
      if(dabs(tau(3,j)).gt.1d-30) go to 2
      dps=(x-tau(2,j))/(1.5d0*tau(4,j))
      dp=dsign(dps*dps,dps)
      dp0=sngl(dp)
      if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 9
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=sngl(dsgn*ps)
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*sngl(.75d0*tau(4,j)/dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 8
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
      go to 8
 2    do 4 jj=1,2
      go to (5,6),jj
 5    arg=9d0*tau(4,j)*tau(4,j)+32d0*tau(3,j)*(x-tau(2,j))
      if(arg.ge.0d0) go to 3
      write(msg,100)arg
 100  format('Bad sqrt argument:',1p,d11.2,'.')
      call warning (msg(1:30))
      stop
 3    dps=-(3d0*tau(4,j)+dsign(dsqrt(dabs(arg)),tau(4,j)))/(8d0*tau(3,j))
      dp=dsign(dps*dps,dps)
      dp0=sngl(dp)
      go to 7
 6    dps=(tau(2,j)-x)/(2d0*tau(3,j)*dps)
      dp=dsign(dps*dps,dps)
 7    if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 4
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=sngl(dsgn*ps)
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*sngl(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 4
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
 4    continue
 9    if(n.gt.le) go to 8
      write(msg,101)phcd(jb),x,dp0,dp,p1,p0
 101  format('Failed to find phase:  ',a,f8.1,4f7.4)
      call warning (msg)
 8    j=i
 1    continue
!
 12   if(x.lt.dbrn(jb,1).or.x.gt.dbrn(jb,2)) go to 10
      if(n.ge.max) go to 13
      j=jndx(jb,1)
      i=jndx(jb,2)
      dp=pt(i)-pt(j)
      dps=dsqrt(dabs(dp))
      n=n+1
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+pt(j)*x)
      dtdd(n)=dsgn*sngl(pt(j))
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-pt(j)*pt(j))))
      dddp(n)=dpn*sngl(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dps,deps))
      ln=index(phcd(jb),' ')-1
      if(ln.le.0) ln=len(phcd(jb))
      phnm(n)=phcd(jb)(1:ln)//'diff'
 10   continue
      return
 13   write(msg,102)max
 102  format('More than',i3,' arrivals found.')
      call warning (msg(1:28))
      return
      end


!***********************************************************************************************************************************      
      subroutine query(ia,log)
!
! $$$$$ calls tnoua $$$$$
!
!   Subroutine query scans character string ia (up to 78 characters) for
!   a question mark or a colon.  It prints the string up to and
!   including the flag character plus two blanks with no newline on the
!   standard output.  If the flag was a question mark, query reads the
!   users response.  If the response is 'y' or 'yes', log is set to
!   true.  If the response is 'n' or 'no', log is set to false.  Any
!   other response causes the question to be repeated.  If the flag was
!   a colon, query simply returns allowing user input on the same line.
!   If there is no question mark or colon, the last non-blank character
!   is treated as if it were a colon.  If the string is null or all
!   blank, query prints an astrisk and returns.  Programmed on 3
!   December 1979 by R. Buland.
!
      implicit none
      save
      
      logical log
      character*(*) ia
      character*81 ib
      character*4 ans
      integer i, k, ifl, nn
      
      nn=len(ia)
      log=.true.
      ifl=1
      k=0
!   Scan ia for flag characters or end-of-string.
      do 1 i=1,nn
      ib(i:i)=ia(i:i)
      if(ib(i:i).eq.':') go to 7
      if(ib(i:i).eq.'?') go to 3
      if(ib(i:i).eq.'\0') go to 5
      if(ib(i:i).ne.' ') k=i
 1    continue
!   If we fell off the end of the string, branch if there were any non-
!   blank characters.
 5    if(k.gt.0) go to 6
!   Handle a null or all blank string.
      i=1
      ib(i:i)='*'
      go to 4
!   Handle a string with no question mark or colon but at least one
!   non-blank character.
 6    i=k
!   Append two blanks and print the string.
 7    i=i+2
      ib(i-1:i-1)=' '
      ib(i:i)=' '
!   Tnoua prints the first i characters of ib without a newline.
 4    call tnoua(ib,i)
      if(ifl.gt.0) return
!   If the string was a yes-no question read the response.
      read 102,ans
 102  format(a4)
      call uctolc(ans,-1)
!   If the response is yes log is already set properly.
      if(ans.eq.'y   '.or.ans.eq.'yes ') return
!   If the response is no set log to false.  Otherwise repeat the
!   question.
      if(ans.ne.'n   '.and.ans.ne.'no  ') go to 4
      log=.false.
      return
 3    ifl=-ifl
      go to 7
      end


!***********************************************************************************************************************************      
      subroutine uctolc(ia,ifl)
!
! $$$$$ calls only library routines $$$$$
!
!   Subroutine uctolc converts alphabetic characters in string ia from
!   upper case to lower case.  If ifl < 0 all characters are converted.
!   Otherwise characters enclosed by single quotes are left unchanged.
!   Programmed on 21 January by R. Buland.  Calling sequence changed
!   on 11 December 1985 by R. Buland.
!
      implicit none
      
      character*(*) ia
      integer i, ib, n, ifl, nfl
      
      data nfl/1/
      
      if(ifl.lt.0) nfl=1
!   Scan the string.
      n=len(ia)
      do 1 i=1,n
      if(ifl.lt.0) go to 2
!   Look for single quotes.
      if(ia(i:i).eq.'''') nfl=-nfl
!   If we are in a quoted string skip the conversion.
      if(nfl.lt.0) go to 1
!   Do the conversion.
 2    ib=ichar(ia(i:i))
      if(ib.lt.65.or.ib.gt.90) go to 1
      ia(i:i)=char(ib+32)
 1    continue
      return
      end


!***********************************************************************************************************************************      
      subroutine r4sort(n,rkey,iptr)
!
! $$$$$ calls no other routine $$$$$
!
!   R4sort sorts the n elements of array rkey so that rkey(i), 
!   i = 1, 2, 3, ..., n are in asending order.  R4sort is a trivial
!   modification of ACM algorithm 347:  "An efficient algorithm for
!   sorting with minimal storage" by R. C. Singleton.  Array rkey is
!   sorted in place in order n*alog2(n) operations.  Coded on
!   8 March 1979 by R. Buland.  Modified to handle real*4 data on
!   27 September 1983 by R. Buland.
!
      implicit none
      save
      
      real rkey(n), r, tmpkey
      integer iptr(n), il(10), iu(10), i, j, k, l, m, n, ij, it, kk, ib
      
!   Note:  il and iu implement a stack containing the upper and
!   lower limits of subsequences to be sorted independently.  A
!   depth of k allows for n<=2**(k+1)-1.

      if(n.le.0) return
      do 1 i=1,n
 1    iptr(i)=i
      if(n.le.1) return
      r=.375
      m=1
      i=1
      j=n
!
!   The first section interchanges low element i, middle element ij,
!   and high element j so they are in order.
!
 5    if(i.ge.j) go to 70
 10   k=i
!   Use a floating point modification, r, of Singleton's bisection
!   strategy (suggested by R. Peto in his verification of the
!   algorithm for the ACM).
      if(r.gt..58984375) go to 11
      r=r+.0390625
      go to 12
 11   r=r-.21875
 12   ij=i+(j-i)*r
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 20
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 20   l=j
      if(rkey(iptr(j)).ge.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(j)
      iptr(j)=it
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 39   tmpkey=rkey(iptr(ij))
      go to 40
!
!   The second section continues this process.  K counts up from i and
!   l down from j.  Each time the k element is bigger than the ij
!   and the l element is less than the ij, then interchange the
!   k and l elements.  This continues until k and l meet.
!
 30   it=iptr(l)
      iptr(l)=iptr(k)
      iptr(k)=it
 40   l=l-1
      if(rkey(iptr(l)).gt.tmpkey) go to 40
 50   k=k+1
      if(rkey(iptr(k)).lt.tmpkey) go to 50
      if(k.le.l) go to 30
!
!   The third section considers the intervals i to l and k to j.  The
!   larger interval is saved on the stack (il and iu) and the smaller
!   is remapped into i and j for another shot at section one.
!
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
!
!   The fourth section pops elements off the stack (into i and j).  If
!   necessary control is transfered back to section one for more
!   interchange sorting.  If not we fall through to section five.  Note
!   that the algorighm exits when the stack is empty.
!
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.11) go to 10
      if(i.eq.1) go to 5
      i=i-1
!
!   The fifth section is the end game.  Final sorting is accomplished
!   (within each subsequence popped off the stack) by rippling out
!   of order elements down to their proper positions.
!
 90   i=i+1
      if(i.eq.j) go to 70
      if(rkey(iptr(i)).le.rkey(iptr(i+1))) go to 90
      k=i
      kk=k+1
      ib=iptr(kk)
 100  iptr(kk)=iptr(k)
      kk=k
      k=k-1
      if(rkey(ib).lt.rkey(iptr(k))) go to 100
      iptr(kk)=ib
      go to 90
      end


!***********************************************************************************************************************************      
      integer function iupcor(phnm,dtdd,xcor,tcor)
      
      implicit none
      save
      include 'ttlim.inc'
      
      character*(*) phnm
      double precision x, dp, dps, ps, cn
      real oldep, dtdd, xcor, tcor, eps
      integer jp, js, ie, i, j, is, jb
      
      common /tabc/ us(2), pt(jout), tau(4,jout), xlim(2,jout), xbrn(jbrn,3), dbrn(jbrn,2), xn, pn, tn, dn, hn, jndx(jbrn,2),&
       idel(jbrn,3), mbr1, mbr2
      double precision us, pt, tau, xlim, xbrn, dbrn 
      real xn, pn, tn, dn, hn
      integer jndx, idel, mbr1, mbr2
       
      common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
       coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
       isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
      double precision zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
      real odep, fcs
      integer nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
      
      common/pcdc/phcd(jbrn)
      character*8 phcd
      
      data oldep,jp,js/-1.,2*0/,cn/57.295779d0/
      
      eps = 2.*epsilon(oldep)
!
      iupcor=1
!     print *,'oldep odep',oldep,odep
      if(abs(oldep-odep) .lt. eps) go to 1
      oldep=odep
!   Find the upgoing P branch.
!     print *,'mbr1 mbr2',mbr1,mbr2
      do 2 jp=mbr1,mbr2
!     print *,'jp phcd xbrn',jp,'  ',phcd(jp),xbrn(jp,1)
      if((phcd(jp).eq.'Pg'.or.phcd(jp).eq.'Pb'.or.phcd(jp).eq.'Pn'.or.phcd(jp).eq.'P').and.xbrn(jp,1).le.0d0) go to 3
 2    continue
      jp=0
!   Find the upgoing S branch.
 3    do 4 js=mbr1,mbr2
!     print *,'js phcd xbrn',js,'  ',phcd(js),xbrn(js,1)
      if((phcd(js).eq.'Sg'.or.phcd(js).eq.'Sb'.or.phcd(js).eq.'Sn'.or.phcd(js).eq.'S').and.xbrn(js,1).le.0d0) go to 1
 4    continue
      js=0
!
!1    print *,'jp js',jp,js
 1    if(phnm.ne.'P'.and.phnm.ne.'p') go to 5
      jb=jp
      if(jb)14,14,6
!
 5    if(phnm.ne.'S'.and.phnm.ne.'s') go to 13
      jb=js
      if(jb)14,14,6
!
 6    is=jndx(jb,1)+1
      ie=jndx(jb,2)
      ps=abs(dtdd)/dn
!     print *,'jb is ie dtdd dn ps',jb,is,ie,dtdd,dn,ps
      if(ps.lt.pt(is-1).or.ps.gt.pt(ie)) go to 13
      do 7 i=is,ie
!     print *,'i pt',i,pt(i)
      if(ps.le.pt(i)) go to 8
 7    continue
      go to 13
!
 8    j=i-1
      dp=pt(ie)-ps
      dps=dsqrt(dabs(dp))
      x=tau(2,j)+2d0*dp*tau(3,j)+1.5d0*dps*tau(4,j)
!     print *,'j pt dp dps x',j,pt(ie),dp,dps,x
      tcor=tn*sngl(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      xcor=sngl(cn*x)
!     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
!
 13   iupcor=-1
 14   xcor=0.
      tcor=0.
!     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
      end
      

!***********************************************************************************************************************************      
      subroutine brnset(nn,pcntl,prflg)
!
!   Brnset takes character array pcntl(nn) as a list of nn tokens to be
!   used to select desired generic branches.  Prflg(3) is the old
!   prnt(2) debug print flags in the first two elements plus a new print
!   flag which controls a branch selected summary from brnset.  Note that
!   the original two flags controlled a list of all tau interpolations
!   and a branch range summary respectively.  The original summary output
!   still goes to logical unit 10 (ttim1.lis) while the new output goes
!   to the standard output (so the caller can see what happened).  Each
!   token of pcntl may be either a generic branch name (e.g., P, PcP,
!   PKP, etc.) or a keyword (defined in the data statement for cmdcd
!   below) which translates to more than one generic branch names.  Note
!   that generic branch names and keywords may be mixed.  The keywords
!   'all' (for all branches) and 'query' (for an interactive token input
!   query mode) are also available.
!
      implicit none
      save
      include 'ttlim.inc'
      
      integer ncmd, lcmd
      parameter(ncmd=4,lcmd=16)
      
      logical prflg(3), fnd, all
      character*(*) pcntl(nn)
      character*8 segcd(jbrn), cmdcd(ncmd), cmdlst(lcmd), phtmp, phlst(jseg)
      integer nsgpt(jbrn), ncmpt(2,ncmd), no, i, kseg, j, l, k, j1, j2, nn 
      
      common /brkc/ zs, pk(jseg), pu(jtsm0,2), pux(jxsm,2), tauu(jtsm,2), xu(jxsm,2), px(jbrn,2), xt(jbrn,2), taut(jout),&
       coef(5,jout), tauc(jtsm), xc(jxsm), tcoef(5,jbrna,2), tp(jbrnu,2), odep, fcs(jseg,3), nin, nph0, int0(2), ki, msrc(2),&
       isrc(2), nseg, nbrn, ku(2), km(2), nafl(jseg,3), indx(jseg,2), kndx(jseg,2), iidx(jseg), jidx(jbrn), kk(jseg)
      double precision zs, pk, pu, pux, tauu, xu, px, xt, taut, coef, tauc, xc, tcoef, tp
      real odep, fcs
      integer nin, nph0, int0, ki, msrc, isrc, nseg, nbrn, ku, km, nafl, indx, kndx, iidx, jidx, kk
      
      common /pcdc/ phcd(jbrn)
      character*8 phcd
      
!   Segmsk is a logical array that actually implements the branch
!   editing in depset and depcor.

      common /prtflc/ segmsk(jseg), prnt(2)
      logical segmsk, prnt
!
!   The keywords do the following:
!      P      gives P-up, P, Pdiff, PKP, and PKiKP
!      P+     gives P-up, P, Pdiff, PKP, PKiKP, PcP, pP, pPdiff, pPKP,
!             pPKiKP, sP, sPdiff, sPKP, and sPKiKP
!      S+     gives S-up, S, Sdiff, SKS, sS, sSdiff, sSKS, pS, pSdiff,
!             and pSKS
!      basic  gives P+ and S+ as well as ScP, SKP, PKKP, SKKP, PP, and
!             P'P'
!   Note that generic S gives S-up, Sdiff, and SKS already and so
!   doesn't require a keyword.
!
      data cmdcd/'P','P+','basic','S+'/
      data cmdlst/'P','PKiKP','PcP','pP','pPKiKP','sP','sPKiKP','ScP','SKP','PKKP','SKKP','PP','S','ScS','sS','pS'/
      data ncmpt/1,2,1,7,1,13,13,16/
!
!   Take care of the print flags.
      prnt(1)=prflg(1)
      prnt(2)=prflg(2)
      if(prnt(1)) prnt(2)=.true.
!   Copy the token list into local storage.
      no=min0(nn,jseg)
      do 23 i=1,no
 23   phlst(i)=pcntl(i)
!   See if we are in query mode.
      if(no.gt.1.or.(phlst(1).ne.'query'.and.phlst(1).ne.'QUERY')) go to 1
!
!   In query mode, get the tokens interactively into local storage.
!
 22   print *,'Enter desired branch control list at the prompts:'
      no=0
 21   call query(' ',fnd)
      if(no.ge.jseg) go to 1
      no=no+1
      read 100,phlst(no)
 100  format(a)
!   Terminate the list of tokens with a blank entry.
      if(phlst(no).ne.' ') go to 21
      no=no-1
      if(no.gt.0) go to 1
!   If the first token is blank, help the user out.
      print *,'You must enter some branch control information!'
      print *,'     possibilities are:'
      print *,'          all'
      print 101,cmdcd
 101  format(11x,a)
      print *,'          or any generic phase name'
      go to 22
!
!   An 'all' keyword is easy as this is already the default.
 1    all=.false.
      if(no.eq.1.and.(phlst(1).eq.'all'.or.phlst(1).eq.'ALL')) all=.true.
      if(all.and..not.prflg(3)) return
!
!   Make one or two generic branch names for each segment.  For example,
!   the P segment will have the names P and PKP, the PcP segment will
!   have the name PcP, etc.
!
      kseg=0
      j=0
!   Loop over the segments.
      do 2 i=1,nseg
      if(.not.all) segmsk(i)=.false.
!   For each segment, loop over associated branches.
 9    j=j+1
      phtmp=phcd(j)
!   Turn the specific branch name into a generic name by stripping out
!   the crustal branch and core phase branch identifiers.
      do 3 l=2,8
 6    if(phtmp(l:l).eq.' ') go to 4
      if(phtmp(l:l).ne.'g'.and.phtmp(l:l).ne.'b'.and.phtmp(l:l).ne.'n') go to 5
      if(l.lt.8) phtmp(l:)=phtmp(l+1:)
      if(l.ge.8) phtmp(l:)=' '
      go to 6
 5    if(l.ge.8) go to 3
      if(phtmp(l:l+1).ne.'ab'.and.phtmp(l:l+1).ne.'ac'.and. phtmp(l:l+1).ne.'df') go to 3
      phtmp(l:)=' '
      go to 4
 3    continue
!4    print *,'j phcd phtmp =',j,' ',phcd(j),' ',phtmp
!
!   Make sure generic names are unique within a segment.
 4    if(kseg.lt.1) go to 7
      if(phtmp.eq.segcd(kseg)) go to 8
 7    kseg=kseg+1
      segcd(kseg)=phtmp
      nsgpt(kseg)=i
!     if(prflg(3)) print *,'kseg nsgpt segcd =',kseg,nsgpt(kseg),' ', segcd(kseg)
 8    if(jidx(j).lt.indx(i,2)) go to 9
 2    continue
      if(all) go to 24
!
!   Interpret the tokens in terms of the generic branch names.
!
      do 10 i=1,no
!   Try for a keyword first.
      do 11 j=1,ncmd
      if(phlst(i).eq.cmdcd(j)) go to 12
 11   continue
!
!   If the token isn't a keyword, see if it is a generic branch name.
      fnd=.false.
      do 14 k=1,kseg
      if(phlst(i).ne.segcd(k)) go to 14
      fnd=.true.
      l=nsgpt(k)
      segmsk(l)=.true.
!     print *,'Brnset:  phase found - i k l segcd =',i,k,l,' ', segcd(k)
 14   continue
!   If no matching entry is found, warn the caller.
      if(.not.fnd) print *,'Brnset:  phase ',phlst(i),' not found.'
      go to 10
!
!   If the token is a keyword, find the matching generic branch names.
 12   j1=ncmpt(1,j)
      j2=ncmpt(2,j)
      do 15 j=j1,j2
      do 15 k=1,kseg
      if(cmdlst(j).ne.segcd(k)) go to 15
      l=nsgpt(k)
      segmsk(l)=.true.
!     print *,'Brnset:  cmdlst found - j k l segcd =',j,k,l,' ', segcd(k)
 15   continue
 10   continue
!
!   Make the caller a list of the generic branch names selected.
!
 24   if(.not.prflg(3)) return
      fnd=.false.
      j2=0
!   Loop over segments.
      do 16 i=1,nseg
      if(.not.segmsk(i)) go to 16
!   If selected, find the associated generic branch names.
      j2=j2+1
      do 17 j1=j2,kseg
      if(nsgpt(j1).eq.i) go to 18
 17   continue
      print *,'Brnset:  Segment pointer (',i,') missing?'
      go to 16
 18   do 19 j2=j1,kseg
      if(nsgpt(j2).ne.i) go to 20
 19   continue
      j2=kseg+1
!   Print the result.
 20   j2=j2-1
      if(.not.fnd .and. verbose_taup) print *,'Brnset:  the following phases have '//'been selected -'
      fnd=.true.
      if (verbose_taup) print 102,i,(segcd(j),j=j1,j2)
 102  format(10x,i5,5(2x,a))
 16   continue
      return
      end


!***********************************************************************************************************************************      
      subroutine warning (msg)
      
      ! Print a message
      
      implicit none
      
      character*(*) msg
      
      write (*,'(1x,a)') msg
 
      return
      end


!***********************************************************************************************************************************      
      subroutine tnoua (ia, n)

      ! Writes the first n characters of string ia to the
      ! standard output without the trailing newline (allowing user input
      ! on the same line).  Programmed on 17 September 1980 by R. Buland.
      ! Modified for the "advance" ioflag by EAB, 12/2/07.
      
      implicit none

      character*(*) ia
      integer n
      
      write (*,'(a)',advance='no') ia(1:n)
      
      return
      end

!***********************************************************************************************************************************      
      subroutine abort1 (msg)
      
      ! Print a message, close all open logical units, and stop
      
      implicit none
      
      character*(*) msg
      
      write (*,'(1x,a)') msg
      stop
      return
      end
