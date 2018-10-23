!c   The "j" parameters (1 st line) are intended to be user settable:
!c        jsrc   Maximum number of discrete model slowness samples above
!c               the maximum source depth of interest.
!c        jseg   Maximum number of different types of travel-times
!c               considered.
!c        jbrn   Maximum number of different travel-time branches to be
!c               searched.
!c        jout   Maximum length of all travel-time branches strung
!c               together.
!c        jtsm   Maximum length of the tau depth increments.
!c        jxsm   Maximum number of x-values needed for the depth
!c               increments.
!c        jbrnu  Maximum length of the up-going branches.
!c        jbrna  Maximum length of branches which may need
!c               re-interpolation.
!c
      integer jsrc, jseg, jbrn, jout, jtsm, jxsm, jbrnu, jbrna
      parameter(jsrc=150,jseg=30,jbrn=100,jout=2500,jtsm=350)
      parameter(jxsm=jbrn,jbrnu=jbrn,jbrna=jbrn)
     
!c   A few derived parameters are also needed.

      integer jrec, jtsm0
      parameter(jrec=jtsm+jxsm,jtsm0=jtsm+1)
      
      common /ttim/ iottim
      integer iottim

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
      
      common /pcdc/ phcd(jbrn)
      character*8 phcd
            
      common /prtflc/ segmsk(jseg), prnt(2)
      logical segmsk, prnt

      common /pdec/ ua(5,2), taua(5,2), deplim, ka
      double precision ua, taua
      real deplim
      integer ka