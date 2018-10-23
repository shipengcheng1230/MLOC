      integer nmax, nmax2
      parameter (nmax = nevmax)
      parameter (nmax2 = nmax*2)

      common /rsv/ rsvx(nevmax), rsvy(nevmax)
      real rsvx, rsvy

      common /tst/ covs(nmax2,nmax2), cvs(nmax,2)
      real covs, cvs
