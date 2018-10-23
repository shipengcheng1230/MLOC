      integer, parameter :: max_layers = 60
      integer, parameter :: np = 40

      real :: v0, h, zmax, zoold, rmax
      integer :: jmod, imoh, icon, ipd, isd
      character(len=4) :: azo
      common /model/  v0(2,max_layers), h(max_layers), rmax, zmax, jmod, azo(max_layers), imoh, icon, ipd, isd, zoold
      
      real :: rr, tt, pa, fa, ttp, ppp, pim, aa, vhq, pa1, piref, v, g, v2, del
      integer :: ion, ib, iqq, iql, jh1
      character(len=8) :: phaseref
      logical :: conv
      common /ref2/ rr(2), tt(2), pa(2), fa(2,max_layers), ttp(np,3,3), ppp(np,3,3), ion(np,3,3), ib, iqq, iql, jh1, pim, aa,&
       vhq(2), pa1(2), piref, v(2,max_layers), g(2,max_layers), v2(2,max_layers), del(3), phaseref, conv

      character(len=8) :: phlist
      common /phasel/ phlist(np)
