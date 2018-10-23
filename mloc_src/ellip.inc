      integer numph, mdel
      parameter (numph = 57)
      parameter (mdel = 50)



      common /tau1/ phcod(numph)
      character*8 phcod
      
      common /tau2/ phnch(numph), np(numph)
      integer phnch, np
            
      common /tau3/ dpth(6), delta(numph,mdel), di1(numph), di2(numph), t0(numph,mdel,6), t1(numph,mdel,6), t2(numph,mdel,6)
      real dpth, delta, di1, di2, t0, t1, t2