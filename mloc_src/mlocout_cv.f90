!***********************************************************************************************************************************
      subroutine mlocout_cv (it)
      
      ! Write spatial components of cluster vectors and covariance matrix
      
      implicit none
      
      include 'mloc.inc'
      
      character*100 outfil
      character(len=132) :: msg
      real covs(mtmax2,mtmax2), al, bl, alpha2, dgkmlo, cv1f, cv2f
      integer it, i, j, iev, jmt, lmt, jev, imt, kmt, nev2, it1
      
      it1 = it + 1
      
      outfil = trim(outfile)//'.cv'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_cv: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      open (io_out,file=outfil,status='new')
      
      write (io_out,'(i4)') nev
      
      ! Confidence ellipses, unscaled and scaled semi-axes
      do iev = 1,nev 
         al = 1./alic(iev)
         bl = 1./blic(iev)
         alpha2 = alphac(iev) + 90.
         if (alpha2 .gt. 360.) alpha2 = alpha2 - 360.
         write (io_out,'(i4,5f10.3,i6,3f10.3)') iev, alphac(iev), alpha2, al, bl, eciev(iev,it), ndatc(iev,it), kcritc(iev),&
          xl1c(iev), xl2c(iev)
      end do
      
      ! Covariance matrix for each event
      do iev = 1,nev 
         write (io_out,'(i4,10f10.3)') iev, ccv(iev,1,1), ccv(iev,2,2), ccv(iev,3,3), ccv(iev,4,4),&
          ccv(iev,1,2), ccv(iev,1,3), ccv(iev,1,4), ccv(iev,2,3), ccv(iev,2,4), ccv(iev,3,4)
      end do
      
      ! Cluster vectors, in km
      do iev = 1,nev 
         cv1f = (latp(iev,it1) - lath(it1))/dgkmla
         cv2f = (lonp(iev,it1) - lonh(it1))/dgkmlo(lath(it1))
         write (io_out,'(i4,2f10.3)') iev, cv1f, cv2f
      end do
      
      ! Full covariance matrix for epicenters
      jmt = 0
      lmt = 0
      do jev = 1,nev
         imt = 0
         kmt = 0
         do iev = 1,nev
            covs(kmt+1,lmt+1) = sngl(vhatc(imt+1,jmt+1))
            covs(kmt+1,lmt+2) = sngl(-vhatc(imt+1,jmt+2)) ! Sign change for geocentric latitude
            covs(kmt+2,lmt+1) = sngl(-vhatc(imt+2,jmt+1)) ! Sign change for geocentric latitude
            covs(kmt+2,lmt+2) = sngl(vhatc(imt+2,jmt+2))
            imt = imt + mtiev(iev)
            kmt = kmt + 2
         end do
         jmt = jmt + mtiev(jev)
         lmt = lmt + 2
      end do
      nev2 = nev*2
      write (io_out,'(10(e12.6,1x))') ((covs(i,j),i=1,nev2),j=1,nev2)
      
      close (io_out)
      
      return
      end
