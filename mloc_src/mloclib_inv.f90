!***********************************************************************************************************************************
subroutine elips (t22a, alphad, al, bl)

! Confidence ellipse

   implicit none

   real, parameter :: rpd=1.7453293e-2

   real :: t22a(2,2), t22b(2,2), alpha, alphad, al, bl, arg1, arg2, c2a, s2a, sca

   call smatinv2 (t22a, t22b, 2, 2, 2, 2)
   arg1 = 2.*t22b(1,2)
   arg2 = t22b(1,1) - t22b(2,2)
   alpha = 0.5*atan2(arg1,arg2)
   alphad = alpha/rpd
   c2a = cos(alpha)*cos(alpha)
   s2a = sin(alpha)*sin(alpha)
   sca = sin(alpha)*cos(alpha)
   al = t22b(1,1)*c2a + 2.*t22b(1,2)*sca + t22b(2,2)*s2a
   bl = t22b(1,1)*s2a - 2.*t22b(1,2)*sca + t22b(2,2)*c2a 

   return
   
end subroutine elips

      
!***********************************************************************************************************************************
subroutine delips (t22a, alphad, al, bl)

! Confidence ellipse, using double precision

   implicit none

   real, parameter :: rpd=1.7453293e-2

   double precision :: t22a(2,2), t22b(2,2), c2a, s2a, sca
   real :: arg1, arg2, alpha, alphad, al, bl

   call dmatinv2 (t22a, t22b, 2, 2, 2, 2)
   arg1 = 2.*sngl(t22b(1,2))
   arg2 = sngl(t22b(1,1) - t22b(2,2))
   alpha = 0.5*atan2(arg1,arg2)
   alphad = alpha/rpd
   c2a = dble(cos(alpha)*cos(alpha))
   s2a = dble(sin(alpha)*sin(alpha))
   sca = dble(sin(alpha)*cos(alpha))
   al = sngl(t22b(1,1)*c2a + 2.0d0*t22b(1,2)*sca + t22b(2,2)*s2a)
   bl = sngl(t22b(1,1)*s2a - 2.0d0*t22b(1,2)*sca + t22b(2,2)*c2a) 

   return
   
end subroutine delips


!***********************************************************************************************************************************
subroutine postq (ain, na, ma, n, m, bout)

! post-multiplies qbhat by matrix ain, without forming qbhat.
! qbhat * ain = bout
! dimension of ain and bout in calling program is (na,ma).
! dimensions used in multiplication are (n,m)
! Uses double precision.

   implicit none

   include 'mloc.inc'

   integer :: na, ma, n, m, j, i, k, ii, jbi
   double precision :: ain(na,ma), bout(na,ma), tempij, qbhatik, boutij 

   do j = 1,m
      do i = 1,n
         jbi = jb(i)
         tempij = dble(wq2(jbi)*sighatj(i))
         boutij = 0.0d0
         do k = 1,ntqi(jbi)
            ii = ntq(k,jbi)
            qbhatik = -1.0d0/(tempij*dble(sighatj(ii)))
            if (ii .eq. i) qbhatik = qbhatik + 1.0d0
            boutij = boutij + (qbhatik*ain(ii,j))
         end do
         bout(i,j) = boutij
      end do
   end do

   return
   
end subroutine postq


!***********************************************************************************************************************************
subroutine postqv (ain, na, n, bout)

! post-multiplies qbhat by vector ain, without forming qbhat.
! qbhat * ain = bout
! dimension of ain and bout in calling program is na.
! vector length used in multiplication is n.
! Uses double precision.

   implicit none

   include 'mloc.inc'

   integer :: na, n, i, jbi, k, ii
   double precision :: ain(na), bout(na), tempi, bouti, qbhatik

   do i = 1,n
      jbi = jb(i)
      tempi = dble(wq2(jbi)*sighatj(i))
      bouti = 0.0d0
      do k = 1,ntqi(jbi)
         ii = ntq(k,jbi)
         qbhatik = -1.0d0/(tempi*dble(sighatj(ii)))
         if (ii .eq. i) qbhatik = qbhatik + 1.0d0
         bouti = bouti + (qbhatik*ain(ii))
      end do
      bout(i) = bouti
   end do

   return
   
end subroutine postqv


!***********************************************************************************************************************************
subroutine postpv (ain, na, n, bout)

! post-multiplies pbhat by vector ain, without forming qbhat.
! pbhat * ain = bout
! dimension of ain and bout in calling program is (na).
! vector length used in multiplication is n.
! Uses double precision.

   implicit none

   include 'mloc.inc'

   integer :: na, n, i, jbi, k, ii
   double precision :: ain(na), bout(na), tempi, bouti, pbhatik

   do i = 1,n
      jbi = jb(i)
      tempi = wq2(jbi)*dble(sighatj(i))
      bouti = 0.0d0
      do k = 1,ntqi(jbi)
         ii = ntq(k,jbi)
         pbhatik = 1.0d0/(tempi*dble(sighatj(ii)))
         bouti = bouti + (pbhatik*ain(ii))
      end do
      bout(i) = bouti
   end do

   return
   
end subroutine postpv

!***********************************************************************************************************************************
subroutine indexx (n, arrin, indx)

! indexes array 'arrin' of length 'n', i.e., outputs the array 'indx'
! such that arrin(indx(j)) is in ascending order for j=1,2,...,n.  the
! input quantities n and arrin are not changed.

   implicit none

   real :: arrin(n), q
   integer :: indx(n), n, j, l, ir, indxt, i

   do j = 1,n
      indx(j) = j
   end do
   if (n .eq. 1) return
   l = n/2+1
   ir = n
10 continue
   if (l .gt. 1) then
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
   else
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir-1
      if (ir .eq. 1) then
         indx(1) = indxt
         return
      end if
   end if
   i = l
   j = l + l
20 if (j .le. ir) then
      if (j .lt. ir) then
         if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
      end if
      if (q .lt. arrin(indx(j))) then
         indx(i) = indx(j)
         i = j
         j = j + j
      else
         j = ir + 1
      end if
      go to 20 
   end if 
   indx(i) = indxt
   go to 10

end subroutine indexx
         
!***********************************************************************************************************************************
subroutine dindexx (n, arrin, indx)

! indexes double precision array 'arrin' of length 'n', i.e., outputs 
! the array 'indx'
! such that arrin(indx(j)) is in ascending order for j=1,2,...,n.  the
! input quantities n and arrin are not changed.

   implicit none

   integer :: indx(n), n, i, j, l, ir, indxt 
   double precision :: arrin(n), q

   do j = 1,n
      indx(j) = j
   end do
   if (n .eq. 1) return
   l = n/2+1
   ir = n
10 continue
   if (l .gt. 1) then
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
   else
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
      if (ir .eq. 1) then
         indx(1) = indxt
         return
      end if
   end if
   i = l
   j = l + l
20 if (j .le. ir) then
      if (j .lt. ir) then
         if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
      end if
      if (q .lt. arrin(indx(j))) then
         indx(i) = indx(j)
         i = j
         j = j + j
      else
         j = ir + 1
      end if
      go to 20 
   end if 
   indx(i) = indxt
   go to 10

end subroutine dindexx
         
!***********************************************************************************************************************************
subroutine rank (n, indx, irank)

! given 'indx' of length 'n' as output from subroutine 'indexx', this
! subroutine returns an array 'irank', the corresponding table of ranks.

   implicit none

   integer :: indx(n), irank(n), n, j

   do j = 1,n
      irank(indx(j)) = j
   end do
   
   return
   
end subroutine rank

!***********************************************************************************************************************************
subroutine dot1 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, single precision.
! a * b = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, l, m, n, i, k, j
   real :: a(ma,na), b(mb,nb), c(mc,nc), bjk

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.
      end do
   end do 

   do k = 1,n
      do j =1,m
         bjk = b(j,k)
         do i = 1,l
            c(i,k) = c(i,k) + a(i,j)*bjk
         end do
      end do
   end do

   return
   
end subroutine dot1


!***********************************************************************************************************************************
subroutine dot2 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, single precision.
! a**tr * b = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, l, m, n, i, j, k
   real :: a(ma,na), b(mb,nb), c(mc,nc), cik

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.
      end do
   end do 

   do k = 1,n
      do i = 1,l
         cik = c(i,k)
         do j = 1,m
            cik = cik + a(j,i)*b(j,k)
         end do
         c(i,k) = cik
      end do
   end do

   return
   
end subroutine dot2


!***********************************************************************************************************************************
subroutine dot3 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, single precision.
! a * b**tr = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, l, m, n, i, j, k
   real :: a(ma,na), b(mb,nb), c(mc,nc), bkj

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.
      end do
   end do 

   do j = 1,m
      do k = 1,n
         bkj = b(k,j)
         do i = 1,l
            c(i,k) = c(i,k) + a(i,j)*bkj
         end do
      end do
   end do

   return
   
end subroutine dot3
      
      
!***********************************************************************************************************************************
subroutine ddot1 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, double precision.
! a * b = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, i, j, k, l, m, n
   double precision :: a(ma,na), b(mb,nb), c(mc,nc), bjk 

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.0d0
      end do
   end do 

   do k = 1,n
      do j = 1,m
         bjk = b(j,k)
         do i = 1,l
            c(i,k) = c(i,k) + a(i,j)*bjk
         end do
      end do
   end do

   return
   
end subroutine ddot1


!***********************************************************************************************************************************
subroutine ddot2 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, double precision.
! a**tr * b = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, i, j, k, l, m, n
   double precision :: a(ma,na), b(mb,nb), c(mc,nc), cik 

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.0d0
      end do
   end do 

   do k = 1,n
      do i = 1,l
         cik = c(i,k)
         do j = 1,m
            cik = cik + a(j,i)*b(j,k)
         end do
         c(i,k) = cik
      end do
   end do

   return
   
end subroutine ddot2

!***********************************************************************************************************************************
subroutine ddot3 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)

! dot product of two matrices, double precision.
! a * b**tr = c
! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).

! no subroutines called.

   implicit none

   integer :: ma, na, mb, nb, mc, nc, i, j, k, l, m, n
   double precision :: a(ma,na), b(mb,nb), c(mc,nc), bkj 

   do k = 1,n
      do i = 1,l
         c(i,k) = 0.0d0
      end do
   end do 

   do j = 1,m
      do k = 1,n
         bkj = b(k,j)
         do i = 1,l
            c(i,k) = c(i,k) + a(i,j)*bkj
         end do
      end do
   end do

   return
   
end subroutine ddot3


!***********************************************************************************************************************************
subroutine smatinv2 (a, b, idim, jdim, ia, ja)

! inverse of square matrix a returned in b.
! single precision.  a is unchanged.

   implicit none

   integer :: ik(500), jk(500), idim, jdim, ia, ja, i, j, k, l
   real :: a(idim,jdim), b(idim,jdim), amax, save_me, eps
   character(len=132) :: msg
   
   eps = 2.*epsilon(amax)

   if (ia .ne. ja) then
     call warnings ('smatinv2: non-square matrix dimensions')
     return
   end if

   do j = 1,ja
      do i = 1,ia
         b(i,j) = a(i,j)
      end do
   end do

   do k = 1,ia

      ! find largest element b(i,j) in rest of matrix.
      amax = 0.0
      do j = k,ia
         do i = k,ia
            if (abs(amax) .lt. abs(b(i,j))) then
               amax = b(i,j)
               ik(k) = i
               jk(k) = j
            end if
         end do
      end do

      ! interchange rows and columns of matrix to put amax in b(k,k)
      ! (ie. do partial pivoting to improve accuracy).
      if (abs(amax) .lt. eps) then
        write (msg,'(a,10e12.3)') 'smatinv2: matrix is singular; ', a
        call oops (trim(msg))
      end if

      i = ik(k)
      if (i .ne. k) then
         do j = 1,ia
            save_me = b(k,j)
            b(k,j) = b(i,j)
            b(i,j) = -save_me
         end do
      end if

      j = jk(k)
      if (j .ne. k) then
         do i = 1,ia
            save_me = b(i,k)
            b(i,k) = b(i,j)
            b(i,j) = -save_me
         end do
      end if

      ! accumulate elements of the inverse matrix.
      save_me = b(k,k)
      do i = 1,ia
         b(i,k) = -b(i,k)/amax
      end do
      b(k,k) = save_me

      do j = 1,ia
         if (k .ne. j) then
            save_me = b(k,j)
            do i = 1,ia
               b(i,j) = b(i,j) + b(i,k)*save_me
            end do
            b(k,j) = save_me
         end if
      end do
      do j = 1,ia
         b(k,j) = b(k,j)/amax
      end do
      b(k,k) = 1.0/amax
   end do 

   ! restore ordering of the matrix.
   do l = 1,ia
      k = ia - l + 1
      j = ik(k)
      if (j .gt. k) then
         do i = 1,ia
            save_me = b(i,k)
            b(i,k) = -b(i,j)
            b(i,j) = save_me
         end do
      end if
      i = jk(k)
      if (i .gt. k) then
         do j = 1,ia
            save_me = b(k,j)
            b(k,j) = -b(i,j)
            b(i,j) = save_me
         end do
      end if
   end do

   return
   
end subroutine smatinv2
      
      
!***********************************************************************************************************************************
subroutine dmatinv2 (a, b, idim, jdim, ia, ja)

! inverse of square matrix a returned in b.
! double precision.  a is unchanged.

   implicit none

   integer :: ik(500), jk(500), idim, jdim, ia, ja, i, j, k, l
   double precision :: a(idim,jdim), b(idim,jdim), amax, save_this, eps
   
   eps = 2.*epsilon(amax)

   if (ia .ne. ja) then
     call warnings ('dmatinv2: non-square matrix dimensions')
     return
   end if

   do j = 1,ja
      do i = 1,ia
         b(i,j) = a(i,j)
      end do
   end do

   do k = 1,ia

      ! find largest element b(i,j) in rest of matrix.
      amax = 0.0d0
      do j = k,ia
         do i = k,ia
            if (abs(amax) .lt. abs(b(i,j))) then
               amax = b(i,j)
               ik(k) = i
               jk(k) = j
            end if
         end do
      end do

      ! interchange rows and columns of matrix to put amax in b(k,k)
      ! (ie. do partial pivoting to improve accuracy).
      if (abs(amax) .lt. eps) then
        call warnings ('dmatinv2: matrix is singular')
        return
      end if

      i = ik(k)
      if (i .ne. k) then
         do j = 1,ia
            save_this = b(k,j)
            b(k,j) = b(i,j)
            b(i,j) = -save_this
         end do
      end if

      j = jk(k)
      if (j .ne. k) then
         do i = 1,ia
            save_this = b(i,k)
            b(i,k) = b(i,j)
            b(i,j) = -save_this
         end do
      end if

      ! accumulate elements of the inverse matrix.
      save_this = b(k,k)
      do i = 1,ia
         b(i,k) = -b(i,k)/amax
      end do
      b(k,k) = save_this

      do j = 1,ia
         if (k .ne. j) then
            save_this = b(k,j)
            do i = 1,ia
               b(i,j) = b(i,j) + b(i,k)*save_this
            end do
            b(k,j) = save_this
         end if
      end do
      do j = 1,ia
         b(k,j) = b(k,j)/amax
      end do
      b(k,k) = 1.0/amax
   end do 

   ! restore ordering of the matrix.
   do l = 1,ia
      k = ia - l + 1
      j = ik(k)
      if (j .gt. k) then
         do i = 1,ia
            save_this = b(i,k)
            b(i,k) = -b(i,j)
            b(i,j) = save_this
         end do
      end if
      i = jk(k)
      if (i .gt. k) then
         do j = 1,ia
            save_this = b(k,j)
            b(k,j) = -b(i,j)
            b(i,j) = save_this
         end do
      end if
   end do

   return
   
end subroutine dmatinv2
