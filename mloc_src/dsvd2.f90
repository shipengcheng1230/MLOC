!***********************************************************************
      subroutine dsvd (ain, m, n, mp, np, u, v, q)

! singular value decomposition.

! further attempts made to optimize 7/8/88 eab

! double precision

! for algol program see wilkinson+reinsch: handbook for automatic
! computation vol 2 - linear algebra, pp. 140-144. translated from
! algol by r.l.parker.

! cleaned up by e. bergman - 

! the matrix ain(m,n) is decomposed. singular values in q, pre-matrix in u,
! post-matrix in v. 

! program altered by p. silver to handle unpacked arrays.
! mp and np are dimensions in main routine. m and n are actual
! dimensions to be used in the subroutine.

! no subroutines called.

! eps = machine accuracy: smallest double precision value which,
!       when added to 1.0, yields a result different from 1.0.
! tol = smallest double precision floating point number considered to be
!       greater than 0.

! the include file is only used to bring in the dimension of the
! error vector 'e' from the main program.

      implicit none
      
      include 'mloc.inc'

      integer m, n, mp, np, i, j, k, l, iback, kback, lback, l1, lplus
      integer i1d, lm1
      double precision ain(mp,np), u(mp,np), v(np,np), q(np), e(mtmax)
      double precision eps, tol, f, g, h, x, y, z, scale, hinv, ginv, c 

      eps=1.0d-15
      tol=1.0d-100

      do j=1,n
         do i=1,m
            u(i,j)=ain(i,j)
         end do
      end do

! householder reduction to bi-diagonal form

      g=0.0d0
      x=0.0d0
      do i=1,n

         e(i)=g
         scale=0.0d0
         l=i+1
         do j=i,m
            scale=u(j,i)*u(j,i)+scale
         end do

         if (scale .ge. tol) then
            f=u(i,i)
            g=-dsign(dsqrt(scale),f)
            h=f*g-scale
            u(i,i)=f-g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  scale=0.0d0
                  do k=i,m
                     scale=u(k,i)*u(k,j)+scale
                  end do
                  f=scale*hinv
                  do k=i,m
                     u(k,j)=u(k,j)+f*u(k,i)
                  end do 
               end do
            end if
         else
            g=0.0d0
         end if

         q(i)=g
         scale=0.0d0
         if (l .le. n) then
            do j=l,n
               scale=u(i,j)*u(i,j)+scale
            end do
         end if

         if (scale .ge. tol) then
            i1d = i+1
            f=u(i,i1d)
            g=-dsign(dsqrt(scale),f)
            h=f*g-scale
            u(i,i+1)=f-g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  e(j)=u(i,j)*hinv
               end do
            end if
            if (l .le. m) then
               do j=l,m
                  scale=0.0d0
                  if (l .le. n) then
                     do k=l,n
                        scale=u(j,k)*u(i,k)+scale
                     end do
                     do k=l,n
                        u(j,k)=u(j,k)+scale*e(k)
                     end do
                  end if
               end do
            end if
         else
            g=0.0d0
         end if

         y=dabs(q(i))+dabs(e(i))
         if (y .gt. x) x=y

      end do

! accumulation of right-hand transforms (v)

      do iback=1,n
         i=n+1-iback
         if (abs(g) .gt. eps) then
            h=u(i,i+1)*g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  v(j,i)=u(i,j)*hinv
               end do
               do j=l,n
                  scale=0.0d0
                  do k=l,n
                     scale=u(i,k)*v(k,j)+scale
                  end do
                  do k=l,n
                     v(k,j)=v(k,j)+scale*v(k,i)
                  end do
               end do
            end if
         end if
         if (l .le. n) then
            do j=l,n
               v(j,i)=0.0d0
               v(i,j)=0.0d0
            end do
         end if
         v(i,i)=1.0d0
         g=e(i)
         l=i
      end do

! accumulation of left-hand transforms

      do iback=1,n
         i=n+1-iback
         l=i+1
         g=q(i)
         if (l .le. n) then
            do j=l,n
               u(i,j)=0.0d0
            end do
         end if
         if (abs(g) .gt. eps) then
            h=u(i,i)*g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  scale=0.0d0
                  do k=l,m
                     scale=u(k,i)*u(k,j)+scale
                  end do
                  f=scale*hinv
                  do k=i,m
                     u(k,j)=u(k,j)+f*u(k,i)
                  end do
               end do
            end if
            ginv=1.0d0/g
            do j=i,m
               u(j,i)=u(j,i)*ginv
            end do
         else
            do j=i,m
               u(j,i)=0.0d0
            end do
         end if
         u(i,i)=u(i,i)+1.0d0
      end do

! diagonalization of bi-diagonal form

      eps=eps*x
      do kback=1,n
         k=n+1-kback

! test f-splitting

 5000    continue
         do lback=1,k
            l=k+1-lback
            if (dabs(e(l)) .le. eps) go to 6500
            lm1 = l-1
            if (dabs(q(lm1)) .le. eps) go to 6000
         end do

! cancellation of e(l) if l .gt. 1

 6000    c=0.0d0
         scale=1.0d0
         l1=l-1
         do i=l,k
            f=scale*e(i)
            e(i)=c*e(i)
            if (dabs(f) .le. eps) go to 6500
            g=q(i)
            q(i)=dsqrt(f*f+g*g)
            h=q(i)
            c=g/h
            scale=-f/h
            do j=1,m
               y=u(j,l1)
               z=u(j,i)
               u(j,l1)=y*c+z*scale
               u(j,i)=-y*scale+z*c
            end do
         end do

! test f-convergence

 6500    z=q(k)
         if (l .ne. k) then
            x=q(l)
            y=q(k-1)
            g=e(k-1)
            h=e(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
            g=sqrt(f*f+1.0d0)
            f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x

! next q-r transformation

            c=1.0d0
            scale=1.0d0
            lplus=l + 1
            do i=lplus,k
               g=e(i)
               y=q(i)
               h=scale*g
               g=c*g
               z=dsqrt(f*f+h*h)
               e(i-1)=z
               c=f/z
               scale=h/z
               f=x*c+g*scale
               g=-x*scale+g*c
               h=y*scale
               y=y*c
               do j=1,n
                  x=v(j,i-1)
                  z=v(j,i)
                  v(j,i-1)=x*c+z*scale
                  v(j,i)=-x*scale+z*c
               end do
               z=dsqrt(f*f+h*h)
               q(i-1)=z
               c=f/z
               scale=h/z
               f=c*g+scale*y
               x=-scale*g+c*y
               do j=1,m
                  y=u(j,i-1)
                  z=u(j,i)
                  u(j,i-1)=y*c+z*scale
                  u(j,i)=-y*scale+z*c
               end do
            end do
            e(l)=0.0d0
            e(k)=f
            q(k)=x
            go to 5000
         end if

! convergence

         if (z .lt. 0.0d0) then
            q(k)=-z
            do j=1,n
               v(j,k)=-v(j,k)
            end do
         end if

      end do

      return
      end

