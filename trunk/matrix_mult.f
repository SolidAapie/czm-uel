      subroutine matrix_mult(a,b,axb,n_)

      implicit none
      integer n_
      real*8 a(n_,n_),b(n_,n_),axb(n_,n_) 
      integer i,j,k

      do i=1,n_
         do j=1,n_
            axb(i,j) = 0
            do k=1,n_
               axb(i,j) = axb(i,j) + a(i,k)*b(k,j)
            enddo
         enddo
      enddo

      return
      end
