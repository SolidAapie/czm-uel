      subroutine compute_dlambda_F_dnu(nu,lambda,dlambda_F_dnu)

      implicit none
      real*8 nu(3),lambda,dlambda_F_dnu(3)
      integer n

c     get rid of division by zero problem, if necessary
      if(lambda.eq.0) then
         do n=1,3
            dlambda_F_dnu(n) = 0
         enddo
      else
         dlambda_F_dnu(1) = nu(1)/lambda
         dlambda_F_dnu(2) = nu(2)/lambda
         if(nu(3).gt.0) then
            dlambda_F_dnu(3) = nu(3)/lambda
         else
            dlambda_F_dnu(3) = 0
         endif
      endif

      
      return
      end
