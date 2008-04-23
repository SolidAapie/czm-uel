      subroutine compute_lambda(nu,lambda)

      implicit none
      real*8 nu(3),lambda
      
      if(nu(3).gt.0) then
         lambda = sqrt(nu(1)**2 + nu(2)**2 + nu(3)**2)
c         lambda = nu(3)
      else
         lambda = sqrt(nu(1)**2 + nu(2)**2)
c         lambda = 0
      endif

      return
      end
