      subroutine INT1_interface_dtdu(u,u_max,dTdu,
     &     properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u,u_max,dTdu,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max
      real*8 u(3),u_max(3),dTdu(3,3),properties(*),
     &     E_elastic(3),delta(3),nu(3),nu_max(3)
c     elcon(0:ncmat_,ntmat_,*),

      real*8 nu_init,sigma_max,lambda

      if(properties(4).eq.1) then
         call interface_dtdu_bilinear(u,u_max,dTdu,
     &        properties)
      elseif(properties(4).eq.2) then
         call interface_dtdu_chaboche(u,u_max,dTdu,
     &        properties)
      else
         write(*,*) "ERROR: invalid CZM: ",properties(4)
         stop
      endif

      return
      end
