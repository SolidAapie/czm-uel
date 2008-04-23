      subroutine INT1_interface_t(u,u_max,T,properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u(3,9),u_max(3,9),T(3,9),E_elastic,delta,nu,nu_max
c      real*8 u,u_max,T,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max
      real*8 u(3),u_max(3),T(3),properties(*),
     &     E_elastic(3),delta(3),nu(3),nu_max(3)
c     elcon(0:ncmat_,ntmat_,*),      

      real*8 nu_init,sigma_max,lambda

      if(properties(4).eq.1) then
         call interface_t_bilinear(u,u_max,T,properties)
      elseif(properties(4).eq.2) then
         call interface_t_chaboche(u,u_max,T,properties)
      endif

      return
      end
