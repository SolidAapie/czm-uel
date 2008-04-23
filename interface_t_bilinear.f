      subroutine interface_t_bilinear(
     &     u,u_max,T,properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u(3,9),u_max(3,9),T(3,9),E_elastic,delta,nu,nu_max
c      real*8 u,u_max,T,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max
      real*8 u(3),u_max(3),T(3),properties(*),
     &     E(3),delta(3),nu(3),nu_max(3),G_c,sigma_max
c     elcon(0:ncmat_,ntmat_,*),
      
      real*8 lambda_init,lambda,lambda_max,lambda_F,F
      integer m



      G_c = properties(1)
      sigma_max = properties(2) 
      lambda_init = .1

                                ! Determine delta,E,nu,nu_max
      do m=1,3
         delta(m) = 2.0*G_c/sigma_max
         E(m) = sigma_max/lambda_init
         nu(m) = u(m)/delta(m)
         nu_max(m) = u_max(m)/delta(m)
      enddo
                                 ! Determine lambda_F
      call compute_lambda(nu,lambda)
      call compute_lambda(nu_max,lambda_max)
      if(lambda.ge.lambda_max) then
         lambda_F = lambda
      else
         lambda_F = lambda_max
      endif
                                ! Determine F(lambda_F)
      call compute_F_lambda_bilinear(lambda_F,lambda_init,F)
                                ! T_m
      do m=1,3
         T(m) = E(m)*nu(m)*F
      enddo
                                ! Check for interpenetration
      if(nu(3).lt.0) then
         T(3) = E(3)*nu(3)
      endif

c      delta(3) = 48/27*G_c/sigma_max
c      nu(3) = u(3)/delta(3)
c      nu_max(3) = u_max(3)/delta(3)
c      nu_init = .1

c      if(nu_max(3).lt.nu_init) then
c         lambda = nu_init
c      else
c         lambda = nu_max(3)
c      endif

c      if(nu(3).lt.0) then
c         T(3) = (1/nu_init - 1)*nu(3)*sigma_max
c      elseif(nu(3).lt.lambda) then
c         T(3) = (1/lambda - 1)*nu(3)*sigma_max
c      elseif(nu(3).lt.1) then
c         T(3) = sigma_max*(1-nu(3))
c      else
c         T(3) = 0
c      endif
c      if(nu_max(3).gt.1) then
c         T(3) = 0
c      endif


      return
      end
