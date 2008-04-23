      subroutine interface_t_chaboche(
     &     u,u_max,T,properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u(3,9),u_max(3,9),T(3,9),E_elastic,delta,nu,nu_max
c      real*8 u,u_max,T,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max,Gc
      real*8 u(3),u_max(3),T(3),properties(*),
     &     E(3),delta(3),nu(3),nu_max(3),G_c,sigma_max
c     elcon(0:ncmat_,ntmat_,*),

      real*8 lambda,lambda_max,lambda_F,F
      integer m
c$$$      delta(3) = elcon(2,1,nmat)
c$$$      sigma_max = elcon(1,1,nmat)


c     Chaboche (normal only)
c     same energy and sigma_max as bilinear
c$$$      Gc = .5*delta*sigma_max*.9
c$$$      E_elastic = sigma_max*27/4
c$$$      delta = Gc*12/E_elastic

c$$$      E_elastic(3) = elcon(1,1,nmat)
c$$$
c$$$      nu(3) = u(3)/delta(3)
c$$$      nu_max(3) = u_max(3)/delta(3)
c$$$      if(nu(3).lt.0) then
c$$$         T(3) =  E_elastic(3)*nu(3)
c$$$      elseif(nu(3).lt.1) then
c$$$         T(3) =  E_elastic(3)*nu(3)*(1-nu(3))**2
c$$$      else
c$$$         T(3) = 0
c$$$      endif
c$$$      if(nu_max(3).gt.1.and.nu(3).gt.0) then
c$$$         T(3) = 0
c$$$      endif
c$$$
c$$$      T(1) = 0
c$$$      T(2) = 0

         G_c = properties(1)
         sigma_max = properties(2)

                                ! Determine delta,E,nu,nu_max
      do m=1,3
c         delta(m) = properties(2)
c         E(m) = properties(1)
         delta(m) = 48.0/27.0*G_c/sigma_max
         E(m) = 27.0/4.0*sigma_max
         nu(m) = u(m)/delta(m)
         nu_max(m) = u_max(m)/delta(m)
      enddo
c      write(*,*) "CZM: s_m,gc = ",sigma_max,G_c
c      write(*,*) "CZM: E,d = ",E(3),delta(3)

                                ! Determine lambda_F
      call compute_lambda(nu,lambda)
      call compute_lambda(nu_max,lambda_max)
      if(lambda.ge.lambda_max) then
         lambda_F = lambda
      else
         lambda_F = lambda_max
      endif
                                ! Determine F(lambda_F)
      call compute_F_lambda_chaboche(lambda_F,F)
                                ! T_m
      do m=1,3
         T(m) = E(m)*nu(m)*F
      enddo
                                ! Check for interpenetration
      if(nu(3).lt.0) then
         T(3) = E(3)*nu(3)
      endif

      return
      end
