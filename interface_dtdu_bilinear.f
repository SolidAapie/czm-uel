      subroutine interface_dtdu_bilinear(u,u_max,dTdu,
     &     properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u,u_max,dTdu,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max
      real*8 u(3),u_max(3),dTdu(3,3),properties(*),
     &     E(3),delta(3),nu(3),nu_max(3),G_c,sigma_max
c     elcon(0:ncmat_,ntmat_,*),

      real*8 lambda_init,lambda,lambda_max,lambda_F,F,dlambda_F_dnu(3),
     &     dFdlambda_F

      integer k,m,n


c     Geubelle (normal only)

      G_c = properties(1)
      sigma_max = properties(2)
      lambda_init = .1

      do m=1,3
         delta(m) = 2.0*G_c/sigma_max
         E(m) = sigma_max/lambda_init
         nu(m) = u(m)/delta(m)
         nu_max(m) = u_max(m)/delta(m)
      enddo
                                ! Determine lambda_F,dlambda_F_dnu
      call compute_lambda(nu,lambda)
      call compute_lambda(nu_max,lambda_max)
      if(lambda.ge.lambda_max) then
         lambda_F = lambda
         call compute_dlambda_F_dnu(nu,lambda,dlambda_F_dnu)
      else
         lambda_F = lambda_max
         do n=1,3
            dlambda_F_dnu(n) = 0
         enddo
      endif
                                !Determine F,dFdlambda_F
      call compute_F_lambda_bilinear(lambda_F,lambda_init,F)
      call compute_dFdlambda_F_bilinear(lambda_F,lambda_init,
     &     dFdlambda_F)
                                ! dTdu
      do m=1,3
         do n=1,3
            dTdu(m,n) = E(m)/delta(n)*nu(m)*dlambda_F_dnu(n)*dFdlambda_F
         enddo
      enddo
                                ! Add diagonal terms
      do m=1,3
         dTdu(m,m) = dTdu(m,m) + E(m)/delta(m)*F
      enddo
                                ! Check for interpenetration
                                ! this does not make the jacobian unsymmetric because d lambda/d nu(3) = 0 when nu(3) < 0
                                ! however, it is unsymmetric if E(m) != E(n)
      if(nu(3).lt.0) then
         dTdu(3,1) = 0
         dTdu(3,2) = 0
         dTdu(3,3) = E(3)/delta(3)
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

c      do k=1,3
c         do m=1,3
c            dTdu(k,m) = 0
c         enddo
c      enddo

c      if(nu(3).lt.0) then
c         dTdu(3,3) = (1/nu_init - 1)*sigma_max/delta(3)
c      elseif(nu(3).lt.lambda) then
c         dTdu(3,3) = (1/lambda - 1)*sigma_max/delta(3)
c      elseif(nu(3).lt.1) then
c         dTdu(3,3) = -sigma_max/delta(3)
c      else
c         dTdu(3,3) = 0
c      endif
c      if(nu_max(3).gt.1) then
c         dTdu(3,3) = 0
c      endif



      return
      end
      
