      subroutine interface_dtdu_chaboche(
     &     u,u_max,dTdu,properties)

      implicit none

c     add area/4 to sum over integration points!
c      real*8 u,u_max,dTdu,elcon(0:ncmat_,ntmat_,*),
c     &     E_elastic,delta,nu,nu_max,Gc
      real*8 u(3),u_max(3),dTdu(3,3),properties(*),
     &     E(3),delta(3),nu(3),nu_max(3),G_c,sigma_max
c     elcon(0:ncmat_,ntmat_,*),

      real*8 nu_init,lambda,lambda_max,lambda_F,F,
     &     dlambda_F_dnu(3),dFdlambda_F

      integer k,m,n

      
c     Chaboche (normal only)
c$$$      Gc = .5*delta*sigma_max*.9
c$$$      E_elastic = sigma_max*27/4
c$$$      delta = Gc*12/E_elastic

c$$$      E_elastic(3) = elcon(1,1,nmat)
c$$$
c$$$      nu(3) = u(3)/delta(3)
c$$$      nu_max(3) = u_max(3)/delta(3)
c$$$
c$$$      do k=1,3
c$$$         do m=1,3
c$$$            dTdu(k,m) = 0
c$$$         enddo
c$$$      enddo
c$$$
c$$$      if(nu(3).lt.0) then
c$$$         dTdu(3,3) = E_elastic(3)/delta(3)
c$$$      elseif(nu(3).lt.1) then
c$$$         dTdu(3,3) = E_elastic(3)/delta(3)*(1-4*nu(3)+3*nu(3)**2)
c$$$      else
c$$$         dTdu(3,3) = 0
c$$$      endif
c$$$      if(nu_max(3).gt.1.and.nu(3).gt.0) then
c$$$         dTdu(3,3) = 0
c$$$      endif


         G_c = properties(1)
         sigma_max = properties(2)

                                ! Determine delta,E,nu,nu_max
      do m=1,3
c         delta(m) = elcon(2,1,nmat)
c         E(m) = elcon(1,1,nmat)
c         delta(m) = properties(2)
c         E(m) = properties(1)
         delta(m) = 48.0/27.0*G_c/sigma_max
         E(m) = 27.0/4.0*sigma_max
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
c$$$$         call compute_dlambda_F_dnu(nu_max,lambda_max,dlambda_F_dnu)
         do n=1,3
            dlambda_F_dnu(n) = 0
         enddo
      endif
                                !Determine F,dFdlambda_F
      call compute_F_lambda_chaboche(lambda_F,F)
      call compute_dFdlambda_F_chaboche(lambda_F,dFdlambda_F)
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
      if(nu(3).lt.0) then
         dTdu(3,1) = 0
         dTdu(3,2) = 0
         dTdu(3,3) = E(3)/delta(3)
      endif


      return
      end
      
