
      subroutine INT1_element_results(xl,vl,calcul_fn,update_state,nope,
     &  properties,fnl,xstate,nstate_)

      implicit none
!
      logical calcul_fn,update_state
!
c      character*6 lakon(*)
!
      integer konl(20)
!
      integer j,k,m,nope,iout,mint_,nstate_
!
      real*8 p(3),q(3)
      real*8 xl(3,20),vl(3,20),
     &  properties(*)
c     &  fn(3,*)
c     &  elcon(0:ncmat_,ntmat_,*),
c     rhcon(0:1,ntmat_,*),

      real*8 xi,et,ze
!
      real*8 xstate(nstate_)

      integer mint, thickness_coord
      real*8 fnl(3,20),
     &     x1,x2,x3,area,
     &     integrand(3,20),u_node(3,20),u_max_node(3,20)

      real*8 node2d3(2,8),node_weight2d3(8),gauss2d2(2,4),weight2d2(4),
     &     gauss2d3(2,9),weight2d3(9),gauss2d4(2,16),weight2d4(16),
     &     gauss2d5(2,25),weight2d5(25),gauss2d6(2,36),weight2d6(36)

c      include "adapt_integrate_vars.f"
      integer minpts,maxpts,iwk      
      real*8 a(2),b_adapt(2),eps,wk(20,21000)
      real*8 result(60),relerr
      integer nfnevl,ifail,nfuncs
      
      logical fixed_order_int
      
      real*8 interface_plane_x,interface_plane_y,interface_z,
     &     vl_element(3,20),xs(3,3),xsi(3,3),
     &     fnl_element(3,20)
      real*8 xi_node(20),et_node(20),ze_node(20)
      integer u_permutation1(20),u_permutation2(20),u_permutation3(20)

      integer unused_node(4,3)

      real*8 nu(3)
      integer coord_permutation(3)
      real*8 u_node_interface(3,20),u_max_node_interface(3,20)

      data unused_node /
     &     9, 11, 13, 15,
     &     17, 18, 19, 20,
     &     10, 12, 14, 16 /

      data u_permutation1 /
     &     2,1,4,3,6,5,8,7,9,12,11,10,13,16,15,14,18,17,20,19 /
      data u_permutation2 /
     &     5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20 /
      data u_permutation3 /
     &     4,3,2,1,8,7,6,5,11,10,9,12,15,14,13,16,20,19,18,17 /

      data xi_node /
     &     -1,1,1,-1, -1,1,1,-1, 0,1,0,-1, 0,1,0,-1, -1,1,1,-1 /
      data et_node /
     &     -1,-1,-1,-1, 1,1,1,1, -1,-1,-1,-1, 1,1,1,1, 0,0,0,0 /
      data ze_node /
     &     1,1,-1,-1, 1,1,-1,-1, 1,0,-1,0, 1,0,-1,0, 1,1,-1,-1 /

      include "integration_points.f"

c     Function to integrate is fnl(3,j) = T(3)*dudq(j)*area/4
c     -Integrator passes (xi,et,ze)
c     -func returns T(3)*dudq(j)*area/4 at (xi,et,ze)
c$$$      do k=1,60
c$$$         write(*,*)"k,(k+2)/3 =",k,k-3*((k+2)/3)+3
c$$$      enddo

c      call INT1_element_info_nope(lakon(i), nope)
c     mint = 9

      if(update_state) then
         iout = 2
      endif
                                ! Get nodal forces
      do j=1,nope
         do k=1,3
c            fnl(k,j) = fn(k,konl(j))
            fnl_element(k,j) = 0
         enddo
      enddo
            ! Transform nodal displacements into element coord system
      do k=1,nope
         xi = xi_node(k)
         et = et_node(k)
         ze = ze_node(k)
         call INT1_coord_transform(xi,et,ze,xl,xs,xsi)
         call INT1_transform_global_2_element(vl,
     &        xsi,vl_element,k)
      enddo
                                ! Determine element orientation

c      if(.false.) then


      x1 = sqrt((xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,2))**2
     &     +(xl(3,1)-xl(3,2))**2)
      x2 = sqrt((xl(1,1)-xl(1,5))**2+(xl(2,1)-xl(2,5))**2
     &     +(xl(3,1)-xl(3,5))**2)
      x3 = sqrt((xl(1,1)-xl(1,4))**2+(xl(2,1)-xl(2,4))**2
     &     +(xl(3,1)-xl(3,4))**2)
c      thickness = min(x1,x2,x3)
c      x2 = x1/2!#################kludge!
      if(x1.lt.x2.and.x1.lt.x3) then
         thickness_coord = 1    !use xi for top and bottom of element
c         xi = (xl(1,1)-xl(1,2))/dabs(xl(1,1)-xl(1,2))
c         xi = dabs(xl(1,1)-xl(1,2))
c         write(*,*) "CZM: xi =",xi
c         if(xi.gt.0) then
c      write(*,*) "CZM: interface_z =1 "
c            interface_z = 1.0
c         else
c      write(*,*) "CZM: interface_z =-1 "
c            interface_z = -1.0
c         endif
c         interface_z = xi
c      write(*,*) "CZM: interface_z = ", interface_z
c         interface_z = (xl(1,1)-xl(1,2))/dabs(xl(1,1)-xl(1,2))
c      write(*,*) "CZM: thickness_coord 1.3"
      elseif(x2.lt.x1.and.x2.lt.x3) then
         thickness_coord = 2    !use et for top and bottom of element
c         interface_z = (xl(3,1)-xl(3,5))/dabs(xl(3,1)-xl(3,5))
      else
         thickness_coord = 3    !use ze for top and bottom of element
c         interface_z = (xl(2,1)-xl(2,4))/dabs(xl(2,1)-xl(2,4))
      endif
c      write(*,*) "CZM: thickness_coord = ",thickness_coord

!     this calculation of the area is linear in the sense that considers the original shape of the element rather than the current shape (xl+u)
!     it also assumes that the original shape of the element is a parallelogram
      do j = 1,3
         p(j) = xl(j,unused_node(2,thickness_coord))-xl(j,unused_node(1,
     &		thickness_coord))
         q(j) = xl(j,unused_node(3,thickness_coord))-xl(j,unused_node(1,
     &	   thickness_coord))
      enddo
      area = sqrt( (p(1)*p(1)+p(2)*p(2)+p(3)*p(3)) 
     &     * (q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) 
     &     - (p(1)*q(1)+p(2)*q(2)+p(3)*q(3))**2 )

c     Find u at the nodes (in the element coords)
      if(thickness_coord.eq.1) then
         do j=1,nope
            do m=1,3
               u_node(m,j) = vl_element(m,j)-
     &              vl_element(m,u_permutation1(j))
            enddo
         enddo               
      elseif(thickness_coord.eq.2) then
         do j=1,nope
            do m=1,3
               u_node(m,j) = vl_element(m,j)-
     &              vl_element(m,u_permutation2(j))
            enddo
         enddo
      else
         do j=1,nope
            do m=1,3
               u_node(m,j) = vl_element(m,j)-
     &              vl_element(m,u_permutation3(j))               
            enddo
         enddo
      endif

c     Find u_max at the nodes
      do j=1,nope
         do m=1,2
            if(dabs(u_node(m,j)).gt.dabs(xstate(3*(j-1)+m+30)).and.
     &           iout.eq.2) then
               u_max_node(m,j) = xstate(3*(j-1)+m+30)
               xstate(3*(j-1)+m+30) = dabs(u_node(m,j))
c ****** and u_max_node = u_node ? no. but, u_max_node = [old] xstate ?
            else
               u_max_node(m,j) = xstate(3*(j-1)+m+30)
            endif
         enddo
         m=3
            if(u_node(m,j).gt.xstate(3*(j-1)+m+30).and.
     &           iout.eq.2) then
               u_max_node(m,j) = xstate(3*(j-1)+m+30)
               xstate(3*(j-1)+m+30) = u_node(m,j)
c ****** and u_max_node = u_node ? no. but, u_max_node = [old] xstate ?
            else
               u_max_node(m,j) = xstate(3*(j-1)+m+30)
            endif
         if(calcul_fn) then
c            write(*,*)"nu_node(",i,",3,",j,") = ",
c     &           u_node(3,j)/elcon(2,1,ielmat(i))
         endif
      enddo

      if(calcul_fn)then
c     Transform u,u_max to interface coordinates.
c     Interface coordinates:
c       u_interface(1) = interface x
c       u_interface(2) = interface y
c       u_interface(3) = interface z
         if(thickness_coord.eq.1) then
            coord_permutation(1) = 2
            coord_permutation(2) = 3
            coord_permutation(3) = 1
         elseif(thickness_coord.eq.2) then
            coord_permutation(1) = 3
            coord_permutation(2) = 1
            coord_permutation(3) = 2
         elseif(thickness_coord.eq.3) then
            coord_permutation(1) = 1
            coord_permutation(2) = 2
            coord_permutation(3) = 3
         endif
         do j=1,nope
            do m=1,3
               u_node_interface(m,j) = u_node(coord_permutation(m),j)
            u_max_node_interface(m,j)=u_max_node(coord_permutation(m),j)
            enddo
         enddo


c     INTEGRATE OVER AREA
         include "integrate_f.f"

                                !Artificial forces for unused nodes
c         k_unused_node = properties(1)*x1*x2*x3
c         m_unused_node = density*x1*x2*x3  !rhcon(1,1,1)*x1*x2*x3
c         do j=1,4
c            do k=1,3
c               fnl_element(k,unused_node(j)) = .5*(vl_element(k,1)+vl_element(k,2))
c            enddo
c         enddo

                    ! Transform nodal forces to global coord system
         do k=1,nope
            xi = xi_node(k)
            et = et_node(k)
            ze = ze_node(k)
            call INT1_coord_transform(xi,et,ze,xl,xs,xsi)
            call INT1_transform_element_2_global(fnl_element,xs,fnl,k)
         enddo

c$$$        do j=1,nope
c$$$           do k=1,2
c$$$              fnl(k,j) = 0
c$$$           enddo
c$$$        enddo

         do j=1,nope
c            write(*,*)"fnl(",i,",3,",j,") = ",fnl(3,j)-fn(3,konl(j))
            do k=1,3
c               fn(k,konl(j)) = fn(k,konl(j)) + fnl(k,j)
c               write(*,*) "CZM: fn uer = ",k,j,fnl(k,j)
c               write(*,*) "CZM: fnl_e = ",k,j,fnl_element(k,j)
c               fnl(k,j) = 42
            enddo
         enddo
      
      endif



c      endif

      if(update_state) then
         do j=1,nope
            do m=1,3
c               nu_max(m) = u_max_node_interface(m,j)/delta(m)
c               nu(m) = u_max_node_interface(m,j)/(48.0/27.0*properties(1)/properties(2))
c               nu(m) = u_max_node(m,j)/(48.0/27.0*properties(1)/properties(2))
c               nu(m) = dabs(nu(m))
               nu(m) = u_max_node(m,j)

c               nu(m) = u_node(m,j)
c               nu(m) = vl_element(m,j)-
c     &              vl_element(m,u_permutation1(j))
            enddo
c            call compute_lambda(nu,xstate(j))
c            xstate(j) = lambda_max
            xstate(j) = sqrt(nu(1)**2 + nu(2)**2 + nu(3)**2)            
         enddo
      endif

      return
      end


