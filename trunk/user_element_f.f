
c      subroutine user_element_f(xi,et,ze,thickness,thickness_coord,nope,
c     &     xl,u_node,u_max_node,elcon,nmat,ncmat_,ntmat_,f)
      subroutine user_element_f(interface_plane_x,interface_plane_y,
     &   interface_z,area,thickness_coord,nope,xl,u_node_interface,
     &     u_max_node_interface,properties,f)
      implicit none
!
      integer thickness_coord,nope
      real*8 interface_plane_x,interface_plane_y,
     &     interface_z,area,properties(*),
     &     xl(3,20),u_node(3,20),u_max_node(3,20),f(3,20)
c     elcon(0:ncmat_,ntmat_,*),

      integer j,k,m,coord_permutation(3)
      real*8 xi,et,ze,xsj,
     &     shp1(4,20),shp2(4,20),u(3),u_max(3),dudq(20),
     &     T(3),u_interface(3),u_max_interface(3),T_interface(3),
     &     u_node_interface(3,20),u_max_node_interface(3,20)
      if(thickness_coord.eq.1) then
c         xi = interface_z
         xi = 1
         et = interface_plane_x
         ze = interface_plane_y
      elseif(thickness_coord.eq.2) then
         xi = interface_plane_y
c         et = interface_z
         et = 1
         ze = interface_plane_x
      elseif(thickness_coord.eq.3) then
         xi = interface_plane_x
         et = interface_plane_y
c         ze = interface_z
         ze = 1
      endif
      call shape20h_centered(xi,et,ze,xl,xsj,shp1)
c     Reverse thickness coordinate
      if(thickness_coord.eq.1) then
         xi = -xi
      elseif(thickness_coord.eq.2) then
         et = -et
      elseif(thickness_coord.eq.3) then
         ze = -ze
      else
         write(*,*) "ERROR in user_element_f(): bad thickness_coord"
      endif
      call shape20h_centered(xi,et,ze,xl,xsj,shp2)
c     Reverse thickness coordinate
      if(thickness_coord.eq.1) then
         xi = -xi
      elseif(thickness_coord.eq.2) then
         et = -et
      elseif(thickness_coord.eq.3) then
         ze = -ze
      endif

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
      
c     Calculate u, u_max, and dudq
      do m=1,3
         u_interface(m) = 0
         u_max_interface(m) = 0
      enddo
      do j=1,nope
         do m=1,3
c            u(m) = u(m) + (shp1(4,j)-shp2(4,j))*u_node(m,j)
c            u_max(m) = u_max(m) + (shp1(4,j)-shp2(4,j))*u_max_node(m,j)
            u_interface(m) = u_interface(m) + 
     &           shp1(4,j)*u_node_interface(m,j)
            u_max_interface(m) = u_max_interface(m) + 
     &           shp1(4,j)*u_max_node_interface(m,j)
c$$$            u(m) = u(m) + shp2(4,j)*u_node(m,j)
c$$$            u_max(m) = u_max(m) + shp2(4,j)*u_max_node(m,j)
         enddo
         dudq(j) = shp1(4,j)-shp2(4,j)
c$$$         dudq(j) = shp2(4,j)-shp1(4,j)
      enddo

c     Get T(u)
c      call INT1_interface_t(u_interface,u_max_interface,T_interface,
      call INT1_interface_t(u_interface,u_max_interface,T_interface,
     &     properties)
      
c     Transform T from interface to element coordinates.
      do j=1,3
         T(coord_permutation(j)) = T_interface(j)
      enddo

      do j=1,nope
c      j = i
                                ! xsj is the jacobian deterimant for the
                                !  volume integral; 2/thickness makes it 
                                !  a surface integral perpendicular to
                                !  the thickness.
         do k=1,3
c     T(3) is not always right. It should be T(thickness_coord)
c     There should also be a T(x_coord) and T(y_coord) (in-plane coords)
c            f(k,j) = T(k)*dudq(j)*xsj/thickness*2d0
            f(k,j) = T(k)*dudq(j)*area/4.0 ! need to divide by 4 bc the area in element coords is 4
c$$$            if (xsj/thickness*2.0.lt.1d-6) then
c$$$               f(k,j) = T(k)*dudq(j)*4.6720125d-7
c$$$            else
c$$$               f(k,j) = T(k)*dudq(j)*5.80644d-5
c$$$            endif
c            if (dabs(T(k)).gt.1e-20.and.dabs(dudq(j)).gt.1e-20) then
c               !write(*,*)"CZM user_element_f: area = ",j,k,f(k,j)/T(k)/dudq(j)
c            endif
         enddo
c         f(j) = T(3)*dudq(j)*.25
      enddo
c      write(*,*) "CZM area = ",xsj/thickness*8d0 ! the sum of the weights for the gauss points = 4
c      write(*,*) "CZM thickness = ", thickness
c      write(*,*) "CZM xsj = ", xsj
      
      return
      end
