
      subroutine user_e_c3d_f(interface_plane_x,interface_plane_y,
     &   interface_z,area,thickness_coord,nope,xl,u_node_interface,
     &     u_max_node_interface,properties,f)
      implicit none
!
      integer thickness_coord,nope
      real*8 interface_plane_x,interface_plane_y,
     &     interface_z,area,properties(*),
     &     xl(3,20),u_node(3,20),u_max_node(3,20),f(60,60)
c     elcon(0:ncmat_,ntmat_,*),

      integer j,m,i1,i2,k,coord_permutation(3)
      real*8 xi,et,ze,xsj,shp1(4,20),shp2(4,20),u(3),u_max(3),dudq(20),
     &     dTdu(3,3),u_interface(3),u_max_interface(3),
     &     dTdu_interface(3,3),
     &     u_node_interface(3,20),u_max_node_interface(3,20)
                                !
      
      if(thickness_coord.eq.1) then
c         xi = interface_z
         xi = 1
         et = interface_plane_x
         ze = interface_plane_y
      elseif(thickness_coord.eq.2) then
         xi = interface_plane_x
c         et = interface_z
         et = 1
         ze = interface_plane_y
      elseif(thickness_coord.eq.3) then
         xi = interface_plane_y
         et = interface_plane_x
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
         write(*,*) "ERROR in user_e_c3d(): bad thickness_coord"
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

c     Calculate u and dudq
      do m=1,3
         u_interface(m) = 0
         u_max_interface(m) = 0
      enddo
      do j=1,nope
         do m=1,3
c            u(k) = u(k) + (shp2(4,j)-shp1(4,j))*voldl(3,j)
            u_interface(m) = u_interface(m) + 
     &           shp1(4,j)*u_node_interface(m,j)
         u_max_interface(m) = u_max_interface(m) + 
     &           shp1(4,j)*u_max_node_interface(m,j)
         enddo
         dudq(j) = shp1(4,j)-shp2(4,j)
      enddo      

c     Calculate dTdu
      call INT1_interface_dtdu(u_interface,u_max_interface,
     &     dTdu_interface,properties)

c     Transform dTdu from interface to element coordinates.
      do k=1,3
         do m=1,3
            dTdu(coord_permutation(k),coord_permutation(m)) = 
     &           dTdu_interface(k,m)
         enddo
      enddo

c      write(*,*) "dudq(1),xsj,thickness = ",dudq(1),xsj,
c     &     thickness
c      do k=1,3
c         do m=1,3
c            write(*,*)"dTdu = ",k,m,dTdu(k,m)
c         enddo
c      enddo

      
      do i1=1,nope
         do i2=1,nope
            do k=1,3
               do m=1,3
c     f(60*(k*i1-1)+m*i2)=dTdu(k,m)*dudq(i1)*dudq(i2)*xsj/thickness*2
c                  f(3*(i1-1)+k,3*(i2-1)+m)=dTdu(k,m)*dudq(i1)*dudq(i2)*
c     &                 4.6720125d-7
c     &                 xsj/thickness*2.0
                  f(3*(i1-1)+k,3*(i2-1)+m)=dTdu(k,m)*dudq(i1)*dudq(i2)*
     &                 area/4.0 ! need to divide by 4 bc the area in element coords is 4
c$$$                  if (xsj/thickness*2.0.lt.1d-6) then
c$$$                     f(3*(i1-1)+k,3*(i2-1)+m)=dTdu(k,m)*dudq(i1)*dudq(i2)*
c$$$     &                    4.6720125d-7
c$$$                  else
c$$$                     f(3*(i1-1)+k,3*(i2-1)+m)=dTdu(k,m)*dudq(i1)*dudq(i2)*
c$$$     &                    5.80644d-5
c$$$                  endif
c                  write(*,*)"CZM user_e_c3d_f: area = ",i1,i2,k,m,f(3*(i1-1)+k,3*(i2-1)+m)/dTdu(k,m)/dudq(i1)/dudq(i2)
c     &                 0.0000595478616
c     &                 .25
         !########
               enddo
            enddo
c           f(60*(3*i1-1)+3*i2)=dTdu(3)*dudq(i1)*dudq(i2)*xsj/thickness*2
c           f(60*(3*i1-1)+3*i2)=dTdu(3)*dudq(i1)*dudq(i2)*.25
         enddo
      enddo

c      j = 0 !hook for breakpoints

      return
      end
      
