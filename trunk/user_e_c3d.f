      subroutine INT1_e_c3d(xl,voldl,nope,s,sm,properties,xstate,nstate_
     &	)

      implicit none
c      character*6 lakonl
      integer 
     &  j,k,m,i1,i2,
     &  nope,mint_,nstate_
      real*8 p(3),q(3)
      real*8 xl(3,20),
     &  s(60,60),properties(*),
     &  density,voldl(3,20),xi,et,ze,sm(60,60)
      real*8 xstate(nstate_)
c     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*)

      integer mint,thickness_coord
      real*8  x1,x2,x3,
     &     integrand(60,60),u_node(3,20),u_max_node(3,20)
!
      real*8 area, integral(3600)
      logical fixed_order_int

      integer minpts,maxpts,iwk      
      real*8 a(2),b_adapt(2),eps,wk(3600,5000)
      real*8 result(3600),relerr
      integer nfnevl,ifail,nfuncs

      real*8 node2d3(2,8),node_weight2d3(8),gauss2d2(2,4),weight2d2(4),
     &     gauss2d3(2,9),weight2d3(9),gauss2d4(2,16),weight2d4(16),
     &     gauss2d5(2,25),weight2d5(25),gauss2d6(2,36),weight2d6(36)

      real*8 interface_plane_x,interface_plane_y,interface_z,
     &     voldl_element(3,20),s_element(60,60),xs(3,3),xsi(3,3)
      real*8 xi_node(20),et_node(20),ze_node(20),xs1(3,3),xsi1(3,3)
     &     ,xs2(3,3),xsi2(3,3),tmp3x3(3,3),s_3x3(3,3),s_element_3x3(3,3)
      integer u_permutation1(20),u_permutation2(20),u_permutation3(20)

      real*8 k_unused_node,m_unused_node
      integer unused_node(4,3)

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

c      call INT1_element_info_nope(lakonl, nope)
c      mint = 9

      density = properties(5)

      do i1=1,3*nope
         do i2=1,3*nope
            integrand(i1,i2) = 0   !initialize, only because abaqus will choke if we don't
         enddo
      enddo

            ! Transform nodal displacements into element coord system
      do k=1,nope
         xi = xi_node(k)
         et = et_node(k)
         ze = ze_node(k)
         call INT1_coord_transform(xi,et,ze,xl,xs,xsi)
         call INT1_transform_global_2_element(voldl,
     &        xsi,voldl_element,k)
      enddo
                                ! Determine elment orientation
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
c         interface_z = (xl(1,1)-xl(1,2))/dabs(xl(1,1)-xl(1,2))
      elseif(x2.lt.x1.and.x2.lt.x3) then
         thickness_coord = 2    !use et for top and bottom of element
c         interface_z = (xl(3,1)-xl(3,5))/dabs(xl(3,1)-xl(3,5))
      else
         thickness_coord = 3    !use ze for top and bottom of element
c         interface_z = (xl(2,1)-xl(2,4))/dabs(xl(2,1)-xl(2,4))
      endif

!     this calculation of the area is linear in the sense that considers the original shape of the element rather than the current shape (xl+u)
!     it also assumes that the original shape of the element is a parallelogram
      do j = 1,3
         p(j) = xl(j,unused_node(2,thickness_coord))-xl(j,unused_node(1,
     &   thickness_coord))
         q(j) = xl(j,unused_node(3,thickness_coord))-xl(j,unused_node(1,
     &	   thickness_coord))
      enddo
      area = sqrt( (p(1)*p(1)+p(2)*p(2)+p(3)*p(3)) 
     &     * (q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) 
     &     - (p(1)*q(1)+p(2)*q(2)+p(3)*q(3))**2 )

c     Find u at the nodes
c$$$      do j=1,nope
c$$$         do m=1,3
c$$$            if((j.lt.5).or.((j.gt.8).and.(j.lt.13))) then
c$$$               u_node(m,j) = voldl(m,j)-voldl(m,j+4)
c$$$            elseif(((j.gt.4).and.(j.lt.9)).or.
c$$$     &              ((j.gt.12).and.(j.lt.17))) then
c$$$               u_node(m,j) = voldl(m,j)-voldl(m,j-4)
c$$$            else
c$$$               u_node(m,j) = 0
c$$$            endif
c$$$         enddo
c$$$      enddo
      if(thickness_coord.eq.1) then
         do j=1,nope
            do m=1,3
               u_node(m,j) = voldl_element(m,j)-
     &              voldl_element(m,u_permutation1(j))
            enddo
         enddo               
      elseif(thickness_coord.eq.2) then
         do j=1,nope
            do m=1,3
               u_node(m,j) = voldl_element(m,j)-
     &              voldl_element(m,u_permutation2(j))
            enddo
         enddo
      else
         do j=1,nope
            do m=1,3
               u_node(m,j) = voldl_element(m,j)-
     &              voldl_element(m,u_permutation3(j))               
            enddo
         enddo
      endif

c      write(*,*) "nelem = ",nelem
c      do j=1,1000
c         write(*,*) "xstate = ",m,j,xstate(1,1,j)
c      enddo

c     Find u_max at the nodes
      do j=1,nope
         do m=1,3
c            if(nelem.eq.889) then
c               write(*,*) "m,j,asdf = ",m,j,3*(j-1)+m+30
c               write(*,*) "u_max_node = ",u_max_node(m,j)
c               write(*,*) "xstate = ",xstate(3*(j-1)+m+30,1,nelem)
c            endif
            u_max_node(m,j) = xstate(3*(j-1)+m+30)
c            u_max_node(m,j) = 0
         enddo
      enddo

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

c     INTEGRATE OVER THE AREA

c$$$      do j=1,3600
c$$$         integral(j) = 0
c$$$      enddo
      do i1=1,3*nope
         do i2=1,3*nope
            s_element(i1,i2) = 0
         enddo
      enddo

      include "integrate_s.f"

c     Update stiffness matrix
c$$$         do i1=1,nope
c$$$            do i2=1,nope
c$$$               
c$$$               s(3*i1,3*i2) = integral(60*(3*i1-1)+3*i2)
c$$$c     Integration area is scaled by area/4 because the Gauss points
c$$$c      are distributed over [-1,1],[-1,1] : Gauss_Area = 4
c$$$c            write(*,*) "s(",nelem,",",3*i1,",",3*i2,") = ",s(3*i1,3*i2)
c$$$            enddo
c$$$         enddo
c     Which nodes are unused depends on the element orientation.
      k_unused_node = properties(1)*x1*x2*x3*10000000
      m_unused_node = density*x1*x2*x3   !rhcon(1,1,1)*x1*x2*x3
      do i1=1,4
         k = 3*unused_node(i1,thickness_coord)-3
         do i2=1,3
            s_element(k+i2,k+i2) = k_unused_node
            sm(k+i2,k+i2) = m_unused_node !give them some mass so they're 
                                !stable in an explicit integration
         enddo
      enddo
c$$$      s_element(3*9-2,3*9-2) = elcon(1,1,ielmat(nelem))*x1*x2*x3
c$$$      s_element((3*9)+1-2,(3*9)+1-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*9)+2-2,(3*9)+2-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*11)+1-2,(3*11)+1-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*11)+2-2,(3*11)+2-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*11)-2,(3*11)-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*13)+1-2,(3*13)+1-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*13)+2-2,(3*13)+2-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*13)-2,3*13-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*15)+1-2,(3*15)+1-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*15)+2-2,(3*15)+2-2) = s_element(3*9-2,3*9-2)
c$$$      s_element((3*15)-2,(3*15)-2) = s_element(3*9-2,3*9-2)
                                ! Transform to global coordinates
      do i1=1,nope
         xi = xi_node(i1)
         et = et_node(i1)
         ze = ze_node(i1)
         call INT1_coord_transform(xi,et,ze,xl,xs1,xsi1)
         do i2=1,nope
            xi = xi_node(i2)
            et = et_node(i2)
            ze = ze_node(i2)
            call INT1_coord_transform(xi,et,ze,xl,xs2,xsi2)
            do k=1,3
               do m=1,3
                  s_element_3x3(k,m) = s_element(3*(i1-1)+k,3*(i2-1)+m)
               enddo
            enddo
c$$$            call matrix_mult(s_element_3x3,xsi2,tmp3x3,3)
c$$$            call matrix_mult(xs1,tmp3x3,s_3x3,3)
            call matrix_mult(s_element_3x3,xs2,tmp3x3,3)
            call matrix_mult(xsi1,tmp3x3,s_3x3,3)
            do k=1,3
               do m=1,3
                  s(3*(i1-1)+k,3*(i2-1)+m) = s_3x3(k,m)
               enddo
            enddo
         enddo
      enddo      

c     KLUDGE: put in a tiny bit of stiffness to prevent singularity
c        -only needed when one side of interface is not tied to another element
c         and it fully debonds.
c      do i1=1,60
c         write(*,*) "CZM: s(",i1,") = ",s(i1,i1)
c         if(abs(s(i1,i1)).lt.1e-6) then
c            s(i1,i1) = 1e-6
c         endif
c      enddo
c$$$      do i1=1,16
c$$$         do i2=1,16
c$$$            do k=1,3
c$$$               do m=1,3
c$$$                  if(k.ne.3.or.m.ne.3) then
c$$$                     s(3*(i1-1)+k,3*(i2-1)+m) = 0
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo

c$$$c     Which nodes are unused depends on the element orientation.
c$$$c      s(9,9) = elcon(1,1,ielmat(nelem))*x1*x2*x3*thickness
c$$$      s(3*9-2,3*9-2) = -elcon(1,1,ielmat(nelem))*x1*x2*x3
c$$$      s((3*9)+1-2,(3*9)+1-2) = s(3*9-2,3*9-2)
c$$$      s((3*9)+2-2,(3*9)+2-2) = s(3*9-2,3*9-2)
c$$$      s((3*11)+1-2,(3*11)+1-2) = s(3*9-2,3*9-2)
c$$$      s((3*11)+2-2,(3*11)+2-2) = s(3*9-2,3*9-2)
c$$$      s((3*11)-2,(3*11)-2) = s(3*9-2,3*9-2)
c$$$      s((3*13)+1-2,(3*13)+1-2) = s(3*9-2,3*9-2)
c$$$      s((3*13)+2-2,(3*13)+2-2) = s(3*9-2,3*9-2)
c$$$      s((3*13)-2,3*13-2) = s(3*9-2,3*9-2)
c$$$      s((3*15)+1-2,(3*15)+1-2) = s(3*9-2,3*9-2)
c$$$      s((3*15)+2-2,(3*15)+2-2) = s(3*9-2,3*9-2)
c$$$      s((3*15)-2,(3*15)-2) = s(3*9-2,3*9-2)

c$$$      do j=49,60                !The unused midside nodes.
c$$$c         s(j,j) = s(3,3)  !fill in diagonal so it's not singular
c$$$         s(j,j) = elcon(1,1,ielmat(nelem))*x1*x2*x3*thickness
c$$$         sm(j,j) = rhcon(1,1,1)*x1*x2*x3 !give them some mass so they're
c$$$                                !stable in an explicit integration
c$$$      enddo

c      write(*,*) "s(",nelem,") = "
c      write(*,*) "s() = "
c      do i1=1,60
c         write(*,*) ">",i1
c         write(*,*) (s(i1,i2),i2=1,60)
c      enddo

! hook for breakpoints
c      E_elastic = 0
!
c      do i1=1,3*nope
c         do i2=1,3*nope
c            s_element(i1,i2) = (s(i1,i2)+s(i2,i1))/2.0
c            s(i1,i2) = 0
c         enddo
c      enddo
c      do i1=1,3*nope
c         do i2=1,3*nope
c            s(i1,i2) = s_element(i1,i2)
c         enddo
c      enddo
c      do i1=1,3*nope
c         s(i1,i1) = 1
c         s(i1,i1) = 1000*s_element(i1,i1)
c      enddo
c      write(*,*) "s = ["
c      do i1=1,3*nope
c         do i2=1,3*nope
c            if(dabs(s(i1,i2)).gt.1e10) then
c               s(i1,i2) = 0
c            endif
c           write(*,*)"i1,i2,=",i1,i2,s(i1,i2)-s(i2,i1),s(i1,i2),s(i2,i1)
c           write(*,*)"i1,i2,=",i1,i2,s(i1,i2)-s(i2,i1)
c           write(*,*)s_element(i1,i2)-s_element(i2,i1),s_element(i1,i2)," ..."
c            write(*,*) s(i1,i2)," ..."
c         enddo
c         write(*,*) ";"
c      enddo
c      write(*,*) "];"

      return
      end


