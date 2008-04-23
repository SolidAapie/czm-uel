c     Determine whether or not we need to integrate adaptively

      fixed_order_int = .false.
      if(properties(3).eq.-1) then
         call adapt_query(u_node_interface,u_max_node_interface,
     &        properties,thickness_coord,fixed_order_int)
      endif
c      write(*,*) "elcon(3,1,nmat) = ",elcon(3,1,ielmat(i))
c      write(*,*) "fixed_order_int = ",fixed_order_int
c      fixed_order_int = .false.
      if(fixed_order_int.or.properties(3).eq.5) then
                                ! Fifth-order integration
c     Loop over integration points
c         write(*,*) "CZM 5th order"
         mint = 9
         do k=1,mint
            interface_plane_x = gauss2d3(1,k)
            interface_plane_y = gauss2d3(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 weight2d3(k)*integrand(m,j)
               enddo
            enddo
         enddo
      elseif(properties(3).eq.3) then
                                ! Third-order integration
c     Loop over integration points
         mint = 4
         do k=1,mint
            interface_plane_x = gauss2d2(1,k)
            interface_plane_y = gauss2d2(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 weight2d2(k)*integrand(m,j)
               enddo
            enddo
         enddo
      elseif(properties(3).eq.7) then
                                ! Seventh-order integration
c     Loop over integration points
c         write(*,*) "CZM 7th order"
         mint = 16
         do k=1,mint
            interface_plane_x = gauss2d4(1,k)
            interface_plane_y = gauss2d4(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 weight2d4(k)*integrand(m,j)
               enddo
            enddo
         enddo
      elseif(properties(3).eq.9) then
                                ! Ninth-order integration
c     Loop over integration points
         mint = 25
         do k=1,mint
            interface_plane_x = gauss2d5(1,k)
            interface_plane_y = gauss2d5(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 weight2d5(k)*integrand(m,j)
               enddo
            enddo
         enddo

      elseif(properties(3).eq.11) then
                                ! Eleventh-order integration
c     Loop over integration points
         mint = 36
         do k=1,mint
            interface_plane_x = gauss2d6(1,k)
            interface_plane_y = gauss2d6(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 weight2d6(k)*integrand(m,j)
               enddo
            enddo
         enddo

      elseif(properties(3).eq.0) then
                                !     Zero-order integration
         mint = 8
         do k=1,mint
            interface_plane_x = node2d3(1,k)
            interface_plane_y = node2d3(2,k)
            
               call user_element_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,nope,
     &              xl,u_node_interface,u_max_node_interface,properties,
     &              integrand)
            
            do j=1,nope
               do m=1,3
                  fnl_element(m,j) = fnl_element(m,j) + 
     &                 node_weight2d3(k)*integrand(m,j)
               enddo
            enddo
         enddo

      else
                                ! Adaptive integration
c#define insert_integrand    include "integrand_f.f"
c#define nfuncs_             20
c#include "adapt_integrate.f"


c$$$c     Loop over integration points
c$$$         do j=1,nope
c$$$c            integral(j) = fnl(3,j)
c$$$            integral(j) = 0
c$$$         enddo
c$$$         do k=1,mint
c$$$            xi = gauss2d3(1,k) 
c$$$            ze = gauss2d3(2,k)
c$$$            
c$$$
c$$$            call user_element_f(xi,et,ze,thickness,thickness_coord,nope,
c$$$     &              xl,u_node,u_max_node,properties,ielmat(i),
c$$$     &              ncmat_,ntmat_,integrand)
c$$$            
c$$$            do j=1,nope
c$$$               integral(j) = integral(j) + weight2d3(k)*integrand(j)
c$$$            enddo
c$$$         enddo


         a(1) = -1
         a(2) = -1
         b_adapt(1) = 1
         b_adapt(2) = 1
c         nfuncs = nope-4
c         nfuncs = 3*(nope-4)
         nfuncs = 3*nope
         minpts = 1
c         maxpts = 10000         !must be >= 17

c         maxpts = 330         !must be >= 17
c         eps = .001

         maxpts = 3300         !must be >= 17
         eps = .0001
         iwk = 21000            !must be >= 7 * (1 + maxpts/17)/2
         
c         write(*,*) "Adaptively integrating f ",i,"..."
         call adapt_integrate_f(a,b_adapt,minpts,maxpts,eps,wk,iwk,
     &        result,relerr,nfnevl,ifail,nfuncs,et,area,
     &        thickness_coord,nope,xl,u_node_interface,
     &        u_max_node_interface,properties,
     &        integrand)
c         write(*,*) "ifail = ",ifail
         do j=1,nope
            do m=1,3
               fnl_element(m,j) = fnl_element(m,j) + result(m+3*(j-1))
            enddo
         enddo
c$$$         write(*,*) "fixed_order_int = ", fixed_order_int
c$$$
c$$$         do j=1,nope-4
c$$$            write(*,*)"i,j",i,j,"   fixed O ",integral(j),"   adapt ",
c$$$     &           result(j),"   re",(integral(j)/result(j)-1)
c$$$          enddo

      endif
