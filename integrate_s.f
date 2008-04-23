c     Determine whether or not we need to integrate adaptively

      fixed_order_int = .false.
      if(properties(3).eq.-1) then
         call adapt_query(u_node_interface,u_max_node_interface,
     &        properties,thickness_coord,fixed_order_int)
      endif
c      write(*,*) "CZM: p(3) = ",properties(3)

c      if(.false.) then



c      write(*,*) "elcon(3,1,nmat) = ",elcon(3,1,ielmat(nelem))
c      write(*,*) "fixed_order_int = ",fixed_order_int
c      fixed_order_int = .false.

      if(fixed_order_int.or.properties(3).eq.5) then
                                ! Fifth-order integration
c         write(*,*) "CZM: 5th O ",properties(3)
         mint = 9
         do k=1,mint
            interface_plane_x = gauss2d3(1,k)
            interface_plane_y = gauss2d3(2,k)

                  call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &                 interface_z,area,thickness_coord,
     &                 nope,xl,u_node_interface,u_max_node_interface,
     &                 properties,integrand)
              
c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
c                  write(*,*) " CZM: selement =",i1,i2,s_element(i1,i2)
c                  write(*,*) " CZM: integrand =",integrand(i1,i2)
c                  write(*,*) " CZM: w =",weight2d3(k)
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 weight2d3(k)*integrand(i1,i2)
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
         enddo
      elseif(properties(3).eq.3) then
                                ! Third-order integration
         mint = 4
         do k=1,mint
            interface_plane_x = gauss2d2(1,k)
            interface_plane_y = gauss2d2(2,k)
                  
                  call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &                 interface_z,area,thickness_coord,
     &                 nope,xl,u_node_interface,u_max_node_interface,
     &                 properties,integrand)
            
c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
c           if(dabs(dabs(weight2d2(k)*integrand(i1,i2))-128).lt.0.1) then
c                  if((i1.eq.1.and.i2.eq.4).or.(i2.eq.1.and.i1.eq.4)) then
c       write(*,*) "i1,i2,k,w*i = ",i1,i2,k,weight2d2(k)*integrand(i1,i2)
c                  endif
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 weight2d2(k)*integrand(i1,i2)
c                  if(dabs(dabs(s_element(i1,i2))-64.0).lt.0.1) then
c                     write(*,*) "i1,i2,k = ",i1,i2,k
c                  endif
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
            
         enddo
c         do i1=1,3*nope
c            do i2=1,3*nope
c               write(*,*)s_element(i1,i2)-s_element(i2,i1)," ..."
c            enddo
c         enddo

      elseif(properties(3).eq.7) then
                                ! Seventh-order integration
         mint = 16
         do k=1,mint
            interface_plane_x = gauss2d4(1,k)
            interface_plane_y = gauss2d4(2,k)
                  
                  call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &                 interface_z,area,thickness_coord,
     &                 nope,xl,u_node_interface,u_max_node_interface,
     &                 properties,integrand)
            
c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 weight2d4(k)*integrand(i1,i2)
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
            
         enddo
      elseif(properties(3).eq.9) then
                                ! Ninth-order integration
         mint = 25
         do k=1,mint
            interface_plane_x = gauss2d5(1,k)
            interface_plane_y = gauss2d5(2,k)
                  
                  call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &                 interface_z,area,thickness_coord,
     &                 nope,xl,u_node_interface,u_max_node_interface,
     &                 properties,integrand)
            
c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 weight2d5(k)*integrand(i1,i2)
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
            
         enddo
      elseif(properties(3).eq.11) then
                                ! Eleventh-order integration
         mint = 36
         do k=1,mint
            interface_plane_x = gauss2d6(1,k)
            interface_plane_y = gauss2d6(2,k)
                  
                  call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &                 interface_z,area,thickness_coord,
     &                 nope,xl,u_node_interface,u_max_node_interface,
     &                 properties,integrand)
            
c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 weight2d6(k)*integrand(i1,i2)
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
            
         enddo
      elseif(properties(3).eq.0) then
                                !     Zero-order integration
c     Sum over Gauss points
         mint = 8
         do k=1,mint
c     Get max u
c     u_max(k) = eei(3,k,nelem)

            interface_plane_x = node2d3(1,k)
            interface_plane_y = node2d3(2,k)
         
c     Get integrand
               call user_e_c3d_f(interface_plane_x,interface_plane_y,
     &              interface_z,area,thickness_coord,
     &              nope,xl,u_node_interface,u_max_node_interface,
     &              properties,integrand)
               
c               write(*,*) "k,integrand = ",k,integrand(1,1) 

c     Update integral
            do i1=1,3*nope
               do i2=1,3*nope
                  
c                  integral(60*(3*i1-1)+3*i2)=integral(60*(3*i1-1)+3*i2)+ 
                  s_element(i1,i2) = s_element(i1,i2)+
     &                 node_weight2d3(k)*integrand(i1,i2)
c     &                 node_weight2d3(k)*integrand(60*(i1-1)+i2)
c     Integration area is scaled by area/4 because the Gauss points
c     are distributed over [-1,1],[-1,1] : Gauss_Area = 4
               enddo
            enddo
            
         enddo

      else
         a(1) = -1
         a(2) = -1
         b_adapt(1) = 1
         b_adapt(2) = 1
         nfuncs = 3600          !could be lower
         minpts = 1
c         maxpts = 10000         !must be >= 17

c         maxpts = 330         !must be >= 17
c         eps = .001

         maxpts = 330         !must be >= 17
         eps = .001
         iwk = 5000            !must be >= 7 * (1 + maxpts/17)/2
         
c         write(*,*) "Adaptively integrating s ..."
         call adapt_integrate_s(a,b_adapt,minpts,maxpts,eps,wk,iwk,
     &        integral,relerr,nfnevl,ifail,nfuncs,et,area,
     &        thickness_coord,nope,xl,u_node_interface,
     &        u_max_node_interface,properties,
     &        result)
c         write(*,*) "called adaptive_integrate_s ..."

         do i1=1,3*nope
            do i2=1,3*nope
              s_element(i1,i2) = integral(60*(i1-1)+i2)
            enddo
         enddo

      endif

c      endif







