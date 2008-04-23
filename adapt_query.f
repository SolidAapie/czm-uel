c     Determine whether or not we need to integrate adaptively

      subroutine adapt_query(
     &     u_node,u_max_node,properties,thickness_coord,
     &     fixed_order_int)

      implicit none

      logical fixed_order_int
      integer thickness_coord
      real*8 u_node(3,20),u_max_node(3,20),properties(*)
c     elcon(0:ncmat_,ntmat_,*)

      integer u_region,j
      real*8 nu,lambda,lambda_init,nu_max

      if(properties(4).eq.1) then
         call adapt_query_bilinear(
     &        u_node,u_max_node,properties,thickness_coord,
     &        fixed_order_int)
      elseif(properties(4).eq.2) then
         call adapt_query_chaboche(
     &        u_node,u_max_node,properties,thickness_coord,
     &        fixed_order_int)
      else
         write(*,*) "ERROR: invalid CZM: ",properties(4)
         stop
      endif

      return
      end
