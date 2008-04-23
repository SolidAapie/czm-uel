c     Determine whether or not we need to integrate adaptively

      subroutine adapt_query_chaboche(
     &     u_node_interface,u_max_node_interface,properties,
     &     thickness_coord,fixed_order_int)

      implicit none

      logical fixed_order_int
      integer thickness_coord,nmat,ncmat_,ntmat_
      real*8 u_node_interface(3,20),u_max_node_interface(3,20),
     &     properties(*)
c     elcon(0:ncmat_,ntmat_,*)

      integer j,M
      real*8 nu(3),lambda,lambda_max,lambda_init,nu_max(3)
      integer nope,lambda_region1,lambda_region,nodes
      real*8 nu1,nu_max1,delta(3),E(3)
      integer node_permutation(8,3)

      data node_permutation /
     &     2, 6, 7, 3, 18, 14, 19, 10,
     &     5, 6, 7, 8, 13, 14, 15, 16,
     &     1, 2, 6, 5, 9, 18, 13, 17 /

c     First, check if all nu(3) > 0 or all nu(3) < 0
c     -There will be a discontinuity in derivatives of lambda at nu(3) = 0 
      fixed_order_int = .true.  !default
      nodes = 8
                                ! Determine delta,E
      do m=1,3
         delta(m) = properties(2)
         E(m) = properties(1)
      enddo
c$$$      call INT1_element_info_nope('INT1  ', nope)
c$$$      nope = nope-4             ! don't count unused midside nodes
c$$$                                !need to do permutation for which nodes to check
      lambda_init = .1
c$$$      if(u_node_interface(3,1).gt.0) then
c$$$         do j=2,nodes
c$$$            if(u_node_interface(3,j).lt.0) then !opposite sign
c$$$               fixed_order_int = .false. 
c$$$               return
c$$$            endif
c$$$         enddo
c$$$      else
c$$$         do j=2,nodes
c$$$            if(u_node_interface(3,j).gt.0) then !opposite sign
c$$$               fixed_order_int = .false.
c$$$               return               
c$$$            endif            
c$$$         enddo
c$$$      endif 

c     Second, check for lambda discontinuities 
c       lambda = lambda_max

                                ! Determine nu,nu_max at first node
      do m=1,3
         nu(m) = u_node_interface(m,node_permutation(1,thickness_coord))
     &        /delta(m)
         nu_max(m) = u_max_node_interface(m,
     &        node_permutation(1,thickness_coord))/delta(m)
      enddo
c$$$      nu = u_node_interface(3,node_permutation(1,thickness_coord))
c$$$     &     /elcon(2,1,nmat)
c$$$      nu_max = u_max_node_interface(3,
c$$$     &     node_permutation(1,thickness_coord))/elcon(2,1,nmat)

                                ! Determine lambda
      call compute_lambda(nu,lambda)
c$$$      if(nu.gt.0) then       ! norm: sqrt(nu1**2 + nu2**2 + nu3**2)
c$$$         lambda = abs(nu)
c$$$      else
c$$$         lambda = 0
c$$$      endif

                                ! Determine lambda_max
      call compute_lambda(nu_max,lambda_max)
c$$$      if(nu_max.gt.0) then       ! norm: sqrt(nu1**2 + nu2**2 + nu3**2)
c$$$         lambda_max = abs(nu_max) 
c$$$      else
c$$$         lambda_max = 0
c$$$      endif
      if(lambda.gt.lambda_max) then
         do j=2,nodes
            do m=1,3
               nu(m) = u_node_interface(m,
     &              node_permutation(j,thickness_coord))/delta(m)
               nu_max(m) = u_max_node_interface(m,
     &              node_permutation(j,thickness_coord))/delta(m)
            enddo
            call compute_lambda(nu,lambda)
            call compute_lambda(nu_max,lambda_max)
c     Don't be too severe if lambda_max is small (lambda_max.gt.lambda_init)
            if(lambda.lt.lambda_max.and.lambda.lt.1
     &           ) then
c     &           .and.lambda.gt.lambda_init) then
c               write(*,*) "     lambda discont. 1"
               fixed_order_int = .false. 
            endif
         enddo
      else
         do j=2,nodes
            do m=1,3
               nu(m) = u_node_interface(m,
     &              node_permutation(j,thickness_coord))/delta(m)
               nu_max(m) = u_max_node_interface(m,
     &              node_permutation(j,thickness_coord))/delta(m)
            enddo
            call compute_lambda(nu,lambda)
            call compute_lambda(nu_max,lambda_max)
c     Don't be too severe if lambda_max is small (lambda_max.gt.lambda_init)
            if(lambda.gt.lambda_max.and.lambda_max.lt.1
     &           ) then
c     &           .and.lambda.gt.lambda_init) then
c               write(*,*) "     lambda discont. 2"
               fixed_order_int = .false. 
            endif
         enddo
      endif


c     Third, check for discontinuities in derivatives of F(lambda)
      do m=1,3
         nu(m) = u_node_interface(m,node_permutation(1,thickness_coord))
     &        /delta(m)
         nu_max(m) = u_max_node_interface(m,
     &        node_permutation(1,thickness_coord))/delta(m)
      enddo
      nu1 = nu(3)
      nu_max1 = nu_max(3)
      call compute_lambda(nu,lambda)
      call compute_lambda(nu_max,lambda_max)
      call chaboche_lambda_region(lambda,lambda_max,lambda_region1)
      do j=2,nodes
         do m=1,3
            nu(m) = u_node_interface(m,
     &           node_permutation(j,thickness_coord))/delta(m)
            nu_max(m) = u_max_node_interface(m,
     &           node_permutation(j,thickness_coord))/delta(m)
         enddo
         call compute_lambda(nu,lambda)
         call compute_lambda(nu_max,lambda_max)
         call chaboche_lambda_region(lambda,lambda_max,lambda_region)
         if(lambda_region.ne.lambda_region1) then
c$$$            write(*,*) "     F      discont.  region1,region = ",
c$$$     &           lambda_region1,lambda_region,"   nu1,nu_max1 = ",
c$$$     &           nu1,nu_max1
            fixed_order_int = .false. !different region; the element 
                                      !contains a discontinuity
            return
         endif         
      enddo



      return
      end
