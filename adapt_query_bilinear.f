c     Determine whether or not we need to integrate adaptively

      subroutine adapt_query_bilinear(
     &     u_node,u_max_node,properties,thickness_coord,
     &     fixed_order_int)

      implicit none

      logical fixed_order_int
      integer thickness_coord
      real*8 u_node(3,20),u_max_node(3,20),properties(*)
c     elcon(0:ncmat_,ntmat_,*)

      integer u_region,j
      real*8 nu,lambda,lambda_init,nu_max

      nu = u_node(3,1)/properties(2)
      nu_max = u_max_node(3,1)/properties(2)
      lambda_init = .1
      lambda = max(nu_max,lambda_init)
c     Should account for discontinuities in u_max, too!!!
      if((nu.lt.0).or.(nu_max.lt.lambda_init.and.nu.lt.lambda)) then
         u_region = 0
      elseif(nu.lt.lambda.and.nu.gt.0.and.nu_max.gt.lambda_init.and.
     &        nu_max.lt.1) then
         u_region = 1
      elseif(nu.lt.1.and.nu.gt.lambda) then
         u_region = 2
      else
         u_region = 3
      endif

c      write(*,*)"node,nu,nu_max ","1",nu,
c     &     u_max_node(3,1)/elcon(2,1,nmat)
      fixed_order_int = .true.  ! default, if possible
      do j = 2,16
         if((j.lt.5).or.((j.gt.8).and.(j.lt.13))) then
            nu = u_node(3,j)/properties(2)
            nu_max = u_max_node(3,j)/properties(2)
            lambda = max(nu_max,lambda_init)
         elseif(((j.gt.4).and.(j.lt.9)).or.
     &           ((j.gt.12).and.(j.lt.17))) then
c     Had to comment out next line because SGI f77 didn't like it.
c     But I don't plan on using bilinear, so right now I just want this to
c     compile without complaint.
c     Whatever compiler ABAQUS uses doesn't like it either.
c            cycle
            write(*,*)"ERROR j = ",j
         else
            write(*,*)"ERROR j = ",j
         endif
         if((nu.lt.0).or.(nu_max.lt.lambda_init.and.nu.lt.lambda)) then
            if(u_region.ne.0) then
               fixed_order_int = .false.
            endif
         elseif(nu.lt.lambda.and.nu.gt.0.and.nu_max.gt.lambda_init.and.
     &           nu_max.lt.1) then
            if(u_region.ne.1) then
               fixed_order_int = .false.
            endif
         elseif(nu.lt.1.and.nu.gt.lambda) then
            if(u_region.ne.2) then
               fixed_order_int = .false.
            endif
         else
            if(u_region.ne.3) then
               fixed_order_int = .false.
            endif
         endif
c         if(fixed_order_int) then
c            
c         else
c            write(*,*)"node,nu,nu_max ",j,nu,
c     &     u_max_node(3,j)/elcon(2,1,nmat)
c         endif
      enddo

      return
      end
