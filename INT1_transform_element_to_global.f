      subroutine INT1_transform_element_2_global(f_element,xs,
     &     f_global,k)
      
      implicit none
      
      real*8 f_element(3,20),f_global(3,20),xs(3,3)
      integer k

      integer i,j

      do i=1,3
         f_global(i,k) = 0
         do j=1,3
            f_global(i,k) = f_global(i,k) + f_element(j,k)*xs(j,i)
         enddo
      enddo

      return
      end
