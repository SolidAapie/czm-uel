      subroutine INT1_transform_global_2_element(u_global,xsi,
     &     u_element,k)
      
      implicit none
      
      real*8 u_element(3,20),u_global(3,20),xsi(3,3)
      integer k

      integer i,j

      do i=1,3
         u_element(i,k) = 0
         do j=1,3
c            write(*,*) "CZM:k,i,j = ",k,i,j
c            write(*,*) "CZM:u_g = ",u_global(j,k)
c            write(*,*) "CZM:xsi = ",xsi(j,i)
            u_element(i,k) = u_element(i,k) + u_global(j,k)*xsi(j,i)
         enddo
      enddo

      return
      end







