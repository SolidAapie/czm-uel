c     Zero mass matrix

      subroutine mass_matrix(amatrix,coords,u,nprops,nnode,nsvars,
     &	properties,svars,ndofel)

c      parameter ( zero = 0.D0)

      integer nope,nstate_

      dimension amatrix(ndofel,ndofel)

      real*8 xl(3,20),vold(3,20),
     &  properties(nprops),
     &  s(ndofel,ndofel),sm(ndofel,ndofel),xstate(nsvars)

c      call INT1_e_c3d(xl,vold,nope,s,sm,properties,xstate,nstate_)
      call INT1_e_c3d(coords,vold,nnode,s,sm,properties,svars,nstate_)

      do i=1,ndofel
         do j=1,ndofel
c      do i=1,4
c         do j=1,3
c            write(*,*) "CZM: i,j = ",i,j
            amatrix(i,j) = 0
         enddo
      enddo

      return
      end
