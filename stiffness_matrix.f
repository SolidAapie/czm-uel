      subroutine stiffness_matrix(amatrx,coords,u,nprops,nnode,nsvars,
     &	properties,svars,ndofel)

      integer nope,nstate_
      
      integer node,direction

c      dimension amatrx(ndofel,ndofel)
      
      real*8 xl(3,20),vl(3,20),
     &     properties(nprops),
     &     s(ndofel,ndofel),sm(ndofel,ndofel),xstate(nsvars)
      

      real*8 svars(nsvars),u(ndofel),amatrx(ndofel,ndofel)
      

      do node = 1,nnode
         do direction = 1,3
            vl(direction,node) = u(3*(node-1) + direction)
         enddo
      enddo      
      
c      call INT1_e_c3d(xl,voldl,nope,s,sm,properties,xstate,nstate_)
      call INT1_e_c3d(coords,vl,nnode,s,sm,properties,svars,nsvars)

      do i1=1,ndofel
         do i2=1,ndofel
c            write(*,*) "CZM: s = ",i1,i2,s(i1,i2)
c            write(*,*) "CZM: amatrx = ",i1,i2,amatrx(i1,i2)
            amatrx(i1,i2) = s(i1,i2)
         enddo
      enddo

      return
      end
