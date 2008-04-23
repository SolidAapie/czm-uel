      subroutine residual_vector(coords,u,nprops,nnode,nsvars,props,
     & svars,ndofel,update_state,fn)
      
      logical calcul_fn,update_state

      integer node,direction

      real*8 xl(3,nnode),vl(3,nnode),
     &  properties(nprops),
     &  fn(3,nnode),xstate(nsvars)

      real*8 svars(nsvars),u(ndofel)
      
      do node = 1,nnode
         do direction = 1,3
            vl(direction,node) = u(3*(node-1) + direction)
c            write(*,*) "CZM: vl = ",direction,node,vl(direction,node)
         enddo
      enddo
      
      
c      svars(1) = 42
c      write(*,*) "CZM: nsvars = ",nsvars
c      call INT1_element_results(xl,vl,calcul_fn,nope,properties,
c     &  fn,xstate,nstate_)
c      write(*,*) "CZM: calling e_r "
      call INT1_element_results(coords,vl,.true.,update_state,nnode,
     &     props,fn,svars,nsvars)
c      write(*,*) "CZM: returned from e_r "

c      do node = 1,nnode
c         do direction = 1,3
c            write(*,*) "CZM: fn rv = ",direction,node,fn(direction,node)
c         enddo
c      enddo
      

      return
      end
