      PROGRAM PROGNAME

      END PROGRAM
      
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,
     &     NDOFEL,NRHS,NSVARS,PROPS,NPROPS,COORDS,
     &     MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &     ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,
     &     DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     &     PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0,
     &     ONE = 1.D0, SEVEN=7.0D0, EIGHT=8.0D0 )
C     
      DIMENSION RHS(MLVARX,*),
     1     AMATRX(NDOFEL,NDOFEL),
     2     SVARS(NSVARS),ENERGY(8),PROPS(*),
     3     COORDS(MCRD,NNODE),U(NDOFEL),
     4     DU(MLVARX,*),V(NDOFEL),A(NDOFEL),
     5     TIME(2),PARAMS(3),JDLTYP(MDLOAD,*),
     6     ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     7     PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     8     JPROPS(*)
C     
c     User element to implement a cohesive zone model (CZM)
C     
c      implicit none

      logical compute_residual_vector, compute_stiffness_matrix,
     &	 compute_mass_matrix
      logical update_state

      real*8 fn(3,nnode)

      integer node,direction   !,i1,i2

c      do i1=1,ndofel
c         do i2=1,ndofel
c            write(*,*) "CZM: s = ",i1,i2,s(i1,i2)
c            amatrx(i1,i2) = 0     !initialize, only because abaqus chokes if we don't
c            s(i1,i2) = 0
c            write(*,*) "CZM: amatrix = ",i1,i2,amatrix(i1,i2)
c         enddo
c      enddo

c      write(*,*) "CZM: Element number ",jelem

c     Error checking based on lflags
      if(.false.) then
      if(lflags(1).gt.12) then
         write(*,*) "ERROR in CZM user element: Procedure type ",
     &        lflags(1)," not supported."
         write(*,*) "     See Abaqus Standard Manual 5.1.2, pp. 377-378
     &        for procedure types."
      endif
      if(lflags(2).eq.0) then
         write(*,*) "ERROR in CZM user element: Small-displacement 
     &        selected. Must be large-displacement"
         stop
      endif
      if(lflags(3).eq.3) then
         write(*,*) "ERROR in CZM user element: Damping matrix
     &        requested"
         stop
      endif
      if(lflags(3).eq.100) then
         write(*,*) "ERROR in CZM user element: Perturbation quantities
     &        requested"
         stop
      endif
      if(lflags(4).eq.1) then
         write(*,*) "ERROR in CZM user element: Linear perturbation
     &        step requested"
         stop
      endif
      endif

c     Find out what to do based on lflags(3)
      if(lflags(3).eq.1) then
         compute_residual_vector = .true.
         compute_stiffness_matrix = .true.
         compute_mass_matrix = .false.
c         update_state = .false.
         update_state = .true.
      endif
      if(lflags(3).eq.2) then
         compute_residual_vector = .false.
         compute_stiffness_matrix = .true.
         compute_mass_matrix = .false.
         update_state = .false.
      endif
      if(lflags(3).eq.4) then
         compute_residual_vector = .false.
         compute_stiffness_matrix = .false.
         compute_mass_matrix = .true.
         update_state = .false.
      endif
      if(lflags(3).eq.5) then
         compute_residual_vector = .true.
         compute_stiffness_matrix = .false.
         compute_mass_matrix = .false.
c         update_state = .true.
         update_state = .false.
      endif
      if(lflags(3).eq.6) then
         compute_residual_vector = .true.
         compute_stiffness_matrix = .false.
         compute_mass_matrix = .true.
         update_state = .false.
      endif
c     for static, always update state (it won't get lflags(3).eq.5, half-step residual)
      if(lflags(1).eq.1.or.lflags(1).eq.2) then
         update_state = .true.
      endif
c      write(*,*) "CZM: lflags(3) = ",lflags(3)
c      write(*,*) "CZM: ndofel = ",ndofel

      if(compute_residual_vector) then
c         write(*,*) "CZM: update_state = ",update_state
         do node = 1,nnode
            do direction = 1,3
               fn(direction,node) = 0   !initialize, only because abaqus chokes if we don't
            enddo
         enddo
         call residual_vector(coords,u,nprops,nnode,nsvars,props,svars,
     &        ndofel,update_state,fn)
c         write(*,*) "CZM: mlvarx = ",mlvarx
        
         do node = 1,nnode
            do direction = 1,3
c               rhs(3*(node-1) + direction,1) = -fn(direction,node)/1.58022963077 !transfer nodal forces to rhs
               rhs(3*(node-1) + direction,1) = -fn(direction,node) !transfer nodal forces to rhs
c               write(*,*) "CZM: fn = ",direction,node,fn(direction,node)
            enddo
         enddo
      endif
      
      if(compute_stiffness_matrix) then
         call stiffness_matrix(amatrx,coords,u,nprops,nnode,nsvars,props
     &   ,svars,ndofel)
      endif

      if(compute_mass_matrix) then
c         call mass_matrix(amatrx,coords,u,nprops,nnode,nsvars,props,svars,ndofel)
      endif

      





      return
      end
      

c     Include all other fortran files

c      include "tmp.f"
      include "./residual_vector.f"
      include "./stiffness_matrix.f"
      include "./mass_matrix.f"
      include "./INT1_coord_transform.f"
      include "./INT1_interface_dtdu.f"
      include "./INT1_interface_t.f"
      include "./INT1_transform_element_to_global.f"
      include "./INT1_transform_global_to_element.f"
c      include "./_adapt_integrate_f.f"
c      include "./_adapt_integrate_s.f"
c      include "./adapt_integrate.f"
      include "./adapt_integrate_f.f"
      include "./adapt_integrate_s.f"
c      include "./adapt_integrate_vars.f"
      include "./adapt_query.f"
      include "./adapt_query_bilinear.f"
      include "./adapt_query_chaboche.f"
      include "./chaboche_lambda_region.f"
      include "./compute_F_lambda_chaboche.f"
      include "./compute_dFdlambda_F_chaboche.f"
      include "./compute_F_lambda_bilinear.f"
      include "./compute_dFdlambda_F_bilinear.f"
      include "./compute_dlambda_F_dnu.f"
      include "./compute_lambda.f"
      include "./element_info.f"
c      include "./integrand_f.f"
c      include "./integrand_s.f"
c      include "./integrand_vars_f.f"
c      include "./integrand_vars_s.f"
c      include "./integrate_f.f"
c      include "./integrate_s.f"
c      include "./integration_points.f"
      include "./interface_dtdu_bilinear.f"
      include "./interface_dtdu_chaboche.f"
      include "./interface_t_bilinear.f"
      include "./interface_t_chaboche.f"
      include "./matrix_mult.f"
      include "./shape20h_centered.f"
      include "./user_element_f.f"
      include "./user_element_results.f"
      include "./user_e_c3d_f.f"
c     user_e_c3d.f will give "Warning: Stack frame size (144118640) larger than system limit (67108864)" unless "limit stacksize 141000" is given to the shell
      include "./user_e_c3d.f"

c      include "/tmp/GMLSPT_verification.f"
