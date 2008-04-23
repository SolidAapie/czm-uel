      subroutine compute_dFdlambda_F_bilinear(lambda_F,lambda_init,
     &     dFdlambda_F)

      implicit none
      real*8 lambda_F,lambda_init,dFdlambda_F

      if(lambda_F.lt.1.and.lambda_F.gt.lambda_init) then
         dFdlambda_F = lambda_init/(lambda_init-1)/lambda_F**2
      else
         dFdlambda_F = 0
      endif
      if(lambda_F.lt.0) then
         write(*,*) "ERROR in compute_dFdlambda_F_bilinear:"
         write(*,*) "  lambda_F = ",lambda_F
         stop
      endif

      
      return
      end
