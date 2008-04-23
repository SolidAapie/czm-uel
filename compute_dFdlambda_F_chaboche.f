      subroutine compute_dFdlambda_F_chaboche(lambda_F,dFdlambda_F)

      implicit none
      real*8 lambda_F,dFdlambda_F

      if(lambda_F.lt.1) then
         dFdlambda_F = 2*(lambda_F-1)
      else
         dFdlambda_F = 0
      endif
      if(lambda_F.lt.0) then
         write(*,*) "ERROR in compute_dFdlambda_F_chaboche:"
         write(*,*) "  lambda_F = ",lambda_F
         stop
      endif

      
      return
      end
