      subroutine compute_F_lambda_chaboche(lambda,F)

      implicit none
      real*8 lambda,F
      
      if(lambda.lt.1) then
         F = (1-lambda)**2
      else
         F = 0
      endif
      if(lambda.lt.0) then
         write(*,*) "ERROR in compute_F_lambda_chaboche:"
         write(*,*) "  lambda = ",lambda
         stop
      endif

      return
      end
