      subroutine compute_F_lambda_bilinear(lambda,lambda_init,F)

      implicit none
      real*8 lambda,lambda_init,F
      
      if(lambda.lt.lambda_init) then
         F = 1
      elseif(lambda.lt.1) then
         F = lambda_init/lambda*(1-lambda)/(1-lambda_init)
      else
         F = 0
      endif
      if(lambda.lt.0) then
         write(*,*) "ERROR in compute_F_lambda_bilinear:"
         write(*,*) "  lambda = ",lambda
         stop
      endif

      return
      end
