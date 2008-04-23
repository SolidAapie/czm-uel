      subroutine chaboche_lambda_region(lambda,lambda_max,lambda_region)

      implicit none
      integer lambda_region
      real*8 lambda,lambda_max,lambda_init

      lambda_init = .1
      if(lambda.gt.1.or.lambda_max.gt.1) then
         lambda_region = 3
      elseif(lambda.ge.lambda_max) then
         lambda_region = 2
      elseif(lambda.ge.0) then
         lambda_region = 1
      else
         write(*,*) "ERROR in chaboche_lambda_region: lambda = ",lambda
         lambda_region = -1
      endif

c     Go easy on initial debonding (where nu(3) >,<,~= 0)
      if(lambda.lt.lambda_init.and.lambda_max.lt.lambda_init) then
         lambda_region = 0
      endif

      return
      end
