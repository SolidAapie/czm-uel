# 1 "_adapt_integrate_f.f"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "_adapt_integrate_f.f"





# 1 "adapt_integrate.f" 1


c#define SAFE_ERROR 1
      subroutine adapt_integrate_f(
     & a,b_adapt,minpts,maxpts,eps,wk,iwk,result,
     & relerr,nfnevl,ifail,nfuncs_,
     & interface_z,area,thickness_coord,nope,xl,u_node_interface,
     & u_max_node_interface,properties,integrand)

      implicit none

      include "integrand_vars_f.f"

      logical ldv

      integer n,ifail,minpts,maxpts,iwk,nfnevl,nfuncs_
      real*8 a(*),b_adapt(*),wk(nfuncs_,*)
      real*8 ctr(15),wth(15),wthl(15),z(15)
      real*8 w(2:15,5),wp(2:15,3)

      real*8 f_return(nfuncs_)

      real*8 eps,relerr,relerr_safe,result(nfuncs_)
      integer ifncls,irgnst,irlcls,isbrgn,isbrgs,j_adapt,idvaxn,idvax0,
     & j1,k_adapt,l,m_adapt,isbtmp,isbtpp,nfuncs_count,mef
      real*8 r1,hf,xl2,xl4,xl5,w2,w4,wp2,wp4,abserr(nfuncs_),
     & abserr_safe(nfuncs_),twondm,
     & rgnvol,sum1(nfuncs_),sum2(nfuncs_),sum3(nfuncs_),
     & sum4(nfuncs_),sum5(nfuncs_),dif,difmax,f2(nfuncs_),
     & f3(nfuncs_),rgncmp(nfuncs_),rgnval(nfuncs_),rgnerr(nfuncs_)


      parameter (n = 2)
      parameter (r1 = 1, hf = r1/2)

      parameter (xl2 = 0.35856 85828 00318 073d0)
      parameter (xl4 = 0.94868 32980 50513 796d0)
      parameter (xl5 = 0.68824 72016 11685 289d0)

      parameter (w2 = 980*r1/6561, w4 = 200*r1/19683)
      parameter (wp2 = 245*r1/486, wp4 = 25*r1/729)

      data (w(j_adapt,1),w(j_adapt,3),j_adapt=2,15)
     1/-0.193872885230909911d+00, 0.518213686937966768d-01,
     2 -0.555606360818980835d+00, 0.314992633236803330d-01,
     3 -0.876695625666819078d+00, 0.111771579535639891d-01,
     4 -0.115714067977442459d+01, -0.914494741655235473d-02,
     5 -0.139694152314179743d+01, -0.294670527866686986d-01,
     6 -0.159609815576893754d+01, -0.497891581567850424d-01,
     7 -0.175461057765584494d+01, -0.701112635269013768d-01,
     8 -0.187247878880251983d+01, -0.904333688970177241d-01,
     9 -0.194970278920896201d+01, -0.110755474267134071d+00,
     a -0.198628257887517146d+01, -0.131077579637250419d+00,
     b -0.198221815780114818d+01, -0.151399685007366752d+00,
     c -0.193750952598689219d+01, -0.171721790377483099d+00,
     d -0.185215668343240347d+01, -0.192043895747599447d+00,
     e -0.172615963013768225d+01, -0.212366001117715794d+00/

      data (w(j_adapt,5),w(j_adapt+1,5),j_adapt=2,14,2)
     1/ 0.871183254585174982d-01, 0.435591627292587508d-01,
     2 0.217795813646293754d-01, 0.108897906823146873d-01,
     3 0.544489534115734364d-02, 0.272244767057867193d-02,
     4 0.136122383528933596d-02, 0.680611917644667955d-03,
     5 0.340305958822333977d-03, 0.170152979411166995d-03,
     6 0.850764897055834977d-04, 0.425382448527917472d-04,
     7 0.212691224263958736d-04, 0.106345612131979372d-04/

      data (wp(j_adapt,1),wp(j_adapt,3),j_adapt=2,15)
     1/-0.133196159122085045d+01, 0.445816186556927292d-01,
     2 -0.229218106995884763d+01, -0.240054869684499309d-01,
     3 -0.311522633744855959d+01, -0.925925925925925875d-01,
     4 -0.380109739368998611d+01, -0.161179698216735251d+00,
     5 -0.434979423868312742d+01, -0.229766803840877915d+00,
     6 -0.476131687242798352d+01, -0.298353909465020564d+00,
     7 -0.503566529492455417d+01, -0.366941015089163228d+00,
     8 -0.517283950617283939d+01, -0.435528120713305891d+00,
     9 -0.517283950617283939d+01, -0.504115226337448555d+00,
     a -0.503566529492455417d+01, -0.572702331961591218d+00,
     b -0.476131687242798352d+01, -0.641289437585733882d+00,
     c -0.434979423868312742d+01, -0.709876543209876532d+00,
     d -0.380109739368998611d+01, -0.778463648834019195d+00,
     e -0.311522633744855959d+01, -0.847050754458161859d+00/

c write(*,*) "in adapt_integrate.f"
c Setup
      do nfuncs_count=1,nfuncs_
         result(nfuncs_count)=0
         abserr(nfuncs_count)=0
         abserr_safe(nfuncs_count)=0
      enddo
      ifail=3
      if(minpts .gt. maxpts) return

      mef = 1 ! max error function; function with
! the current maximum relative error (controls subdivision)
      ifncls=0 !number of function evaluations
      ldv=.false.
      twondm=2**n !constant = 4 subregion volume,
! to be scaled later with wth
      irgnst=2*n+3 !constant = 7 abserr,result,idvax0,
!ctr,wth for each subregion
      irlcls=2**n+2*n*(n+1)+1 !constant = 17 points per subregion
      isbrgn=irgnst !index in wk of current subregion ?
      isbrgs=irgnst !last index in wk ?
      if(maxpts .lt. irlcls) return
      do j_adapt = 1,n
         ctr(j_adapt)=(b_adapt(j_adapt)+a(j_adapt))*hf !center of
! subregion to integrate
         wth(j_adapt)=(b_adapt(j_adapt)-a(j_adapt))*hf !half width
! of subregion to integrate
      enddo

c Begin adaptive integration loop
   20 rgnvol=twondm
c write(*,*) "mef = ",mef
      do j_adapt = 1,n
         rgnvol=rgnvol*wth(j_adapt)
         z(j_adapt)=ctr(j_adapt)
      enddo
c < Sample Integrand >
      include "integrand_f.f"
      do nfuncs_count=1,nfuncs_
         sum1(nfuncs_count)=f_return(nfuncs_count)
      enddo

      difmax=0
      do nfuncs_count=1,nfuncs_
         sum2(nfuncs_count)=0
         sum3(nfuncs_count)=0
      enddo
      do j_adapt = 1,n
         z(j_adapt)=ctr(j_adapt)-xl2*wth(j_adapt)
         include "integrand_f.f"
         do nfuncs_count=1,nfuncs_
            f2(nfuncs_count)=f_return(nfuncs_count)
         enddo
         z(j_adapt)=ctr(j_adapt)+xl2*wth(j_adapt)
         include "integrand_f.f"
         do nfuncs_count=1,nfuncs_
            f2(nfuncs_count)=f2(nfuncs_count)+f_return(nfuncs_count)
         enddo
         wthl(j_adapt)=xl4*wth(j_adapt)
         z(j_adapt)=ctr(j_adapt)-wthl(j_adapt)
         include "integrand_f.f"
         do nfuncs_count=1,nfuncs_
            f3(nfuncs_count)=f_return(nfuncs_count)
         enddo
         z(j_adapt)=ctr(j_adapt)+wthl(j_adapt)
         include "integrand_f.f"
         do nfuncs_count=1,nfuncs_
            f3(nfuncs_count)=f3(nfuncs_count)+f_return(nfuncs_count)
         enddo
         do nfuncs_count=1,nfuncs_
            sum2(nfuncs_count)=sum2(nfuncs_count)+f2(nfuncs_count)
            sum3(nfuncs_count)=sum3(nfuncs_count)+f3(nfuncs_count)
         enddo
         dif=abs(7*f2(mef)-f3(mef)-12*sum1(mef))
         difmax=max(dif,difmax)
         if(difmax .eq. dif) then
            idvaxn=j_adapt
         endif
         z(j_adapt)=ctr(j_adapt)
      enddo

      do nfuncs_count=1,nfuncs_
         sum4(nfuncs_count)=0
      enddo
      do j_adapt = 2,n
         j1=j_adapt-1
         do k_adapt = j_adapt,n
            do l = 1,2
               wthl(j1)=-wthl(j1)
               z(j1)=ctr(j1)+wthl(j1)
               do m_adapt = 1,2
                  wthl(k_adapt)=-wthl(k_adapt)
                  z(k_adapt)=ctr(k_adapt)+wthl(k_adapt)
                  include "integrand_f.f"
                  do nfuncs_count=1,nfuncs_
                     sum4(nfuncs_count)=sum4(nfuncs_count)+
     & f_return(nfuncs_count)
                  enddo
               enddo
            enddo
            z(k_adapt)=ctr(k_adapt)
         enddo
         z(j1)=ctr(j1)
      enddo

      do nfuncs_count=1,nfuncs_
         sum5(nfuncs_count)=0
      enddo
      do j_adapt = 1,n
         wthl(j_adapt)=-xl5*wth(j_adapt)
         z(j_adapt)=ctr(j_adapt)+wthl(j_adapt)
      enddo

   90 continue
      include "integrand_f.f"
      do nfuncs_count=1,nfuncs_
         sum5(nfuncs_count)=sum5(nfuncs_count)+f_return(nfuncs_count)
      enddo
      do j_adapt = 1,n
         wthl(j_adapt)=-wthl(j_adapt)
         z(j_adapt)=ctr(j_adapt)+wthl(j_adapt)
         if(wthl(j_adapt) .gt. 0) then
            go to 90
         endif
c continue
      enddo
c </ Sample Integrand >

      do nfuncs_count=1,nfuncs_
         rgncmp(nfuncs_count)=rgnvol*(wp(n,1)*sum1(nfuncs_count)+
     & wp2*sum2(nfuncs_count)+wp(n,3)*sum3(nfuncs_count)+
     & wp4*sum4(nfuncs_count))
         rgnval(nfuncs_count)=w(n,1)*sum1(nfuncs_count)+
     & w2*sum2(nfuncs_count)+w(n,3)*sum3(nfuncs_count)+
     & w4*sum4(nfuncs_count)+w(n,5)*sum5(nfuncs_count)
         rgnval(nfuncs_count)=rgnvol*rgnval(nfuncs_count)




         rgnerr(nfuncs_count)=rgnval(nfuncs_count)-
     & rgncmp(nfuncs_count)

         result(nfuncs_count)=result(nfuncs_count)+
     & rgnval(nfuncs_count)
         abserr(nfuncs_count)=abserr(nfuncs_count)+rgnerr(nfuncs_count)
         abserr_safe(nfuncs_count)=abserr_safe(nfuncs_count)+
     & abs(rgnval(nfuncs_count)-rgncmp(nfuncs_count))
      enddo
      ifncls=ifncls+irlcls

c Divide?
      if(ldv) then !integrate first half (or very first
! subregion)
 110
     & isbtmp=2*isbrgn
         if(isbtmp .gt. isbrgs) then !isbtmp is beyond previous max
! index
            go to 160
         endif
         if(isbtmp .lt. isbrgs) then !isbtmp is within already used
! memory
            isbtpp=isbtmp+irgnst
                                !isbtmp set to index of greater
! abserr of isbtmp and next abserr




            if(abs(wk(mef,isbtmp)) .lt. abs(wk(mef,isbtpp))) then
! abserr of isbtmp lt abserr of next abserr

               isbtmp=isbtpp
            endif
         endif




         if(abs(rgnerr(mef)) .ge. abs(wk(mef,isbtmp))) then

            go to 160
         endif
! copy isbtmp to isbrgn
         do k_adapt = 0,irgnst-1
            do nfuncs_count=1,nfuncs_
               wk(nfuncs_count,isbrgn-k_adapt)=
     & wk(nfuncs_count,isbtmp-k_adapt)
            enddo
         enddo
         isbrgn=isbtmp
         go to 110
      else !already integrated first half
  140 isbtmp=(isbrgn/(2*irgnst))*irgnst
         ! while...
         if(isbtmp .ge. irgnst .and.



     & abs(rgnerr(mef)) .gt. abs(wk(mef,isbtmp))) then

            do k_adapt = 0,irgnst-1
               do nfuncs_count=1,nfuncs_
                  wk(nfuncs_count,isbrgn-k_adapt)=
     & wk(nfuncs_count,isbtmp-k_adapt)
               enddo
            enddo
            isbrgn=isbtmp
            go to 140
         endif
      endif

  160 do nfuncs_count=1,nfuncs_



         wk(nfuncs_count,isbrgn)=rgnerr(nfuncs_count)

         wk(nfuncs_count,isbrgn-1)=rgnval(nfuncs_count)
         wk(nfuncs_count,isbrgn-2)=idvaxn
      enddo
      do j_adapt = 1,n
         isbtmp=isbrgn-2*j_adapt-2
         do nfuncs_count=1,nfuncs_
            wk(nfuncs_count,isbtmp+1)=ctr(j_adapt)
            wk(nfuncs_count,isbtmp)=wth(j_adapt)
         enddo
      enddo
      if(ldv) then !integrate other half of newly
! divided subregion
         ldv=.false.
         ctr(idvax0)=ctr(idvax0)+2*wth(idvax0)
         isbrgs=isbrgs+irgnst
         isbrgn=isbrgs
         go to 20
      endif
c pick the max relative error of all the integrated functions
      relerr = 0
      do nfuncs_count=1,nfuncs_
         if(abs(abserr(nfuncs_count)/result(nfuncs_count)).gt.
     & abs(relerr))then
            relerr = abs(abserr(nfuncs_count)/result(nfuncs_count))
            mef = nfuncs_count
         endif
         if(abs(abserr_safe(nfuncs_count)/result(nfuncs_count)).gt.
     & abs(relerr_safe))then
            relerr_safe = abs(abserr_safe(nfuncs_count)/
     & result(nfuncs_count))
         endif
      enddo
c Check for termination conditions

c IWK is too small for the specified number MAXPTS?
      if(isbrgs+irgnst .gt. iwk) then
         ifail=2
      endif
c MAXPTS is too small for the specified accuracy EPS?
      if(ifncls+2*irlcls .gt. maxpts) then
         ifail=1
      endif
c Normal exit, relerr < eps?
      if(((abs(relerr) .lt. eps).or.(abs(relerr_safe).lt.eps)) .and.
     & ifncls .ge. minpts) then
         ifail=0
      endif
      if(ifail .eq. 3) then

c No termination; Prep next cycle
!
         ldv=.true.
         isbrgn=irgnst
         do nfuncs_count=1,nfuncs_
            abserr(nfuncs_count)=abserr(nfuncs_count)-
     & wk(nfuncs_count,isbrgn)
            abserr_safe(nfuncs_count)=abserr_safe(nfuncs_count)-
     & abs(wk(nfuncs_count,isbrgn))
            result(nfuncs_count)=result(nfuncs_count)-
     & wk(nfuncs_count,isbrgn-1)
         enddo
         idvax0=wk(mef,isbrgn-2)
         do j_adapt = 1,n
            isbtmp=isbrgn-2*j_adapt-2
            ctr(j_adapt)=wk(mef,isbtmp+1)
            wth(j_adapt)=wk(mef,isbtmp)
         enddo
         wth(idvax0)=hf*wth(idvax0)
         ctr(idvax0)=ctr(idvax0)-wth(idvax0)
         go to 20
      endif
c Terminate

      write(*,*) "func evals   ",ifncls
      nfnevl=ifncls
      return
      end
# 7 "_adapt_integrate_f.f" 2


c$$$ subroutine adapt_integrate_f(a,b_adapt,minpts,maxpts,eps,wk,iwk,
c$$$ & result,relerr,nfnevl,ifail,nfuncs_,et,thickness,
c$$$ & thickness_coord,nope,xl,u_node,u_max_node,elcon,nmat,
c$$$ & ncmat_,ntmat_,integrand)
c$$$
c$$$ implicit none
c$$$
c$$$ integer ifail,minpts,maxpts,iwk,nfnevl,nfuncs_
c$$$ real*8 a(*),b_adapt(*),wk(nfuncs_,*)
c$$$ real*8 eps,relerr,result(nfuncs_)
c$$$
c$$$ integer thickness_coord,nope,nmat,ncmat_,ntmat_
c$$$ real*8 et,thickness,elcon(0:ncmat_,ntmat_,*),
c$$$ & xl(3,20),u_node(3,20),u_max_node(3,20),integrand(20)
c$$$ integer nfuncs_count,l,m
c$$$ real*8 xi,ze
c$$$
c$$$c Prep
c$$$ do nfuncs_count=1,nfuncs_
c$$$ result(nfuncs_count) = 0
c$$$ enddo
c$$$ maxpts = 100
c$$$c loop over fixed points
c$$$ do l=1,maxpts
c$$$ do m=1,maxpts
c$$$ xi = 2d0*(l-1)/maxpts - 1 + 1d0/maxpts
c$$$ ze = 2d0*(m-1)/maxpts - 1 + 1d0/maxpts
c$$$ call user_element_f(xi,et,ze,thickness,thickness_coord,
c$$$ & nope,xl,u_node,u_max_node,elcon,nmat,
c$$$ & ncmat_,ntmat_,integrand)
c$$$ do nfuncs_count=1,nfuncs_
c$$$ result(nfuncs_count) = result(nfuncs_count) +
c$$$ & 4d0*integrand(nfuncs_count)/(maxpts*maxpts)
c$$$c write(*,*) "integrand(",l,",",m,",",nfuncs_count,") = ",
c$$$c & integrand(nfuncs_count)
c$$$ enddo
c$$$c write(*,*)"xi,ze = ",xi,ze,"    l,m = ",l,m,
c$$$c & "   f = ",result(1)
c$$$c write(*,*)"integrand(1) = ",integrand(1),"  result(1) = ",
c$$$c & result(1)
c$$$ enddo
c$$$ enddo
c$$$
c$$$ write(*,*)"non-adaptive: "
c$$$ do nfuncs_count=1,nfuncs_
c$$$ write(*,*) "result(",nfuncs_count,") = ",result(nfuncs_count)
c$$$ enddo
c$$$ return
c$$$ end
