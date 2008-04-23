#define insert_integrand    include "integrand_f.f"
#define dadmul              adapt_integrate_f
#define integrand_vars1      interface_z,area,thickness_coord,nope,xl,u_node_interface
#define integrand_vars2      u_max_node_interface,elcon,nmat,ncmat_,ntmat_,integrand
#define declare_integrand_vars    include "integrand_vars_f.f"
#include  "adapt_integrate.f"


c$$$      subroutine adapt_integrate_f(a,b_adapt,minpts,maxpts,eps,wk,iwk,
c$$$     &        result,relerr,nfnevl,ifail,nfuncs_,et,thickness,
c$$$     &        thickness_coord,nope,xl,u_node,u_max_node,elcon,nmat,
c$$$     &        ncmat_,ntmat_,integrand)
c$$$
c$$$      implicit none
c$$$
c$$$      integer ifail,minpts,maxpts,iwk,nfnevl,nfuncs_
c$$$      real*8 a(*),b_adapt(*),wk(nfuncs_,*)
c$$$      real*8 eps,relerr,result(nfuncs_)
c$$$
c$$$      integer thickness_coord,nope,nmat,ncmat_,ntmat_
c$$$      real*8 et,thickness,elcon(0:ncmat_,ntmat_,*),
c$$$     &     xl(3,20),u_node(3,20),u_max_node(3,20),integrand(20)
c$$$      integer nfuncs_count,l,m
c$$$      real*8 xi,ze
c$$$
c$$$c     Prep
c$$$      do nfuncs_count=1,nfuncs_
c$$$         result(nfuncs_count) = 0
c$$$      enddo
c$$$      maxpts = 100
c$$$c     loop over fixed points
c$$$      do l=1,maxpts
c$$$         do m=1,maxpts
c$$$            xi = 2d0*(l-1)/maxpts - 1 + 1d0/maxpts
c$$$            ze = 2d0*(m-1)/maxpts - 1 + 1d0/maxpts
c$$$            call user_element_f(xi,et,ze,thickness,thickness_coord,
c$$$     &           nope,xl,u_node,u_max_node,elcon,nmat,
c$$$     &           ncmat_,ntmat_,integrand)
c$$$            do nfuncs_count=1,nfuncs_
c$$$               result(nfuncs_count) = result(nfuncs_count) + 
c$$$     &              4d0*integrand(nfuncs_count)/(maxpts*maxpts)
c$$$c               write(*,*) "integrand(",l,",",m,",",nfuncs_count,") = ",
c$$$c     &              integrand(nfuncs_count)
c$$$            enddo
c$$$c            write(*,*)"xi,ze = ",xi,ze,"    l,m = ",l,m,
c$$$c     &           "   f = ",result(1)
c$$$c            write(*,*)"integrand(1) = ",integrand(1),"  result(1) = ",
c$$$c     &           result(1)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      write(*,*)"non-adaptive: "
c$$$      do nfuncs_count=1,nfuncs_
c$$$         write(*,*) "result(",nfuncs_count,") = ",result(nfuncs_count)
c$$$      enddo
c$$$      return
c$$$      end
