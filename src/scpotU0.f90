!*******************************************************************
!   RHEED multi slice method :  subroutines for bulk & surf
!    v.3    90/4/12   v.4:2014/4    T.Hanada
!   contains asfcrr asfcef scpot
!*******************************************************************
! # of gaussians for real atomic scattering factor : nabr
!**********************************************************
!       energy correction of atomic scattering factor
!**********************************************************
        subroutine asfcrr(be,wn,nelm,nabr,a,b,c,dth,dtk,dtz,elm,tem,smesh)
        implicit none
        integer :: nelm,nabr
        real(8) :: a(nabr,nelm),b(nabr,nelm),c(nabr,nelm)
        real(8) :: dtz(nelm),dth(nelm),dtk(nelm),elm(nelm),be,wn,tem,smesh

        real(8) :: s,rel,rat,vib,dwz
        integer :: i,j
        real(8), parameter :: pi=atan(1d0)*4d0, c2m=511.001d0, ek=.262466d0

        wn=sqrt(1d3*be*ek*(1d0+0.5d0*be/c2m))
          s=4d0*pi/smesh
          rel=(1d0+be/c2m)*s
          rat=sqrt(100d0/be*(1d0+50d0/c2m)/(1d0+be*0.5d0/c2m))*s
          vib=11407d0*tem
        do i=1,nelm
! real part
          if (tem > 1d-5) then       ! Debye T
            dwz=vib/(dtz(i)*dtz(i)*elm(i))                 ! B
            dth(i)=sqrt((vib+vib)/elm(i))/(4d0*pi*dth(i))  ! sqrt(<u^2>)
            dtk(i)=sqrt((vib+vib)/elm(i))/(4d0*pi*dtk(i))  ! sqrt(<u^2>)
          else if (tem < -1d-5) then ! sqrt(<u^2>)
            dwz=(8d0*pi*pi)*dtz(i)*dtz(i)                  ! B
          else                       ! B = 8 pi^2 <u^2>
            dwz=dtz(i)                         ! B
            dth(i)=sqrt(dth(i)*0.5d0)/(pi+pi)  ! sqrt(<u^2>)
            dtk(i)=sqrt(dtk(i)*0.5d0)/(pi+pi)  ! sqrt(<u^2>)
          endif
! dth(i) & dtk(i) will be used only in the subroutine asfcef.

          do j=1,nabr                   ! a <= 0 for positron
            c(j,i)=b(j,i)/(-16d0*pi*pi)
            a(j,i)=a(j,i)*sqrt(4d0*pi/(b(j,i)+dwz))*rel
            b(j,i)=-4d0*pi*pi/(b(j,i)+dwz)
          end do
        end do
        return
        end
!**********************************************************
!       a.s.f. coefficient
!**********************************************************
        subroutine asfcef(nv,nelm,igh,igk,ghx,ghy,gky &
                         ,nabr,a,c,asf,dth,dtk,Ghk)
        implicit none
        integer :: igh(nv),igk(nv),nv,nelm,nabr, i,j,k
        real(8) :: dth(nelm),dtk(nelm),ghx,ghy,gky,Ghk(nv)
        real(8) :: a(nabr,nelm),c(nabr,nelm),asf(nv,nabr,nelm)
! let exp(uf) -> 0
        real(8), parameter :: uf=-70d0
        real(8) :: xh,yh,yk,s,xj,yj,sj,ex
!
        do k=1,nv
          xh=ghx*igh(k)
          yh=ghy*igh(k)
          yk=gky*igk(k)
          s=xh*xh+(yh+yk)*(yh+yk)
          Ghk(k)=sqrt(s)                 ! for ion
          do j=1,nelm
            xj=dth(j)*xh
            yj=dth(j)*yh+dtk(j)*yk
            sj=-0.5*(xj*xj+yj*yj)
! real part
            do i=1,nabr                  ! a <= 0 for positron
              ex=c(i,j)*s+sj             ! Debye-Waller
              if (ex < uf) then
                asf(k,i,j)=0d0
              else
                asf(k,i,j)=a(i,j)*exp(ex)
              endif
            end do
          end do
        end do
        return
        end
!**********************************************************
!       scattering potential
!**********************************************************
! mv is necessary for inclusion of the bulk top layer
!     to evaluate the surface potential.
        subroutine scpot(iclr,mv,nv,ns,v,vi,natm,nelm,ielm,z,zo,dz,nabr &
              ,asf,b,sap,rmesh,nsg,ma,mb,na,nb,gh,gk,Ghk,ocr,x,y,dx,dy)
        implicit none
        integer :: iclr,mv,nv,ns,natm,nelm,nabr,nsg,ma,mb,na,nb
        complex(8) :: v(mv,ns),vi(mv,ns)
        real(8) :: asf(mv,nabr,nelm),b(nabr,nelm),sap(nelm),gh(nv),gk(nv),Ghk(nv)
        real(8) :: ocr(natm),x(natm),y(natm),z(natm),zo,dz,dx,dy,rmesh
        integer :: ielm(natm)

        complex(8) :: st(nv)
        integer :: i,j,l,k,ie
! let exp(uf) -> 0
        real(8), parameter :: uf=-70d0
        real(8) :: ocrmesh,zi,z1,z2,ex,aex
!
        if (iclr /= 0) then
           v(1:mv,1:ns)=(0d0,0d0)
          vi(1:mv,1:ns)=(0d0,0d0)
        endif
!
      do j=1,natm
        call strfac(nsg,nv,ma,mb,na,nb,gh,gk,x(j),y(j),dx,dy,st)
        ocrmesh=rmesh*ocr(j)
        ie=ielm(j)
        zi=dz*0.5d0+zo
        do l=1,ns
          z1=zi-z(j)
          z2=z1*z1
! real part : v
          do i=1,nabr
            ex=b(i,ie)*z2
            if (ex > uf) then
              ex=exp(ex)*ocrmesh
              do k=1,nv
                aex=asf(k,i,ie)*ex
                v(k,l)=v(k,l)+dcmplx(aex)*st(k)
! imaginary part : i * | sap(ie) * (real part)|
                                           ! asf <= 0 for positron
                  vi(k,l)=vi(k,l)+dcmplx(0d0,abs(sap(ie)*aex))*st(k)
              end do
            endif
          end do
          zi=zi+dz
        end do ! l=1,ns
      end do ! j=1,natm
      return
      end
