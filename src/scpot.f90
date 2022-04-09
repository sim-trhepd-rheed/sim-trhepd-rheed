!*******************************************************************
!   RHEED multi slice method :  subroutines for bulk & surf
!    v.3    90/4/12   v.4:2014/4   v.5:2022/1    T.Hanada
!   contains asfcrr asfcef scpot
!*******************************************************************
!**********************************************************
!       energy correction of atomic scattering factor
!       domain independent
!**********************************************************
subroutine asfcrr(inegpos,be,wn,nelm,iz,da1,nab,a,b,c,Bz,smesh)
        implicit none
        integer :: iz(nelm),inegpos,nelm,nab
        real(8) :: a(nab,nelm),b(nab,nelm),c(nab,nelm)
        real(8) :: da1(nelm),Bz(nelm),be,wn,smesh
        real(8) :: rel
        integer :: i,j
        real(8), parameter :: pi=atan(1d0)*4d0, c2m=511.001d0, ek=.262466d0

        wn=sqrt(1d3*be*ek*(1d0+0.5d0*be/c2m))
        rel=(1d0+be/c2m)*4d0*pi/smesh
! Bz = 8 pi^2 <uz^2>  ! a < 0 for positron
        do i=1,nelm
          call asfparam(iabs(iz(i)),a(1,i),b(1,i))  ! electron
          if (a(1,i) < 0d0) then
            write (*,*) ' |iz| must be less than 99 !'
            stop
          endif
          if (iz(i) > 0) then
            a(1,i)=a(1,i)-da1(i)
          else
            a(1,i)=a(1,i)*abs(da1(i))
          endif
          if (inegpos > 0) a(1:nab,i)=-a(1:nab,i)   ! positron

          do j=1,nab
            c(j,i)=b(j,i)/(16d0*pi*pi)
            a(j,i)=a(j,i)*sqrt(4d0*pi/(b(j,i)+Bz(i)))*rel
            b(j,i)=4d0*pi*pi/(b(j,i)+Bz(i))
          end do
        end do
end subroutine asfcrr
!**********************************************************
!       a.s.f. coefficient
!**********************************************************
subroutine asfcef(nv,nelm,igh,igk,ghx,ghy,gky &
                         ,nab,a,c,asf,Bh,Bk)
        implicit none
        integer :: igh(nv),igk(nv),nv,nelm,nab, i,j,k
        real(8) :: Bh(nelm),Bk(nelm),ghx,ghy,gky
        real(8) :: a(nab,nelm),c(nab,nelm),asf(nv,nab,nelm)
! let exp(-uf) -> 0
        real(8), parameter :: pi2=atan(1d0)*8d0, uf=70d0
        real(8) :: xh,yh,yk,s,xi,yi,si,ex
        real(8) :: uh(nelm),uk(nelm)

! Bh, Bk = 8 pi^2 <u^2>
        do i=1,nelm
          uh(i)=sqrt(Bh(i)*0.5d0)/pi2  ! sqrt(<u^2>)
          uk(i)=sqrt(Bk(i)*0.5d0)/pi2  ! sqrt(<u^2>)
        end do

        do k=1,nv
          xh=ghx*igh(k)
          yh=ghy*igh(k)
          yk=gky*igk(k)
          s=xh*xh+(yh+yk)*(yh+yk)
          do i=1,nelm
            xi=uh(i)*xh
            yi=uh(i)*yh+uk(i)*yk
            si=0.5*(xi*xi+yi*yi)        ! Debye-Waller
! real part
            do j=1,nab
              ex=c(j,i)*s+si
              if (ex > uf) then
                asf(k,j,i)=0d0
              else
                asf(k,j,i)=a(j,i)*exp(-ex) ! a < 0 for positron
              endif
            end do
          end do
        end do
end subroutine asfcef
!**********************************************************
!       scattering potential
!**********************************************************
! mv is necessary for inclusion of the bulk top layer
!     to evaluate the surface potential.
subroutine scpot(iclr,mv,nv,ns,v,vi,natm,nelm,ielm,z,zo,dz &
          ,nab,asf,b,sap,nsg,ma,mb,na,nb,gh,gk,ocr,x,y,dx,dy)
        implicit none
        integer :: iclr,mv,nv,ns,natm,nelm,nab,nsg,ma,mb,na,nb
        complex(8) :: v(mv,ns),vi(mv,ns)
        real(8) :: asf(mv,nab,nelm),b(nab,nelm),sap(nelm),gh(nv),gk(nv)
        real(8) :: ocr(natm),x(natm),y(natm),z(natm),zo,dz,dx,dy
        integer :: ielm(natm)

        complex(8) :: st(nv)
        integer :: i,j,l,k,ie
! let exp(uf) -> 0
        real(8), parameter :: uf=70d0
        real(8) :: rmesh,ocrmesh,zi,z1,z2,ex,aex
!
        if (iclr /= 0) then
           v(1:mv,1:ns)=(0d0,0d0)
          vi(1:mv,1:ns)=(0d0,0d0)
        endif
!
      rmesh=1d0/abs(ma*nb-mb*na)
      do j=1,natm
        call strfac(nsg,nv,ma,mb,na,nb,gh,gk,x(j),y(j),dx,dy,st)
        ocrmesh=rmesh*ocr(j)
        ie=ielm(j)
        zi=dz*0.5d0+zo
        do l=1,ns
          z1=zi-z(j)
          z2=z1*z1
! real: v
          do i=1,nab
            ex=b(i,ie)*z2
            if (ex < uf) then
              ex=exp(-ex)*ocrmesh
              do k=1,nv
! asf < 0 for positron, >0 for electron
                aex=asf(k,i,ie)*ex
                v(k,l)=v(k,l)+dcmplx(aex)*st(k)
! imaginary (absorption): vi
! positive |sap*asf| must be used for both positron and electron
                if (k == 1) then
                  vi(k,l)=vi(k,l)+dcmplx(0d0,abs(sap(ie)*aex))*st(k)
                else if (sap(ie) > 0d0) then
                  vi(k,l)=vi(k,l)+dcmplx(0d0,abs(sap(ie)*aex))*st(k)
! set sap(ie) < 0 if you want to use only 00(k=1) component of vi.
                endif
              end do
            endif
          end do
          zi=zi+dz
        end do ! l=1,ns
      end do ! j=1,natm
end subroutine scpot
