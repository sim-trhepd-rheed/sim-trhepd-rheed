!*******************************************************************
!   rheed multi slice method : surface
!   subroutine srfghk,asfcrr,asfcef,scpot
!   U0:2019/5
!    T.Hanada
!*******************************************************************
! max # of potential fourier components : nvm
!*******************************************************************
        subroutine surfioU
        implicit none
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:),allocatable :: rdom,wdom,gh,gk,Ghk
        real(8), dimension(:),allocatable :: dth,dtk,dtz,elm,sap
        real(8), dimension(:,:),allocatable :: a,b,c
        real(8), dimension(:,:,:),allocatable :: asf
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: ielm,nb,iord,ih,ik,nbg
        integer, dimension(:),allocatable :: igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: nvm, ndom
        real(8) :: be,wn,azi,daz,gi,dg,tem,dz,epsb,aa,bb,gam,cc,dxb,dyb,zoffset
        real(8) :: smesh,dthick,dxs,dys,ghx,ghy,gky,da1,zz,amass
        integer :: inegpos,nh,nk,naz,ng,ml,nelmb,nsgb,mp,natmb
        integer :: nelms,nelm,nsgs,msa,msb,nsa,nsb,natms,natm,ns, ngr,nvb,nv
        integer :: nbt,ibt,nabr,iabr,iabp,iabe
        integer :: iz,idom, i,j
        real(8), parameter :: pi=atan(1d0)*4d0, deg=180d0/pi

        real(8) :: zout
        real(8), parameter :: eV=-1.0545888d-34*1.0545888d-34*1d20 &
                                /(2d0*9.1095345d-31*1.6021892d-19)
!----------input(1)----------
! inegpos <=0: electron, >=1: positron
        read (1) ndom
        read (1) inegpos,nh,nk,iabr,iabp,iabe,nabr
allocate (nb(ndom)); allocate (rdom(ndom)); allocate (wdom(ndom))
        read (1) (nb(i),rdom(i),i=1,ndom)
        nbt=0
        do i=1,ndom
          nbt=nbt+nb(i)
        end do
allocate (ih(nbt)); allocate (ik(nbt))
        read (1) (ih(j),ik(j),j=1,nbt)
        read (1) be,azi,daz,naz,gi,dg,ng
        read (1) tem,dz,ml,epsb
! bulk atomic & structural parameters
        read (1) nsgb,mp,aa,bb,gam,cc,dxb,dyb,zoffset
        read (1) nelmb,natmb
!----------input(2)----------
! surface atomic parameters
        read (2,*) nelms
        nelm=nelmb+nelms

allocate (dth(nelm)); allocate (dtk(nelm)); allocate (dtz(nelm))
allocate (elm(nelm)); allocate (sap(nelm))
allocate (a(nabr,nelm)); allocate (b(nabr,nelm)); allocate (c(nabr,nelm))
        do i=1,nelmb
          read (1) a(1:nabr,i),b(1:nabr,i)
          read (1) dth(i),dtk(i),dtz(i),elm(i),sap(i)
        end do

        do i=nelmb+1,nelm
          read (2,*) iz,da1,sap(i)
          elm(i)=amass(iz)

          call asfparam(iz,iabr,nabr,a(1,i),b(1,i))
            if (nabr <= 0) return
          a(1,i)=a(1,i)-da1                           ! electron
          if (inegpos > 0) a(1:nabr,i)=-a(1:nabr,i)   ! positron

          read (2,*) dth(i),dtk(i),dtz(i)
          dth(i)=abs(dth(i)); dtk(i)=abs(dtk(i)); dtz(i)=abs(dtz(i))
        end do
! surface structural parameters
        read (2,*) nsgs,msa,msb,nsa,nsb,dthick,dxs,dys
        read (2,*) natms
        natm=natmb+natms
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        do i=1,natmb
          read (1) ielm(i),ocr(i),x(i),y(i),z(i) ! zoffset has been added to z(i)
        end do
        do i=natmb+1,natm
          read (2,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (ielm(i) >= -nelmb .and. ielm(i) <= -1) then
            ielm(i)=-ielm(i)
          else if (ielm(i) >= 1 .and. ielm(i) <= nelms) then
            ielm(i)=ielm(i)+nelmb
          else
            return
          endif
          z(i)=z(i)+zoffset
        end do

        if (ndom > 1) then
          read (2,*) (wdom(i),i=1,ndom)
        else
          wdom(1)=1d0
        endif
!----------surface slices----------
        if (nsgs < 1) then
          ns=0 ! to see bulk intensity
        else
          zz=0d0
          do i=natmb+1,natm
            if (z(i) > zz) zz=z(i)
          end do
          ns=int((zz+dthick+cc)/dz)+1
        endif
!-----------energy correction of atomic scattering factor--------
        smesh=abs(msa*nsb-msb*nsa)
        call asfcrr(be,wn,nelm,nabr,a,b,c,dth,dtk,dtz,elm &
                   ,tem,aa*bb*sin(gam)*smesh)
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!----------domain----------
      ibt=1; idom=1
!----------input(1)----------
! information of beams in bulk
allocate (iord(nb(idom)))
        read (1) (iord(i),i=1,nb(idom))
        read (1) ngr,nvb
          nvm=1+nb(idom)*(nb(idom)-1)/2
          if (nvb > nvm) return
allocate (nbg(ngr))
allocate (igh(nvm)); allocate (igk(nvm))
        read (1) (nbg(i),i=1,ngr)
        read (1) (igh(i),igk(i),i=1,nvb)
!-----------scattering vector & atomic scattering factor--------
allocate (iv(nb(idom),nb(idom)))
        call srfghk(nvm,nv,nvb,igh,igk,nb(idom),ih(ibt),ik(ibt),iv)
allocate (asf(nv,nabr,nelm)); allocate (Ghk(nv))
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,nabr,a,c,asf,dth,dtk,Ghk)
!-----------elements of scattering matrix : v/vi----------
allocate (gh(nv)); allocate (gk(nv))
        do i=1,nv
          gh(i)=(-(pi+pi)/nh)*dble(igh(i))
          gk(i)=(-(pi+pi)/nk)*dble(igk(i))
        end do
allocate (v(nv,ns)); allocate (vi(nv,ns))
! surface reconstructed layer
        call scpot(1,nv,nv,ns,v,vi,natms,nelm,ielm(natmb+1),z(natmb+1),-cc,dz &
             ,nabr,asf,b,sap,1d0,nsgs,msa,msb,nsa,nsb,gh,gk,Ghk &
             ,ocr(natmb+1),x(natmb+1),y(natmb+1),dxs+dxb,dys+dyb)
! topmost bulk unit layer (in the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm,z,0d0,dz &
          ,nabr,asf,b,sap,smesh,nsgb,1,0,0,1,gh,gk,Ghk,ocr,x,y,dxb,dyb)
! second bulk unit layer (below the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm,z,cc,dz &
          ,nabr,asf,b,sap,smesh,nsgb,1,0,0,1,gh,gk,Ghk,ocr,x,y,0d0,0d0)
!-----------output U00----------
      zout=0.5d0*dz-zoffset-cc ! one bulk unit is added below surface layer
      do i=1,ns
        write (4,'(E12.4,2(",",E12.4))') zout,dble(v(1,i))*eV,imag(vi(1,i))*eV
        zout=zout+dz
      end do
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (Ghk); deallocate (asf)
deallocate (iv); deallocate (igk); deallocate (igh)
deallocate (nbg); deallocate (iord)
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (c); deallocate (b); deallocate (a)
deallocate (sap); deallocate (elm)
deallocate (dtz); deallocate (dtk); deallocate (dth)
deallocate (ik); deallocate (ih)
deallocate (wdom); deallocate (rdom); deallocate (nb)
      return
      end
