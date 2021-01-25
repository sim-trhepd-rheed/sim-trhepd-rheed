!*******************************************************************
!   rheed multi slice method : bulk
!   subroutine scanrg,blkibg,blkghk,asfcrr,asfcef,scpot
!   U0:2019/5
!    T.Hanada
!*******************************************************************
! max # of potential fourier components : nvm
! max # of unit layers : mlm
!*******************************************************************
        subroutine bulkioU(inegpos,nubulk)
        implicit none
        integer, parameter :: mlm=1000, mab=4
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:,:,:),allocatable :: asf
        real(8), dimension(:),allocatable :: rdom,gh,gk,Ghk
        real(8), dimension(:),allocatable :: dth,dtk,dtz,elm,sap
        real(8), dimension(:,:),allocatable :: a,b,c
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: ielm,nb,iord,ih,ik,nbg
        integer, dimension(:),allocatable :: igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: nbgm,nvm
        real(8) :: be,wn,azi,azf,daz,gi,gf,dg,tem,dz,epsb,aa,bb,gam,cc,dx,dy
        real(8) :: ghx,ghy,gky,da1, amass,zmin,zmax,zoffset
        integer :: inegpos,nh,nk,ndom,ml,nelm,nsg,mp,natm,naz,ng,ns,ngr,nv
        integer :: nbt,ibt,nabr,iabr,iabp,iabe
        integer :: idom,iz, i,j,k
        real(8), parameter :: eps=1d-20, pi=atan(1d0)*4d0, rad=pi/180d0

        integer :: nubulk
        real(8) :: zout,suma,U0,gh0,gk0
        complex(8) :: st
        real(8), parameter :: eV=-1.0545888d-34*1.0545888d-34*1d20 &
                                /(2d0*9.1095345d-31*1.6021892d-19)
!----------input(3)----------
! beam parameters
        read (3,*) nh,nk,ndom
        if (nh < 1 .or. nk < 1) return
        if (ndom < 1) return
allocate (nb(ndom)); allocate (rdom(ndom))
        read (3,*) (nb(i),i=1,ndom)
        nbt=0
        do i=1,ndom
          if (nb(i) < 1) return
          nbt=nbt+nb(i)
        end do
        read (3,*) (rdom(i),i=1,ndom)
allocate (ih(nbt)); allocate (ik(nbt))
        ibt=0
        do i=1,ndom
          read (3,*) (ih(j),ik(j),j=ibt+1,ibt+nb(i))
          ibt=ibt+nb(i)
        end do
!------ first domain only
        ndom=1; nbt=nb(1)
!
        read (3,*) be,azi,azf,daz,gi,gf,dg
        read (3,*) dz,ml
        tem=0d0; epsb=1d-10
!        if (epsb < eps) epsb=eps

! atomic parameters
        iabr=0; iabp=1; iabe=0
        read (3,*) nelm
        if (nelm < 1) return
allocate (dth(nelm)); allocate (dtk(nelm)); allocate (dtz(nelm))
allocate (elm(nelm)); allocate (sap(nelm))
allocate (a(mab,nelm)); allocate (b(mab,nelm)); allocate (c(mab,nelm))

        do i=1,nelm
          read (3,*) iz,da1,sap(i)
          elm(i)=amass(iz)

          call asfparam(iz,iabr,nabr,a(1,i),b(1,i))
            if (nabr <= 0) return
          a(1,i)=a(1,i)-da1                           ! electron
          if (inegpos > 0) a(1:nabr,i)=-a(1:nabr,i)   ! positron
          if (nabr > mab) return

          read (3,*) dth(i),dtk(i),dtz(i)
          dth(i)=abs(dth(i)); dtk(i)=abs(dtk(i)); dtz(i)=abs(dtz(i))
       end do

! structural parameters
        read (3,*) nsg,aa,bb,gam,cc,dx,dy
        mp=0
        if (ndom > 1) then
          if (abs(aa-bb) > 1d-4) return
          if (abs(gam-90d0) < 1d-4) then
            do i=1,ndom
              if (abs(rdom(i)/90d0 - nint(rdom(i)/90d0)) > 1d-4) return
            end do
          else if (abs(gam-120d0) < 1d-4) then
            do i=1,ndom
              if (abs(rdom(i)/60d0 - nint(rdom(i)/60d0)) > 1d-4) return
            end do
          else
            return
          endif
        endif

! atomic structural parameters
        read (3,*) natm
        if (natm < 1) return
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        zmin=1d10; zmax=-1d10
        do i=1,natm
          read (3,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (ielm(i) < 1 .or. ielm(i) > nelm) return
          if (z(i) < zmin) zmin=z(i)
          if (z(i) > zmax) zmax=z(i)
        end do
          zoffset=0.5d0*(cc-zmax-zmin)
          do i=1,natm
            z(i)=z(i)+zoffset
          end do
!----------scan range----------
        gam=gam*rad
          ns=int(cc/dz)+1
          dz=cc/ns
        if (ns < 2) return
!----------mean inner potential in bulk (eV)----------
        gh0=0d0; gk0=0d0
        U0=0d0
        do i=1,natm
          suma=0d0
          do j=1,nabr
            suma=suma+a(j,ielm(i))
          end do
          call strfac(nsg,1,1,0,0,1,gh0,gk0,x(i),y(i),0d0,0d0,st)
! real(st) is number of symmetrically equivalent atoms 
          U0=U0+suma*ocr(i)*dble(st)
        end do
        U0=U0/(aa*bb*sin(gam)*cc)*4d0*pi*eV
        write (*,'(A,F10.5)') ' mean inner potential in bulk (eV) =',U0
        write (4,'(A,F10.5)') '# mean inner potential in bulk (eV) =',U0
        write (4,'(A)') '#z:Angstrom,real(U00):eV,imag(U00):eV'
!----------output(1)----------
        write (1) ndom
        write (1) inegpos,nh,nk,iabr,iabp,iabe,nabr
        write (1) (nb(i),rdom(i),i=1,ndom)
        write (1) (ih(i),ik(i),i=1,nbt)
        write (1) be,azi,daz,naz,gi,dg,ng
        write (1) tem,dz,ml,epsb
        write (1) nsg,mp,aa,bb,gam,cc,dx,dy,zoffset
        write (1) nelm,natm

        do i=1,nelm
          write (1) a(1:nabr,i),b(1:nabr,i)
          write (1) dth(i),dtk(i),dtz(i),elm(i),sap(i)
        end do
        do i=1,natm
          write (1) ielm(i),ocr(i),x(i),y(i),z(i)
        end do
!-----------energy correction of atomic scattering factor--------
        call asfcrr(be,wn,nelm,nabr,a,b,c,dth,dtk,dtz,elm &
                   ,tem,aa*bb*sin(gam))
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!----------domain----------
      ibt=1;      idom=1
!----------groups of interacting beams in bulk layer----------
! order of ih,ik will be sorted
allocate (iord(nb(idom))); allocate (nbg(nb(idom)))
        call blkibg(nb(idom),ih(ibt),ik(ibt),nh,nk, iord,ngr,nbg)
        nbgm=0; nvm=1; j=0
        do i=1,ngr
          if (nbg(i) > nbgm) nbgm=nbg(i)
          nvm=nvm+nbg(i)*(nbg(i)-1)/2
          j=j+nbg(i)
        end do
        if (j /= nb(idom)) then
          write (*,*) 'grouping error ',j,nb(idom)
          return
        endif
!-----------scattering vector & atomic scattering factor--------
allocate (iv(nbgm,nb(idom)))
allocate (igh(nvm)); allocate (igk(nvm))
        call blkghk(ngr,nbg,nvm,nv,igh,igk,nbgm,nb(idom),ih(ibt),ik(ibt),iv)
allocate (asf(nv,nabr,nelm)); allocate (Ghk(nv))
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,nabr,a,c,asf,dth,dtk,Ghk)
!----------output(1)----------
        write (1) (iord(i),i=1,nb(idom))
        write (1) ngr,nv
        write (1) (nbg(i),i=1,ngr)
        write (1) (igh(i),igk(i),i=1,nv)
!-----------elements of scattering matrix : v/vi----------
allocate (gh(nv)); allocate (gk(nv))
        do i=1,nv
          gh(i)=(-(pi+pi)/nh)*dble(igh(i))
          gk(i)=(-(pi+pi)/nk)*dble(igk(i))
        end do
allocate (v(nv,ns)); allocate (vi(nv,ns))
! bulk unit (in the 'ns' slices)
        call scpot(1,nv,nv,ns,v,vi,natm,nelm,ielm,z,0d0,dz &
             ,nabr,asf,b,sap,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,0d0,0d0)
! lower unit (below the 'ns' slices)
        call scpot(0,nv,nv,ns,v,vi,natm,nelm,ielm,z,cc,dz &
             ,nabr,asf,b,sap,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,-dx,-dy)
! upper unit (above the 'ns' slices)
        call scpot(0,nv,nv,ns,v,vi,natm,nelm,ielm,z,-cc,dz &
             ,nabr,asf,b,sap,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,dx,dy)
!-----------output U00----------
        zout=0.5d0*dz-zoffset-cc-nubulk*cc
! one bulk unit is added below surface layer
        do j=1,nubulk
        do i=1,ns
          write (4,'(E12.4,2(",",E12.4))') zout,dble(v(1,i))*eV,imag(vi(1,i))*eV
          zout=zout+dz
        end do
        end do
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (Ghk); deallocate (asf)
deallocate (igk); deallocate (igh)
deallocate (iv); deallocate (nbg); deallocate (iord)
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (c); deallocate (b); deallocate (a)
deallocate (sap); deallocate (elm)
deallocate (dtz); deallocate (dtk); deallocate (dth)
deallocate (ik); deallocate (ih)
deallocate (rdom); deallocate (nb)
      return
      end
