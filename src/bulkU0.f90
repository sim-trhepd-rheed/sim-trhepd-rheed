!!!*******************************************************************
!   rheed multi slice method : bulk
!   subroutine scanrg,blkibg,blkghk,asfcrr,asfcef,scpot
!   U0:2019/5  U1:2021/8 modifications from bulkio.f90 are indicated by !!!-lines
!    T.Hanada
!*******************************************************************
subroutine bulkioU(inegpos,nubulk)
        implicit none
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:),allocatable :: rdom,gh,gk
        real(8), dimension(:),allocatable :: Bh,Bk,Bz,da1,sap
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: iz,ielm,nb,jorg,ih,ik,nbg
        integer, dimension(:),allocatable :: igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: nbgm,nvm
        real(8) :: be,wn,azi,azf,daz,gi,gf,dg,dz,epsb,aa,bb,gam,cc,dx,dy
        real(8) :: ghx,ghy,gky
        integer :: inegpos,nh,nk,ndom,ml,nelm,nsg,natm,naz,ng,ns,ngr,nv
        integer :: nbt,ibt,idiag
        integer :: iprn, idom, i,j,k
        real(8), parameter :: pi=atan(1d0)*4d0, rad=pi/180d0
!!!
        integer :: nubulk
        real(8) :: zout,U0,gh0,gk0,asf0
        complex(8) :: st
        real(8), parameter :: eV=-1.0545888d-34*1.0545888d-34*1d20 &
                                /(2d0*9.1095345d-31*1.6021892d-19)
!!!
!----------input(3)----------
! beam parameters
        read (3,*) nh,nk,ndom
        if (nh < 1 .or. nk < 1) then
          write (*,*) ' bulkU0.f90 line 36: bad nh,nk ',nh,nk
          stop
        endif
        if (ndom < 1)  then
          write (*,*) ' bulkU0.f90 line 36: bad ndom ',ndom
          stop
        endif
allocate (nb(ndom)); allocate (rdom(ndom))
        read (3,*) (nb(i),i=1,ndom)
        nbt=0
        do i=1,ndom
          if (nb(i) < 1) then
            write (*,*) ' bulkU0.f90 line 46: bad nb ',nb(i)
            stop
          endif
          nbt=nbt+nb(i)
        end do
        read (3,*) (rdom(i),i=1,ndom)
allocate (ih(nbt)); allocate (ik(nbt))
        ibt=0
        do i=1,ndom
          read (3,*) (ih(j),ik(j),j=ibt+1,ibt+nb(i))
          ibt=ibt+nb(i)
        end do
        read (3,*) be,azi,azf,daz,gi,gf,dg
        read (3,*) dz,ml
        epsb=1d-10

! atomic parameters
        read (3,*) nelm
        if (nelm < 1) then
          write (*,*) ' bulkU0.f90 line 67: bad nelm ',nelm
          stop
        endif
allocate (iz(nelm)); allocate (da1(nelm)); allocate (sap(nelm))
allocate (Bh(nelm)); allocate (Bk(nelm)); allocate (Bz(nelm)); 

        do i=1,nelm
          read (3,*) iz(i),da1(i),sap(i)
          if (iz(i) == 0 .or. iabs(iz(i)) > 98) then
            write (*,*) ' bulkU0.f90 line 76: max. atomic number iz is 98.',iz(i)
            stop
          endif
          read (3,*) Bh(i),Bk(i),Bz(i)
       end do

! structural parameters
        read (3,*) nsg,aa,bb,gam,cc,dx,dy
!!!        if (ndom > 1) then
!!! deleted
!!!        endif

! atomic structural parameters
        read (3,*) natm
        if (natm < 1) then
          write (*,*) ' bulkU0.f90 line 87: bad natm ',natm
          stop
        endif
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        do i=1,natm
          read (3,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (ielm(i) < 1 .or.  ielm(i) > nelm) then
            write (*,*) ' bulkU0.f90 line 99: bad ielm ',ielm(i)
            stop
          endif
        end do
!----------scan range----------
!!!        call scanrg(azi,azf,daz,naz,gi,gf,dg,ng)
        rdom(1:ndom)=rdom(1:ndom)*rad
        gam=gam*rad
        ns=int(cc/dz)+1
        if (ns < 2) then
          write (*,*) ' bulkio.f90 line 108: too small ns ',ns
          stop
        endif
        dz=cc/ns
!----------output(1)----------
        write (1) ndom
        write (1) inegpos,nh,nk,idiag
        write (1) (nb(i),rdom(i),i=1,ndom)
        write (1) (ih(i),ik(i),i=1,nbt)
        write (1) be,azi,daz,naz,gi,dg,ng
        write (1) dz,ml,epsb
        write (1) nsg,aa,bb,gam,cc,dx,dy
        write (1) nelm,natm

        do i=1,nelm
          write (1) iz(i),da1(i),sap(i)
          write (1) Bh(i),Bk(i),Bz(i)
        end do
        do i=1,natm
          write (1) ielm(i),ocr(i),x(i),y(i),z(i)
        end do
!-----------atomic scattering factor--------
! domain independent
        call asfcrr(inegpos,be,wn,nelm,iz,da1,Bz,aa*bb*sin(gam))
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!----------domain----------
      ibt=1;      idom=1
!!!      ibt=1
!!!      do idom=1,ndom
!----------groups of interacting beams in bulk layer----------
! order of ih,ik will be sorted
allocate (jorg(nb(idom))); allocate (nbg(nb(idom)))
        call blkibg(nb(idom),nh,nk,ih(ibt),ik(ibt), jorg,ngr,nbg)
        nbgm=maxval(nbg(1:ngr))
        nvm=1; j=0
        do i=1,ngr
          nvm=nvm+nbg(i)*(nbg(i)-1)/2
          j=j+nbg(i)
        end do
        if (j /= nb(idom)) then
          write (*,*) 'grouping error ',j,nb(idom)
          stop
        endif
!-----------scattering vector & atomic scattering factor--------
allocate (iv(nbgm,nb(idom)))
allocate (igh(nvm)); allocate (igk(nvm))
        call blkghk(ngr,nbg,nbgm,nb(idom),ih(ibt),ik(ibt),jorg,nvm, nv,igh,igk,iv)
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,Bh,Bk)
!----------output(1)----------
        write (1) jorg(1:nb(idom))
        write (1) ngr,nv
        write (1) nbg(1:ngr)
        write (1) igh(1:nv)
        write (1) igk(1:nv)
!-----------elements of scattering matrix : v/vi----------
allocate (gh(nv)); allocate (gk(nv))
        do i=1,nv
          gh(i)=(-(pi+pi)/nh)*dble(igh(i))
          gk(i)=(-(pi+pi)/nk)*dble(igk(i))
        end do
allocate (v(nv,ns)); allocate (vi(nv,ns))
! bulk unit (in the 'ns' slices)
        call scpot(1,nv,nv,ns,v,vi,natm,nelm,ielm,z,0d0,dz &
             ,sap,nsg,1,0,0,1,gh,gk,ocr,x,y,0d0,0d0)
! lower unit (below the 'ns' slices)
        call scpot(0,nv,nv,ns,v,vi,natm,nelm,ielm,z,cc,dz &
             ,sap,nsg,1,0,0,1,gh,gk,ocr,x,y,-dx,-dy)
! upper unit (above the 'ns' slices)
        call scpot(-1,nv,nv,ns,v,vi,natm,nelm,ielm,z,-cc,dz &
             ,sap,nsg,1,0,0,1,gh,gk,ocr,x,y,dx,dy)
!!!----------mean inner potential in bulk (eV)----------
        gh0=0d0; gk0=0d0
        U0=0d0
        do i=1,natm
          call strfac(nsg,1,1,0,0,1,gh0,gk0,x(i),y(i),0d0,0d0,st)
! real(st) is number of symmetrically equivalent atoms 
          U0=U0+ocr(i)*dble(st)*asf0(iz(ielm(i)),da1(ielm(i)))
        end do
        U0=U0/(aa*bb*sin(gam)*cc)*4d0*pi*eV
        if (inegpos > 0) U0=-U0     ! positron
        write (*,'(A,F0.5)') ' mean inner potential in bulk (eV) = ',U0
        write (4,'(A,F0.5)') '# mean inner potential in bulk (eV) = ',U0
        write (4,'(A)') '#z:Angstrom,real(U00):eV,imag(U00):eV'
!!!
!!!-----------output U00----------
        zout=0.5d0*dz-cc-nubulk*cc
! one bulk unit is added below surface layer
        do j=1,nubulk
        do i=1,ns
          write (4,'(ES12.4,2(",",ES12.4))') zout,dble(v(1,i))*eV,imag(vi(1,i))*eV
          zout=zout+dz
        end do
        end do
!!!
!-----------incident beam rocking-----------
!!!       call blkref(nv,nbgm,nb(idom),ns,v,vi,iv,dz,epsb,ngr,nbg &
!!!                  ,ih(ibt),ik(ibt),jorg,nh,nk,ml,dx,dy,wn,ghx,ghy,gky &
!!!                  ,azi+rdom(idom),daz,naz,gi,dg,ng,idiag,iprn)
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (igk); deallocate (igh)
deallocate (iv); deallocate (nbg); deallocate (jorg)
!!!        ibt=ibt+nb(idom)
!!!      end do ! idom=1,ndom
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (Bz); deallocate (Bk); deallocate (Bh)
deallocate (sap); deallocate (da1); deallocate (iz)
deallocate (ik); deallocate (ih)
deallocate (rdom); deallocate (nb)

end subroutine bulkioU
!**********************************************************
!       atomic scattering factor (s=0)
!**********************************************************
double precision function asf0(iz,da1)
        implicit none
        integer, parameter :: nab=4
        integer :: inegpos,iz
        real(8) :: da1,be
        real(8) :: a(nab),b(nab) !,rel
        integer :: j
!        real(8), parameter :: c2m=511.001d0

!        rel=1d0+be/c2m
        call asfparam(iabs(iz),a,b)  ! electron
        if (a(1) < 0d0) then
          write (*,*) ' |iz| must be less than 99 !'
          stop
        endif
        if (iz > 0) then
          a(1)=a(1)-da1
        else
          a(1)=a(1)*abs(da1)
        endif

        asf0=sum(a(1:nab))
end function asf0
