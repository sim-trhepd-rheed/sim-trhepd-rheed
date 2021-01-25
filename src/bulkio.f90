!*******************************************************************
!   rheed multi slice method : bulk
!   subroutine scanrg,blkibg,blkghk,asfcrr,asfcef,scpot,blkref
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2
!    T.Hanada
!*******************************************************************
! max # of potential fourier components : nvm
! max # of unit layers : mlm
!*******************************************************************
        subroutine bulkio(inegpos,idiag,iprn)
        implicit none
        integer, parameter :: mlm=1000, mab=6
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:,:,:),allocatable :: asf
        real(8), dimension(:),allocatable :: rdom,gh,gk,Ghk
        real(8), dimension(:),allocatable :: dth,dtk,dtz,elm,dchg
        real(8), dimension(:,:),allocatable :: a,b,c
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: ielm,nb,iord,ih,ik,nbg
        integer, dimension(:),allocatable :: igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: nbgm,nvm
        real(8) :: be,wn,azi,azf,daz,gi,gf,dg,tem,dz,epsb,aa,bb,gam,cc,dx,dy
        real(8) :: ghx,ghy,gky,da1,sap,sae, amass,zmin,zmax,zoffset
        integer :: inegpos,isym,nh,nk,ndom,ml,nelm,nsg,mp,natm,naz,ng,ns,ngr,nv
        integer :: nbt,ibt,nabr,nabi,iabr,iabp,iabe,idiag
        integer :: iprn, idom,iz,ichg, i,j,k
        real(8), parameter :: eps=1d-20, pi=atan(1d0)*4d0, rad=pi/180d0
!----------input(3)----------
! beam parameters
        isym=0
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
        read (3,*) be,azi,azf,daz,gi,gf,dg
        read (3,*) dz,ml
        tem=0d0; epsb=1d-10
        if (ml < 1 .or. ml > mlm) ml=mlm
!        if (epsb < eps) epsb=eps

! atomic parameters
        iabr=0; iabp=1; iabe=0
        read (3,*) nelm
        if (nelm < 1) return
allocate (dth(nelm)); allocate (dtk(nelm)); allocate (dtz(nelm))
allocate (elm(nelm)); allocate (dchg(nelm))
allocate (a(mab,nelm)); allocate (b(mab,nelm)); allocate (c(mab,nelm))

        do i=1,nelm
          read (3,*) iz,da1,sap;  ichg=0; sae=0d0
          dchg(i)=dble(ichg)
          elm(i)=amass(iz)

          call asfparam(iz,iabr,nabr,a(1,i),b(1,i))
            if (nabr <= 0) return
          a(1,i)=a(1,i)-da1                           ! electron
          if (inegpos > 0) a(1:nabr,i)=-a(1:nabr,i)   ! positron

          nabi=1 
          if (nabr+nabi > mab) return
          a(nabr+1,i)=sap; b(nabr+1,i)=0d0 ! V(imag)=a*V(real), b is dummy

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
        call scanrg(azi,azf,daz,naz,gi,gf,dg,ng)
        rdom(1:ndom)=rdom(1:ndom)*rad
        gam=gam*rad
          ns=int(cc/dz)+1
          dz=cc/ns
        if (ns < 2) return
!----------output(1)----------
        write (1) ndom
        write (1) inegpos,isym,nh,nk,iabr,iabp,iabe,nabr,nabi,idiag
        write (1) (nb(i),rdom(i),i=1,ndom)
        write (1) (ih(i),ik(i),i=1,nbt)
        write (1) be,azi,daz,naz,gi,dg,ng
        write (1) tem,dz,ml,epsb
        write (1) nsg,mp,aa,bb,gam,cc,dx,dy,zoffset
        write (1) nelm,natm

        do i=1,nelm
          write (1) a(1:nabr+nabi,i),b(1:nabr+nabi,i)
          write (1) dth(i),dtk(i),dtz(i),elm(i),dchg(i)
        end do
        do i=1,natm
          write (1) ielm(i),ocr(i),x(i),y(i),z(i)
        end do
!-----------energy correction of atomic scattering factor--------
        call asfcrr(be,wn,nelm,mab,nabr,nabi,a,b,c,dchg,dth,dtk,dtz,elm &
                   ,tem,aa*bb*sin(gam))
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!----------domain----------
      ibt=1
      do idom=1,ndom
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
allocate (asf(nv,mab,nelm)); allocate (Ghk(nv))
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,mab,nabr,nabi,a,c,asf,dth,dtk,Ghk)
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
             ,mab,nabr,nabi,asf,b,dchg,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,0d0,0d0)
! lower unit (below the 'ns' slices)
        call scpot(0,nv,nv,ns,v,vi,natm,nelm,ielm,z,cc,dz &
             ,mab,nabr,nabi,asf,b,dchg,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,-dx,-dy)
! upper unit (above the 'ns' slices)
        call scpot(0,nv,nv,ns,v,vi,natm,nelm,ielm,z,-cc,dz &
             ,mab,nabr,nabi,asf,b,dchg,1d0,nsg,1,0,0,1,gh,gk,Ghk,ocr,x,y,dx,dy)
!-----------incident beam rocking-----------
       call blkref(nv,nbgm,nb(idom),ns,v,vi,iv,dz,epsb,ngr,nbg,ih(ibt),ik(ibt) &
                  ,nh,nk,ml,dx,dy,wn,ghx,ghy,gky &
                  ,azi+rdom(idom),daz,naz,gi,dg,ng,idiag,iprn)
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (Ghk); deallocate (asf)
deallocate (igk); deallocate (igh)
deallocate (iv); deallocate (nbg); deallocate (iord)
        ibt=ibt+nb(idom)
      end do ! idom=1,ndom
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (c); deallocate (b); deallocate (a)
deallocate (dchg); deallocate (elm)
deallocate (dtz); deallocate (dtk); deallocate (dth)
deallocate (ik); deallocate (ih)
deallocate (rdom); deallocate (nb)
      return
      end
