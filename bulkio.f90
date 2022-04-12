!*******************************************************************
!   rheed multi slice method : bulk
!   subroutine scanrg,blkibg,blkghk,asfcrr,asfcef,scpot,blkref
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2   v.5:2022/1
!    T.Hanada
!*******************************************************************
subroutine bulkio(inegpos,idiag,iprn)
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
        real(8), parameter :: pi2=atan(1d0)*8d0, rad=pi2/360d0
!----------input(3)----------
! beam parameters
        read (3,*) nh,nk,ndom
        if (nh < 1 .or. nk < 1) then
          write (*,*) ' bulkio.f90 line 29: bad nh,nk ',nh,nk
          stop
        endif
        if (ndom < 1)  then
          write (*,*) ' bulkio.f90 line 29: bad ndom ',ndom
          stop
        endif
allocate (nb(ndom)); allocate (rdom(ndom))
        read (3,*) (nb(i),i=1,ndom)
        nbt=0
        do i=1,ndom
          if (nb(i) < 1) then
            write (*,*) ' bulkio.f90 line 39: bad nb ',nb(i)
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
          write (*,*) ' bulkio.f90 line 60: bad nelm ',nelm
          stop
        endif
allocate (iz(nelm)); allocate (da1(nelm)); allocate (sap(nelm))
allocate (Bh(nelm)); allocate (Bk(nelm)); allocate (Bz(nelm)); 
        do i=1,nelm
          read (3,*) iz(i),da1(i),sap(i)
          if (iz(i) == 0 .or. iabs(iz(i)) > 98) then
            write (*,*) ' bulkio.f90 line 68: max. atomic number iz is 98.',iz(i)
            stop
          endif
          read (3,*) Bh(i),Bk(i),Bz(i)
        end do

! structural parameters
        read (3,*) nsg,aa,bb,gam,cc,dx,dy
        if (ndom > 1) then
          if (abs(aa-bb) > 1d-4) then
            write (*,*) ' bulkio.f90 line 73: aa = bb required for multidomain ',aa,bb
            stop
          endif
          if (abs(gam-90d0) < 1d-4) then
            do i=1,ndom
              if (abs(rdom(i)/90d0 - nint(rdom(i)/90d0)) > 1d-4) then
                write (*,*) ' bulkio.f90 line 48: rdom must be an integer multiple of 90 deg ',rdom(i)
                stop
              endif
            end do
          else if (abs(gam-120d0) < 1d-4) then
            do i=1,ndom
              if (abs(rdom(i)/60d0 - nint(rdom(i)/60d0)) > 1d-4) then
                write (*,*) ' bulkio.f90 line 48: rdom must be an integer multiple of 60 deg ',rdom(i)
                stop
              endif
            end do
          else
            write (*,*) ' bulkio.f90 line 73: gam must be 90 or 60 deg for multidomain',gam
            stop
          endif
        endif

! atomic structural parameters
        read (3,*) natm
        if (natm < 1) then
          write (*,*) ' bulkio.f90 line 90: bad natm ',natm
          stop
        endif
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        do i=1,natm
          read (3,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (ielm(i) < 1 .or.  ielm(i) > nelm) then
            write (*,*) ' bulkio.f90 line 98: bad ielm ',ielm(i)
            stop
          endif
        end do
!----------scan range----------
        call scanrg(azi,azf,daz,naz,gi,gf,dg,ng)
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
        ghx=pi2/(aa*nh)
        ghy=-pi2/(aa*tan(gam)*nh)
        gky=pi2/(bb*sin(gam)*nk)
!----------domain----------
      ibt=1
      do idom=1,ndom
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
          gh(i)=(-pi2/nh)*dble(igh(i))
          gk(i)=(-pi2/nk)*dble(igk(i))
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
!-----------incident beam rocking-----------
       call blkref(nv,nbgm,nb(idom),ns,v,vi,iv,dz,epsb,ngr,nbg &
                  ,ih(ibt),ik(ibt),jorg,nh,nk,ml,dx,dy,wn,ghx,ghy,gky &
                  ,azi+rdom(idom),daz,naz,gi,dg,ng,idiag,iprn)
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (igk); deallocate (igh)
deallocate (iv); deallocate (nbg); deallocate (jorg)
        ibt=ibt+nb(idom)
      end do ! idom=1,ndom

deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (Bz); deallocate (Bk); deallocate (Bh)
deallocate (sap); deallocate (da1); deallocate (iz)
deallocate (ik); deallocate (ih)
deallocate (rdom); deallocate (nb)

end subroutine bulkio
