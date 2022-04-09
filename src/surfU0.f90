!!!*******************************************************************
!   rheed multi slice method : surface
!   subroutine srfghk,asfcrr,asfcef,scpot
!   U0:2019/5  U1:2021/8 modifications from surfio.f90 are indicated by !!!-lines
!    T.Hanada
!*******************************************************************
subroutine surfioU
        implicit none
        integer, parameter :: nab=4
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:,:,:),allocatable :: asf
        real(8), dimension(:,:),allocatable :: asf_a,asf_b,asf_c
        real(8), dimension(:),allocatable :: rdom,wdom,gh,gk
        real(8), dimension(:),allocatable :: Bh,Bk,Bz,da1,sap
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: iz,ielm,nb,ih,ik
        integer, dimension(:),allocatable :: jorg,nbg,igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: inegpos,idiag,iprn,ndom,nh,nk,naz,ng,ml
        integer :: nelmb,nelms,nelm,natmb,natms,natm,nsgb,nsgs,msa,msb,nsa,nsb
        real(8) :: be,azi,daz,gi,dg,dz,epsb,aa,bb,gam,cc,dxb,dyb,dxs,dys,dthick
        real(8) :: wn,ghx,ghy,gky
        integer :: ns, ngr,nvb,nv,nvm,nbt,ibt
        integer :: idom, i,j
        real(8), parameter :: pi=atan(1d0)*4d0, deg=180d0/pi
!!!
        real(8) :: zout
        real(8), parameter :: eV=-1.0545888d-34*1.0545888d-34*1d20 &
                                /(2d0*9.1095345d-31*1.6021892d-19)
!!!
!----------input(1)----------
! inegpos <=0: electron, >=1: positron
        read (1) ndom
        read (1) inegpos,nh,nk,idiag
allocate (nb(ndom)); allocate (rdom(ndom)); allocate (wdom(ndom))
        read (1) (nb(i),rdom(i),i=1,ndom)
        nbt=0
        do i=1,ndom
          nbt=nbt+nb(i)
        end do
allocate (ih(nbt)); allocate (ik(nbt))
        read (1) (ih(j),ik(j),j=1,nbt)
        read (1) be,azi,daz,naz,gi,dg,ng
        read (1) dz,ml,epsb
! bulk atomic & structural parameters
        read (1) nsgb,aa,bb,gam,cc,dxb,dyb
        read (1) nelmb,natmb
!----------input(2)----------
! surface atomic parameters
        read (2,*) nelms
        nelm=nelms+nelmb

allocate (iz(nelm)); allocate (da1(nelm)); allocate (sap(nelm))
allocate (Bh(nelm)); allocate (Bk(nelm)); allocate (Bz(nelm))
        do i=1,nelms
          read (2,*) iz(i),da1(i),sap(i)
          if (iz(i) == 0 .or. iabs(iz(i)) > 98) then
            write (*,*) ' surfU0.f90 line 57: max. atomic number iz is 98.',iz(i)
            stop
          endif
          read (2,*) Bh(i),Bk(i),Bz(i)
        end do

        do i=nelms+1,nelm
          read (1) iz(i),da1(i),sap(i)
          read (1) Bh(i),Bk(i),Bz(i)
        end do

! surface structural parameters
        read (2,*) nsgs,msa,msb,nsa,nsb,dthick,dxs,dys
        read (2,*) natms
        natm=natms+natmb
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        do i=1,natms
          read (2,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (ielm(i) < 1 .or.  ielm(i) > nelm) then
            write (*,*) ' element #',ielm(i),' is out of range !'
            stop
          endif
        end do

        do i=natms+1,natm
          read (1) ielm(i),ocr(i),x(i),y(i),z(i)
          ielm(i)=ielm(i)+nelms
        end do

        read (2,*,IOSTAT = i) wdom
        if (i /= 0) wdom=1d0
!----------surface slices----------
        if (nsgs < 1) then
          ns=0 ! to see bulk intensity
        else
          ns=int( ( maxval(z(1:natms))+dthick+cc )/dz )+1
        endif
!!!
!-----------atomic scattering factor--------
! domain independent
allocate (asf_a(nab,nelm)); allocate (asf_b(nab,nelm)); allocate (asf_c(nab,nelm))
        call asfcrr(inegpos,be,wn,nelm,iz,da1,nab,asf_a,asf_b,asf_c,Bz,aa*bb*sin(gam))
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!!!
      ibt=1; idom=1
!!!      do idom=1,ndom
!----------input(1)----------
! information of beams in bulk
allocate (jorg(nb(idom)))
        read (1) jorg(1:nb(idom))
        read (1) ngr,nvb
          nvm=1+nb(idom)*(nb(idom)-1)/2
          if (nvb > nvm) then
            write (*,*) ' surfU0.f90 line110:',nvb,' > ',nvm
            stop
          endif
allocate (nbg(ngr))
allocate (igh(nvm)); allocate (igk(nvm))
        read (1) nbg(1:ngr)
        read (1) igh(1:nvb)
        read (1) igk(1:nvb)
!-----------scattering vector & atomic scattering factor--------
allocate (iv(nb(idom),nb(idom)))
        call srfghk(nvm,nvb,nb(idom),ih(ibt),ik(ibt), nv,igh,igk,iv)
allocate (asf(nv,nab,nelm))
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,nab,asf_a,asf_c,asf,Bh,Bk)
!!!        if (iprn >= 0) write (4,'(A,I0)') 'nv = ',nv
!-----------elements of scattering matrix : v/vi----------
allocate (gh(nv)); allocate (gk(nv))
        do i=1,nv
          gh(i)=(-(pi+pi)/nh)*dble(igh(i))
          gk(i)=(-(pi+pi)/nk)*dble(igk(i))
        end do
allocate (v(nv,ns)); allocate (vi(nv,ns))
! surface reconstructed layer
        call scpot(1,nv,nv,ns,v,vi,natms,nelm,ielm,z &
             ,-cc,dz,nab,asf,asf_b,sap,nsgs,msa,msb,nsa,nsb,gh,gk &
             ,ocr,x,y,dxs+dxb,dys+dyb)
! topmost bulk unit layer (in the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm(natms+1),z(natms+1) &
          ,0d0,dz,nab,asf,asf_b,sap,nsgb,1,0,0,1,gh,gk &
          ,ocr(natms+1),x(natms+1),y(natms+1),dxb,dyb)
! second bulk unit layer (below the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm(natms+1),z(natms+1) &
          ,cc,dz,nab,asf,asf_b,sap,nsgb,1,0,0,1,gh,gk &
          ,ocr(natms+1),x(natms+1),y(natms+1),0d0,0d0)
!!!-----------output U00----------
      zout=0.5d0*dz-cc ! one bulk unit is added below surface layer
      do i=1,ns
        write (4,'(ES12.4,2(",",ES12.4))') zout,dble(v(1,i))*eV,imag(vi(1,i))*eV
        zout=zout+dz
      end do
!!!
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (asf)
deallocate (iv); deallocate (igk); deallocate (igh)
deallocate (nbg); deallocate (jorg)
!!!        ibt=ibt+nb(idom)
!!!      end do ! idom=1,ndom

deallocate (asf_c); deallocate (asf_b); deallocate (asf_a)
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (Bz); deallocate (Bk); deallocate (Bh)
deallocate (sap); deallocate (da1); deallocate (iz)
deallocate (ik); deallocate (ih)
deallocate (wdom); deallocate (rdom); deallocate (nb)

end subroutine surfioU
