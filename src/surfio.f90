!*******************************************************************
!   rheed multi slice method : surface
!   subroutine srfghk,asfcrr,asfcef,scpot,srfref
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2
!    T.Hanada
!*******************************************************************
! max # of potential fourier components : nvm
!*******************************************************************
        subroutine surfio(iprn)
        implicit none
        integer, parameter :: mab=6
        complex(8), dimension(:,:),allocatable :: v,vi
        real(8), dimension(:),allocatable :: rdom,wdom,gh,gk,Ghk
        real(8), dimension(:),allocatable :: dth,dtk,dtz,elm,dchg
        real(8), dimension(:,:),allocatable :: a,b,c
        real(8), dimension(:,:,:),allocatable :: asf
        real(8), dimension(:),allocatable :: ocr,x,y,z
        integer, dimension(:),allocatable :: ielm,nb,iord,ih,ik,nbg
        integer, dimension(:),allocatable :: igh,igk
        integer, dimension(:,:),allocatable :: iv

        integer :: nvm, ndom,iprn
        real(8) :: be,wn,azi,daz,gi,dg,tem,dz,epsb,aa,bb,gam,cc,dxb,dyb,zoffset
        real(8) :: smesh,dthick,dxs,dys,ghx,ghy,gky,da1,sap,sae,zz,amass
        integer :: inegpos,isym,nh,nk,naz,ng,ml,nelmb,nsgb,mp,natmb
        integer :: imd,nelms,nelm,nsgs,msa,msb,nsa,nsb,natms,natm,ns, ngr,nvb,nv
        integer :: nbt,ibt,nabr,nabi,iabr,iabp,iabe,iabrs,iabps,iabes,idiag
        integer :: iz,ichg,idom, i,j
        real(8), parameter :: pi=atan(1d0)*4d0, deg=180d0/pi
!----------input(1)----------
! inegpos <=0: electron, >=1: positron
        read (1) ndom
        read (1) inegpos,isym,nh,nk,iabr,iabp,iabe,nabr,nabi,idiag
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
        imd=0
! surface atomic parameters
        read (2,*) nelms
        nelm=nelmb+nelms

allocate (dth(nelm)); allocate (dtk(nelm)); allocate (dtz(nelm))
allocate (elm(nelm)); allocate (dchg(nelm))
allocate (a(mab,nelm)); allocate (b(mab,nelm)); allocate (c(mab,nelm))
        do i=1,nelmb
          read (1) a(1:nabr+nabi,i),b(1:nabr+nabi,i)
          read (1) dth(i),dtk(i),dtz(i),elm(i),dchg(i)
        end do

        do i=nelmb+1,nelm
          read (2,*) iz,da1,sap;  ichg=0; sae=0d0
          dchg(i)=dble(ichg)
          elm(i)=amass(iz)

          call asfparam(iz,iabr,nabr,a(1,i),b(1,i))
            if (nabr <= 0) return
          a(1,i)=a(1,i)-da1                           ! electron
          if (inegpos > 0) a(1:nabr,i)=-a(1:nabr,i)   ! positron

            nabi=1 
          if (nabr+nabi > mab) return
          a(nabr+1,i)=sap; b(nabr+1,i)=0d0 ! V(imag)=a*V(real), b is dummy

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
!----------output(4)----------
        write (4,'(A)') 'inegpos,idiag, nh,nk,ndom, naz,ng' !isym,imd
        write (4,'(9I4)') inegpos,idiag, nh,nk,ndom, naz,ng !isym,imd
        write (4,'(A)') '(nb(i),rdom(i)*deg,wdom(i),i=1,ndom)'
        write (4,'(6(I4,F10.3,E12.3))') (nb(i),rdom(i)*deg,wdom(i),i=1,ndom)
        write (4,'(A)') '(ih(j),ik(j),j=1,nb(i))'
        ibt=0
        do i=1,ndom
          write (4,'(20I4)') (ih(j),ik(j),j=ibt+1,ibt+nb(i))
          ibt=ibt+nb(i)
        end do
        write (4,'(A)') 'be,tem,azi*deg,daz*deg,gi*deg,dg*deg,epsb'
        write (4,'(6F10.3,E12.4)') be,tem,azi*deg,daz*deg,gi*deg,dg*deg,epsb
! atomic parameters
        write (4,'(A)') '----- atomic parameters -----'
        write (4,'(A)')  'nelm, iabr,iabp,iabe, nabr,nabi'
        write (4,'(6I4)') nelm, iabr,iabp,iabe, nabr,nabi
          write (4,'(A)') 'ielm'
          write (4,'(A)') '(a(j,i),b(j,i),j=1,nabr)'
          write (4,'(A)') '(a(j,i),b(j,i),j=nabr+1,nabr+nabi)'
          write (4,'(A)') 'dth(i),dtk(i),dtz(i),elm(i),dchg(i)'
        do i=nelmb+1,nelm
          write (4,'(I4)') i-nelmb
          write (4,'(20F10.5)') (a(j,i),b(j,i),j=1,nabr)
          write (4,'(20F10.5)') (a(j,i),b(j,i),j=nabr+1,nabr+nabi)
          write (4,'(6F10.5)') dth(i),dtk(i),dtz(i),elm(i),dchg(i)
        end do
        do i=1,nelmb
          write (4,'(I4)') -i
          write (4,'(20F10.5)') (a(j,i),b(j,i),j=1,nabr)
          write (4,'(20F10.5)') (a(j,i),b(j,i),j=nabr+1,nabr+nabi)
          write (4,'(6F10.5)') dth(i),dtk(i),dtz(i),elm(i),dchg(i)
        end do
! structural parameters
        write (4,'(A)') '----- structural parameters -----'
        write (4,'(A)') 'aa,bb,gam*deg,cc,dxb,dyb,dxs,dys,zoffset'
        write (4,'(9F9.4)') aa,bb,gam*deg,cc,dxb,dyb,dxs,dys,zoffset
        write (4,'(A)') 'nsgb,ml,nsgs,msa,msb,nsa,nsb,natmb,natms,ns,dz'
        write (4,'(10I4,F8.5)') nsgb,ml,nsgs,msa,msb,nsa,nsb,natmb,natms,ns,dz
          write (4,'(A)') 'ielm(i),ocr(i),x(i),y(i),z(i)'
        do i=natmb+1,natm
          write (4,'(I4,4F10.5)') ielm(i)-nelmb,ocr(i),x(i),y(i),z(i)-zoffset
        end do
        do i=1,natmb
          write (4,'(I4,4F10.5)') -ielm(i),ocr(i),x(i),y(i),z(i)-zoffset
        end do
!-----------energy correction of atomic scattering factor--------
        smesh=abs(msa*nsb-msb*nsa)
        call asfcrr(be,wn,nelm,mab,nabr,nabi,a,b,c,dchg,dth,dtk,dtz,elm &
                   ,tem,aa*bb*sin(gam)*smesh)
        ghx=(pi+pi)/(aa*nh)
        ghy=-(pi+pi)/(aa*tan(gam)*nh)
        gky=(pi+pi)/(bb*sin(gam)*nk)
!----------domain----------
        if (ndom == 1) then
          write (3,'("#azimuths,g-angles,beams")')
          write (3,'(3I4)') naz,ng,nb(1)
          write (3,'("#ih,ik")')
          write (3,'(400I4)') (ih(j),ik(j),j=1,nb(1))
        else
          write (3,'("#azimuths,g-angles,domains,nh,nk,gam")')
          write (3,'(5I4,F10.3)') naz,ng,ndom,nh,nk,gam*deg
          write (3,'("#beams,domain angle,domain weight")')
          do i=1,ndom
            write (3,'(I4,F10.3,E12.3)') nb(i),rdom(i)*deg,wdom(i)
          end do
          write (3,'("#ih,ik")')
          ibt=0
          do i=1,ndom
            write (3,'(400I4)') (ih(j),ik(j),j=ibt+1,ibt+nb(i))
            ibt=ibt+nb(i)
          end do
        endif

      ibt=1
      do idom=1,ndom
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
allocate (asf(nv,mab,nelm)); allocate (Ghk(nv))
        call asfcef(nv,nelm,igh,igk,ghx,ghy,gky,mab,nabr,nabi,a,c,asf,dth,dtk,Ghk)
        write (4,'(A,I4)') 'nv',nv
!-----------elements of scattering matrix : v/vi----------
allocate (gh(nv)); allocate (gk(nv))
        do i=1,nv
          gh(i)=(-(pi+pi)/nh)*dble(igh(i))
          gk(i)=(-(pi+pi)/nk)*dble(igk(i))
        end do
allocate (v(nv,ns)); allocate (vi(nv,ns))
! surface reconstructed layer
        call scpot(1,nv,nv,ns,v,vi,natms,nelm,ielm(natmb+1),z(natmb+1),-cc,dz &
             ,mab,nabr,nabi,asf,b,dchg,1d0,nsgs,msa,msb,nsa,nsb,gh,gk,Ghk &
             ,ocr(natmb+1),x(natmb+1),y(natmb+1),dxs+dxb,dys+dyb)
! topmost bulk unit layer (in the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm,z,0d0,dz,mab &
          ,nabr,nabi,asf,b,dchg,smesh,nsgb,1,0,0,1,gh,gk,Ghk,ocr,x,y,dxb,dyb)
! second bulk unit layer (below the 'ns' slices)
        call scpot(0,nv,nvb,ns,v,vi,natmb,nelm,ielm,z,cc,dz,mab &
          ,nabr,nabi,asf,b,dchg,smesh,nsgb,1,0,0,1,gh,gk,Ghk,ocr,x,y,0d0,0d0)
!-----------incident beam rocking-----------
        if (iprn == 777) ns=0 ! to see bulk intensity
        call srfref(nv,nb(idom),ns,v,vi,iv,dz,ngr,nbg,iord,ih(ibt),ik(ibt) &
                ,wn,ghx,ghy,gky,azi+rdom(idom),daz,naz,gi,dg,ng,idiag,iprn)
deallocate (vi); deallocate (v)
deallocate (gk); deallocate (gh)
deallocate (Ghk); deallocate (asf)
deallocate (iv); deallocate (igk); deallocate (igh)
deallocate (nbg); deallocate (iord)
        ibt=ibt+nb(idom)
      end do ! idom=1,ndom
deallocate (z); deallocate (y); deallocate (x)
deallocate (ocr); deallocate (ielm)
deallocate (c); deallocate (b); deallocate (a)
deallocate (dchg); deallocate (elm)
deallocate (dtz); deallocate (dtk); deallocate (dth)
deallocate (ik); deallocate (ih)
deallocate (wdom); deallocate (rdom); deallocate (nb)
      return
      end
