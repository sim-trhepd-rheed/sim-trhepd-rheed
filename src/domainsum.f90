!*******************************************************************
!   rheed multi slice method 
!   weighted sum of intensities from multi-domains
!     AA = BB, GAM = 90 or 120
!   v.1:2014/4  v.1b:2017/3  T.Hanada
!*******************************************************************
        subroutine domainsum(nwdom)
        implicit none
        integer :: nwdom
        real(8), dimension(:),allocatable :: f2,angle,rdom,wdom
        real(8), dimension(:,:,:),allocatable :: s2
        integer, dimension(:),allocatable :: nb,ih,ik
        real(8) :: dr,cdr,sdr,gx,gy,gxr,gyr,gam
        integer :: naz,ng,ndom,mb,nh,nk,nbt,ibt,ibtr,j
        integer :: nrep1,nrep2,irep1,irep2,idom,idomr,ib,ibr
        real(8), parameter :: rad=atan(1d0)/45d0
        character str*80

      read (2,'(A)') str     ! "#azimuths,g-angles,domains,nh,nk,gam"
      if (str /= '#azimuths,g-angles,domains,nh,nk,gam') then
        write (*,'(A)') ' single domain'
        return
      endif
      read (2,*) naz,ng,ndom,nh,nk,gam
      if (ng > naz) then
        nrep1=naz; nrep2=ng
      else
        nrep1=ng;  nrep2=naz
      endif

allocate (nb(ndom)); allocate (rdom(ndom)); allocate (wdom(ndom))
      read (2,'(A)') str     ! "#beams,domain angle,domain weight"
      do idom=1,ndom
        read (2,*) nb(idom),rdom(idom),wdom(idom)
      end do
!----- change weights -----
      if (nwdom == 0) then
        write (*,'(A,I2,A)') ' intensity weights of ',ndom,' domains ? :'
        read (*,*) (wdom(idom),idom=1,ndom)
      endif

      mb=0; idomr=1
      nbt=0
      do idom=1,ndom
        nbt=nbt+nb(idom)
        if (nb(idom) > mb) then
          mb=nb(idom); idomr=idom
        endif
      end do

allocate (ih(nbt)); allocate (ik(nbt))
      read (2,'(A)') str     ! "#ih,ik"
      ibt=0
      do idom=1,ndom
        read (2,*) (ih(j),ik(j),j=ibt+1,ibt+nb(idom))
        if (idom == idomr) ibtr=ibt
        ibt=ibt+nb(idom)
      end do

! clear
allocate (s2(mb,nrep2,nrep1)); allocate (f2(mb)); allocate (angle(nrep2))
      s2(1:mb,1:nrep2,1:nrep1)=0d0
! sum
      ibt=0
      do idom=1,ndom
        if (idom /= idomr) then
          dr=rdom(idomr)-rdom(idom)
          cdr=cos(dr*rad);  sdr=sin(dr*rad)
        endif
        do irep1=1,nrep1
          do irep2=1,nrep2
            read (2,*) angle(irep2),f2(1:nb(idom))
            if (idom == idomr) then
              do ib=1,nb(idom)
                s2(ib,irep2,irep1)=s2(ib,irep2,irep1)+f2(ib)*wdom(idom)
              end do
            else
              do ib=1,nb(idom)
                gx=dble(ih(ib+ibt))/dble(nh)
                gy=dble(ik(ib+ibt))/dble(nk)
                if (abs(gam-120d0) < 1d-4) then ! GAM = 120 deg
                  gx=(gx+gx+gy)/sqrt(3d0)
                endif
                gxr=cdr*gx-sdr*gy
                gyr=sdr*gx+cdr*gy
                do ibr=1,mb
                  gx=dble(ih(ibr+ibtr))/dble(nh)
                  gy=dble(ik(ibr+ibtr))/dble(nk)
                  if (abs(gam-120d0) < 1d-4) then ! GAM = 120 deg
                    gx=(gx+gx+gy)/sqrt(3d0)
                  endif
                  if (abs(gxr-gx) < 1d-4 .and. abs(gyr-gy) < 1d-4) then
                     s2(ibr,irep2,irep1)=s2(ibr,irep2,irep1)+f2(ib)*wdom(idom)
                     exit
                  endif
                end do
              end do
            endif
          end do ! irep2=1,nrep2
          read (2,'(A)') str     ! empty line
        end do ! irep1=1,nrep1
        ibt=ibt+nb(idom)
      end do ! idom=1,ndom

! output (single domain format)
      write (3,'("#azimuths,g-angles,beams")')
      write (3,'(3I4)') naz,ng,mb
      write (3,'("#ih,ik")')
      write (3,'(400I4)') (ih(ib),ik(ib),ib=1+ibtr,mb+ibtr)
      do irep1=1,nrep1
        do irep2=1,nrep2
          write (3,'(E12.4,200(",",E12.4))') &
                angle(irep2),(s2(ib,irep2,irep1),ib=1,mb)
        end do
        write (3,*)
      end do
deallocate (s2); deallocate (f2); deallocate (angle)
deallocate (ih); deallocate (ik)
deallocate (nb); deallocate (rdom); deallocate (wdom)
      return
      end
