!*******************************************************************
!   rheed multi slice method :  subroutines for bulk
!    v.3    90/4/25           T.Hanada
!   contains scanrg blkibg blkghk
!*******************************************************************
!**********************************************************
!       scan range
!**********************************************************
        subroutine scanrg(azi,azf,daz,naz,gi,gf,dg,ng)
        implicit none
        integer :: naz,ng
        real(8) :: azi,azf,daz,gi,gf,dg
        integer, parameter :: nazm=200, ngm=200
        real(8), parameter :: rad=atan(1d0)/45d0 ! rad=3.141592654d0/180d0
! max # of azimuth steps :nazm
! max # of glancing angle steps : ngm
        if (abs(daz) < 1d-4) then
          naz=1
        else
          naz=int((azf-azi)/daz+0.1)+1
          if (naz < 1) then
            naz=1
          else if (naz > nazm) then
            naz=nazm
          endif
        endif
        azi=azi*rad
        daz=daz*rad
        if (abs(dg) < 1d-4) then
          ng=1
        else
          ng=int((gf-gi)/dg+.1)+1
          if (ng < 1) then
            ng=1
          else if (ng > ngm) then
            ng=ngm
          endif
        endif
        gi=gi*rad
        dg=dg*rad
        return
        end
!**********************************************************
!       groups of interacting beams in bulk layer
!**********************************************************
      subroutine blkibg(nb,ih,ik,nh,nk, iord,ngr,nbg)
      implicit none
      integer :: nb,nh,nk,ngrm,ngr
      integer :: ih(nb),ik(nb),iord(nb),nbg(nb)
      integer :: ih0,ik0,i,j,k,l,m,n,jj

      do i=1,nb
        iord(i)=i
      end do

      if (nh == 1 .and. nk == 1) then
        ngr=1
        nbg(1)=nb
      else
        j=1
        ngr=0

        do
          i=j
          j=i+1
          ih0=ih(i)
          ik0=ik(i)
          ngr=ngr+1
          do k=i+1,nb
            if (mod(ih(k)-ih0,nh) == 0 .and. mod(ik(k)-ik0,nk) == 0) then
              if (k > j) then
                l=ih(k)
                n=ik(k)
                m=iord(k)
                do jj=k,j+1,-1
                  ih(jj)=ih(jj-1)
                  ik(jj)=ik(jj-1)
                  iord(jj)=iord(jj-1)
                end do
                ih(j)=l
                ik(j)=n
                iord(j)=m
              endif
              j=j+1
            endif
          end do
          nbg(ngr)=j-i
          if (j > nb) return
        end do
      endif
      return
      end
!**********************************************************
!       bulk ghk
!**********************************************************
        subroutine blkghk(ngr,nbg,nvm,nv,igh,igk,nbgm,nb,ih,ik,iv)
        implicit none
        integer :: ngr,nvm,nv,nbgm,nb
        integer :: nbg(ngr),igh(nvm),igk(nvm),ih(nb),ik(nb),iv(nbgm,nb)
        integer :: ibas,igr,iend,ih0,ik0,k,l,m,iflag
        igh(1)=0
        igk(1)=0
        nv=1
        ibas=0
        do igr=1,ngr
          iend=ibas+nbg(igr)
          do k=1+ibas,iend-1
            do l=k+1,iend
              ih0=ih(l)-ih(k)
              ik0=ik(l)-ik(k)
!----------
              iflag=0
              do m=2,nv
                if (ih0 == igh(m) .and. ik0 == igk(m)) then
                  iv(l-ibas,k)=m; iflag=1
                  exit
                endif
              end do

              if (iflag == 0) then
                nv=nv+1
                if (nv > nvm) then
                  write (*,*) ' blkghk: nv=',nv,' > nvm=',nvm
                  stop
                endif
                igh(nv)=ih0
                igk(nv)=ik0
                iv(l-ibas,k)=nv
              endif
            end do
          end do
          ibas=iend
        end do
        return
        end
