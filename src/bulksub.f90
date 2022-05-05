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
        real(8), parameter :: rad=atan(1d0)/45d0
        if (abs(daz) < 1d-4) then
          naz=1
        else
          naz=max(int((azf-azi)/daz+0.1d0)+1,1)
        endif
        azi=azi*rad
        daz=daz*rad
        if (abs(dg) < 1d-4) then
          ng=1
        else
          ng=max(int((gf-gi)/dg+0.1d0)+1,1)
        endif
        gi=gi*rad
        dg=dg*rad
end subroutine scanrg
!**********************************************************
!       groups of interacting beams in bulk layer
!**********************************************************
subroutine blkibg(nb,nh,nk,ih,ik, jorg,ngr,nbg)
      implicit none
      integer :: nb,nh,nk,ngr
      integer :: ih(nb),ik(nb),jorg(nb),nbg(nb), idone(nb)
      integer :: ih0,ik0,init,jinit,jnext,inext,iflag, i,j

      ngr=1
      if (nh == 1 .and. nk == 1) then
        nbg(1)=nb
        do j=1,nb
          jorg(j)=j
        end do
      else
        init=1
        jinit=1
        idone(2:nb)=0
        do ! ngr loop
          jorg(jinit)=init
          jnext=jinit+1
          idone(init)=1
          ih0=ih(init)
          ik0=ik(init)
          iflag=0
          do i=init+1,nb
            if (idone(i) == 0) then
              if (mod(ih(i)-ih0,nh) == 0 .and. mod(ik(i)-ik0,nk) == 0) then
                idone(i)=1
                jorg(jnext)=i
                jnext=jnext+1
              else if (iflag == 0) then
                inext=i
                iflag=1
              endif
            endif
          end do
          nbg(ngr)=jnext-jinit
          if (jnext > nb) return ! normally finish
          ngr=ngr+1
          jinit=jnext
          init=inext
        end do
      endif
end subroutine blkibg
!**********************************************************
!       bulk ghk
!**********************************************************
subroutine blkghk(ngr,nbg,nbgm,nb,ih,ik,jorg,nvm, nv,igh,igk,iv)
        implicit none
        integer :: ngr,nvm,nv,nbgm,nb
        integer :: nbg(ngr),igh(nvm),igk(nvm),ih(nb),ik(nb),jorg(nb),iv(nbgm,nb)
        integer :: ibas,igr,iend,ih0,ik0,k,l,m,iflag
        igh(1)=0
        igk(1)=0
        nv=1
        ibas=0
        do igr=1,ngr
          iend=ibas+nbg(igr)
          do k=1+ibas,iend-1
            do l=k+1,iend
              ih0=ih(jorg(l))-ih(jorg(k))
              ik0=ik(jorg(l))-ik(jorg(k))
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
end subroutine blkghk
