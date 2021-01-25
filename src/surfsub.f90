!*******************************************************************
!   rheed multi slice method :  subroutines for surf
!    v.3    90/4/12         T.Hanada
!   contains srfghk
!*******************************************************************
!**********************************************************
!       surface ghk
!**********************************************************
        subroutine srfghk(nvm,nv,nvb,igh,igk,nb,ih,ik,iv)
        implicit none
        integer :: nvm,nv,nvb,nb
        integer :: igh(nvm),igk(nvm),ih(nb),ik(nb),iv(nb,nb)
        integer :: ih0,ik0,k,l,m,iflag
        nv=nvb
        do k=1,nb-1
          do l=k+1,nb
            ih0=ih(l)-ih(k)
            ik0=ik(l)-ik(k)
!---------
            iflag=0
            do m=2,nv
              if (ih0 == igh(m) .and. ik0 == igk(m)) then
                iv(l,k)=m; iflag=1
                exit
              endif
            end do

            if (iflag == 0) then
              nv=nv+1
              if (nv > nvm) then
                write (*,*) ' srfghk: nv=',nv,' > nvm=',nvm
                stop
              endif
              igh(nv)=ih0
              igk(nv)=ik0
              iv(l,k)=nv
            endif
          end do
        end do
        return
        end
