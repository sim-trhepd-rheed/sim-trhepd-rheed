!*******************************************************************
!   RHEED multi slice method
!      transfer matrix of 1 slice --> t1
!   v.1:84/10   v.2:86/11   v.3:90/4    T.Hanada
!*******************************************************************
      subroutine trmat(nv,nbm,nb,gma,v,vi,iv,dz,t1,t2,t3,idiag)
      implicit none
      integer :: iv(nbm,nb),nv,nbm,nb,idiag
      complex(8) :: v(nv),vi(nv),gma(nb)
      complex(8) :: t1(nbm+nbm,nbm+nbm),t2(nbm+nbm,nbm+nbm),t3(nbm+nbm,nbm+nbm)
      real(8) :: dz

      complex(8) :: ph,ph2,tau,rou,x1,x2,x3, x4
      real(8) :: eig(nb), s,r
      integer :: nbm2, i,j,icom, i2,j2, it,k

      nbm2=nbm+nbm
        if (nb == 1) then
          t2(2,1)=v(1)+vi(1)+gma(1)*gma(1)
          t2(1,1)=(1d0,0d0)
        else
! h=t1, v=t2, e=t2(nb+1,nb)
          do i=1,nb
            t1(i,i)=v(1)+vi(1)+gma(i)*gma(i)
            do j=i+1,nb
              t1(j,i)=v(iv(j,i))+vi(iv(j,i))
              t1(i,j)=dconjg(v(iv(j,i))-vi(iv(j,i)))
            end do
          end do
! nb >= 2
          call zgeev('n','v',nb,t1,nbm2,t2(nb+1,nb),t2(1,nb+1) &
                  ,nbm2,t2,nbm2,t3,nbm2*nbm2,t1(1,nb+1),icom)
          if (icom /= 0) write (*,*) ' zgeev failed ',icom
        endif
!----------transfer matrix of 1 slice----------
        do j=1,nb
          j2=j+nb
          x2=sqrt(t2(j2,nb))
          ph=exp(dcmplx(0d0,dz)*x2)
          ph2=(1d0,0d0)/ph
          do i=1,nb
            i2=i+nb
            x1=gma(i)
            x3=t2(i,j)
            tau=(x1+x2)*x3
            rou=(x1-x2)*x3
            t1(i,j)=tau*ph
            t1(i,j2)=rou*ph2
            t1(i2,j)=rou*ph
            t1(i2,j2)=tau*ph2
            t3(i,j)=tau
            t3(i,j2)=rou
            t3(i2,j)=rou
            t3(i2,j2)=tau
          end do
        end do
        call gcmi2(nbm2,nb+nb,t1,t3)
      return
      end
