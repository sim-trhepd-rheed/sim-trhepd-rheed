!*******************************************************************
!   rheed multi slice method
!   subroutine trmat gcmi2 gcim2
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4   v.4b:2017/2
!    T.Hanada
!*******************************************************************
       subroutine blkref(nv,nbgm,nb,ns,v,vi,iv,dz    &
                  ,epsb,ngr,nbg,ih,ik,nh,nk,ml &
                  ,dx,dy,wn,ghx,ghy,gky,azi,daz,naz,gi,dg,ng,idiag,iprn)
        implicit none
        complex(8) :: v(nv,ns),vi(nv,ns)
        real(8) :: dz,epsb,dx,dy,wn,ghx,ghy,gky,azi,daz,gi,dg
        integer :: iv(nbgm,nb),ih(nb),ik(nb),nbg(ngr)
        integer :: nv,nbgm,nb,ns,ngr,nh,nk,ml,naz,ng,idiag,iprn

        complex(8) :: t1(nbgm+nbgm,nbgm+nbgm),t2(nbgm+nbgm,nbgm+nbgm)
        complex(8) :: t3(nbgm+nbgm,nbgm+nbgm),t4(nbgm+nbgm,nbgm+nbgm),gma(nbgm),ph
        real(8) :: r20(nbgm)
        real(8) :: az,caz,saz,ga,cga,wnx,wny,wgx,wgy,s,pih,pik,rmax,r2,rat
        integer :: nbgm2,nrep1,nrep2,irep1,irep2,ibas,igr,nbf,nbf2,i,i2,j,j2,k,l
        real(8), parameter :: pi2=atan(1d0)*8d0

        nbgm2=nbgm+nbgm
        if (ng > naz) then
          nrep1=naz
          nrep2=ng
        else
          nrep1=ng
          nrep2=naz
        endif
!----------azimuth,glancing angle scan----------
      do irep1=1,nrep1
        if (ng > naz) then
          az=azi+(irep1-1)*daz
          caz=cos(az)
          saz=sin(az)
        else
          cga=cos(gi+(irep1-1)*dg)
        endif
      do irep2=1,nrep2
        if (ng > naz) then
          cga=cos(gi+(irep2-1)*dg)
        else
          az=azi+(irep2-1)*daz
          caz=cos(az)
          saz=sin(az)
        endif
        wnx=wn*cga*caz
        wny=wn*cga*saz
!----------interacting beams----------
        ibas=0
      do igr=1,ngr
          nbf=nbg(igr)
          nbf2=nbf+nbf
          if (iprn == 1) write (*,'(A,5I5)') &
             ' irep1,irep2,ngr,igr,nbg=',irep1,irep2,ngr,igr,nbf
!----------z component of scattering vector----------
          do i=1,nbf
            wgx=ghx*ih(i+ibas)+wnx
            wgy=ghy*ih(i+ibas)+gky*ik(i+ibas)+wny
            s=wn*wn-wgx*wgx-wgy*wgy
            if (s >= 0d0) then
              gma(i)=dcmplx(sqrt(s))
            else
              gma(i)=dcmplx(0d0,sqrt(-s))
            endif
          end do
!----------transfer matrix of 1 slice (t1)-----------
        do l=1,ns
          call trmat(nv,nbgm,nbf,gma &
               ,v(1,l),vi(1,l),iv(1,ibas+1),dz,t1,t2,t3,idiag)
!------------products of 1 slice transfer matrix (t4 <- t3=t1*t4)
          if (ns == 1) then
            t3(1:nbf2,1:nbf2)=t1(1:nbf2,1:nbf2)
          else if (l == 1) then
            t4(1:nbf2,1:nbf2)=t1(1:nbf2,1:nbf2)
          else
            do i=1,nbf2
              do j=1,nbf2
                t3(i,j)=sum(t1(i,1:nbf2)*t4(1:nbf2,j))
              end do
            end do
            if (l < ns) t4(1:nbf2,1:nbf2)=t3(1:nbf2,1:nbf2)
          endif
        end do
!-----------transfer matrix of unit layer (t3)-----------
        if (abs(dx) > 1d-4 .or. abs(dy) > 1d-4) then
          pih=dx*pi2/nh
          pik=dy*pi2/nk
          do j=1,nbf
            s=pih*ih(j+ibas)+pik*ik(j+ibas)
            ph=dcmplx(cos(s),sin(s))
            j2=j+nbf
            t3(1:nbf2,j) =t3(1:nbf2,j)*ph
            t3(1:nbf2,j2)=t3(1:nbf2,j2)*ph
          end do
        endif
!----------diffraction from bulk layer----------
        do l=1,ml
          do j=1,nbf
            j2=j+nbf
            do i=1,nbf
              t1(i,j)=t3(i,j2)
              t2(i,j)=t3(i+nbf,j2)
            end do
            if (l > 1) then
              do k=1,nbf
                do i=1,nbf
                  t1(i,j)=t1(i,j)+t3(i,k)*t4(k,j)
                  t2(i,j)=t2(i,j)+t3(i+nbf,k)*t4(k,j)
                end do
              end do
            endif
          end do
          call gcmi2(nbgm2,nbf,t1,t2)
          t4(1:nbf,1:nbf)=t1(1:nbf,1:nbf)
!----------judgement----------
          if (l == 20) then
            do i=1,nbf
              r20(i)=dble(t4(i,i)*dconjg(t4(i,i)))
            end do
          else if (l > 20) then
            rmax=0d0
            do i=1,nbf
              r2=dble(t4(i,i)*dconjg(t4(i,i)))
              rat=abs(r2-r20(i))
              if (rat > rmax) rmax=rat
              r20(i)=r2
            end do
            if (rmax < epsb) exit
          endif
        end do ! l=1,ml
!----------output----------
        if (iprn == 1) write (*,*) l,' layers'
        write (1) ((t4(i,j),i=1,nbf),j=1,nbf)
        ibas=ibas+nbf
      end do ! igr=1,ngr
      end do ! irep2=1,nrep2
      end do ! irep1=1,nrep1
      return
      end
