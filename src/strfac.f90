!**********************************************************
!   rheed multislice method :  structure factor
!    v.3    90/4/12   v.4:2014/4    T.Hanada
!**********************************************************
        subroutine strfac(nsg,nv,ma,mb,na,nb,gh,gk,x,y,dx,dy,st)
        implicit none
        integer :: nsg,nv,ma,mb,na,nb, ndeg,m,i
        complex(8) :: st(nv)
        real(8) :: gh(nv),gk(nv),x,y,dx,dy
        real(8) :: reduce01, s,xs(12),ys(12),xb(12),yb(12),gr

! (x,y) 1x1 lattice --> (xs,ys) surface lattice
        s=dble(ma*nb-mb*na)
        xs(1)=reduce01((nb*x-na*y)/s)
        ys(1)=reduce01((ma*y-mb*x)/s)

        if (nsg < 2 .or. nsg > 17) then ! p1
          nsg=1; ndeg=1
        else if (nsg == 2) then
          call p211(xs,ys,ndeg)
        else if (nsg == 3) then
          call p1m1(xs,ys,ndeg)
        else if (nsg == 4) then
          call p1g1(xs,ys,ndeg)
        else if (nsg == 5) then
          call c1m1(xs,ys,ndeg)
        else if (nsg == 6) then
          call p2mm(xs,ys,ndeg)
        else if (nsg == 7) then
          call p2mg(xs,ys,ndeg)
        else if (nsg == 8) then
          call p2gg(xs,ys,ndeg)
        else if (nsg == 9) then
          call c2mm(xs,ys,ndeg)
        else if (nsg == 10) then
          call p4(xs,ys,ndeg)
        else if (nsg == 11) then
          call p4mm(xs,ys,ndeg)
        else if (nsg == 12) then
          call p4gm(xs,ys,ndeg)
        else if (nsg == 13) then
          call p3(xs,ys,ndeg)
        else if (nsg == 14) then
          call p3m1(xs,ys,ndeg)
        else if (nsg == 15) then
          call p31m(xs,ys,ndeg)
        else if (nsg == 16) then
          call p6(xs,ys,ndeg)
        else
          call p6mm(xs,ys,ndeg)
        endif

! (xs,ys) surface lattice --> (xb,yb) 1x1 lattice
        do i=1,ndeg
          xs(i)=reduce01(xs(i))
          ys(i)=reduce01(ys(i))
          xb(i)=ma*xs(i)+na*ys(i)
          yb(i)=mb*xs(i)+nb*ys(i)
!          write (4,'(4F10.5)') xb(i),yb(i), xs(i),ys(i)
        end do
!          write (4,*) ' '

        do m=1,nv
          st(m)=(0d0,0d0)
          do i=1,ndeg
            gr=gh(m)*(dx+xb(i))+gk(m)*(dy+yb(i))
            st(m)=st(m)+dcmplx(cos(gr),sin(gr))
          end do
        end do
        return
        end
!**********************************************************
!       subtract integer part to let -1d-4 <= reduce01 < 1-1d-4
!**********************************************************
        real(8) function reduce01(x)
        real(8) :: x
        reduce01=x-floor(x+1d-4)
        return
        end
!**********************************************************
!       p211, nsg=2
!**********************************************************
        subroutine p211(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(2),ys(2)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) .and. &
             (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) ) then
          ndeg=1 ! twofold axis
        else
          ndeg=2
          xs(2)=-x; ys(2)=-y
        endif
        return
        end
!**********************************************************
!       p1m1, nsg=3
!**********************************************************
        subroutine p1m1(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(2),ys(2)

        x=xs(1); y=ys(1)
        if (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) then
          ndeg=1 ! mirror plane
        else
          ndeg=2
          xs(2)=x; ys(2)=-y
        endif
        return
        end
!**********************************************************
!       p1g1, nsg=4
!**********************************************************
        subroutine p1g1(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(2),ys(2)

        x=xs(1); y=ys(1)
        ndeg=2
        xs(2)=x+0.5d0; ys(2)=-y
        return
        end
!**********************************************************
!       c1m1, nsg=5
!**********************************************************
        subroutine c1m1(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(4),ys(4)

        x=xs(1); y=ys(1)
        if (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) then
          ndeg=2 ! mirror plane
          xs(2)=x+0.5d0; ys(2)=0.5d0-y
        else
          ndeg=4
          xs(2)=x+0.5d0; ys(2)=0.5d0-y
          xs(3)=x;       ys(3)=-y
          xs(4)=x+0.5d0; ys(4)=0.5d0+y
        endif
        return
        end
!**********************************************************
!       p2mm, nsg=6
!**********************************************************
        subroutine p2mm(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(4),ys(4)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) .and. &
             (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) ) then
          ndeg=1 ! twofold axis
        else if (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) then
          ndeg=2 ! mirror plane along y
          xs(2)=x; ys(2)=-y
        else if (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) then
          ndeg=2 ! mirror plane along x
          xs(2)=-x; ys(2)=y
        else
          ndeg=4
          xs(2)=x;  ys(2)=-y
          xs(3)=-x; ys(3)=y
          xs(4)=-x; ys(4)=-y
        endif
        return
        end
!**********************************************************
!       p2mg, nsg=7
!**********************************************************
        subroutine p2mg(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(4),ys(4)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) .and. &
             (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis
          xs(2)=x;  ys(2)=0.5d0-y
        else if (abs(y-0.25d0) < 2d-4 .or. abs(y-0.75d0) < 2d-4) then
          ndeg=2 ! mirror plane
          xs(2)=-x; ys(2)=-y
        else
          ndeg=4
          xs(2)=x;  ys(2)=0.5d0-y
          xs(3)=-x; ys(3)=-y
          xs(4)=-x; ys(4)=0.5d0+y
        endif
        return
        end
!**********************************************************
!       p2gg, nsg=8
!**********************************************************
        subroutine p2gg(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(4),ys(4)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) .and. &
             (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis
          xs(2)=0.5d0-x; ys(2)=0.5d0+y
        else
          ndeg=4
          xs(2)=0.5d0-x; ys(2)=0.5d0+y
          xs(3)=0.5d0+x; ys(3)=0.5d0-y
          xs(4)=-x;      ys(4)=-y
        endif
        return
        end
!**********************************************************
!       c2mm, nsg=9
!**********************************************************
        subroutine c2mm(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(8),ys(8)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) .and. &
             (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis at m-cross
          xs(2)=0.5d0-x; ys(2)=0.5d0+y
        else if ( (abs(x-0.25d0) < 2d-4 .or. abs(x-0.75d0) < 2d-4) .and. &
                  (abs(y-0.25d0) < 2d-4 .or. abs(y-0.75d0) < 2d-4) ) then
          ndeg=4 ! twofold axis at g-cross
          xs(2)=x;       ys(2)=-y
          xs(3)=-x;      ys(3)=y
          xs(4)=-x;      ys(4)=-y
        else if (abs(x) < 2d-4 .or. abs(x-0.5d0) < 2d-4) then
          ndeg=4 ! mirror plane along y
          xs(2)=x;       ys(2)=-y
          xs(3)=0.5d0-x; ys(3)=0.5d0+y
          xs(4)=0.5d0-x; ys(4)=0.5d0-y
        else if (abs(y) < 2d-4 .or. abs(y-0.5d0) < 2d-4) then
          ndeg=4 ! mirror plane along x
          xs(2)=-x;      ys(2)=y
          xs(3)=0.5d0+x; ys(3)=0.5d0-y
          xs(4)=0.5d0-x; ys(4)=0.5d0-y
        else
          ndeg=8
          xs(2)=x;       ys(2)=-y
          xs(3)=-x;      ys(3)=y
          xs(4)=-x;      ys(4)=-y
          xs(5)=0.5d0+x; ys(5)=0.5d0+y
          xs(6)=0.5d0+x; ys(6)=0.5d0-y
          xs(7)=0.5d0-x; ys(7)=0.5d0+y
          xs(8)=0.5d0-x; ys(8)=0.5d0-y
        endif
        return
        end
!**********************************************************
!       p4, nsg=10
!**********************************************************
        subroutine p4(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(4),ys(4)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .and. abs(y) < 2d-4) .or. &
             (abs(x-0.5d0) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=1 ! fourfold axis
        else if ( (abs(x-0.5d0) < 2d-4 .and. abs(y) < 2d-4) .or. &
                  (abs(x) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis
          xs(2)=y;       ys(2)=-x
        else
          ndeg=4
          xs(2)=-y; ys(2)=x
          xs(3)=-x; ys(3)=-y
          xs(4)=y;  ys(4)=-x
        endif
        return
        end
!**********************************************************
!       p4mm, nsg=11
!**********************************************************
        subroutine p4mm(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(8),ys(8)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .and. abs(y) < 2d-4) .or. &
             (abs(x-0.5d0) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=1 ! fourfold axis
        else if ( (abs(x-0.5d0) < 2d-4 .and. abs(y) < 2d-4) .or. &
                  (abs(x) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis
          xs(2)=y;  ys(2)=x
        else if (abs(y-x) < 2d-4 .or. abs(x+y-1d0) < 2d-4 .or. &
                 abs(x) < 2d-4   .or. abs(x-0.5d0) < 2d-4 .or. &
                 abs(y) < 2d-4   .or. abs(y-0.5d0) < 2d-4) then
          ndeg=4 ! mirror plane
          xs(2)=-y; ys(2)=x
          xs(3)=-x; ys(3)=-y
          xs(4)=y;  ys(4)=-x
        else
          ndeg=8
          xs(2)=y;  ys(2)=x
          xs(3)=-y; ys(3)=x
          xs(4)=-x; ys(4)=y
          xs(5)=-x; ys(5)=-y
          xs(6)=-y; ys(6)=-x
          xs(7)=y;  ys(7)=-x
          xs(8)=x;  ys(8)=-y
        endif
        return
        end
!**********************************************************
!       p4gm, nsg=12
!**********************************************************
        subroutine p4gm(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(8),ys(8)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .and. abs(y) < 2d-4) .or. &
             (abs(x-0.5d0) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! fourfold axis
          xs(2)=0.5d0-y; ys(2)=0.5d0-x
        else if ( (abs(x-0.5d0) < 2d-4 .and. abs(y) < 2d-4) .or. &
                  (abs(x) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=2 ! twofold axis
          xs(2)=y;  ys(2)=-x
        else if (abs(y-x+0.5d0) < 2d-4 .or. abs(y-x-0.5d0) < 2d-4 .or. &
                 abs(x+y-0.5d0) < 2d-4 .or. abs(x+y-1.5d0) < 2d-4) then
          ndeg=4 ! mirror plane
          xs(2)=-y; ys(2)=x
          xs(3)=-x; ys(3)=-y
          xs(4)=y;  ys(4)=-x
        else
          ndeg=8
          xs(2)=-y; ys(2)=x
          xs(3)=-x; ys(3)=-y
          xs(4)=y;  ys(4)=-x
          xs(5)=0.5d0-y; ys(5)=0.5d0-x
          xs(6)=0.5d0+x; ys(6)=0.5d0-y
          xs(7)=0.5d0+y; ys(7)=0.5d0+x
          xs(8)=0.5d0-x; ys(8)=0.5d0+y
        endif
        return
        end
!**********************************************************
!       p3, nsg=13
!**********************************************************
        subroutine p3(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(3),ys(3)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .and. abs(y) < 2d-4) .or. &
             (abs(x-1d0/3d0) < 2d-4 .and. abs(y-2d0/3d0) < 2d-4) .or. &
             (abs(x-2d0/3d0) < 2d-4 .and. abs(y-1d0/3d0) < 2d-4) ) then
          ndeg=1 ! threefold axis
        else
          ndeg=3
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
        endif
        return
        end
!**********************************************************
!       p3m1, nsg=14
!**********************************************************
        subroutine p3m1(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(6),ys(6)

        x=xs(1); y=ys(1)
        if ( (abs(x) < 2d-4 .and. abs(y) < 2d-4) .or. &
             (abs(x-1d0/3d0) < 2d-4 .and. abs(y-2d0/3d0) < 2d-4) .or. &
             (abs(x-2d0/3d0) < 2d-4 .and. abs(y-1d0/3d0) < 2d-4) ) then
          ndeg=1 ! threefold axis
        else if ( abs(x+x-y) < 2d-4 .or. abs(y+y-x) < 2d-4 &
             .or. abs(x+x-y-1d0) < 2d-4 .or. abs(y+y-x-1d0) < 2d-4 &
             .or. abs(x+y-1d0) < 2d-4 ) then
          ndeg=3 ! mirror plane
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
        else
          ndeg=6
          xs(2)=x;   ys(2)=x-y
          xs(3)=-y;  ys(3)=x-y
          xs(4)=y-x; ys(4)=y
          xs(5)=y-x; ys(5)=-x
          xs(6)=-y;  ys(6)=-x
        endif
        return
        end
!**********************************************************
!       p31m, nsg=15
!**********************************************************
        subroutine p31m(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(6),ys(6)

        x=xs(1); y=ys(1)
        if (abs(x) < 2d-4 .and. abs(y) < 2d-4) then
          ndeg=1 ! threefold axis at (0,0)
        else if ( (abs(x-1d0/3d0) < 2d-4 .and. abs(y-2d0/3d0) < 2d-4) &
             .or. (abs(x-2d0/3d0) < 2d-4 .and. abs(y-1d0/3d0) < 2d-4) ) then
          ndeg=2 ! threefold axis
          xs(2)=y;   ys(2)=x
        else if (abs(x) < 2d-4 .or. abs(y) < 2d-4 .or. abs(x-y) < 2d-4) then
          ndeg=3 ! mirror plane
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
        else
          ndeg=6
          xs(2)=y;   ys(2)=x
          xs(3)=-y;  ys(3)=x-y
          xs(4)=-x;  ys(4)=y-x
          xs(5)=y-x; ys(5)=-x
          xs(6)=x-y; ys(6)=-y
        endif
        return
        end
!**********************************************************
!       p6, nsg=16
!**********************************************************
        subroutine p6(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(6),ys(6)

        x=xs(1); y=ys(1)
        if (abs(x) < 2d-4 .and. abs(y) < 2d-4) then
          ndeg=1 ! sixfold axis
        else if ( (abs(x-1d0/3d0) < 2d-4 .and. abs(y-2d0/3d0) < 2d-4) &
             .or. (abs(x-2d0/3d0) < 2d-4 .and. abs(y-1d0/3d0) < 2d-4) ) then
          ndeg=2 ! threefold axis
          xs(2)=y;   ys(2)=x
        else if ( (abs(x-0.5d0) < 2d-4 .and. abs(y) < 2d-4) .or. &
                  (abs(x) < 2d-4 .and. abs(y-0.5d0) < 2d-4) .or. &
                  (abs(x-0.5d0) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=3 ! twofold axis
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
        else
          ndeg=6
          xs(2)=x-y; ys(2)=x
          xs(3)=-y;  ys(3)=x-y
          xs(4)=-x;  ys(4)=-y
          xs(5)=y-x; ys(5)=-x
          xs(6)=y;   ys(6)=y-x
        endif
        return
        end
!**********************************************************
!       p6mm, nsg=17
!**********************************************************
        subroutine p6mm(xs,ys,ndeg)
        implicit none
        integer :: ndeg
        real(8) :: x,y, xs(12),ys(12)

        x=xs(1); y=ys(1)
        if (abs(x) < 2d-4 .and. abs(y) < 2d-4) then
          ndeg=1 ! sixfold axis
        else if ( (abs(x-1d0/3d0) < 2d-4 .and. abs(y-2d0/3d0) < 2d-4) &
             .or. (abs(x-2d0/3d0) < 2d-4 .and. abs(y-1d0/3d0) < 2d-4) ) then
          ndeg=2 ! threefold axis
          xs(2)=y;   ys(2)=x
        else if ( (abs(x-0.5d0) < 2d-4 .and. abs(y) < 2d-4) .or. &
                  (abs(x) < 2d-4 .and. abs(y-0.5d0) < 2d-4) .or. &
                  (abs(x-0.5d0) < 2d-4 .and. abs(y-0.5d0) < 2d-4) ) then
          ndeg=3 ! twofold axis
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
        else if (abs(x) < 2d-4 .or. abs(y) < 2d-4 .or. abs(x-y) < 2d-4) then
          ndeg=6 ! mirror plane
          xs(2)=-y;  ys(2)=x-y
          xs(3)=y-x; ys(3)=-x
          xs(4)=-x;  ys(4)=-y
          xs(5)=y;   ys(5)=y-x
          xs(6)=x-y; ys(6)=x
        else
          ndeg=12
          xs(2)=x;    ys(2)=x-y
          xs(3)=x-y;  ys(3)=x
          xs(4)=y;    ys(4)=x
          xs(5)=-y;   ys(5)=x-y
          xs(6)=y-x;  ys(6)=y
          xs(7)=-x;   ys(7)=-y
          xs(8)=-x;   ys(8)=y-x
          xs(9)=y-x;  ys(9)=-x
          xs(10)=-y;  ys(10)=-x
          xs(11)=y;   ys(11)=y-x
          xs(12)=x-y; ys(12)=-y
        endif
        return
        end
