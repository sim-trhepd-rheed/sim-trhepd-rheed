!*******************************************************************
!   rheed multi slice method
!   bulk+surf input file --> xyz file
!   v.1:2014/4    T.Hanada
!*******************************************************************
        implicit none
        real(8), dimension(:),allocatable :: ocr,x,y,z,ocrs,xs,ys,zs
        real(8), dimension(:),allocatable :: oocr,xx,yy,zz
        integer, dimension(:),allocatable :: ielm,ielms,iz,izs,nb
        real(8) :: xb(12),yb(12),aa,bb,gam,cc,dx,dy,dthick,dxs,dys,xbl,ybl
        real(8) :: dummy,zmin,zmax,zoffset,bbx,bby,pxa,pxb,pya,pyb, socr,r2
        integer :: nh,nk,ndom,nelm,nsg,natm,nelms,nsgs,msa,msb,nsa,nsb,natms
        integer :: iabr,iabp
        integer :: nall,iall,natoms,natomb,ndeg, ix,iy,nx,ny,l,nbl, idotb,idots,i,j
        character sname*20,bname*20,fname*48,as*2
!        write (*,'(A)') ' sequence-filename (console=return) ? '
!        read (*,'(A)') fname
!        if (fname /= ' ') open (5,file=fname)
!----------file open-----------
      do
!----------input bulk data (3)----------
        write (*,'(A)') ' bulk-input-filename (end=e) ? :'
        read (*,'(A)') bname
        write (*,'(" ",A)') bname
        if (bname == 'e' .or. bname == 'E') stop
        open (3,file=bname,status='old')

        read (3,*) nh,nk,ndom
allocate (nb(ndom))
        read (3,*) (nb(i),i=1,ndom)
        read (3,*) (dummy,i=1,ndom) ! rdom
        do i=1,ndom
          read (3,*) (ix,iy,j=1,nb(i)) ! ih,ik
        end do
        read (3,*) dummy !be,azi,azf,daz,gi,gf,dg
        read (3,*) dummy !tem,dz,ml,epsb
! atomic parameters
        read (3,*) nelm
allocate (iz(nelm))
        do i=1,nelm
          read (3,*) iz(i) !,da1,sap;  sae=0d0
          read (3,*) dummy !dtz(i),dth(i),dtk(i)
       end do
! structural parameters
        read (3,*) nsg,aa,bb,gam,cc,dx,dy
        read (3,*) natm
allocate (ielm(natm)); allocate (ocr(natm))
allocate (x(natm)); allocate (y(natm)); allocate (z(natm))
        zmin=1d10; zmax=-1d10
        do i=1,natm
          read (3,*) ielm(i),ocr(i),x(i),y(i),z(i)
          if (z(i) < zmin) zmin=z(i)
          if (z(i) > zmax) zmax=z(i)
        end do
        zoffset=0.5d0*(cc-zmax-zmin)
        do i=1,natm
          z(i)=z(i)+zoffset
        end do
        close (3)
!----------input surface data (2)----------
        write (*,'(A)') ' surface-input-filename ? :'
        read (*,'(A)') sname
        write (*,'(" ",A)') sname
        open (2,file=sname,status='old')

! surface atomic parameters
        read (2,*) nelms
allocate (izs(nelms))
        do i=1,nelms
          read (2,*) izs(i) !,da1,sap;  sae=0d0
          read (2,*) dummy !dtz(i),dth(i),dtk(i)
        end do
! surface structural parameters
        read (2,*) nsgs,msa,msb,nsa,nsb,dthick,dxs,dys
        read (2,*) natms
allocate (ielms(natms)); allocate (ocrs(natms))
allocate (xs(natms)); allocate (ys(natms)); allocate (zs(natms))
        do i=1,natms
          read (2,*) ielms(i),ocrs(i),xs(i),ys(i),zs(i)
          zs(i)=zs(i)+zoffset
        end do
        close (2)
!----------output xyz data (1)----------
        idotb=scan(bname,".",BACK=.true.)
        if (idotb > 0) then
          idotb=idotb-1
        else
          idotb=LEN_TRIM(bname)
        endif
        idots=scan(sname,".",BACK=.true.)
        if (idots > 0) then
          idots=idots-1
        else
          idots=LEN_TRIM(sname)
        endif
        fname=sname(:idots)//'-'//bname(:idotb)//'.xyz'

        write (*,'(" ",A)') fname
        open (1,file=fname)
!---- count atoms
        natoms=0
        do i=1,natms
          call extendxyz(nsgs,msa,msb,nsa,nsb,xs(i),ys(i),ndeg,xb,yb)
          natoms=natoms+ndeg
        end do

        natomb=0
        do i=1,natm
          call extendxyz(nsg,1,0,0,1,x(i),y(i),ndeg,xb,yb)
          natomb=natomb+ndeg
        end do
!---- output xyz
        write (*,'(A)') ' # of a-units, b-units & bulk layers ? :'
        read (*,*) nx,ny,nbl
        write (*,*) nx,ny,nbl

        nall=(natoms+natomb*nbl)*nx*ny
        write (1,'(I5)') nall
        write (1,'(A," / ",A)') sname,bname

        bbx=bb*cos(gam*atan(1d0)/45d0)
        bby=bb*sin(gam*atan(1d0)/45d0)

allocate (xx(nall)); allocate (yy(nall)); allocate (zz(nall))
allocate (oocr(nall))
        iall=0
        do i=1,natms
          call extendxyz(nsgs,msa,msb,nsa,nsb,xs(i),ys(i),ndeg,xb,yb)
          call atomicsymbol(izs(ielms(i)),as)
          do ix=0,nx-1
            pxa=dble(msa*ix); pxb=dble(msb*ix)
            do iy=0,ny-1
              pya=pxa+dble(nsa*iy); pyb=pxb+dble(nsb*iy)
              do j=1,ndeg
                iall=iall+1
                xx(iall)=(xb(j)+pya)*aa+(yb(j)+pyb)*bbx
                yy(iall)=(yb(j)+pyb)*bby
                zz(iall)=zs(i)
                oocr(iall)=ocr(i)
                write (1,'(A,3F12.5)') as,xx(iall),yy(iall),zz(iall)
              end do
            end do
          end do
        end do

        do l=1,nbl
          do i=1,natm
            xbl=x(i)-dxs-dx*(l-1)
            ybl=y(i)-dys-dy*(l-1)
            call extendxyz(nsg,1,0,0,1,xbl,ybl,ndeg,xb,yb)
            call atomicsymbol(iz(ielm(i)),as)
            do ix=0,nx-1
              pxa=dble(ix)
              do iy=0,ny-1
                pyb=dble(iy)
                do j=1,ndeg
                  iall=iall+1
                  xx(iall)=(xb(j)+pxa)*aa+(yb(j)+pyb)*bbx
                  yy(iall)=(yb(j)+pyb)*bby
                  zz(iall)=z(i)-cc*l
                  oocr(iall)=ocr(i)
                  write (1,'(A,3F12.5)') as,xx(iall),yy(iall),zz(iall)
                end do
              end do
            end do
          end do
        end do
        close (1)
!---- find overlap
        do i=1,nall
          socr=oocr(i)
          do j=i+1,nall
            r2=(xx(i)-xx(j))**2+(yy(i)-yy(j))**2+(zz(i)-zz(j))**2
            if (r2 < 0.16d0) then ! distance < 0.4 Angstrom
              socr=socr+oocr(j)
            endif
          end do
          if (socr-1d0 > 1d-5) then
            write (*,*) ' sum of ocr > 1 at ',xx(i),yy(i),zz(i)
          endif
        end do

deallocate (nb)
deallocate (iz)
deallocate (ielm); deallocate (ocr)
deallocate (x); deallocate (y); deallocate (z)
deallocate (izs)
deallocate (ielms); deallocate (ocrs)
deallocate (xs); deallocate (ys); deallocate (zs)
deallocate (xx); deallocate (yy); deallocate (zz); deallocate (oocr)
      end do
      end
!**********************************************************
        subroutine extendxyz(nsg,ma,mb,na,nb,x,y,ndeg,xb,yb)
        implicit none
        integer :: nsg,ma,mb,na,nb,ndeg, i
        real(8) :: x,y, reduce01, s,xs(12),ys(12),xb(12),yb(12)

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
        end do
        return
        end
