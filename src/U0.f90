!*******************************************************************
!   rheed multi slice method : U0
!   v.1:2019/5    T.Hanada
!*******************************************************************
        program U0
        implicit none
        integer :: inegpos, idotb,idots, nubulk
        character fname*60,bname*24,sname*24,ep*1 !,dia*1 ! ,yn*1
!---------------------
!       write (*,'(A)') ' 0:electron 1:positron ? '
!       read (*,*) inegpos
!       if (inegpos <= 0) then
!         inegpos=0; ep='E'
!       else
!         inegpos=1; ep='P'
!       endif
!       write (*,*) ep
        ep='P' ! ep='E' for electron
        write (*,'(A)') ' # of bulk-unit repetition ? '
        read (*,*) nubulk
        if (nubulk < 1) nubulk=1
!----------file open-----------
!     do
!       write (*,'(A)') ' bulk-input-filename (end=e) ? :'
!       read (*,'(A)') bname
!       write (*,'(" ",A)') bname
!       if (bname == 'e' .or. bname == 'E') stop
        bname = 'bulk.txt'
         
        idotb=scan(bname,".",BACK=.true.)
        if (idotb > 0) then
          idotb=idotb-1
        else
          idotb=LEN_TRIM(bname)
        endif

!       write (*,'(A)') ' surface-structure-filename (end=e) ? :'
!       read (*,'(A)') sname
!       write (*,'(" ",A)') sname
!       if (sname == 'e' .or. sname == 'E') exit
        sname = 'surf.txt'
!
        idots=scan(sname,".",BACK=.true.)
        if (idots > 0) then
          idots=idots-1
        else
          idots=LEN_TRIM(sname)
        endif

!       write (*,'(A)') ' output-filename :'
!        read (*,'(A)') fname
          fname=sname(:idots)//'-'//bname(:idotb)//'-U0.txt'
        write (*,'(" ",A)') fname
!----------main routine----------
        open (4,file=fname)

        open (3,file=bname,status='old')
        open (1,file='temporary.b',form='unformatted')
        call bulkioU(inegpos,nubulk-1)
! one bulk unit is added below surface layer
        close (3)

        rewind 1
        open (2,file=sname,status='old')
        call surfioU
        close (1)
        close (2)

        close (4)
!     end do
      end
