!*******************************************************************
!   rheed multi slice method : surface
!   subroutine surfio
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4    T.Hanada
!*******************************************************************
        program surf
        integer :: iprn,idotb,idots,ndom
        character bname*20,sname*20,fname*44,dname*46,yn*1,ep*1
        integer :: date_time(8)
        character*10 cdate,ctime,czone
!        write (*,'(A)') ' sequence-filename (console=return) ? '
!        read (*,'(A)') fname
!        if (fname /= ' ') open (5,file=fname)
!---------------------
      ep='P' ! 'P' for positoron, 'E' for electron
!
      iprn=0
      call date_and_time(cdate,ctime,czone, date_time)
      open (4,file='SURF'//cdate(:8)//'-'//ctime(:6)//'log.txt')
!----------beam loop-----------
    do
!      write (*,'(A)') ' bulk-filename (end=e) ? :'
!      read (*,'(A)') bname
      bname = 'bulkP.b' ! 'bulkP.b' for positron, 'bulkE.b' for electron
      if (ep == 'E') bname = 'bulkE.b'
!
!     write (*,'(" ",A)') bname
!     if (bname == 'e' .or. bname == 'E') then
!       close (4)
!       stop
!     endif
!
      open (1,file=bname,status='old',form='unformatted')
      read (1) ndom

      idotb=scan(bname,".",BACK=.true.)
      if (idotb > 0) then
        idotb=idotb-1
      else
        idotb=LEN_TRIM(bname)
      endif
!----------structure loop----------
      do
!        write (*,'(A)') ' structure-filename (end=e) ? :'
!        read (*,'(A)') sname
        sname = 'surf.txt'
!       write (*,'(" ",A)') sname
        if (sname == 'e' .or. sname == 'E') then
          close (1)
          exit
        endif
        open (2,file=sname,status='old')
        idots=scan(sname,".",BACK=.true.)
        if (idots > 0) then
          idots=idots-1
        else
          idots=LEN_TRIM(sname)
        endif

!        write (*,'(A)') ' output-filename :'
!        read (*,'(A)') fname
        if (ndom > 1) then
          fname=sname(:idots)//'-'//bname(:idotb)//'.md'
        else
          fname=sname(:idots)//'-'//bname(:idotb)//'.s'
        endif
!       write (*,'(" ",A)') fname
        open (3,file=fname)

        write (4,'(A)') '----------------------------------'
        write (4,'(2A)') 'BULK output: ',bname
        write (4,'(2A)') 'SURF input : ',sname
        write (4,'(2A)') 'SURF output: ',fname
!----------main routine----------
        rewind 1
        call surfio(iprn)
        close (2)
        close (3)
!----------domain sum----------
        if (ndom > 1) then
          open (2,file=fname)
          dname=sname(:idots)//'-'//bname(:idotb)//'.s'
          write (*,'(" ",A)') dname
          open (3,file=dname)
          write (4,'(2A)') 'domainsum output: ',dname
          call domainsum(ndom)
        endif
        goto 1000
!        stop!#####altered#####
      end do
    end do
    1000 continue
    end
