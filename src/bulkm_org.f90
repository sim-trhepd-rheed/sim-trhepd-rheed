!*******************************************************************
!   rheed multi slice method : bulk
!   subroutine bulkio
!   v.1:84/10   v.2:86/11   v.3:90/4   v.4:2014/4    T.Hanada
!*******************************************************************
program bulk
        implicit none
        integer :: inegpos, idiag, iprn, idot
        character fname*20,bname*24,ep*1 !,dia*1 ! ,yn*1
        real(4) :: etime,etime1,etime2(2)
!---------------------
        write (*,'(A)') ' 0:electron 1:positron ? '
        read (*,*) inegpos
        if (inegpos <= 0) then
          inegpos=0; ep='E'
        else
          inegpos=1; ep='P'
        endif
        write (*,*) ep
        idiag=3
        iprn=0
!----------file open-----------
      etime1=etime(etime2)
      do
        write (*,'(A)') ' input-filename (end=e) ? :'
        read (*,'(A)') fname
        write (*,'(" ",A)') fname
        if (fname == 'e' .or. fname == 'E') stop
        open (3,file=fname,status='old')

        write (*,'(A)') ' output-filename :'
!        read (*,'(A)') bname
        idot=scan(fname,".",BACK=.true.)
        if (idot > 0) then
          idot=idot-1
        else
          idot=LEN_TRIM(fname)
        endif
        bname=fname(:idot)//ep//'.b'
        write (*,'(" ",A)') bname
        open (1,file=bname,form='unformatted')
!----------main routine----------
        call bulkio(inegpos,idiag,iprn)
        close (1)
        close (3)
        etime1=etime(etime2)
        print *, 'elapsed(sec):', etime1, ', user:', etime2(1), ', sys:', etime2(2)
      end do
end program bulk
