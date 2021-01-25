!*******************************************************************
!   rheed multi slice method 
!   weighted sum of intensities from multi-domains
!     AA = BB, GAM = 90 or 120
!   v.1:2014/4  v.1b:2017/3  T.Hanada
!*******************************************************************
        implicit none
        character fname*48
!        write (*,'(A)') ' sequence-filename (console=return) ? '
!        read (*,'(A)') fname
!        if (fname /= ' ') open (5,file=fname)
!----------
    do
      write (*,'(A)') ' .md input-filename (end=e) ? :'
      read (*,'(A)') fname
      write (*,'(" ",A)') fname
      if (fname == 'e' .or. fname == 'E') then
        stop
      endif
      open (2,file=fname)

      write (*,'(A)') ' .s output-filename (end=e) ? :'
      read (*,'(A)') fname
      write (*,'(" ",A)') fname
      if (fname == 'e' .or. fname == 'E') then
        stop
      endif
      open (3,file=fname)

      call domainsum(0)
      close (2)
      close (3)
    end do
    end
