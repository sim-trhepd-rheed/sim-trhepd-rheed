!***********************************************************
!       complex matrix  a*inv(b) --> a
!     v.1:84/10  v.4 2014/2/26  T.Hanada
!***********************************************************
      subroutine gcmi2(n1,n,a,b)
      implicit none
      integer :: n1,n, i,j,k,m,j1
      complex(8) :: a(n1,n),b(n1,n),c
      real(8) :: s,ab
!---------- scaling
      do j=1,n
        s=-1d0
        do i=1,n
          ab=abs(dble(b(i,j)))+abs(dimag(b(i,j)))
          if (ab > s) then
            s=ab
            m=i
          endif
        end do
        c=(1d0,0d0)/b(m,j)
        a(1:n,j)=a(1:n,j)*c
        b(1:n,j)=b(1:n,j)*c
      end do
!---------- gauss elimination
!----- partial pivotting
      do j=1,n
        s=-1d0
        do k=j,n
          ab=abs(dble(b(j,k)))+abs(dimag(b(j,k)))
          if (ab > s) then
            s=ab
            m=k
          endif
        end do
        if (m /= j) then
          do i=j,n
            c=b(i,j)
            b(i,j)=b(i,m)
            b(i,m)=c
          end do
          do i=1,n
            c=a(i,j)
            a(i,j)=a(i,m)
            a(i,m)=c
          end do
        endif
!----- elimination
        j1=j+1
        c=(1d0,0d0)/b(j,j)
        a(1:n,j)=a(1:n,j)*c
        b(j1:n,j)=b(j1:n,j)*c
        do k=j1,n
          c=-b(j,k)
          a(1:n,k)=a(1:n,k)+a(1:n,j)*c
          b(j1:n,k)=b(j1:n,k)+b(j1:n,j)*c
        end do
      end do
!---------- inverse
      do j=n-1,1,-1
        do k=j+1,n
          a(1:n,j)=a(1:n,j)-a(1:n,k)*b(k,j)
        end do
      end do
      return
      end
