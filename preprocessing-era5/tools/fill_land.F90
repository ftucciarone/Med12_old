subroutine extrap (a, land, il, jl)

      implicit none

      real, dimension(il,jl),intent(inout) :: a, land
      integer, intent(in) :: il, jl

      real, parameter :: crit=1.e-5
      integer, parameter :: maxscn = 500, gtype= 1

      real, parameter :: c0=0.0, p25=0.25

      real, dimension(il,jl) :: res, sor
      logical :: done

      integer :: i,j,n
      real :: relc,resmax,absres

!
!-----------------------------------------------------------------------
!     set the relaxation coefficient to zero over ocean or air
!     relc is somewhat arbitrary
!-----------------------------------------------------------------------
!
      relc = 0.6
      do j=1,jl
        do i=1,il
          if (land(i,j) .lt. 0.5) then
            sor(i,j) = relc
          else
            sor(i,j) = c0
          endif
        enddo
      enddo

!
!-----------------------------------------------------------------------
!     iterate until errors are acceptable.
!-----------------------------------------------------------------------
!     
      n = 0
100   continue
        resmax = c0
        done   = .true.
        n    = n + 1
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = p25*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) &
     &                 - a(i,j)
          enddo
        enddo
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = res(i,j)*sor(i,j)
            a(i,j) = a(i,j) + res(i,j)
            absres = abs(res(i,j))
            if (absres .gt. crit) done = .false.
            resmax = max(absres,resmax)
          enddo
        enddo
!
!-----------------------------------------------------------------------
!       set conditions at edge of grid
!-----------------------------------------------------------------------
!
        if (gtype .eq. 1) then
!
!         use cyclic or no flux conditions on ocean grids
!
          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)

          enddo
        elseif (gtype .eq. 2) then

!         always put cyclic conditions on atmosphere grids

          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
          enddo
        endif

!       no flux condition at northern and southern boundaries

        do i=1,il
          a(i,1)  = a(i,2)
          a(i,jl) = a(i,jl-1)
          enddo

      if (.not. done .and. n .le. maxscn) go to 100
   
      return
      
end subroutine extrap
