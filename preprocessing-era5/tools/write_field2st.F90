subroutine Write_Field2st(cfile,cvar1,i,j,nt,zff,nc,id)

Use netcdf

implicit none

integer,intent(in) :: i,j,nt
integer,intent(inout) :: nc,id
character(len=*)  :: cfile
character(len=*)  :: cvar1
real,intent(out)   :: zff(i,j)

integer :: stat

if ( nc .eq. -999 .and. id .eq. -999) THEN
  stat = nf90_open(Cfile,NF90_WRITE,nc)
  if (stat /= nf90_noerr) call handleerr(stat)
  stat = nf90_inq_varid (nc, cvar1, id)
  if (stat /= nf90_noerr) call handleerr(stat)
endif
if ( nc .ne. -999 .and. id .eq. -999) THEN
  stat = nf90_close(nc)
  if (stat /= nf90_noerr) call handleerr(stat)
  return
endif

stat = nf90_put_var (nc, id, zff,start=(/1,1,nt/),&
& count=(/i,j,1/))
if (stat /= nf90_noerr) call handleerr(stat)

End subroutine Write_Field2st
