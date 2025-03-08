subroutine Write_Field2dt(cfile,cvar1,i1,j1,nt,zff)

Use netcdf

implicit none

integer,intent(in) :: i1,j1,nt
character(len=*)  :: cfile
character(len=*)  :: cvar1
real,intent(in)   :: zff(i1,j1,nt)

integer :: stat,nc,id

stat = nf90_open(Cfile,NF90_WRITE,nc)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_inq_varid (nc, cvar1, id)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_put_var (nc, id, zff)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_close(nc)

End subroutine Write_Field2dt
