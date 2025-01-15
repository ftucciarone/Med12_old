subroutine Read_Field2dt(cfile,cvar1,i,j,nt,zff)

Use netcdf

implicit none

integer,intent(in) :: i,j,nt
character(len=*)  :: cfile
character(len=*)  :: cvar1
real,intent(out)   :: zff(i,j,nt)

integer :: stat,nc,id

stat = nf90_open(Cfile,NF90_NOWRITE,nc)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_inq_varid (nc, cvar1, id)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_get_var (nc, id, zff,start=(/1,1,1/),&
& count=(/i,j,nt/))
if (stat /= nf90_noerr) call handleerr(stat)

End subroutine Read_Field2dt
