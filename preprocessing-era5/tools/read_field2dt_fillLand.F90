subroutine Read_Field2dt(cfile,cvar1,i1,i2,j1,j2,nt,zff)
Use netcdf
implicit none
integer,intent(in) :: i1,i2,j1,j2,nt
character(len=*)  :: cfile
character(len=*)  :: cvar1
real,intent(out)   :: zff(i2-i1+1,j2-j1+1,nt)
integer :: stat,nc,id,ia1,ia2,ia3
ia1=i2-i1+1
ia2=j2-j1+1
stat = nf90_open(Cfile,NF90_NOWRITE,nc)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_inq_varid (nc, cvar1, id)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_get_var (nc, id, zff,start=(/i1,j1,1/),&
& count=(/ia1,ia2,nt/))
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_close (nc)
if (stat /= nf90_noerr) call handleerr(stat)
End subroutine Read_Field2dt
