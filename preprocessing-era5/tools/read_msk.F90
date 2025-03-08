subroutine Read_Msk(cfile,cvar1,ii,jj,kk,msk)

Use netcdf

implicit none

integer,intent(in) :: ii,jj,kk
character(len=*)  :: cvar1,cfile
real,intent(out)   :: msk(ii,jj,kk)
integer :: stat,nc,id

stat = nf90_open(Cfile,NF90_NOWRITE,nc)
if (stat /= nf90_noerr) call handleerr(stat)

stat = nf90_inq_varid (nc, cvar1, id)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_get_var (nc, id, msk)
if (stat /= nf90_noerr) call handleerr(stat)

stat = nf90_close(nc)
if (stat /= nf90_noerr) call handleerr(stat)

End subroutine Read_Msk
