Subroutine handleerr(ist)
Use Netcdf
implicit none
integer, intent(in) :: ist
if(ist /= nf90_noerr) then
   print *, 'Error: ', trim(nf90_strerror(ist))
   call abort
endif
End subroutine handleerr
