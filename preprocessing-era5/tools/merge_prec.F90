PROGRAM merge

Use netcdf
IMPLICIT NONE

integer :: stat,nc1,id1,nc2,id2,nco,ido

REAL, ALLOCATABLE, DIMENSION(:,:) :: F, F2

INTEGER :: NX, NY, NYY, ND
INTEGER :: JD, jk1, jko
REAL :: zcst

CHARACTER(LEN=100) :: CF1, CFO
CHARACTER(LEN=100) :: CV1, CVO
CHARACTER(LEN=4) :: CAUX

CALL GETARG(1,CF1)
CALL GETARG(2,CV1)
CALL GETARG(3,CFO)
CALL GETARG(4,CVO)
CALL GETARG(5,CAUX)  ; READ(CAUX,*) NX
CALL GETARG(6,CAUX)  ; READ(CAUX,*) NY
CALL GETARG(7,CAUX)  ; READ(CAUX,*) NYY

IF( IARGC() .EQ. 8 ) THEN
   CALL GETARG(8,CAUX)  ; READ(CAUX,*) ND
ELSEIF( MOD( NYY,4) .EQ. 0 ) THEN
   ND = 366
ELSE
   ND = 365
ENDIF

WRITE(*,*) ' DIMENSIONS ARE ', NX, NY, ND

ALLOCATE( F( NX, NY) )
ALLOCATE( F2( NX, NY) )

stat = nf90_open(Cf1,NF90_NOWRITE,nc1)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_inq_varid (nc1, cv1, id1)
if (stat /= nf90_noerr) call handleerr(stat)

stat = nf90_open(Cfo,NF90_WRITE,nco)
if (stat /= nf90_noerr) call handleerr(stat)
stat = nf90_inq_varid (nco, cvo, ido)
if (stat /= nf90_noerr) call handleerr(stat)

jko=0
jk1=0
zcst=1000./86400.

DO JD=1,ND

  WRITE(*,*) ' DOING DAY ', JD
  jk1=jk1+1
  stat = nf90_get_var (nc1, id1, f,start=(/1,1,jk1/),&
  & count=(/nx,ny,1/))
  if (stat /= nf90_noerr) call handleerr(stat)
  jk1=jk1+1
  stat = nf90_get_var (nc1, id1, f2,start=(/1,1,jk1/),&
  & count=(/nx,ny,1/))
  if (stat /= nf90_noerr) call handleerr(stat)

  jko=jko+1
  stat = nf90_put_var (nco, ido, zcst*(f+f2),start=(/1,1,jko/),&
  & count=(/nx,ny,1/))
  if (stat /= nf90_noerr) call handleerr(stat)

ENDDO

  stat = nf90_close (nc1)
  if (stat /= nf90_noerr) call handleerr(stat)
  stat = nf90_close (nco)
  if (stat /= nf90_noerr) call handleerr(stat)

END PROGRAM MERGE
