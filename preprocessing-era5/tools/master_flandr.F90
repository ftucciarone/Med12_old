PROGRAM MASTER

IMPLICIT NONE

INTEGER :: im, jm, km
INTEGER :: jt
INTEGER :: ist,ien,jst,jen,kst,ken,ntim,ndim
REAL, DIMENSION(:,:), ALLOCATABLE :: msk
REAL, DIMENSION(:,:), ALLOCATABLE :: mskt
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZF
LOGICAL :: ll_mv

CHARACTER(LEN=99) :: cmfile,cdfile,ctmask,&
& cdx,cdy,cdz,cfin,cvar,caux,cmvar

CALL GETARG(1,cfin)
CALL GETARG(2,cvar)
CALL GETARG(3,caux)
READ(caux,*) im
CALL GETARG(4,caux)
READ(caux,*) jm
CALL GETARG(5,caux)
READ(caux,*) ntim
CALL GETARG(6,cmvar)

ll_mv=.false.
IF(TRIM(cmvar).EQ.'MISSING_VALUE') ll_mv=.true.
WRITE(*,*) ' MISSING VALUE :',ll_mv

WRITE(*,*) ' DIMENSION     :',im,jm,ntim

ist=1 ; ien=im
jst=1 ; jen=jm

WRITE(*,*) ' Allocating field'
ALLOCATE(ZF(im,jm,ntim))
ALLOCATE(MSK(im,jm))
ALLOCATE(MSKT(im,jm))

WRITE(*,*) ' Reading mask'
if(.not. ll_mv) Call Read_Msk('lsm.nc',cmvar,im,jm,1,msk)

WRITE(*,*) ' Number of ocean points : ',100*SUM(msk)/(im*jm)

WRITE(*,*) ' Reading field'
Call Read_Field2dt(cfin,cvar,ist,ien,jst,jen,ntim,ZF)

!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(jt,MSKT)
!$OMP DO SCHEDULE(DYNAMIC,1)
DO jt=1,ntim
    WRITE(*,*) ' Doing time ',jt
    IF(ll_mv) then
       MSKT=1.
       WHERE(.NOT. ABS( ZF(:,:,jt) ) .LT. 998. ) 
         MSKT=0. 
         ZF(:,:,jt)=0.
       ENDWHERE
    ELSE
       MSKT = MSK
    ENDIF
    Call EXTRAP(ZF(:,:,jt),MSKT(:,:),IM,JM)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

Call Write_Field2dt(cfin,cvar,im,jm,ntim,ZF)

DEALLOCATE(ZF,MSK,MSKT)

WRITE(*,*) ' DONE '

END PROGRAM MASTER
