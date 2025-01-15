PROGRAM td2qs

IMPLICIT NONE

REAL, PARAMETER :: c1es=610.78, &
                   c2es=17.269, &
                   c2is=21.875, &
                   c3es=273.16, &
                   c4es=35.86,  &
                   c4is=7.66

REAL, ALLOCATABLE, DIMENSION(:,:) :: TD, PS, QS

INTEGER :: NX, NY, NT, NYY
INTEGER :: JX, JY, JT
INTEGER :: NC1, ID1, NC2, ID2, NCO, IDO

CHARACTER(LEN=100) :: CFTD, CFPS, CFQS
CHARACTER(LEN=100) :: CVTD, CVPS, CVQS, CAUX

LOGICAL :: LL_CELSIUS, LL_HPA

REAL :: QSW, QSI

LL_CELSIUS = .FALSE.
LL_HPA = .FALSE.

CALL GETARG(1,CFTD)
CALL GETARG(2,CVTD)
CALL GETARG(3,CFPS)
CALL GETARG(4,CVPS)
CALL GETARG(5,CFQS)
CALL GETARG(6,CVQS)
CALL GETARG(7,CAUX)  ; READ(CAUX,*) NX
CALL GETARG(8,CAUX)  ; READ(CAUX,*) NY
CALL GETARG(9,CAUX)  ; READ(CAUX,*) NYY
CALL GETARG(10,CAUX) ; IF(CAUX(1:1).EQ.'T') LL_CELSIUS = .TRUE.
CALL GETARG(11,CAUX) ; IF(CAUX(1:1).EQ.'T') LL_HPA = .TRUE.

IF( MOD(NYY,4) .EQ. 0 .AND. NYY .NE. 1900) THEN
    NT=366*4
ELSE
    NT=365*4
ENDIF

WRITE(*,*) ' DIMENSIONS ARE ', NX, NY, NT
WRITE(*,*) ' DEWPOINT IN CELSIUS ', LL_CELSIUS
WRITE(*,*) ' PRESSURE IN HPA     ', LL_HPA

ALLOCATE( TD( NX, NY),&
          PS( NX, NY),&
          QS( NX, NY) )

IF( LL_CELSIUS ) TD = TD + 273.15
IF( LL_HPA     ) PS = PS * 100.

nc1 = -999 ; id1 = -999
nc2 = -999 ; id2 = -999
nco = -999 ; ido = -999

DO JT = 1, NT
  WRITE(*,*) ' DOING ', JT, ' over ', NT
  CALL Read_Field2st(cftd,cvtd,nx,ny,jt,td,nc1,id1)
  CALL Read_Field2st(cfps,cvps,nx,ny,jt,ps,nc2,id2)
  DO JY = 1,NY
    DO JX = 1,NX
      QSW = 0.622 * C1ES * EXP ( C2ES * (TD(JX,JY)-C3ES) / (TD(JX,JY)-C4ES) ) / (PS(JX,JY))
      QSI = 0.622 * C1ES * EXP ( C2IS * (TD(JX,JY)-C3ES) / (TD(JX,JY)-C4IS) ) / (PS(JX,JY))
      QS(JX,JY) = MIN( QSW, QSI )
    ENDDO
  ENDDO
  CALL Write_Field2st(cfqs,cvqs,nx,ny,jt,qs,nco,ido)
ENDDO

CALL Read_Field2st(cftd,cvtd,nx,ny,jt,td,nc1,-999)
CALL Read_Field2st(cfps,cvps,nx,ny,jt,ps,nc2,-999)
CALL Write_Field2st(cfqs,cvqs,nx,ny,jt,qs,nco,-999)

END PROGRAM td2qs
